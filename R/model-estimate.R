# Estimating the candidates --------------------------------------------------
#
# The only stage that reads patient data and the only one that needs `nlmixr2`.
# The candidate set is exactly the models `.pk_single_dose()` can evaluate in
# closed form: a candidate the fitter could estimate and the generator could not
# simulate would be a model that fits and then generates nothing, so the two
# lists are one list.

# Starting values read off the data rather than guessed. A population fit
# started far from the answer either converges slowly or reports the starting
# values back, and every one of these is a textbook non-compartmental reading:
# volume from the peak, clearance from the area under the curve, absorption from
# where the peak falls.
.model_initial_estimates <- function(observations, structural, pk_endpoint) {
  rows <- observations[observations$endpoint == pk_endpoint, , drop = FALSE]
  by_subject <- split(rows, rows$subject)
  dose <- stats::median(vapply(by_subject, function(part) {
    part$first_dose_amt[1L]
  }, numeric(1)), na.rm = TRUE)
  if (!is.finite(dose) || dose <= 0) dose <- 1

  peak <- stats::median(vapply(by_subject, function(part) max(part$dv),
                               numeric(1)), na.rm = TRUE)
  auc <- stats::median(vapply(by_subject, function(part) {
    part <- part[order(part$actual_tad), , drop = FALSE]
    if (nrow(part) < 2L) return(NA_real_)
    sum(diff(part$actual_tad) *
          (utils::head(part$dv, -1L) + utils::tail(part$dv, -1L)) / 2)
  }, numeric(1)), na.rm = TRUE)
  tmax <- stats::median(vapply(by_subject, function(part) {
    part$tad[which.max(part$dv)]
  }, numeric(1)), na.rm = TRUE)

  v <- if (is.finite(peak) && peak > 0) dose / peak else 1
  cl <- if (is.finite(auc) && auc > 0) dose / auc else v / 10
  ka <- if (is.finite(tmax) && tmax > 0) 4 / tmax else 1
  out <- c(cl = cl, v = v)
  if (grepl("oral", structural)) out <- c(out, ka = ka)
  if (grepl("^2cmt", structural)) out <- c(out, q = cl, v2 = v * 2)
  out
}

# The estimation method is `"focei"` rather than `"saem"`, which reverses what
# the design assumed, and the reason is the selection criterion. Choosing among
# candidates on AIC needs every candidate to have one. SAEM's log-likelihood is
# a Gaussian-quadrature step run after the fit, and on cohorts the size of a
# phase 1 study it returns a non-finite value: `theo_sd` fits perfectly well
# under SAEM -- clearance 2.75, volume 32.3, absorption 1.51, which are the
# textbook values -- and reports `AIC = Inf`, so a search over two candidates
# has nothing to compare. Under FOCEi the same fit reports AIC -42.1. SAEM
# remains available through `estimation` for a study large enough to give it a
# likelihood, and a candidate whose AIC is not finite is recorded as not
# converged whichever method produced it.

# One `nlmixr2` model function per candidate, written as text because that is
# what the shape of these functions is: a fixed block of parameter declarations
# and a fixed block of assignments, keyed by which parameters the structural
# model needs. `linCmt()` reads the parameter names and picks the same solution
# `.pk_single_dose()` evaluates, which is what keeps the two lists one list.
.model_nlmixr_function <- function(structural, start, error, error_start) {
  parameters <- names(start)
  ini <- c(
    sprintf("    t%s <- log(%.10g)", parameters, start),
    sprintf("    eta.%s ~ %.10g", parameters, ifelse(parameters == "cl", 0.1,
                                                     0.1)),
    sprintf("    %s.err <- %.10g", error, error_start)
  )
  assignments <- sprintf("    %s <- exp(t%s + eta.%s)", parameters, parameters,
                         parameters)
  # `linCmt()` names the central volume `v` and the peripheral one `vp`.
  assignments <- sub("^    v2 <- ", "    vp <- ", assignments)
  predicted <- switch(error,
                      prop = "    linCmt() ~ prop(prop.err)",
                      add = "    linCmt() ~ add(add.err)")
  text <- paste(c(
    "function() {", "  ini({", ini, "  })", "  model({", assignments,
    predicted, "  })", "}"
  ), collapse = "\n")
  eval(parse(text = text))
}

# The dataset `nlmixr2` reads: the concentration endpoint and the dosing
# records, on recorded times.
#
# Recorded, not planned. Fitting against the planned schedule where a patient's
# dose was reduced would push the drop in concentration that followed into
# clearance, and the model would report a population that eliminates the drug
# faster than the real one.
.model_estimation_data <- function(source, roles, pk_endpoint) {
  time <- suppressWarnings(as.numeric(source[[roles$time]]))
  dv <- suppressWarnings(as.numeric(source[[roles$dv]]))
  endpoint <- .endpoint(source, roles)
  dosed <- .dose_rows(source, roles)
  observed <- .observation_rows(source, roles, require_present = TRUE) &
    endpoint == pk_endpoint
  keep <- (dosed | observed) & is.finite(time)

  amount <- if (is.null(roles$amt)) rep(0, nrow(source)) else
    suppressWarnings(as.numeric(source[[roles$amt]]))
  amount[!is.finite(amount)] <- 0

  out <- data.frame(
    ID = as.character(source[[roles$id]])[keep],
    TIME = time[keep],
    DV = ifelse(observed[keep], dv[keep], NA_real_),
    AMT = ifelse(dosed[keep], amount[keep], 0),
    EVID = ifelse(dosed[keep], 1L, 0L),
    stringsAsFactors = FALSE
  )
  if (!is.null(roles$rate)) {
    rate <- suppressWarnings(as.numeric(source[[roles$rate]]))
    rate[!is.finite(rate)] <- 0
    out$RATE <- ifelse(out$EVID != 0L, rate[keep], 0)
  }
  out <- out[order(out$ID, out$TIME, out$EVID == 0L), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# The search. Every candidate is fitted, the ones that converge are compared on
# AIC, and the ones that do not stay in the table carrying their reason -- a
# search that came down to one survivor should not look like a search that had
# one candidate.
.model_fit_candidates <- function(data, candidates, observations, pk_endpoint,
                                  error, estimation, quiet) {
  fits <- list()
  rows <- list()
  for (candidate in candidates) {
    start <- .model_initial_estimates(observations, candidate, pk_endpoint)
    error_start <- if (identical(error, "prop")) 0.2 else
      stats::sd(data$DV, na.rm = TRUE) / 5
    spec <- .model_nlmixr_function(candidate, start, error, error_start)
    fit <- try(suppressWarnings(suppressMessages(
      nlmixr2est::nlmixr(spec, data, est = estimation,
                         control = list(print = 0L))
    )), silent = TRUE)
    converged <- !inherits(fit, "try-error") &&
      is.finite(suppressWarnings(stats::AIC(fit)))
    rows[[length(rows) + 1L]] <- data.frame(
      model = candidate, converged = converged,
      aic = if (converged) as.numeric(stats::AIC(fit)) else NA_real_,
      note = if (converged) "" else
        if (inherits(fit, "try-error")) .model_first_line(fit) else
          "no finite objective function",
      stringsAsFactors = FALSE
    )
    if (converged) fits[[candidate]] <- fit
    if (!quiet) {
      message(sprintf("  %-14s %s", candidate,
                      if (converged) sprintf("AIC %.1f", stats::AIC(fit))
                      else "did not converge"))
    }
  }
  table <- do.call(rbind, rows)
  table <- table[order(table$aic, na.last = TRUE), , drop = FALSE]
  rownames(table) <- NULL
  if (!any(table$converged)) {
    stop("No candidate model converged, so there is nothing to generate from. ",
         "Candidates tried: ", paste(table$model, collapse = ", "), ". ",
         "`synpmx_avatar()` and `synpmx_pca()` need no identifiable structure ",
         "and will run on this study.", call. = FALSE)
  }
  list(table = table, selected = table$model[which.min(table$aic)],
       fits = fits)
}

.model_first_line <- function(x) {
  message <- conditionMessage(attr(x, "condition"))
  trimws(strsplit(message, "\n", fixed = TRUE)[[1L]][1L])
}

# Reading a converged fit back into the numbers the generator needs, and no
# others. Empirical Bayes estimates are read to report a correlation and are
# then discarded, because they are a description of each real patient.
.model_read_fit <- function(fit, structural, error) {
  fixed <- stats::setNames(
    exp(as.numeric(fit$parFixedDf[paste0("t", .required_pk_params[[structural]]),
                                  "Estimate"])),
    .required_pk_params[[structural]]
  )
  omega <- fit$omega
  rownames(omega) <- sub("^eta\\.", "", rownames(omega))
  colnames(omega) <- rownames(omega)
  keep <- rownames(omega) %in% names(fixed)
  omega <- omega[keep, keep, drop = FALSE]
  residual_value <- as.numeric(
    fit$parFixedDf[paste0(error, ".err"), "Estimate"]
  )
  list(
    fixed = fixed, omega = omega,
    residual = if (identical(error, "prop")) {
      list(kind = "proportional", cv = residual_value)
    } else {
      list(kind = "additive", sd = residual_value)
    },
    # `$eta` carries the subject identifier in its first column. Only the
    # random effects themselves are read, and only to report a correlation.
    etas = {
      frame <- as.data.frame(fit$eta)
      frame[, grepl("^eta\\.", names(frame)), drop = FALSE]
    }
  )
}

# The PD shapes --------------------------------------------------------------
#
# Three time courses with no exposure dependence, fitted by least squares on the
# pooled observations and selected on AIC. Least squares rather than `nlmixr2`
# deliberately: these are three-parameter curves on one endpoint, the fit is
# well conditioned, and routing them through a population fitter would put a
# compiler in the path of every PD endpoint for no gain in what the generator
# then draws.
#
# A PD endpoint driven by concentration is reproduced as a time course that
# happens to resemble the average subject's response. A dataset whose point is
# the exposure-response relationship is not served by this generator, and this
# is where that is true.
.model_fit_pd <- function(observations, endpoint, shapes = NULL) {
  rows <- observations[observations$endpoint == endpoint &
                         is.finite(observations$dv), , drop = FALSE]
  if (nrow(rows) < 4L) return(NULL)
  time <- rows$tad
  value <- rows$dv

  candidates <- list()
  constant <- stats::lm(value ~ 1)
  candidates$constant <- list(pd = "constant",
                              typical = c(baseline = unname(stats::coef(constant)[1L])),
                              aic = stats::AIC(constant))
  linear <- stats::lm(value ~ time)
  coefficients <- stats::coef(linear)
  candidates$linear <- list(
    pd = "linear",
    typical = c(baseline = unname(coefficients[1L]),
                slope = unname(coefficients[2L])),
    aic = stats::AIC(linear)
  )
  exponential <- try(stats::nls(
    value ~ plateau + (baseline - plateau) * exp(-rate * pmax(time, 0)),
    start = list(plateau = stats::median(value[time > stats::median(time)]),
                 baseline = stats::median(value[time <= stats::median(time)]),
                 rate = 1 / max(stats::median(time), 1e-6))
  ), silent = TRUE)
  if (!inherits(exponential, "try-error")) {
    candidates$exponential <- list(
      pd = "exponential", typical = stats::coef(exponential),
      aic = stats::AIC(exponential)
    )
  }
  if (!is.null(shapes) && endpoint %in% names(shapes)) {
    chosen <- candidates[[shapes[[endpoint]]]]
    if (is.null(chosen)) {
      stop("`pd` names the shape `", shapes[[endpoint]], "` for endpoint `",
           endpoint, "`, which could not be fitted to it.", call. = FALSE)
    }
  } else {
    chosen <- candidates[[which.min(vapply(candidates, function(c) c$aic,
                                           numeric(1)))]]
  }

  # Between-subject variability on the baseline, read from the spread of each
  # subject's earliest observation. That is what the generator draws on, and it
  # is a variance rather than any subject's own value.
  first <- vapply(split(rows, rows$subject), function(part) {
    part$dv[which.min(part$tad)]
  }, numeric(1))
  first <- first[is.finite(first) & first > 0]
  chosen$baseline_cv <- if (length(first) > 1L) stats::sd(log(first)) else 0
  chosen$residual <- list(kind = "additive",
                          sd = stats::sd(value - .pd_profile(
                            list(pd = chosen$pd), time, numeric(), numeric(),
                            chosen$typical)))
  chosen$candidates <- data.frame(
    shape = names(candidates),
    aic = vapply(candidates, function(c) c$aic, numeric(1)),
    row.names = NULL, stringsAsFactors = FALSE
  )
  chosen
}

# Allometric scaling on clearance and volume, and nothing else, under
# `covariate_effects = "auto"`. The exponents are the standard 0.75 and 1 rather
# than estimated ones, and the effect is kept only where it improves AIC.
#
# The cost of this default is explicit. A covariate that influences the real
# profiles and is not in the model is generated independently of them, so the
# synthetic data carries no relationship between the two. `synpmx_avatar()`
# preserves those relationships without modelling them, because a blended
# subject's covariates and profile come from the same donors. `model_report()`
# reports the correlation between each declared covariate and the individual
# random effects, which is where an unmodelled relationship shows up.
.model_weight_covariate <- function(source, roles) {
  for (covariate in roles$covariates) {
    values <- suppressWarnings(as.numeric(source[[covariate]]))
    values <- values[is.finite(values)]
    if (!length(values) || any(values <= 0)) next
    if (grepl("^(wt|weight|bw|bodyweight|body_weight)$", covariate,
              ignore.case = TRUE)) {
      return(list(covariate = covariate, reference = stats::median(values)))
    }
  }
  NULL
}

.model_fit_allometry <- function(data, source, roles, structural, start, error,
                                 error_start, estimation, weight, base_aic) {
  column <- weight$covariate
  by_subject <- vapply(split(source[[column]],
                             as.character(source[[roles$id]])),
                       function(x) suppressWarnings(as.numeric(x[1L])),
                       numeric(1))
  data[[column]] <- unname(by_subject[data$ID])
  if (any(!is.finite(data[[column]]))) return(NULL)

  parameters <- names(start)
  exponent <- c(cl = 0.75, v = 1, q = 0.75, v2 = 1)
  ini <- c(sprintf("    t%s <- log(%.10g)", parameters, start),
           sprintf("    eta.%s ~ 0.1", parameters),
           sprintf("    %s.err <- %.10g", error, error_start))
  assignments <- vapply(parameters, function(parameter) {
    scaling <- if (parameter %in% names(exponent)) {
      sprintf(" * (%s / %.10g)^%.2f", column, weight$reference,
              exponent[[parameter]])
    } else ""
    sprintf("    %s <- exp(t%s + eta.%s)%s", parameter, parameter, parameter,
            scaling)
  }, character(1))
  assignments <- sub("^    v2 <- ", "    vp <- ", assignments)
  text <- paste(c("function() {", "  ini({", ini, "  })", "  model({",
                  assignments,
                  sprintf("    linCmt() ~ %s(%s.err)", error, error),
                  "  })", "}"), collapse = "\n")
  spec <- eval(parse(text = text))
  fit <- try(suppressWarnings(suppressMessages(
    nlmixr2est::nlmixr(spec, data, est = estimation, control = list(print = 0L))
  )), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  aic <- suppressWarnings(stats::AIC(fit))
  if (!is.finite(aic) || aic >= base_aic) return(NULL)
  list(fit = fit, aic = aic, effects = stats::setNames(
    lapply(intersect(parameters, names(exponent)), function(parameter) {
      list(covariate = column, reference = weight$reference,
           exponent = unname(exponent[[parameter]]))
    }), intersect(parameters, names(exponent))
  ))
}

#' Estimate a population model from a trial
#'
#' The only stage that reads patient data, and the only one that needs
#' `nlmixr2`. It works out which endpoint is the drug concentration and what
#' design produced it, fits the candidate models that design admits, picks one
#' on AIC, and returns that fit alongside the dosing and visit models the
#' generated subjects are built from. No patient row survives it.
#'
#' **The fitted parameters are not estimates to report.** They exist to make
#' simulated profiles look like the source study. The candidate set is five
#' linear models and the covariate model is allometric scaling or nothing, which
#' is too little to answer a scientific question, and the object prints that
#' warning with itself because its contents look exactly like the output of a
#' real population analysis.
#'
#' `nominal_time` is required, for two reasons. The dosing and visit models sit
#' on the nominal grid, and a grid inferred from recorded times is a statement
#' about the protocol only the caller can make. Estimation, separately, reads
#' the recorded times and the recorded dosing history, because a population fit
#' is a statement about the dose that was actually given.
#'
#' No formal privacy guarantee is offered. No patient's measured value reaches
#' the output, which is the claim `synpmx_pca()` makes and is stronger than
#' `synpmx_avatar()`'s, but the fixed effects and the covariance matrix are
#' functions of the individuals in the source and neither is noised. The cohort
#' floor is the whole defence and it is a threshold rather than an accounting.
#'
#' @param data Source PMX event data.
#' @param roles Explicit column roles from [pmx_roles()], including
#'   `nominal_time`.
#' @param pk One of the five built-in structural models, forcing it and skipping
#'   the search. `NULL` searches the candidates the design admits.
#' @param pd Named character vector of PD shapes per endpoint, skipping that
#'   search. One of `"constant"`, `"linear"` or `"exponential"` each.
#' @param endpoint_roles Named character vector naming which endpoint is the
#'   drug concentration, as `c(pk = "cp")`, overriding the inference.
#' @param covariate_effects `"auto"` fits allometric scaling on clearance and
#'   volume where a weight-like covariate is declared and keeps it where it
#'   improves AIC. `"none"` fits nothing.
#' @param min_subjects Cohort floor. Below it the covariance matrix describes
#'   the subjects it was fitted to rather than a population.
#' @param min_arm_patients Minimum patients in every arm, as
#'   [synpmx_pca_summarize()] uses.
#' @param min_time_bins Minimum distinct nominal times after a dose across the
#'   cohort. Below it no linear model is identifiable.
#' @param estimation Passed to `nlmixr2`. `"focei"` by default because the
#'   selection criterion is AIC and `"saem"` does not reliably produce one at
#'   these cohort sizes.
#' @param seed Seed for the one random step, which is imputing censored values
#'   before the fit.
#' @param quiet Suppress the per-candidate progress messages.
#'
#' @return A `pmx_fitted_model`.
#' @seealso [synpmx_model_generate()], [synpmx_model()], [model_report()],
#'   [model_candidates()], [model_parameters()].
#' @export
synpmx_model_estimate <- function(data, roles, pk = NULL, pd = NULL,
                                  endpoint_roles = NULL,
                                  covariate_effects = "auto",
                                  min_subjects = 20L, min_arm_patients = 3L,
                                  min_time_bins = 6L, estimation = "focei",
                                  seed = NULL, quiet = FALSE) {
  if (!inherits(roles, "pmx_roles")) {
    stop("`roles` must come from `pmx_roles()`.", call. = FALSE)
  }
  # Arguments are checked before the fitter is required, so that a
  # mistyped model name reads as a mistyped model name rather than as a
  # missing suggested package.
  if (!is.null(pk)) pk <- match.arg(pk, .pk_models)
  if (!identical(covariate_effects, "auto") &&
      !identical(covariate_effects, "none")) {
    stop("`covariate_effects` must be \"auto\" or \"none\".", call. = FALSE)
  }
  if (!requireNamespace("nlmixr2est", quietly = TRUE)) {
    stop("`synpmx_model_estimate()` needs the nlmixr2 package, which is in ",
         "Suggests. Install it, or use `synpmx_avatar()` or `synpmx_pca()`, ",
         "which fit no structural model.", call. = FALSE)
  }
  data <- as.data.frame(data)
  source <- data[, intersect(.retained_role_columns(roles), names(data)),
                 drop = FALSE]

  n_source <- length(.unique_in_order(source[[roles$id]]))
  .model_require_subjects(n_source, min_subjects)
  .model_require_nominal_time(source, roles)
  .model_require_time_coverage(source, roles, min_time_bins)

  censoring_source <- source
  source <- if (is.null(seed)) .impute_censored(source, roles) else
    .with_local_seed(seed, .impute_censored(source, roles))

  observations <- .model_observations(source, roles)
  classified <- .model_classify_endpoints(source, roles, observations,
                                          endpoint_roles)
  design <- .model_detect_design(source, roles, observations, classified$pk)
  if (!is.null(pk)) {
    design$candidates <- pk
    design$reason <- paste0("declared through `pk = \"", pk, "\"`")
  }

  subject_group <- .model_subject_arms(source, roles)
  .require_arms(subject_group, min_arm_patients, "synpmx_model_estimate()")

  # The apparatus, on the nominal grid, exactly as `synpmx_pca_summarize()`
  # builds it. `.model_cells()` is this generator's adapter over the same
  # `.arm_models()` the other one calls.
  fittable <- c(classified$pk, classified$pd)
  planned <- source
  planned[[roles$time]] <-
    suppressWarnings(as.numeric(source[[roles$nominal_time]]))
  cells <- .model_cells(source, roles, c(fittable, classified$discrete),
                        min_arm_patients)
  arm_models <- .arm_models(planned, roles, cells, subject_group,
                            min_arm_patients)

  estimation_data <- .model_estimation_data(source, roles, classified$pk)
  values <- estimation_data$DV[!is.na(estimation_data$DV)]
  error <- if (any(values <= 0)) "add" else "prop"
  if (!quiet) {
    message("Fitting ", length(design$candidates), " candidate(s) for `",
            classified$pk, "` (", design$route, "):")
  }
  search <- .model_fit_candidates(estimation_data, design$candidates,
                                  observations, classified$pk, error,
                                  estimation, quiet)
  selected <- search$selected
  parameters <- .model_read_fit(search$fits[[selected]], selected, error)

  effects <- list()
  if (identical(covariate_effects, "auto")) {
    weight <- .model_weight_covariate(source, roles)
    if (!is.null(weight)) {
      base_aic <- search$table$aic[search$table$model == selected]
      allometric <- .model_fit_allometry(
        estimation_data, source, roles, selected,
        .model_initial_estimates(observations, selected, classified$pk), error,
        if (identical(error, "prop")) 0.2 else stats::sd(values) / 5,
        estimation, weight, base_aic
      )
      if (!is.null(allometric)) {
        parameters <- .model_read_fit(allometric$fit, selected, error)
        effects <- allometric$effects
        search$table$note[search$table$model == selected] <-
          sprintf("allometric on %s improved AIC to %.1f", weight$covariate,
                  allometric$aic)
      }
    }
  }

  pd_fits <- stats::setNames(lapply(classified$pd, function(endpoint) {
    .model_fit_pd(observations, endpoint, pd)
  }), classified$pd)
  pd_fits <- pd_fits[!vapply(pd_fits, is.null, logical(1))]

  correlations <- .model_covariate_correlations(source, roles, subject_group,
                                                parameters$etas)
  parameters$etas <- NULL

  .pmx_fitted_model(
    structural = selected, candidates = search$table, parameters = parameters,
    endpoints = list(pk = classified$pk, pd = names(pd_fits),
                     discrete = classified$discrete, signals = classified$signals,
                     decided_by = classified$decided_by),
    arms = list(arms = arm_models$arms, sizes = arm_models$sizes),
    dosing = arm_models$dosing, visits = arm_models$visits,
    schema = .source_schema(censoring_source, roles, fittable, subject_group),
    roles = roles,
    settings = list(min_subjects = min_subjects,
                    min_arm_patients = min_arm_patients,
                    min_time_bins = min_time_bins, estimation = estimation,
                    covariate_effects = covariate_effects, error = error),
    n_source = n_source,
    cells = cells, pd = pd_fits, covariate_effects = effects,
    covariates = .covariate_model(source, roles, subject_group),
    discrete = .discrete_model(source, roles, cells, subject_group),
    design = design, correlations = correlations
  )
}

.model_subject_arms <- function(source, roles) {
  subjects <- .unique_in_order(source[[roles$id]])
  strata_key <- as.character(.subject_strata(source, roles))
  vapply(subjects, function(subject) {
    rows <- which(!is.na(source[[roles$id]]) & source[[roles$id]] == subject)
    strata_key[rows[1L]]
  }, character(1))
}

# Where an unmodelled covariate relationship shows up. The random effects are
# per-subject quantities and are not stored; only these correlations are, which
# is a number per covariate and parameter rather than a description of anybody.
.model_covariate_correlations <- function(source, roles, subject_group, etas) {
  if (is.null(etas) || !length(roles$covariates)) return(NULL)
  subjects <- .unique_in_order(source[[roles$id]])
  first_row <- vapply(subjects, function(subject) {
    which(!is.na(source[[roles$id]]) & source[[roles$id]] == subject)[1L]
  }, integer(1))
  rows <- list()
  for (covariate in roles$covariates) {
    values <- suppressWarnings(as.numeric(source[[covariate]][first_row]))
    if (sum(is.finite(values)) < 3L) next
    for (parameter in names(etas)) {
      eta <- etas[[parameter]]
      if (length(eta) != length(values)) next
      correlation <- suppressWarnings(
        stats::cor(values, eta, use = "complete.obs")
      )
      if (!is.finite(correlation)) next
      rows[[length(rows) + 1L]] <- data.frame(
        covariate = covariate, parameter = sub("^eta\\.", "", parameter),
        correlation = correlation, stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Estimate a population model and generate a synthetic dataset from it
#'
#' A single call for [synpmx_model_estimate()] followed by
#' [synpmx_model_generate()]. Use the two separately to look at the fit before
#' generating from it; it is on the result either way, as the `pmx_fitted_model`
#' attribute.
#'
#' @param data Source PMX event data.
#' @param roles Explicit column roles from [pmx_roles()], including
#'   `nominal_time`.
#' @param n_subjects Number of synthetic subjects. Defaults to the source count.
#' @param seed Seed, used for both stages.
#' @param ... Passed to [synpmx_model_estimate()].
#'
#' @return A data frame in the source's shape, carrying the fitted model as an
#'   attribute.
#' @seealso [synpmx_model_estimate()], [synpmx_model_generate()],
#'   [synpmx_pca()], [synpmx_avatar()].
#' @export
synpmx_model <- function(data, roles, n_subjects = NULL, seed = NULL, ...) {
  synpmx_model_generate(
    synpmx_model_estimate(data, roles, seed = seed, ...),
    n_subjects = n_subjects, seed = seed
  )
}
