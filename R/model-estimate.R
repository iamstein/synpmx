# Estimating the candidates --------------------------------------------------
#
# The only stage that reads patient data and the only one that needs `nlmixr2`.
# The candidate set is exactly the models `.pk_single_dose()` can evaluate in
# closed form: a candidate the fitter could estimate and the generator could not
# simulate would be a model that fits and then generates nothing, so the two
# lists are one list.

# Non-compartmental analysis, per subject, within one dose interval.
#
# One interval, and that restriction is the whole of what makes this correct on
# a multiple-dose study. Sorting a subject's samples by time after dose pools
# every cycle they were in: on a twelve-week daily regimen the resulting
# sequence interleaves a sample drawn 0.1 h into cycle 85 with one drawn 0.2 h
# into cycle 1, which is not a concentration-time profile and has no area and no
# terminal slope. Read that way, `case1_pkpd` reports a half-life of 371 hours
# and the fit is started from a number that describes nothing.
#
# The terminal slope is the quantity everything else leans on. It is fitted by
# log-linear regression over the points after the peak -- at least three of
# them, all positive -- which is the standard reading.
#
# `extrapolate` says whether anything follows this interval. Where nothing does,
# the tail is real and the area runs to infinity, so dose over it is clearance.
# Where another dose follows, the area is taken over the interval as it stands:
# at steady state dose over the interval area is clearance too, and before
# steady state it is an overestimate that is still the right order of magnitude.
.model_nca_subject <- function(part, extrapolate = TRUE) {
  part <- part[order(part$actual_tad), , drop = FALSE]
  part <- part[is.finite(part$actual_tad) & is.finite(part$dv), , drop = FALSE]
  empty <- list(cmax = NA_real_, tmax = NA_real_, auc = NA_real_,
                aumc = NA_real_, lambda_z = NA_real_)
  if (nrow(part) < 3L) return(empty)
  time <- part$actual_tad
  value <- part$dv
  peak <- which.max(value)

  width <- diff(time)
  auc <- sum(width * (utils::head(value, -1L) + utils::tail(value, -1L)) / 2)
  aumc <- sum(width * (utils::head(time * value, -1L) +
                         utils::tail(time * value, -1L)) / 2)

  # The terminal phase: everything from the peak on, positive, and at least
  # three points. A slope that comes out non-negative is not a terminal phase
  # and is discarded rather than used.
  terminal <- seq(peak, length(value))
  terminal <- terminal[value[terminal] > 0]
  lambda_z <- NA_real_
  if (length(terminal) >= 3L) {
    slope <- stats::coef(stats::lm(log(value[terminal]) ~ time[terminal]))[2L]
    if (is.finite(slope) && slope < 0) lambda_z <- -as.numeric(slope)
  }
  if (extrapolate && is.finite(lambda_z)) {
    last_value <- value[length(value)]
    last_time <- time[length(time)]
    auc <- auc + last_value / lambda_z
    aumc <- aumc + last_time * last_value / lambda_z +
      last_value / lambda_z^2
  }
  list(cmax = value[peak], tmax = time[peak], auc = auc, aumc = aumc,
       lambda_z = lambda_z)
}

# Absorption from the peak position rather than from a rule of thumb. For a
# one-compartment oral model the peak sits at log(ka/ke)/(ka - ke), which is one
# equation in one unknown once the terminal slope has given `ke`. Solving it is
# worth the ten lines: `4 / tmax` was the previous answer and it is off by
# whatever the elimination rate happens to be.
.model_ka_from_tmax <- function(tmax, ke) {
  if (!is.finite(tmax) || tmax <= 0 || !is.finite(ke) || ke <= 0) return(1)
  gap <- function(ka) log(ka / ke) / (ka - ke) - tmax
  lower <- ke * 1.01
  upper <- ke * 1000
  if (!is.finite(gap(lower)) || !is.finite(gap(upper)) ||
      gap(lower) * gap(upper) > 0) {
    return(max(4 / tmax, ke * 1.5))
  }
  as.numeric(stats::uniroot(gap, c(lower, upper))$root)
}

# Starting values read off the data rather than guessed. A population fit
# started far from the answer either converges slowly, reports the starting
# values back, or -- as the two-compartment candidates did before this was
# written -- grinds against a flat likelihood for a minute where the same fit
# from good starts takes seconds.
#
# Every value here is a textbook non-compartmental reading. Clearance is dose
# over the area extrapolated to infinity, not over the trapezoid alone, which on
# a study sampled to four half-lives understates the area by about a tenth and
# overstates clearance by the same. Volume is clearance over the terminal slope,
# which is the quantity the terminal phase actually identifies -- dose over the
# peak, the previous answer, is not a volume of any kind for an oral dose.
.model_initial_estimates <- function(observations, structural, pk_endpoint) {
  rows <- observations[observations$endpoint == pk_endpoint, , drop = FALSE]
  by_subject <- split(rows, rows$subject)
  dose <- stats::median(vapply(by_subject, function(part) {
    part$first_dose_amt[1L]
  }, numeric(1)), na.rm = TRUE)
  if (!is.finite(dose) || dose <= 0) dose <- 1

  # Each subject read inside the cohort's richest dose interval, which is the
  # same interval every other signal is read in.
  interval <- .model_richest_interval(rows)
  nca <- lapply(by_subject, function(part) {
    within <- part[!is.na(part$interval) & part$interval == interval, ,
                   drop = FALSE]
    if (nrow(within) < 3L) return(.model_nca_subject(part[0L, ], TRUE))
    # Nothing follows the subject's last interval, so its tail is real.
    last <- max(part$interval, na.rm = TRUE)
    .model_nca_subject(within, extrapolate = identical(interval, last))
  })
  middle <- function(name) {
    stats::median(vapply(nca, function(x) x[[name]], numeric(1)), na.rm = TRUE)
  }
  auc <- middle("auc")
  aumc <- middle("aumc")
  cmax <- middle("cmax")
  tmax <- middle("tmax")
  lambda_z <- middle("lambda_z")

  cl <- if (is.finite(auc) && auc > 0) dose / auc else 1
  volume_terminal <- if (is.finite(lambda_z) && lambda_z > 0) cl / lambda_z else
    if (is.finite(cmax) && cmax > 0) dose / cmax else 10 * cl
  ka <- .model_ka_from_tmax(tmax, if (is.finite(lambda_z)) lambda_z else
    cl / volume_terminal)

  out <- c(cl = cl, v = volume_terminal)
  if (grepl("oral", structural)) out <- c(out, ka = ka)
  if (grepl("^2cmt", structural)) {
    # Steady-state volume from the mean residence time, which is what the first
    # moment of the curve is for. For an oral dose the residence time includes
    # the time spent absorbing, so the absorption mean is taken back off.
    mrt <- if (is.finite(aumc) && is.finite(auc) && auc > 0) aumc / auc else NA
    if (grepl("oral", structural) && is.finite(mrt) && is.finite(ka) && ka > 0) {
      mrt <- mrt - 1 / ka
    }
    vss <- if (is.finite(mrt) && mrt > 0) cl * mrt else volume_terminal / 2
    # The central compartment is the smaller half and the peripheral the rest,
    # with intercompartmental clearance started at elimination clearance. These
    # are starting values for a search, not a decomposition of the curve: what
    # matters is that they are the right order of magnitude and that the central
    # volume is below the terminal one, which the previous `v * 2` was not.
    central <- max(min(vss / 2, volume_terminal * 0.75), 1e-6)
    out[["v"]] <- central
    out <- c(out, q = cl, v2 = max(vss - central, central))
  }
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
#
# `weight` folds allometric scaling into the same model rather than fitting a
# second one to compare against. The exponents are the standard 0.75 and 1 and
# are not estimated: this is a shape that makes a synthetic cohort's spread look
# right, not a covariate analysis, and testing it against AIC would double the
# cost of the only fit this function performs.
.model_allometric_exponents <- c(cl = 0.75, v = 1, q = 0.75, v2 = 1)

.model_nlmixr_function <- function(structural, start, error, error_start,
                                   weight = NULL) {
  parameters <- names(start)
  ini <- c(
    sprintf("    t%s <- log(%.10g)", parameters, start),
    sprintf("    eta.%s ~ 0.1", parameters),
    sprintf("    %s.err <- %.10g", error, error_start)
  )
  assignments <- vapply(parameters, function(parameter) {
    scaling <- if (!is.null(weight) &&
                   parameter %in% names(.model_allometric_exponents)) {
      sprintf(" * (%s / %.10g)^%.2f", weight$covariate, weight$reference,
              .model_allometric_exponents[[parameter]])
    } else ""
    sprintf("    %s <- exp(t%s + eta.%s)%s", parameter, parameter, parameter,
            scaling)
  }, character(1))
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
                                  error, estimation, quiet, weight = NULL) {
  fits <- list()
  rows <- list()
  for (candidate in candidates) {
    start <- .model_initial_estimates(observations, candidate, pk_endpoint)
    error_start <- if (identical(error, "prop")) 0.2 else
      stats::sd(data$DV, na.rm = TRUE) / 5
    spec <- .model_nlmixr_function(candidate, start, error, error_start,
                                   weight)
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
                         is.finite(observations$dv) &
                         is.finite(observations$aligned), , drop = FALSE]
  if (nrow(rows) < 4L) return(NULL)
  # Study time from the first dose, which is the axis generation evaluates the
  # shape on. Fitting against time after dose instead makes every sample in a
  # daily regimen land at nearly the same place.
  time <- rows$aligned
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
  # The earliest observation each subject has, whether or not a dose preceded
  # it. A placebo arm has no dose records at all, so every one of its rows is
  # outside any dose interval; asking for the earliest by time after dose
  # returns nothing and takes the whole fit down with it.
  first <- vapply(split(rows, rows$subject), function(part) {
    if (!any(is.finite(part$aligned))) return(NA_real_)
    part$dv[which.min(replace(part$aligned, !is.finite(part$aligned), Inf))]
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

# The covariate the allometric scaling is applied to, or nothing. Weight-like by
# name, positive, and numeric -- there is no way to recognise a body weight from
# its values alone, and guessing from a distribution would be worse than asking.
#
# Applied where it is found rather than tested for improvement. Testing costs a
# second fit, which is the whole budget of the default path, and allometry on a
# weight is a shape this generator asserts rather than a hypothesis it examines.
# `covariate_effects = "none"` switches it off.
#
# The cost of the default is explicit. A covariate that influences the real
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
    # A baseline weight is rarely called `WT`. `WEIGHTB`, `WTBL` and `BWT` are
    # all the same column with a suffix or a prefix saying it was measured at
    # baseline, and refusing them means silently fitting no covariate at all on
    # a study that declared one. Anything not on this list is not scaled, which
    # is the safe direction: `HEIGHT` and `AGE` must not match.
    if (grepl("^(b|base|baseline)?_?(wt|wgt|weight|bw|bodywt|body_?weight)_?(b|bl|base|baseline|0)?$",
              covariate, ignore.case = TRUE)) {
      return(list(covariate = covariate, reference = stats::median(values)))
    }
  }
  NULL
}

# The per-subject weight, on the estimation dataset.
.model_attach_weight <- function(data, source, roles, weight) {
  if (is.null(weight)) return(data)
  by_subject <- vapply(split(source[[weight$covariate]],
                             as.character(source[[roles$id]])),
                       function(x) suppressWarnings(as.numeric(x[1L])),
                       numeric(1))
  data[[weight$covariate]] <- unname(by_subject[data$ID])
  if (any(!is.finite(data[[weight$covariate]]))) return(NULL)
  data
}

# What the assay limit cost, per endpoint.
#
# Values below the limit are imputed before anything is fitted -- a uniform draw
# inside the censoring region rather than a fixed LLOQ/2, which would replace one
# artificial spike with another -- and the boundary is put back at emit. That is
# the intended behaviour and not a shortcut around the M3 likelihood: the same
# imputation is what lets the apparatus, the PD shapes and the covariate model
# read a latent value rather than a stack of identical boundary substitutions.
#
# It is an assumption all the same, and its weight is the share of the endpoint
# that carries it, so that share is measured here and reported with the fit. A
# study whose concentrations are 46% below the limit is a study whose fitted
# parameters are substantially a statement about the draw.
.model_censoring_summary <- function(source, roles, endpoints) {
  if (is.null(roles$cens)) return(NULL)
  observed <- .observation_rows(source, roles, require_present = TRUE)
  endpoint <- .endpoint(source, roles)
  cens <- suppressWarnings(as.numeric(as.character(source[[roles$cens]])))
  dv <- suppressWarnings(as.numeric(source[[roles$dv]]))
  rows <- lapply(endpoints, function(name) {
    at <- observed & endpoint == name
    if (!any(at)) return(NULL)
    censored <- at & is.finite(cens) & cens != 0
    data.frame(
      endpoint = name, observations = sum(at), imputed = sum(censored),
      fraction = sum(censored) / sum(at),
      limit = if (any(censored)) stats::median(dv[censored]) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
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
  if (!is.null(pk)) {
    unknown <- setdiff(pk, .pk_models)
    if (length(unknown)) {
      stop("`pk` names model(s) outside the closed-form set: ",
           paste(unknown, collapse = ", "), ". Available: ",
           paste(.pk_models, collapse = ", "), ".", call. = FALSE)
    }
  }
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
    design$reason <- if (length(pk) == 1L) {
      paste0("declared through `pk = \"", pk, "\"`")
    } else {
      paste0("searched over the ", length(pk), " models named in `pk`")
    }
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

  # Allometric scaling is folded into the fit rather than compared against one
  # without it, so the default path performs exactly one fit.
  weight <- if (identical(covariate_effects, "auto")) {
    .model_weight_covariate(source, roles)
  } else NULL
  with_weight <- .model_attach_weight(estimation_data, source, roles, weight)
  if (is.null(with_weight)) weight <- NULL else estimation_data <- with_weight

  if (!quiet) {
    message("Fitting ", length(design$candidates), " model(s) for `",
            classified$pk, "` (", design$route, ")",
            if (!is.null(weight)) paste0(", allometric on ", weight$covariate),
            ":")
  }
  search <- .model_fit_candidates(estimation_data, design$candidates,
                                  observations, classified$pk, error,
                                  estimation, quiet, weight)
  selected <- search$selected
  parameters <- .model_read_fit(search$fits[[selected]], selected, error)

  effects <- if (is.null(weight)) list() else stats::setNames(
    lapply(intersect(names(parameters$fixed),
                     names(.model_allometric_exponents)), function(parameter) {
      list(covariate = weight$covariate, reference = weight$reference,
           exponent = unname(.model_allometric_exponents[[parameter]]))
    }), intersect(names(parameters$fixed), names(.model_allometric_exponents))
  )

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
    design = design, correlations = correlations,
    censoring = .model_censoring_summary(censoring_source, roles, fittable)
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
