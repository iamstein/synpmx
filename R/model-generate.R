# Generating from a fitted model ---------------------------------------------
#
# This stage reads no patient data. Its inputs are the fitted model and a
# subject count, and it runs on base R: `.pk_profile()` already superposes
# linear doses, takes a per-subject parameter vector and handles an infusion
# duration, so there is no solver and no compiler here.
#
# The order of the two draws is fixed. The schedule is drawn first and the
# profile is computed from it, so a subject who steps down a dose level has a
# lower exposure from that cycle on and one who skips a cycle has the trough
# that implies. Drawing a profile and then a schedule would reproduce the
# principal-component generator's disconnect between the two with extra steps,
# and closing that disconnect is the main thing this generator offers over it.

# The grid attendance is measured on, built from the nominal times the source
# actually holds. This is the adapter `.pca_cells()` is for the other generator:
# the shared arm models take the cells rather than a fitted object, so the two
# representations meet here and nowhere else.
#
# A cell is kept where at least `floor` distinct patients have an observation
# there. A nominal time one patient attended is that patient, and generating
# from it would put them back.
.model_cells <- function(source, roles, endpoints, floor) {
  nominal <- suppressWarnings(as.numeric(source[[roles$nominal_time]]))
  planned <- source
  planned[[roles$time]] <- nominal
  aligned <- .aligned_time(planned, roles)
  observed <- .observation_rows(source, roles, require_present = TRUE) &
    is.finite(aligned)
  endpoint <- .endpoint(source, roles)
  id <- as.character(source[[roles$id]])

  rows <- list()
  for (name in endpoints) {
    selected <- observed & endpoint == name
    if (!any(selected)) next
    for (time in sort(unique(aligned[selected]))) {
      at <- selected & abs(aligned - time) < sqrt(.Machine$double.eps)
      if (length(unique(id[at])) < floor) next
      rows[[length(rows) + 1L]] <- data.frame(
        name = sprintf("%s@%.10g", name, time), endpoint = name, time = time,
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) {
    stop("No nominal grid cell is held by at least ", floor, " patients, so ",
         "there is no visit model to build. Every generated observation would ",
         "sit at a time one real patient attended.", call. = FALSE)
  }
  out <- do.call(rbind, rows)
  out$index <- seq_len(nrow(out))
  out[, c("index", "name", "endpoint", "time")]
}

# Baseline covariates, per arm. A mean and a standard deviation on the log scale
# for a positive continuous covariate, on the natural scale otherwise, and level
# frequencies for a categorical one. Deliberately independent draws: this
# generator models the relationship between a covariate and a profile only where
# a covariate effect was fitted, and `model_report()` reports the correlations
# that says nothing about.
.covariate_model <- function(source, roles, subject_group) {
  subjects <- .unique_in_order(source[[roles$id]])
  first_row <- vapply(subjects, function(subject) {
    which(!is.na(source[[roles$id]]) & source[[roles$id]] == subject)[1L]
  }, integer(1))
  out <- list()
  for (arm in unique(subject_group)) {
    members <- first_row[subject_group == arm]
    out[[arm]] <- lapply(stats::setNames(roles$covariates, roles$covariates),
                         function(column) {
      values <- source[[column]][members]
      values <- values[!is.na(values)]
      if (!length(values)) return(list(kind = "missing"))
      if (is.numeric(values) && !is.factor(values)) {
        if (all(values > 0)) {
          logged <- log(values)
          return(list(kind = "lognormal", meanlog = mean(logged),
                      sdlog = stats::sd(logged) %|na|% 0,
                      median = stats::median(values)))
        }
        return(list(kind = "normal", mean = mean(values),
                    sd = stats::sd(values) %|na|% 0,
                    median = stats::median(values)))
      }
      counts <- table(as.character(values))
      list(kind = "categorical", levels = names(counts),
           probability = as.numeric(counts) / sum(counts))
    })
  }
  out
}

`%|na|%` <- function(x, y) if (is.na(x)) y else x

.draw_covariates <- function(model, n) {
  lapply(model, function(spec) {
    switch(
      spec$kind,
      lognormal = exp(stats::rnorm(n, spec$meanlog, spec$sdlog)),
      normal = stats::rnorm(n, spec$mean, spec$sd),
      # Indexed rather than sampled by value: `sample(x, n)` on a length-one
      # numeric `x` samples `seq_len(x)` instead of `x`, so a covariate whose
      # source held one level would draw from that level's own value.
      categorical = spec$levels[sample.int(length(spec$levels), n,
                                           replace = TRUE,
                                           prob = spec$probability)],
      rep(NA, n)
    )
  })
}

# Per-visit marginals for the endpoints that are not time courses. A binary or
# ordinal endpoint is drawn from the level frequencies its arm holds at that
# nominal time, which is the same thing the visit model does for attendance.
.discrete_model <- function(source, roles, cells, subject_group) {
  nominal <- suppressWarnings(as.numeric(source[[roles$nominal_time]]))
  planned <- source
  planned[[roles$time]] <- nominal
  aligned <- .aligned_time(planned, roles)
  observed <- .observation_rows(source, roles, require_present = TRUE)
  endpoint <- .endpoint(source, roles)
  subjects <- .unique_in_order(source[[roles$id]])
  arm_of <- stats::setNames(subject_group, as.character(subjects))
  row_arm <- arm_of[as.character(source[[roles$id]])]

  out <- list()
  for (arm in unique(subject_group)) {
    out[[arm]] <- lapply(seq_len(nrow(cells)), function(i) {
      at <- observed & row_arm == arm & endpoint == cells$endpoint[i] &
        abs(aligned - cells$time[i]) < sqrt(.Machine$double.eps)
      values <- source[[roles$dv]][at]
      values <- values[!is.na(values)]
      if (!length(values)) return(NULL)
      counts <- table(as.character(values))
      list(levels = as.numeric(names(counts)),
           probability = as.numeric(counts) / sum(counts))
    })
  }
  out
}

# Between-subject random effects. Drawn from the covariance matrix rather than
# read off any subject: an empirical Bayes estimate is a per-subject quantity,
# and a model that carried one would be writing out a description of a real
# patient.
.draw_random_effects <- function(omega, n) {
  names <- rownames(omega)
  decomposed <- eigen(omega, symmetric = TRUE)
  root <- decomposed$vectors %*%
    diag(sqrt(pmax(decomposed$values, 0)), nrow(omega))
  draws <- matrix(stats::rnorm(n * nrow(omega)), nrow = n) %*% t(root)
  colnames(draws) <- names
  draws
}

# The typical parameters, with this subject's covariates and random effects on
# them. Allometric scaling is the only covariate model fitted by default, and
# the exponents are the standard 0.75 and 1 rather than estimated ones: this is
# a shape that makes a synthetic cohort's spread look right, not a covariate
# analysis.
.subject_parameters <- function(fixed, effects, etas, covariates) {
  p <- fixed
  for (parameter in names(effects)) {
    effect <- effects[[parameter]]
    value <- covariates[[effect$covariate]]
    if (is.null(value) || !is.finite(value) || value <= 0) next
    p[[parameter]] <- p[[parameter]] *
      (value / effect$reference)^effect$exponent
  }
  for (parameter in names(etas)) {
    if (parameter %in% names(p)) p[[parameter]] <- p[[parameter]] *
      exp(etas[[parameter]])
  }
  p
}

# A concentration cannot be negative, and a proportional error on a value near
# the assay limit will produce one often enough to matter: on `warfarin` it puts
# the minimum at -1.4 against a source minimum of 0. Floored rather than
# redrawn, because redrawing until the value is positive is a truncated
# distribution nobody declared.
# The proportional multiplier is lognormal rather than `1 + N(0, cv)`.
#
# The two agree while the coefficient of variation is small, and part company
# exactly where a misfitted structural model puts it. `1 + N(0, cv)` is
# non-positive with probability `pnorm(-1 / cv)`: 2.3% of draws at a CV of 0.5
# and 7.7% at 0.7, each of them a concentration clamped to zero and, on a log
# axis, a profile falling off the bottom of the figure. A CV that high is a
# statement that the structural model does not describe the data, and it should
# come out as a wide band rather than as a scatter of zeros.
#
# `exp(N(0, sqrt(log(1 + cv^2))))` has the same coefficient of variation and a
# median of one, so the spread the fit estimated is preserved and no draw ever
# reaches zero.
.add_residual_error <- function(value, residual, floor = NULL) {
  out <- switch(
    residual$kind,
    proportional = value * exp(stats::rnorm(
      length(value), 0, sqrt(log(1 + residual$cv^2))
    )),
    additive = value + stats::rnorm(length(value), 0, residual$sd),
    value
  )
  if (is.null(floor)) out else pmax(out, floor)
}

#' Generate a synthetic PMX dataset from a fitted model
#'
#' Draws new subjects from a [synpmx_model_estimate()] fit. This stage reads no
#' patient data: its arguments are the model and a subject count, so everything
#' about the source that reaches the output has already passed through the fit.
#'
#' Per subject: an arm is assigned keeping the source arm shares, covariates are
#' drawn from the arm's covariate model, random effects from the between-subject
#' covariance matrix, and the dose schedule from the arm's dosing model. The
#' concentration is then evaluated at the visits drawn from the arm's visit
#' model, against the schedule that was drawn, so a reduced or skipped dose
#' reaches the concentrations rather than appearing only in the dosing records.
#'
#' @param fitted_model A `pmx_fitted_model` from [synpmx_model_estimate()].
#' @param n_subjects Number of synthetic subjects. Defaults to the source count.
#' @param seed Generation seed.
#'
#' @return A data frame in the source's shape, carrying the fitted model as a
#'   `pmx_fitted_model` attribute.
#' @seealso [synpmx_model_estimate()], [synpmx_model()], [model_report()].
#' @export
synpmx_model_generate <- function(fitted_model, n_subjects = NULL,
                                  seed = NULL) {
  if (!inherits(fitted_model, "pmx_fitted_model")) {
    stop("`fitted_model` must come from `synpmx_model_estimate()`.",
         call. = FALSE)
  }
  n_subjects <- as.integer(n_subjects %||% fitted_model$n_source)
  if (!is.finite(n_subjects) || n_subjects < 1L) {
    stop("`n_subjects` must be one positive integer.", call. = FALSE)
  }
  out <- if (is.null(seed)) {
    .model_generate(fitted_model, n_subjects)
  } else {
    .with_local_seed(seed, .model_generate(fitted_model, n_subjects))
  }
  attr(out, "pmx_fitted_model") <- fitted_model
  attr(out, "pmx_source") <- "model"
  out
}

.model_generate <- function(fit, n_subjects) {
  roles <- fit$roles
  schema <- fit$schema
  cells <- fit$cells

  assignment <- .assign_arms(fit$arms$arms, fit$arms$sizes, n_subjects)
  # An environment rather than a counter, because the emit loop below writes to
  # it from inside two nested `for`s.
  floored <- new.env(parent = emptyenv())
  floored$seen <- 0L
  floored$raised <- 0L

  etas <- .draw_random_effects(fit$parameters$omega, n_subjects)
  # Between-subject variability on a PD baseline is its own draw. It is not in
  # the PK covariance matrix, because the PD shapes are fitted separately and a
  # baseline is not a parameter of the concentration-time curve.
  pd_etas <- vapply(fit$pd, function(shape) {
    stats::rnorm(n_subjects, 0, shape$baseline_cv %||% 0)
  }, numeric(n_subjects))
  if (!is.matrix(pd_etas)) {
    pd_etas <- matrix(pd_etas, nrow = n_subjects,
                      dimnames = list(NULL, names(fit$pd)))
  }
  covariates <- stats::setNames(lapply(fit$arms$arms, function(arm) {
    .draw_covariates(fit$covariates[[arm]], sum(assignment == arm))
  }), fit$arms$arms)
  taken <- stats::setNames(integer(length(fit$arms$arms)), fit$arms$arms)

  pieces <- vector("list", n_subjects)
  subject_covariates <- vector("list", n_subjects)
  doses <- numeric(n_subjects)
  for (i in seq_len(n_subjects)) {
    arm <- assignment[[i]]
    taken[arm] <- taken[arm] + 1L
    mine <- lapply(covariates[[arm]], function(column) column[taken[arm]])
    subject_covariates[[i]] <- mine

    schedule <- .draw_schedule(fit$dosing[[arm]])
    doses[i] <- if (nrow(schedule)) sum(schedule$amt) else 0
    p <- .subject_parameters(fit$parameters$fixed, fit$covariate_effects,
                             etas[i, ], mine)

    rows <- list()
    if (nrow(schedule)) {
      rows[[length(rows) + 1L]] <- data.frame(
        TIME = schedule$time, DV = NA_real_, AMT = schedule$amt, EVID = 1L,
        CMT = if (is.null(schema$cmt_dose)) NA else schema$cmt_dose,
        DVID = NA_character_, stringsAsFactors = FALSE
      )
    }
    visits <- fit$visits[[arm]]
    attended <- stats::runif(nrow(cells)) < visits$probability
    for (index in which(attended)) {
      endpoint_name <- cells$endpoint[index]
      time <- cells$time[index]
      value <- if (identical(endpoint_name, fit$endpoints$pk)) {
        if (!nrow(schedule)) next
        concentration <- .pk_profile(list(pk = fit$structural), time,
                                     schedule$amt, schedule$time, p,
                                     fit$parameters$duration %||% 0)
        .add_residual_error(concentration, fit$parameters$residual, floor = 0)
      } else if (endpoint_name %in% names(fit$pd)) {
        shape <- fit$pd[[endpoint_name]]
        baseline <- shape$typical[["baseline"]] *
          exp(pd_etas[i, endpoint_name])
        .add_residual_error(
          .pd_profile(list(pd = shape$pd), time, schedule$amt, schedule$time,
                      params = replace(shape$typical, "baseline", baseline)),
          shape$residual
        )
      } else {
        marginal <- fit$discrete[[arm]][[index]]
        if (is.null(marginal)) next
        # Indexed for the same reason as the covariate draw above: a visit
        # every patient recorded the same level at holds one level, and
        # `sample()` would read that level as a count of levels.
        marginal$levels[[sample.int(length(marginal$levels), 1L,
                                    prob = marginal$probability)]]
      }
      if (!is.finite(value)) next
      # The floor the study's own smallest reported value implies, where it
      # declared no censoring column of its own. A value below it is one the
      # assay could not have returned, so it is reported at the floor rather
      # than at whatever the residual draw produced.
      floor_value <- fit$quantification_floor[[endpoint_name]]
      if (!is.null(floor_value)) {
        floored$seen <- floored$seen + 1L
        if (value < floor_value) {
          floored$raised <- floored$raised + 1L
          value <- floor_value
        }
      }
      value <- .snap_endpoint_values(value,
                                     schema$endpoint_specs[[endpoint_name]])
      rows[[length(rows) + 1L]] <- data.frame(
        TIME = time, DV = value, AMT = 0, EVID = 0L,
        CMT = if (is.null(schema$cmt_obs[[endpoint_name]])) NA else
          schema$cmt_obs[[endpoint_name]],
        DVID = endpoint_name, stringsAsFactors = FALSE
      )
    }
    if (!length(rows)) next
    piece <- do.call(rbind, rows)
    piece <- piece[order(piece$TIME, piece$EVID == 0L), , drop = FALSE]
    piece$.subject <- i
    piece$.arm <- arm
    pieces[[i]] <- piece
  }

  .model_emit(pieces, fit, n_subjects, doses, subject_covariates, floored)
}

# The finished table, in the source's shape. Everything here is bookkeeping the
# schema already decided: which columns exist, what class the identifiers are,
# which compartment each endpoint sits in, and where the assay limit is.
.model_emit <- function(pieces, fit, n_subjects, doses, subject_covariates,
                       floored = NULL) {
  roles <- fit$roles
  schema <- fit$schema
  frame <- do.call(rbind, pieces[!vapply(pieces, is.null, logical(1))])
  if (is.null(frame) || !nrow(frame)) {
    stop("Every generated subject came out empty. The visit model retained no ",
         "cell any subject attended.", call. = FALSE)
  }
  rownames(frame) <- NULL

  out <- data.frame(row.names = seq_len(nrow(frame)))
  ids <- .new_subject_ids(schema, n_subjects)
  out[[roles$id]] <- ids[frame$.subject]
  out[[roles$time]] <- frame$TIME
  out[[roles$nominal_time]] <- frame$TIME
  out[[roles$dv]] <- frame$DV
  if (!is.null(roles$amt)) out[[roles$amt]] <- frame$AMT
  out[[roles$evid]] <- frame$EVID
  if (!is.null(roles$cmt)) out[[roles$cmt]] <- frame$CMT
  if (!is.null(roles$dvid)) {
    for (column in roles$dvid) out[[column]] <- frame$DVID
  }
  if (!is.null(roles$mdv)) out[[roles$mdv]] <- as.integer(frame$EVID != 0L)
  if (!is.null(roles$cens)) out[[roles$cens]] <- 0L
  if (!is.null(roles$limit)) out[[roles$limit]] <- NA_real_
  if (!is.null(roles$rate)) {
    out[[roles$rate]] <- ifelse(frame$EVID != 0L,
                                fit$parameters$rate %||% 0, 0)
  }
  if (!is.null(roles$occasion)) {
    out[[roles$occasion]] <- unlist(lapply(
      split(frame, frame$.subject),
      function(part) pmax(1L, cumsum(part$EVID != 0L))
    ), use.names = FALSE)
  }
  if (!is.null(roles$tad)) out[[roles$tad]] <- .derived_tad(out, roles)
  if (!is.null(roles$assigned_dose)) {
    out[[roles$assigned_dose]] <- doses[frame$.subject]
  }
  for (column in schema$carried) {
    lookup <- stats::setNames(
      lapply(fit$arms$arms, function(arm) schema$arm_values[[arm]][[column]]),
      fit$arms$arms
    )
    out[[column]] <- .restore_column(
      unlist(lookup[frame$.arm], use.names = FALSE), schema$prototypes[[column]]
    )
  }
  for (covariate in roles$covariates) {
    # Kept in the type the draw produced: a numeric covariate stays numeric so
    # that `.restore_column()` can round it back to the source's own class.
    values <- unlist(lapply(subject_covariates, function(mine) {
      value <- mine[[covariate]]
      if (is.null(value) || !length(value)) NA else value
    }), use.names = FALSE)
    out[[covariate]] <- .restore_column(values[frame$.subject],
                                        schema$prototypes[[covariate]])
  }
  # The evaluated value is the latent one: what the subject would have measured
  # with no assay limit. DV, CENS and LIMIT are reconstructed from it together,
  # after any discrete endpoint has been snapped, so the three always agree.
  if (!is.null(roles$cens)) {
    for (endpoint_name in names(schema$censoring)) {
      censoring <- schema$censoring[[endpoint_name]]
      if (is.null(censoring)) next
      rows <- which(frame$EVID == 0L & frame$DVID == endpoint_name)
      if (!length(rows)) next
      out <- .censor_latent(out, rows, out[[roles$dv]][rows], roles,
                            public = censoring)
    }
  }

  out <- out[, intersect(schema$columns, names(out)), drop = FALSE]
  out <- out[order(out[[roles$id]], out[[roles$time]],
                   out[[roles$evid]] == 0L), , drop = FALSE]
  rownames(out) <- NULL

  # The floor keeps a value the assay could not have returned out of the
  # output. It is not a way of making a bad fit look reasonable, and a floor
  # that catches a large share of the output is doing exactly that, so the
  # share is reported rather than left in the figure for someone to notice.
  if (!is.null(floored) && floored$seen > 0L && floored$raised > 0L) {
    share <- floored$raised / floored$seen
    attr(out, "pmx_floored") <- c(raised = floored$raised, seen = floored$seen)
    if (share > 0.05) {
      warning(sprintf(
        paste0("%.0f%% of generated observations (%d of %d) fell below the ",
               "smallest value the study reported and were raised to half of ",
               "it. A floor catching this much is a fitted model that does ",
               "not describe the low end of the data, not an assay limit; ",
               "read `model_report()` before using this dataset."),
        100 * share, floored$raised, floored$seen), call. = FALSE)
    }
  }
  out
}
