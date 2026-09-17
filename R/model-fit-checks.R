# pmxmodel-algorithm, Step 3: acceptance before selection (SIM-087).
# These are deliberately coarse screens for a generator, not inferential
# standard-error or identifiability criteria. No source rows survive them.
.model_check_parameters <- function(parameters) {
  p <- parameters$fixed
  o <- parameters$omega
  r <- unlist(parameters$residual[c("cv", "sd")], use.names = FALSE)
  if (!length(p) || any(!is.finite(p) | p <= 0))
    return("fixed effects are not finite and positive")
  if (length(r) != 1L || !is.finite(r) || r < 0)
    return("residual error is not finite and non-negative")
  if (!is.matrix(o) || nrow(o) != ncol(o) ||
      any(!is.finite(o)) || !isTRUE(all.equal(o, t(o))))
    return("between-subject covariance is not finite and symmetric")
  if (nrow(o) && (is.null(rownames(o)) ||
                 any(!rownames(o) %in% names(p)) ||
                 min(eigen(o, symmetric = TRUE, only.values = TRUE)$values) <
                   -1e-10 * max(1, max(abs(o)))))
    return("between-subject covariance is not usable for population draws")
  NULL
}

.model_check_convergence <- function(fit, estimation) {
  if (!identical(estimation, "saem")) {
    code <- fit$convergence
    bad_message <- any(grepl(
      "false convergence|iteration limit|evaluation limit|failed|failure",
      fit$message %||% "", ignore.case = TRUE))
    if (length(code) != 1L || !is.finite(code) || code != 0 || bad_message)
      return(paste0("optimizer did not report successful convergence",
                     if (length(fit$message)) paste0(": ", fit$message)))
    return(NULL)
  }
  # SAEM has no deterministic optimizer termination code. Compare two blocks
  # in the final fifth of its parameter history. A large continuing drift is
  # evidence against accepting the run; a flat history is not proof of an
  # identified model. The 50% threshold is a gross drift screen.
  h <- fit$parHistData
  if (!is.data.frame(h) || nrow(h) < 20L)
    return("SAEM has no usable parameter trajectory")
  if ("type" %in% names(h)) h <- h[h$type == "Unscaled", , drop = FALSE]
  columns <- setdiff(names(h)[vapply(h, is.numeric, logical(1))],
                     c("iter", "objf"))
  if (nrow(h) < 20L || !length(columns))
    return("SAEM has no usable parameter trajectory")
  tail <- utils::tail(as.matrix(h[, columns, drop = FALSE]),
                       max(10L, floor(nrow(h) / 5)))
  if (any(!is.finite(tail))) return("SAEM trajectory contains non-finite values")
  half <- floor(nrow(tail) / 2)
  a <- colMeans(tail[seq_len(half), , drop = FALSE])
  b <- colMeans(tail[-seq_len(half), , drop = FALSE])
  if (any(abs(b - a) / pmax(1, abs(a)) > 0.5))
    return("SAEM parameter trajectory is still drifting")
  NULL
}

# Monolix-style censoring: DV is the stated boundary and LIMIT, when
# present, closes the interval on the other side. Exact readings have equal
# bounds. Physical concentrations below a left limit are bounded below by zero.
.model_observation_bounds <- function(data) {
  cens <- data$CENS %||% rep(0, nrow(data))
  other <- data$LIMIT %||% rep(NA_real_, nrow(data))
  lower <- upper <- data$DV
  left <- cens == 1
  right <- cens == -1
  lower[left] <- ifelse(is.finite(other[left]), pmax(0, other[left]), 0)
  upper[right] <- ifelse(is.finite(other[right]), other[right], Inf)
  list(lower = lower, upper = upper, observed = cens == 0)
}

# The original, expanded estimation records supply doses and sample times.
# Four new subjects per source design, with independent population effects;
# never the source subjects' fitted effects. Local seeding makes the acceptance
# decision reproducible without consuming the caller's random-number stream.
.model_check_generation <- function(parameters, structural, data, weight = NULL) {
  if (!any(data$EVID == 0L & is.finite(data$DV)))
    return(list(failure = "no observed design available for generation checks", note = ""))
  .with_local_seed(104729L, {
    subjects <- split(data, factor(data$ID, levels = unique(data$ID)))
    etas <- .draw_random_effects(parameters$omega, 4L * length(subjects))
    source <- generated <- typical <- source_lower <- source_upper <- numeric()
    measured <- extreme <- logical()
    omitted_upper <- FALSE
    index <- 0L
    for (part in subjects) {
      obs <- part[part$EVID == 0L & is.finite(part$DV), , drop = FALSE]
      dose <- part[part$EVID != 0L & part$AMT > 0, , drop = FALSE]
      if (!nrow(obs)) next
      bounds <- .model_observation_bounds(obs)
      if (any(bounds$lower > bounds$upper))
        return(list(failure = "censoring bounds are inconsistent", note = ""))
      source_lower <- c(source_lower, bounds$lower)
      source_upper <- c(source_upper, bounds$upper)
      measured <- c(measured, bounds$observed)
      ceiling <- max(bounds$upper)
      if (!is.finite(ceiling)) omitted_upper <- TRUE
      rates <- if ("RATE" %in% names(dose)) dose$RATE else rep(0, nrow(dose))
      duration <- ifelse(rates > 0, dose$AMT / rates, 0)
      routes <- if (grepl("mixed$", structural)) {
        ifelse(dose$CMT == 1L, "extravascular",
                "iv")
      } else NULL
      base <- parameters$fixed
      if (!is.null(weight)) {
        value <- part[[weight$covariate]][1L]
        for (nm in intersect(names(base), names(.model_allometric_exponents)))
          base[nm] <- base[nm] *
            (value / weight$reference)^.model_allometric_exponents[nm]
      }
      profile <- function(p) .pk_profile(list(pk = structural), obs$TIME,
        dose$AMT, dose$TIME, p, duration, routes)
      center <- profile(base)
      if (any(!is.finite(center) | center < -1e-8))
        return(list(failure = "typical profile is non-finite or negative", note = ""))
      typical <- c(typical, pmax(center, 0))
      source <- c(source, obs$DV)
      for (j in seq_len(4L)) {
        index <- index + 1L
        p <- .subject_parameters(base, list(), etas[index, ], list())
        if (any(!is.finite(p) | p <= 0))
          return(list(failure = "population draws produce invalid parameters", note = ""))
        pred <- profile(p)
        if (any(!is.finite(pred) | pred < -1e-8))
          return(list(failure = "population draws produce non-finite or negative profiles", note = ""))
        y <- .add_residual_error(pmax(pred, 0), parameters$residual, floor = 0)
        if (any(!is.finite(y)))
          return(list(failure = "residual draws produce non-finite concentrations", note = ""))
        generated <- c(generated, y)
        if (is.finite(ceiling) && ceiling > 0)
          extreme <- c(extreme, max(y) > 100 * ceiling)
      }
    }
    # Bounds on each reading induce bounds on the source median. This keeps
    # partial left censoring from hiding a collapsed population, without
    # pretending the assay limits were observed concentrations (SIM-087).
    censored <- any(!measured)
    lower <- stats::median(source_lower)
    upper <- stats::median(source_upper)
    median_generated <- stats::median(generated)
    level <- if (!censored && stats::median(source) > 0)
      median_generated / stats::median(source) else NA_real_
    spread <- if (!censored && stats::IQR(source) > 0)
      stats::IQR(generated) / stats::IQR(source) else NA_real_
    positive <- measured & source > 0
    profile_error <- if (any(positive))
      stats::median(abs(log(pmax(typical[positive], .Machine$double.xmin) /
                             source[positive]))) else NA_real_
    failure <- NULL
    notes <- character()
    if (length(extreme) && mean(extreme) > 0.1)
      failure <- "over 10% of generated profiles exceed 100 times their source-design maximum"
    if ((lower > 0 && median_generated / lower < 0.01) ||
        (is.finite(upper) && upper > 0 && median_generated / upper > 100))
      failure <- "generated median is more than 100-fold outside the source median bounds"
    if (!censored) {
      if (!is.na(level) && (level < 0.5 || level > 2))
        notes <- c(notes, sprintf("generated/source median ratio %.3g", level))
      if (!is.na(spread) && (spread < 0.5 || spread > 2))
        notes <- c(notes, sprintf("generated/source interquartile-range ratio %.3g", spread))
    } else {
      notes <- "median checked against censoring bounds; spread comparison omitted for censored observations"
      if ((lower > 0 && median_generated / lower < 0.5) ||
          (is.finite(upper) && upper > 0 && median_generated / upper > 2))
        notes <- c(notes, "generated median is over twofold outside the source median bounds")
    }
    if (!is.na(profile_error) && profile_error > log(4))
      notes <- c(notes, "typical profile differs by over fourfold at the median positive uncensored observation")
    if (omitted_upper) notes <- c(notes,
      "upper-tail comparison omitted for designs whose right censoring leaves the source maximum unknown")
    list(failure = failure, note = paste(notes, collapse = "; "),
         median_ratio = level, spread_ratio = spread,
         source_median_bounds = c(lower = lower, upper = upper),
         generated_median = median_generated,
         extreme_fraction = if (length(extreme)) mean(extreme) else NA_real_)

  })
}

.model_assess_fit <- function(fit, candidate, start, error, estimation,
                              data, weight = NULL) {
  status <- .model_check_convergence(fit, estimation)
  converged <- is.null(status)
  # SIM-090: false convergence alone is a diagnostic, not a veto on usable
  # generation. Keep converged=FALSE in the report and run every other check.
  message <- paste(fit$message %||% "", collapse = "; ")
  false_only <- !identical(estimation, "saem") &&
    length(fit$convergence) == 1L && is.finite(fit$convergence) &&
    fit$convergence %in% c(0, 1, 8) &&
    grepl("false convergence", message, ignore.case = TRUE) &&
    !grepl("iteration limit|evaluation limit|failed|failure", message,
           ignore.case = TRUE)
  failure <- if (false_only) NULL else status
  parameters <- try(.model_read_fit(fit, candidate, error), silent = TRUE)
  if (inherits(parameters, "try-error")) {
    failure <- c(failure, "fitted parameters could not be read")
  } else {
    failure <- c(failure, .model_check_parameters(parameters))
  }
  notes <- if (false_only) paste(
    "Warning: the optimizer stalled before convergence was confirmed.",
    "Check that the synthetic data reasonably reproduces the source data's patterns and variability."
  ) else character()
  if (!length(failure)) {
    movement <- .model_fit_movement(parameters$fixed, parameters$omega, start)
    if (!movement$moved) notes <- c(notes, "all parameters remain within 1% of their starts")
    if (identical(movement$omega_moved, FALSE))
      notes <- c(notes, "between-subject variances remain within 1% of their starts")
    generation <- try(.model_check_generation(parameters, candidate, data, weight),
                        silent = TRUE)
    if (inherits(generation, "try-error")) {
      failure <- "generation check could not run"
    } else {
      failure <- generation$failure
      notes <- c(notes, generation$note)
    }
  }
  list(converged = converged, accepted = !length(failure),
       note = paste(c(failure, notes[nzchar(notes)]), collapse = "; "))
}
