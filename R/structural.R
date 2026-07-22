# Public structural models --------------------------------------------------
#
# The structural model is a public input: it exists independently of the
# confidential dataset, so it consumes no privacy budget. See
# design/PROTOTYPE_SPEC.md section 6.
#
# Built-in analytic solutions cover the common linear cases and need no
# compiler. An rxode2 model may be supplied instead for anything else.

.pk_models <- c("1cmt_iv", "1cmt_oral", "1cmt_infusion", "2cmt_iv",
                "2cmt_oral")
# PD is a simple time course with no exposure dependence. This is both adequate
# for mock data and far better conditioned to calibrate than an exposure-driven
# model, whose effect is a small deviation on a large baseline. See
# design/PROTOTYPE_SPEC.md section 6.
.pd_models <- c("none", "constant", "linear", "exponential")

.required_pk_params <- list(
  `1cmt_iv`       = c("cl", "v"),
  `1cmt_oral`     = c("cl", "v", "ka"),
  `1cmt_infusion` = c("cl", "v"),
  `2cmt_iv`       = c("cl", "v", "q", "v2"),
  `2cmt_oral`     = c("cl", "v", "q", "v2", "ka")
)

# Macro-constants for a two-compartment model. `v` is the central volume.
.two_cmt_rates <- function(p) {
  k10 <- p[["cl"]] / p[["v"]]
  k12 <- p[["q"]] / p[["v"]]
  k21 <- p[["q"]] / p[["v2"]]
  sum_k <- k10 + k12 + k21
  root <- sqrt(max(sum_k^2 - 4 * k10 * k21, 0))
  list(alpha = (sum_k + root) / 2, beta = (sum_k - root) / 2, k21 = k21)
}

.required_pd_params <- list(
  none                 = character(),
  constant             = "baseline",
  linear               = c("baseline", "slope"),
  exponential          = c("baseline", "plateau", "rate")
)

# Concentration from one dose, evaluated at times measured from that dose.
# Superposition is valid because every built-in model is linear in dose.
.pk_single_dose <- function(model, time, dose, p, duration = 0) {
  time <- pmax(as.numeric(time), 0)
  ke <- p[["cl"]] / p[["v"]]
  # `[[` on a named numeric vector errors for a missing name, so test first.
  f <- if ("f" %in% names(p)) p[["f"]] else 1
  switch(
    model,
    `1cmt_iv` = f * dose / p[["v"]] * exp(-ke * time),
    `1cmt_oral` = {
      ka <- p[["ka"]]
      # ka == ke is a removable singularity; use the limiting form.
      if (abs(ka - ke) < 1e-8) {
        f * dose * ka * time / p[["v"]] * exp(-ke * time)
      } else {
        f * dose * ka / (p[["v"]] * (ka - ke)) *
          (exp(-ke * time) - exp(-ka * time))
      }
    },
    `1cmt_infusion` = {
      if (!is.finite(duration) || duration <= 0) {
        return(f * dose / p[["v"]] * exp(-ke * time))
      }
      rate <- f * dose / duration
      during <- time <= duration
      out <- numeric(length(time))
      out[during] <- rate / p[["cl"]] * (1 - exp(-ke * time[during]))
      tail_start <- rate / p[["cl"]] * (1 - exp(-ke * duration))
      out[!during] <- tail_start * exp(-ke * (time[!during] - duration))
      out
    },
    `2cmt_iv` = {
      r <- .two_cmt_rates(p)
      gap <- r$alpha - r$beta
      if (gap < 1e-10) return(f * dose / p[["v"]] * exp(-r$beta * time))
      a <- f * dose / p[["v"]] * (r$alpha - r$k21) / gap
      b <- f * dose / p[["v"]] * (r$k21 - r$beta) / gap
      a * exp(-r$alpha * time) + b * exp(-r$beta * time)
    },
    `2cmt_oral` = {
      r <- .two_cmt_rates(p)
      ka <- p[["ka"]]
      # Nudge off any coincidence of rate constants; the closed form has
      # removable singularities there and mock data does not need the limits.
      eps <- 1e-6
      if (abs(ka - r$alpha) < eps) ka <- ka * (1 + eps)
      if (abs(ka - r$beta) < eps) ka <- ka * (1 - eps)
      if (abs(r$alpha - r$beta) < eps) r$alpha <- r$alpha * (1 + eps)
      coef <- f * dose * ka / p[["v"]]
      coef * (
        (r$k21 - r$alpha) / ((ka - r$alpha) * (r$beta - r$alpha)) *
          exp(-r$alpha * time) +
        (r$k21 - r$beta) / ((ka - r$beta) * (r$alpha - r$beta)) *
          exp(-r$beta * time) +
        (r$k21 - ka) / ((r$alpha - ka) * (r$beta - ka)) * exp(-ka * time)
      )
    },
    stop("Unknown PK model `", model, "`.", call. = FALSE)
  )
}

#' Concentration-time profile for a public structural model
#'
#' @param model A `pmx_structural_model`.
#' @param time Numeric times, measured from the first dose.
#' @param doses Numeric dose amounts.
#' @param dose_times Times at which those doses are given.
#' @param params Named parameter vector overriding the model's typical values.
#' @param duration Infusion duration, recycled over doses.
#'
#' @return Numeric concentrations, one per `time`.
#' @keywords internal
.pk_profile <- function(model, time, doses, dose_times, params = NULL,
                        duration = 0) {
  p <- params %||% model$typical
  time <- as.numeric(time)
  out <- numeric(length(time))
  duration <- rep_len(duration, length(doses))
  for (i in seq_along(doses)) {
    active <- time >= dose_times[i]
    if (!any(active)) next
    out[active] <- out[active] + .pk_single_dose(
      model$pk, time[active] - dose_times[i], doses[i], p, duration[i]
    )
  }
  out
}

# PD is a simple time course with no exposure dependence. This is deliberate:
# such a shape is adequate for exercising longitudinal analysis code, and it
# calibrates through a well-conditioned level correction where an exposure-driven
# deviation statistic does not. See design/PROTOTYPE_SPEC.md section 6. `doses`
# and `dose_times` are accepted for a common signature with `.pk_profile()` but
# are unused.
.pd_profile <- function(model, time, doses, dose_times, params = NULL,
                        duration = 0) {
  p <- params %||% model$typical
  if (model$pd == "none") return(rep(NA_real_, length(time)))
  time_positive <- pmax(as.numeric(time), 0)
  switch(
    model$pd,
    constant = rep(p[["baseline"]], length(time)),
    linear = p[["baseline"]] + p[["slope"]] * time_positive,
    # Covers both decay and growth: the sign is set by plateau vs baseline.
    exponential = p[["plateau"]] +
      (p[["baseline"]] - p[["plateau"]]) * exp(-p[["rate"]] * time_positive),
    stop("Unknown PD model `", model$pd, "`.", call. = FALSE)
  )
}

#' Declare a public structural model
#'
#' The structural model and its typical parameter values are public inputs:
#' they must be established without inspecting the confidential dataset. For a
#' first-in-human compound they normally come from preclinical allometric
#' scaling, which is also what selected the starting dose.
#'
#' @param pk One of `"1cmt_iv"`, `"1cmt_oral"`, `"1cmt_infusion"`,
#'   `"2cmt_iv"`, `"2cmt_oral"`.
#' @param typical Named numeric vector of typical parameter values, interpreted
#'   as the median of a lognormal population. Requires `cl` and `v` (the central
#'   volume), plus `ka` for oral models and `q` and `v2` for two-compartment
#'   models. Optional `f` defaults to 1. PD shapes additionally require
#'   `baseline`, plus `slope` for `"linear"` and `plateau` and `rate` for
#'   `"exponential"`.
#' @param pd The PD time course, with no exposure dependence. One of `"none"`
#'   (PK only), `"constant"`, `"linear"` (needs `slope`), or `"exponential"`
#'   (needs `plateau` and `rate`, covering both decay and growth). A simple
#'   shape is adequate for exercising longitudinal code and calibrates through a
#'   well-conditioned level correction; see `design/PROTOTYPE_SPEC.md`.
#' @param source Required provenance string recording where the model and its
#'   typical values came from. Recorded in the release ledger.
#' @param rx Optional `rxode2` model used in place of the built-in analytic
#'   solution. Requires the `rxode2` package.
#' @param iiv Named vector of between-subject variability, as CV on the log
#'   scale. A public assumption; it consumes no privacy budget.
#' @param residual_cv Proportional residual error, as a CV.
#'
#' @return A `pmx_structural_model`.
#' @export
pmx_structural_model <- function(pk, typical, pd = "none", source,
                                 rx = NULL, iiv = c(cl = 0.3, v = 0.2),
                                 residual_cv = 0.15) {
  pk <- match.arg(pk, .pk_models)
  pd <- match.arg(pd, .pd_models)
  if (missing(source) || !is.character(source) || length(source) != 1L ||
      !nzchar(trimws(source))) {
    stop(
      "`source` is required and must record where this model came from. ",
      "A structural model without data-independent provenance cannot be ",
      "treated as a public input.", call. = FALSE
    )
  }
  if (!is.numeric(typical) || is.null(names(typical)) ||
      anyNA(typical) || any(!is.finite(typical)) || any(typical <= 0)) {
    stop("`typical` must be a named vector of finite positive values.",
         call. = FALSE)
  }
  needed <- c(.required_pk_params[[pk]], .required_pd_params[[pd]])
  missing_params <- setdiff(needed, names(typical))
  if (length(missing_params)) {
    stop("`typical` is missing required parameters: ",
         paste(missing_params, collapse = ", "), ".", call. = FALSE)
  }
  if (!is.null(rx) && !requireNamespace("rxode2", quietly = TRUE)) {
    stop("`rx` was supplied but the rxode2 package is not installed.",
         call. = FALSE)
  }
  if (!is.numeric(residual_cv) || length(residual_cv) != 1L ||
      !is.finite(residual_cv) || residual_cv < 0) {
    stop("`residual_cv` must be one finite nonnegative number.", call. = FALSE)
  }
  structure(list(
    pk = pk, pd = pd, typical = typical, source = source, rx = rx,
    iiv = iiv, residual_cv = residual_cv,
    endpoints = c("cp", if (pd != "none") "pd")
  ), class = "pmx_structural_model")
}

#' @export
print.pmx_structural_model <- function(x, ...) {
  cat("Public structural model\n",
      "  PK: ", x$pk, "\n",
      "  PD: ", x$pd, "\n",
      "  typical: ",
      paste(names(x$typical), signif(x$typical, 3), sep = "=", collapse = ", "),
      "\n  source: ", x$source, "\n", sep = "")
  invisible(x)
}

#' Declare a public trial design
#'
#' Every field is a design fact from the protocol and consumes no privacy
#' budget. See `design/DATA_ELICITATION.md` for which parts of a realized
#' design need a provenance note.
#'
#' @param dose_levels Dose amounts, one per cohort.
#' @param cohort_sizes Planned subjects per cohort, recycled over `dose_levels`.
#' @param sampling Nominal sampling times after each dose, from the protocol.
#' @param n_doses Number of doses per subject.
#' @param dose_interval Time between doses when `n_doses > 1`.
#' @param duration Infusion duration; zero for bolus or oral.
#' @param visit_window Fractional jitter applied to nominal times.
#' @param source Required provenance string.
#'
#' @return A `pmx_trial_design`.
#' @export
pmx_trial_design <- function(dose_levels, cohort_sizes, sampling, n_doses = 1L,
                             dose_interval = 24, duration = 0,
                             visit_window = 0.05, source) {
  if (missing(source) || !is.character(source) || length(source) != 1L ||
      !nzchar(trimws(source))) {
    stop("`source` is required and must record the protocol this came from.",
         call. = FALSE)
  }
  if (!is.numeric(dose_levels) || !length(dose_levels) ||
      any(!is.finite(dose_levels)) || any(dose_levels <= 0)) {
    stop("`dose_levels` must be finite positive numbers.", call. = FALSE)
  }
  if (!is.numeric(sampling) || length(sampling) < 2L ||
      any(!is.finite(sampling)) || any(sampling < 0) ||
      is.unsorted(sampling, strictly = TRUE)) {
    stop("`sampling` must be at least two increasing nonnegative times.",
         call. = FALSE)
  }
  cohort_sizes <- rep_len(as.integer(cohort_sizes), length(dose_levels))
  if (any(!is.finite(cohort_sizes)) || any(cohort_sizes < 1L)) {
    stop("`cohort_sizes` must be positive integers.", call. = FALSE)
  }
  structure(list(
    dose_levels = as.numeric(dose_levels), cohort_sizes = cohort_sizes,
    sampling = as.numeric(sampling),
    n_doses = max(1L, as.integer(n_doses)),
    dose_interval = as.numeric(dose_interval),
    duration = as.numeric(duration),
    visit_window = as.numeric(visit_window),
    source = source
  ), class = "pmx_trial_design")
}

#' @export
print.pmx_trial_design <- function(x, ...) {
  cat("Public trial design\n",
      "  doses: ", paste(x$dose_levels, collapse = ", "),
      "  (n = ", paste(x$cohort_sizes, collapse = ", "), ")\n",
      "  ", x$n_doses, " dose(s), interval ", x$dose_interval, "\n",
      "  sampling: ", paste(signif(x$sampling, 3), collapse = ", "), "\n",
      "  source: ", x$source, "\n", sep = "")
  invisible(x)
}

# Dose times for one subject under a design.
.design_dose_times <- function(design) {
  (seq_len(design$n_doses) - 1L) * design$dose_interval
}

# Nominal observation times: the protocol grid repeated after each dose, then
# deduplicated. This is a public schedule, not an inferred one.
.design_observation_times <- function(design) {
  times <- unlist(lapply(.design_dose_times(design),
                         function(d) d + design$sampling), use.names = FALSE)
  sort(unique(round(times, 8)))
}
