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
# for synthetic data and far better conditioned to calibrate than an exposure-driven
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
      # removable singularities there and synthetic data does not need the limits.
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
#' @param rx Reserved for an `rxode2` model. **Not yet implemented**: the value
#'   is stored on the returned object but the generator always evaluates the
#'   built-in analytic solution, so supplying it warns. See `REV-020` in
#'   `design/REVIEW_BACKLOG.md`.
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
  if (!is.null(rx)) {
    if (!requireNamespace("rxode2", quietly = TRUE)) {
      stop("`rx` was supplied but the rxode2 package is not installed.",
           call. = FALSE)
    }
    # The generator never reads `model$rx`; every profile comes from the
    # built-in analytic solution. Silently returning an analytic curve for a
    # user-supplied ODE model would be a fidelity claim the package cannot
    # keep. See REV-020 in design/REVIEW_BACKLOG.md.
    warning(
      "`rx` is not yet used: profiles are always evaluated from the built-in ",
      "analytic `", pk, "` solution. Supplying an rxode2 model does not ",
      "change the generated data.", call. = FALSE
    )
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
#' Two dosing patterns are supported. A parallel design gives each cohort one
#' dose level (`dose_levels`) repeated `n_doses` times. A within-subject
#' escalation gives every subject the same increasing sequence of doses
#' (`dose_escalation`), one per occasion; this is prespecified design, not an
#' outcome, when the escalation follows a fixed protocol schedule.
#'
#' @param dose_levels Dose amounts, one per cohort. Omit when using
#'   `dose_escalation`.
#' @param cohort_sizes Planned subjects per cohort, recycled over `dose_levels`.
#'   Defaults to equal cohorts.
#' @param sampling Nominal sampling times after each dose, from the protocol.
#' @param n_doses Number of doses per subject in a parallel design.
#' @param dose_interval Time between doses when doses are equally spaced.
#' @param dose_escalation Per-occasion dose amounts for a within-subject
#'   escalation, for example `c(10, 30, 100)`. Applied to every subject. For a
#'   trial that also escalates between cohorts, supply a list of sequences, one
#'   per cohort, for example `list(c(1, 2, 4), c(2, 4, 8), c(4, 8, 16))`, with
#'   `cohort_sizes` giving the subjects in each. Every sequence must have the
#'   same length, since the dosing schedule is shared.
#' @param dose_times Explicit dose times, for example `c(0, 7, 14)`. Defaults to
#'   equally spaced times at `dose_interval`.
#' @param duration Infusion duration; zero for bolus or oral.
#' @param visit_window Fractional jitter applied to nominal times.
#' @param source Required provenance string.
#'
#' @return A `pmx_trial_design`.
#' @export
pmx_trial_design <- function(dose_levels = NULL, cohort_sizes = NULL,
                             sampling, n_doses = 1L, dose_interval = 24,
                             dose_escalation = NULL, dose_times = NULL,
                             duration = 0, visit_window = 0.05, source) {
  if (missing(source) || !is.character(source) || length(source) != 1L ||
      !nzchar(trimws(source))) {
    stop("`source` is required and must record the protocol this came from.",
         call. = FALSE)
  }
  escalating <- !is.null(dose_escalation)
  escalation_list <- NULL
  if (escalating) {
    if (!is.null(dose_levels)) {
      stop("Supply either `dose_levels` (parallel) or `dose_escalation` ",
           "(within-subject), not both.", call. = FALSE)
    }
    # One sequence applied to all subjects, or a list of sequences, one per
    # cohort (e.g. list(c(1, 2, 4), c(2, 4, 8)) for two escalation arms).
    sequences <- if (is.list(dose_escalation)) {
      dose_escalation
    } else {
      list(dose_escalation)
    }
    ok <- vapply(sequences, function(s) {
      is.numeric(s) && length(s) >= 1L && !anyNA(s) && all(is.finite(s)) &&
        all(s > 0)
    }, logical(1))
    if (!length(sequences) || !all(ok)) {
      stop("`dose_escalation` must be a vector of positive dose amounts, or a ",
           "list of such vectors, one per cohort.", call. = FALSE)
    }
    if (length(unique(lengths(sequences))) != 1L) {
      stop("Every cohort's escalation sequence must have the same number of ",
           "doses, because the dosing schedule is shared.", call. = FALSE)
    }
    escalation_list <- lapply(sequences, as.numeric)
    n_doses <- length(escalation_list[[1L]])
    # One cohort per sequence; its starting dose drives cohort assignment.
    dose_levels <- vapply(escalation_list, `[`, numeric(1), 1L)
    if (is.null(cohort_sizes)) cohort_sizes <- 1L
  } else {
    if (!is.numeric(dose_levels) || !length(dose_levels) ||
        any(!is.finite(dose_levels)) || any(dose_levels <= 0)) {
      stop("`dose_levels` must be finite positive numbers.", call. = FALSE)
    }
    if (is.null(cohort_sizes)) cohort_sizes <- 1L
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
  n_doses <- max(1L, as.integer(n_doses))
  if (!is.null(dose_times)) {
    if (!is.numeric(dose_times) || length(dose_times) != n_doses ||
        any(!is.finite(dose_times)) || any(dose_times < 0) ||
        is.unsorted(dose_times, strictly = TRUE)) {
      stop("`dose_times` must give ", n_doses,
           " increasing nonnegative times.", call. = FALSE)
    }
    dose_times <- as.numeric(dose_times)
  }
  structure(list(
    dose_levels = as.numeric(dose_levels), cohort_sizes = cohort_sizes,
    escalation = escalation_list,
    sampling = as.numeric(sampling),
    n_doses = n_doses,
    dose_interval = as.numeric(dose_interval),
    dose_times = dose_times,
    duration = as.numeric(duration),
    visit_window = as.numeric(visit_window),
    source = source
  ), class = "pmx_trial_design")
}

#' @export
print.pmx_trial_design <- function(x, ...) {
  cat("Public trial design\n", sep = "")
  if (!is.null(x$escalation)) {
    if (length(x$escalation) == 1L) {
      cat("  within-subject escalation: ",
          paste(x$escalation[[1L]], collapse = " -> "), "\n", sep = "")
    } else {
      cat("  within-subject escalation, ", length(x$escalation), " cohorts:\n",
          sep = "")
      for (c in seq_along(x$escalation)) {
        cat("    cohort ", c, " (n = ", x$cohort_sizes[c], "): ",
            paste(x$escalation[[c]], collapse = " -> "), "\n", sep = "")
      }
    }
  } else {
    cat("  doses: ", paste(x$dose_levels, collapse = ", "),
        "  (n = ", paste(x$cohort_sizes, collapse = ", "), ")\n", sep = "")
  }
  cat("  dose times: ", paste(signif(.design_dose_times(x), 3),
                              collapse = ", "), "\n",
      "  sampling: ", paste(signif(x$sampling, 3), collapse = ", "), "\n",
      "  source: ", x$source, "\n", sep = "")
  invisible(x)
}

# Dose times for one subject under a design.
.design_dose_times <- function(design) {
  design$dose_times %||% ((seq_len(design$n_doses) - 1L) * design$dose_interval)
}

# Per-occasion dose amounts for one subject in a given cohort. Equal doses for a
# parallel cohort; the cohort's own increasing sequence under escalation.
.design_dose_amounts <- function(design, cohort = 1L) {
  if (!is.null(design$escalation)) {
    design$escalation[[cohort]]
  } else {
    rep(design$dose_levels[cohort], design$n_doses)
  }
}

# Nominal observation times: the protocol grid repeated after each dose, then
# deduplicated. This is a public schedule, not an inferred one.
.design_observation_times <- function(design) {
  times <- unlist(lapply(.design_dose_times(design),
                         function(d) d + design$sampling), use.names = FALSE)
  sort(unique(round(times, 8)))
}
