# Public structural models --------------------------------------------------
#
# The structural model is a public input: it exists independently of the
# confidential dataset, so it consumes no privacy budget. See
# design/PROTOTYPE_SPEC.md section 6.
#
# Built-in analytic solutions cover the common linear cases and need no
# compiler. An rxode2 model may be supplied instead for anything else.

.pk_models <- c("1cmt_iv", "1cmt_oral", "1cmt_infusion")
.pd_models <- c("none", "direct_emax", "idr_inhibit_loss", "idr_stimulate_loss")

.required_pk_params <- list(
  `1cmt_iv`       = c("cl", "v"),
  `1cmt_oral`     = c("cl", "v", "ka"),
  `1cmt_infusion` = c("cl", "v")
)

.required_pd_params <- list(
  none                 = character(),
  direct_emax          = c("baseline", "emax", "ec50"),
  idr_inhibit_loss     = c("baseline", "emax", "ec50", "kout"),
  idr_stimulate_loss   = c("baseline", "emax", "ec50", "kout")
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

# Turnover PD needs integration. A fixed-step RK4 on a grid dense relative to
# the response half-life is ample for mock data and avoids a solver dependency.
.pd_profile <- function(model, time, doses, dose_times, params = NULL,
                        duration = 0, steps_per_unit = 4) {
  p <- params %||% model$typical
  if (model$pd == "none") return(rep(NA_real_, length(time)))
  conc_at <- function(tt) {
    .pk_profile(model, tt, doses, dose_times, p, duration)
  }
  drive <- function(cp) {
    p[["emax"]] * cp / (p[["ec50"]] + cp)
  }
  if (model$pd == "direct_emax") {
    return(p[["baseline"]] * (1 - drive(conc_at(time))))
  }
  horizon <- max(c(time, dose_times), na.rm = TRUE)
  n_steps <- max(50L, as.integer(ceiling(horizon * steps_per_unit)))
  grid <- seq(0, horizon, length.out = n_steps + 1L)
  h <- grid[2L] - grid[1L]
  kout <- p[["kout"]]
  kin <- p[["baseline"]] * kout
  sign <- if (model$pd == "idr_inhibit_loss") -1 else 1
  deriv <- function(tt, r) kin - kout * (1 + sign * drive(conc_at(tt))) * r
  r <- numeric(length(grid))
  r[1L] <- p[["baseline"]]
  for (i in seq_len(n_steps)) {
    tt <- grid[i]; y <- r[i]
    k1 <- deriv(tt, y)
    k2 <- deriv(tt + h / 2, y + h / 2 * k1)
    k3 <- deriv(tt + h / 2, y + h / 2 * k2)
    k4 <- deriv(tt + h, y + h * k3)
    r[i + 1L] <- max(y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4), 0)
  }
  stats::approx(grid, r, xout = pmin(pmax(time, 0), horizon), rule = 2)$y
}

#' Declare a public structural model
#'
#' The structural model and its typical parameter values are public inputs:
#' they must be established without inspecting the confidential dataset. For a
#' first-in-human compound they normally come from preclinical allometric
#' scaling, which is also what selected the starting dose.
#'
#' @param pk One of `"1cmt_iv"`, `"1cmt_oral"`, `"1cmt_infusion"`.
#' @param typical Named numeric vector of typical parameter values. Requires
#'   `cl` and `v`, plus `ka` for the oral model. Optional `f` defaults to 1.
#'   PD models additionally require `baseline`, `emax`, `ec50`, and `kout`.
#' @param pd One of `"none"`, `"direct_emax"`, `"idr_inhibit_loss"`,
#'   `"idr_stimulate_loss"`.
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
