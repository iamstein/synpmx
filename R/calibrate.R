# Public priors, pre-flight feasibility, and the calibrated release ----------
#
# The release is a multiplicative correction to the structural model's
# prediction, not an absolute parameter. The prior on "how wrong is the
# prediction" is far tighter than any prior on the parameter itself, and it is
# free. See design/PROTOTYPE_SPEC.md section 6.

#' Declare one public prior range
#'
#' The range is the dominant sensitivity term in the whole design, so its
#' provenance matters more than any other input. It must be chosen without
#' inspecting the confidential data.
#'
#' @param range Two increasing positive multipliers bracketing the correction
#'   factor, for example `c(1/4, 4)` for a prediction believed accurate to
#'   about four-fold.
#' @param source Required provenance string.
#'
#' @return A `pmx_prior`.
#' @export
pmx_prior <- function(range, source) {
  if (missing(source) || !is.character(source) || length(source) != 1L ||
      !nzchar(trimws(source))) {
    stop(
      "`source` is required. A prior range without data-independent ",
      "provenance voids the privacy guarantee, and nothing downstream ",
      "detects it.", call. = FALSE
    )
  }
  if (!is.numeric(range) || length(range) != 2L || anyNA(range) ||
      any(!is.finite(range)) || any(range <= 0) || range[1L] >= range[2L]) {
    stop("`range` must be two increasing positive multipliers.", call. = FALSE)
  }
  structure(list(range = as.numeric(range), source = source,
                 span = log(range[2L]) - log(range[1L])),
            class = "pmx_prior")
}

#' Collect public priors for the released corrections
#'
#' @param ... Named [pmx_prior()] objects. Recognized names are `pk` and `pd`.
#'
#' @return A `pmx_priors` object.
#' @export
pmx_priors <- function(...) {
  priors <- list(...)
  if (!length(priors) || is.null(names(priors)) || any(!nzchar(names(priors)))) {
    stop("`pmx_priors()` needs named `pmx_prior()` objects.", call. = FALSE)
  }
  if (!all(vapply(priors, inherits, logical(1), "pmx_prior"))) {
    stop("Every element must come from `pmx_prior()`.", call. = FALSE)
  }
  unknown <- setdiff(names(priors), c("pk", "pd"))
  if (length(unknown)) {
    stop("Unknown prior name(s): ", paste(unknown, collapse = ", "), ".",
         call. = FALSE)
  }
  structure(priors, class = "pmx_priors")
}

#' @export
print.pmx_prior <- function(x, ...) {
  cat(sprintf("prior [%.4g, %.4g] (%.3g-fold, span %.2f log units)\n  source: %s\n",
              x$range[1L], x$range[2L], x$range[2L] / x$range[1L], x$span,
              x$source))
  invisible(x)
}

#' Check whether a private release is worth its budget, before spending it
#'
#' Reports `f = d / (epsilon * N)`, the fraction of each prior's width that
#' survives as noise, and the resulting fold-error. This consumes no privacy
#' budget and reads no data: it depends only on the configuration.
#'
#' The fold-error is `exp(f * span)`, capped at the prior's half-width because
#' clipping prevents a release from landing outside the prior. The uncapped
#' form is accurate for `f` below roughly 0.25 and increasingly pessimistic
#' above it; see `design/FEASIBILITY.md` section 8.
#'
#' @param priors A [pmx_priors()] object.
#' @param epsilon The privacy budget under consideration.
#' @param n_subjects Number of subjects in the fit.
#'
#' @return A `pmx_preflight` report. See `design/PRIVACY_BACKGROUND.md`.
#' @export
pmx_preflight <- function(priors, epsilon, n_subjects) {
  if (!inherits(priors, "pmx_priors")) {
    stop("`priors` must come from `pmx_priors()`.", call. = FALSE)
  }
  if (!is.numeric(epsilon) || length(epsilon) != 1L || !is.finite(epsilon) ||
      epsilon <= 0) {
    stop("`epsilon` must be one finite positive value.", call. = FALSE)
  }
  if (!is.numeric(n_subjects) || length(n_subjects) != 1L ||
      !is.finite(n_subjects) || n_subjects < 1) {
    stop("`n_subjects` must be one positive number.", call. = FALSE)
  }
  d <- length(priors) + 1L          # + the subject count release
  f <- d / (epsilon * n_subjects)
  rows <- lapply(names(priors), function(name) {
    data.frame(
      quantity = name,
      prior_fold = priors[[name]]$range[2L] / priors[[name]]$range[1L],
      f = f,
      # Clipping to the prior bounds the damage: a release cannot land outside
      # the prior range, so the error saturates near the prior's half-width
      # however large f becomes. Without the cap the formula reports absurd
      # values in exactly the regime where the release is worthless anyway.
      expected_fold_error = min(exp(f * priors[[name]]$span),
                                exp(priors[[name]]$span / 2)),
      stringsAsFactors = FALSE
    )
  })
  table <- do.call(rbind, rows)
  verdict <- if (f >= 1) {
    "worthless"
  } else if (f > 0.5) {
    "marginal"
  } else if (f > 0.1) {
    "worthwhile"
  } else {
    "consider a smaller epsilon"
  }
  structure(list(
    d = d, epsilon = epsilon, n_subjects = n_subjects, f = f,
    verdict = verdict, table = table
  ), class = "pmx_preflight")
}

#' @export
print.pmx_preflight <- function(x, ...) {
  cat(sprintf("Pre-flight: d = %d, epsilon = %.4g, N = %.0f  ->  f = %.3f\n",
              x$d, x$epsilon, x$n_subjects, x$f))
  print(x$table, row.names = FALSE)
  cat("\nVerdict: ", x$verdict, "\n", sep = "")
  cat(switch(
    x$verdict,
    worthless = paste(
      "The noise is as wide as the prior. This release would tell you nothing",
      "you did not already assume.\nUse prior-mode generation instead, or",
      "raise epsilon only if governance allows.\n"),
    marginal = "The release beats the prior, but not by much.\n",
    worthwhile = "The release meaningfully narrows the prior.\n",
    paste("The prior contributes almost nothing. Consider spending less",
          "epsilon rather than banking accuracy you do not need.\n")
  ))
  invisible(x)
}

# Per-subject non-compartmental summaries -----------------------------------
#
# Each subject's value depends only on that subject's own rows. That is what
# bounds the sensitivity, and it is why an NLME fit cannot be used here:
# shrinkage would couple every subject to every other.

.trapezoid <- function(time, value) {
  keep <- is.finite(time) & is.finite(value)
  time <- time[keep]; value <- value[keep]
  if (length(time) < 2L) return(NA_real_)
  o <- order(time); time <- time[o]; value <- value[o]
  sum(diff(time) * (utils::head(value, -1L) + utils::tail(value, -1L)) / 2)
}

# Ratio of predicted to observed exposure for one subject, on the subject's own
# observation times. Working with a ratio on a shared grid avoids needing F or
# an extrapolation to infinity, both of which add assumptions without adding
# information.
.subject_corrections <- function(data, roles, model, design) {
  id <- data[[roles$id]]
  subjects <- .unique_in_order(id[!is.na(id)])
  dvid <- if (is.null(roles$dvid)) NULL else data[[roles$dvid]]
  out <- lapply(subjects, function(s) {
    rows <- which(!is.na(id) & id == s)
    piece <- data[rows, , drop = FALSE]
    evid <- suppressWarnings(as.numeric(piece[[roles$evid]]))
    amt <- if (is.null(roles$amt)) rep(NA_real_, nrow(piece)) else
      suppressWarnings(as.numeric(piece[[roles$amt]]))
    time <- suppressWarnings(as.numeric(piece[[roles$time]]))
    dose_rows <- which(is.finite(evid) & evid != 0 & is.finite(amt) & amt > 0)
    if (!length(dose_rows)) return(NULL)
    doses <- amt[dose_rows]; dose_times <- time[dose_rows]

    endpoint <- if (is.null(roles$dvid)) rep("cp", nrow(piece)) else
      as.character(piece[[roles$dvid]])
    obs <- is.finite(evid) & evid == 0
    dv <- suppressWarnings(as.numeric(piece[[roles$dv]]))

    pk_rows <- which(obs & endpoint == "cp" & is.finite(dv))
    pk <- NA_real_
    if (length(pk_rows) >= 3L) {
      auc_obs <- .trapezoid(time[pk_rows], dv[pk_rows])
      pred <- .pk_profile(model, time[pk_rows], doses, dose_times,
                          duration = design$duration)
      auc_pred <- .trapezoid(time[pk_rows], pred)
      if (is.finite(auc_obs) && auc_obs > 0 &&
          is.finite(auc_pred) && auc_pred > 0) {
        # CL is inversely proportional to AUC, so this ratio is the CL
        # correction directly.
        pk <- auc_pred / auc_obs
      }
    }

    pd <- NA_real_
    if (model$pd != "none") {
      pd_rows <- which(obs & endpoint == "pd" & is.finite(dv))
      if (length(pd_rows) >= 3L) {
        base <- model$typical[["baseline"]]
        eff_obs <- max(abs(dv[pd_rows] - base))
        pred <- .pd_profile(model, time[pd_rows], doses, dose_times,
                            duration = design$duration)
        eff_pred <- max(abs(pred - base))
        if (is.finite(eff_obs) && eff_obs > 0 &&
            is.finite(eff_pred) && eff_pred > 0) {
          pd <- eff_obs / eff_pred
        }
      }
    }
    c(pk = pk, pd = pd)
  })
  out <- out[!vapply(out, is.null, logical(1))]
  if (!length(out)) {
    stop("No subject yielded a usable non-compartmental summary. Check the ",
         "role declarations and that dose rows carry a positive amount.",
         call. = FALSE)
  }
  do.call(rbind, out)
}

#' Calibrate a public structural model to confidential data
#'
#' The only stage that reads source data. Each subject is reduced to bounded
#' multiplicative corrections against the structural model's own prediction,
#' clipped to public prior ranges, and released through a validated
#' differential-privacy backend.
#'
#' @param data Confidential PMX event data.
#' @param roles Column roles from [pmx_roles()].
#' @param model A public [pmx_structural_model()].
#' @param design A public [pmx_trial_design()].
#' @param priors Public [pmx_priors()] for each released correction.
#' @param epsilon Requested subject-level privacy budget.
#' @param backend `"opendp"`, or `"public"` for an explicitly public fixture.
#' @param public_source Logical assertion that the input is already public.
#'
#' @return A `pmx_calibrated_model`, carrying corrected typical parameters,
#'   accounting, provenance, and a release ledger. It contains no raw records.
#' @export
fit_calibrated_pmx <- function(data, roles, model, design, priors, epsilon,
                               backend = "opendp", public_source = FALSE) {
  if (!is.data.frame(data) || !nrow(data)) {
    stop("`data` must be a nonempty data frame.", call. = FALSE)
  }
  .assert_roles(data, roles)
  if (!inherits(model, "pmx_structural_model")) {
    stop("`model` must come from `pmx_structural_model()`.", call. = FALSE)
  }
  if (!inherits(design, "pmx_trial_design")) {
    stop("`design` must come from `pmx_trial_design()`.", call. = FALSE)
  }
  if (!inherits(priors, "pmx_priors")) {
    stop("`priors` must come from `pmx_priors()`.", call. = FALSE)
  }
  if (!is.numeric(epsilon) || length(epsilon) != 1L || !is.finite(epsilon) ||
      epsilon <= 0) {
    stop("`epsilon` must be supplied explicitly as one finite positive value.",
         call. = FALSE)
  }
  if (model$pd == "none" && !is.null(priors$pd)) {
    stop("A `pd` prior was supplied but the model has no PD component.",
         call. = FALSE)
  }
  resolved <- .resolve_backend(backend, public_source)

  corrections <- .subject_corrections(data, roles, model, design)
  n_subjects <- nrow(corrections)
  accountant <- .new_accountant(epsilon, 0, resolved)

  released <- names(priors)
  d <- length(released) + 1L
  per_query <- epsilon / d

  count <- .private_release(accountant, "subject_count", n_subjects,
                            sensitivity = 1, epsilon = per_query)
  private_count <- max(as.numeric(count), 1)

  corrected <- model$typical
  results <- list()
  for (name in released) {
    prior <- priors[[name]]
    bounds <- log(prior$range)
    values <- corrections[, name]
    values <- values[is.finite(values) & values > 0]
    # Clipping to the public prior is what creates the sensitivity bound.
    unit <- if (length(values)) .to_unit(log(values), bounds) else numeric()
    total <- .private_release(
      accountant, paste0(name, "_correction"), sum(unit),
      sensitivity = 1, epsilon = per_query
    )
    mean_unit <- min(max(as.numeric(total) / private_count, 0), 1)
    factor <- exp(.from_unit(mean_unit, bounds))
    # A correction pressed against its clipping boundary means the prior was
    # wrong and the release is censored, not that the data said this.
    edge <- min(mean_unit, 1 - mean_unit) < 0.02
    results[[name]] <- list(
      factor = factor, at_prior_boundary = edge,
      prior = prior, n_used = length(values)
    )
  }

  if (!is.null(results$pk)) {
    corrected[["cl"]] <- unname(corrected[["cl"]] * results$pk$factor)
  }
  if (!is.null(results$pd)) {
    corrected[["emax"]] <- unname(min(
      corrected[["emax"]] * results$pd$factor, 0.99
    ))
  }

  accounting <- .finalize_accounting(accountant)
  formal_dp <- isTRUE(resolved$validated) && isTRUE(resolved$production)

  warnings <- character()
  f <- d / (epsilon * n_subjects)
  if (f >= 1) {
    warnings <- c(warnings, paste0(
      "f = ", signif(f, 3), ": the noise is as wide as the prior, so this ",
      "release conveys nothing beyond it. Prior-mode generation would give ",
      "the same output at no privacy cost."
    ))
  }
  for (name in names(results)) {
    if (results[[name]]$at_prior_boundary) {
      warnings <- c(warnings, paste0(
        "The ", name, " correction is pressed against its prior boundary. ",
        "The prior is probably wrong and the release is censored; the ",
        "generated data reflects the boundary, not the study."
      ))
    }
  }

  out <- structure(list(
    version = 3L,
    engine = "calibrated_structural_generator",
    model = model, design = design, priors = priors,
    corrections = results,
    corrected_typical = corrected,
    private_subject_count = private_count,
    preflight = pmx_preflight(priors, epsilon, n_subjects),
    privacy = list(
      formal_dp = formal_dp,
      unit = "one subject's complete bounded longitudinal contribution",
      adjacency = "add-or-remove one complete subject",
      epsilon = as.numeric(epsilon), delta = 0,
      backend = list(name = resolved$name, version = resolved$version,
                     mechanism = resolved$mechanism,
                     validated = isTRUE(resolved$validated),
                     production = isTRUE(resolved$production)),
      accounting = accounting,
      proof_assumptions = c(
        "The structural model, its typical parameters, the trial design, and every prior range were established independently of the confidential data.",
        "Each subject's correction is computed only from that subject's own rows, so clipping to the public prior bounds the per-subject L1 sensitivity at one.",
        "OpenDP's Laplace measurement and privacy map are correct for the stated sensitivities.",
        "Basic sequential composition covers every released source-dependent computation in this fit.",
        "Realized trial-design quantities recorded as public were separately disclosed or are treated as public by assertion.",
        "Generation after fitting consults only this released model and public inputs."
      )
    ),
    provenance = data.frame(
      input = c("structural model", "trial design",
                paste(names(priors), "prior")),
      source = c(model$source, design$source,
                 vapply(priors, `[[`, character(1), "source")),
      stringsAsFactors = FALSE
    ),
    ledger = list(
      release_id = .release_id(),
      created_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"),
      requested_epsilon = epsilon,
      realized_epsilon = accounting$realized_epsilon,
      backend = resolved$name
    ),
    warnings = unique(warnings)
  ), class = "pmx_calibrated_model")

  if (length(out$warnings)) {
    warning(paste(out$warnings, collapse = "\n"), call. = FALSE)
  }
  out
}

#' @export
print.pmx_calibrated_model <- function(x, ...) {
  cat("Calibrated structural model (v3)\n")
  cat("  released subject count: ", round(x$private_subject_count, 1), "\n",
      sep = "")
  for (name in names(x$corrections)) {
    cat(sprintf("  %s correction: %.3gx%s\n", name,
                x$corrections[[name]]$factor,
                if (x$corrections[[name]]$at_prior_boundary)
                  "  [AT PRIOR BOUNDARY]" else ""))
  }
  cat("  corrected typical: ",
      paste(names(x$corrected_typical), signif(x$corrected_typical, 3),
            sep = "=", collapse = ", "), "\n", sep = "")
  cat("  epsilon: ", x$privacy$accounting$realized_epsilon,
      "  (formal DP: ", x$privacy$formal_dp, ")\n", sep = "")
  cat("  f = ", signif(x$preflight$f, 3), " (", x$preflight$verdict, ")\n",
      sep = "")
  invisible(x)
}
