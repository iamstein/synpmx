# Generation from a structural model ----------------------------------------
#
# Shared by prior mode and calibrated mode. The only difference between them is
# where the typical parameters came from; everything below is post-processing
# and consumes no privacy budget.

# `typical` is the median (equivalently the geometric mean), following the usual
# population-PK convention for a lognormal parameter. This also has to match
# what the calibration estimates: the released correction is a mean on the log
# scale, so it targets the geometric mean. Centering the arithmetic mean here
# instead would leave a systematic exp(sigma^2 / 2) gap between what is fitted
# and what is generated, which does not shrink with N or epsilon.
.draw_subject_params <- function(typical, iiv, n) {
  out <- matrix(rep(typical, each = n), nrow = n,
                dimnames = list(NULL, names(typical)))
  for (name in intersect(names(iiv), colnames(out))) {
    cv <- iiv[[name]]
    if (!is.finite(cv) || cv <= 0) next
    out[, name] <- out[, name] * stats::rlnorm(n, 0, sqrt(log(1 + cv^2)))
  }
  out
}

# Jitter is applied to time *since the qualifying dose*, not to absolute time.
# That keeps a predose sample exactly predose, keeps every sample within its own
# occasion, and guarantees a nonnegative TAD. Jittering absolute time does none
# of those things and reproduces defect SIM-005.
.structural_clock <- function(nominal, dose_times, window) {
  occ <- pmax(findInterval(nominal, dose_times), 1L)
  anchor <- dose_times[occ]
  tad <- nominal - anchor
  if (is.finite(window) && window > 0) {
    tad <- tad * (1 + stats::runif(length(tad), -window, window))
  }
  # Never before the qualifying dose, never at or past the next one.
  next_dose <- c(dose_times[-1L], Inf)[occ]
  headroom <- (next_dose - anchor) * (1 - sqrt(.Machine$double.eps))
  tad <- pmin(pmax(tad, 0), headroom)
  list(actual = anchor + tad, tad = tad, occasion = occ)
}

.assign_cohorts <- function(design, n_subjects) {
  weights <- design$cohort_sizes / sum(design$cohort_sizes)
  counts <- as.integer(round(weights * n_subjects))
  short <- n_subjects - sum(counts)
  if (short != 0L) {
    order_idx <- order(weights, decreasing = short > 0)
    for (i in seq_len(abs(short))) {
      j <- order_idx[(i - 1L) %% length(counts) + 1L]
      counts[j] <- counts[j] + sign(short)
    }
  }
  counts <- pmax(counts, 0L)
  rep(seq_along(design$dose_levels), times = counts)[seq_len(n_subjects)]
}

#' Generate a mock PMX event table from a structural model
#'
#' Works in two modes. Supplied a [pmx_structural_model()] it generates purely
#' from public inputs, reads no confidential data, and makes no privacy claim.
#' Supplied a [fit_calibrated_pmx()] result it uses the privately corrected
#' parameters; that is post-processing and consumes no further budget.
#'
#' @param x A `pmx_structural_model` or a `pmx_calibrated_model`.
#' @param design A [pmx_trial_design()]. Taken from `x` when it is a calibrated
#'   model.
#' @param n_subjects Number of subjects. Defaults to the planned cohort total,
#'   or to the released private count for a calibrated model.
#' @param seed Ordinary generation seed. Unrelated to privacy noise.
#' @param dropout Fraction of subjects who discontinue early. A public
#'   assumption from the protocol.
#' @param lloq Lower limit of quantification. Observations below it are flagged
#'   `CENS = 1` with `DV` at the limit, following the Monolix convention.
#'
#' @return A data frame in PMX event-table form.
#' @export
pmx_generate <- function(x, design = NULL, n_subjects = NULL, seed = NULL,
                         dropout = 0, lloq = NULL) {
  if (inherits(x, "pmx_calibrated_model")) {
    model <- x$model
    design <- design %||% x$design
    typical <- x$corrected_typical
    n_subjects <- n_subjects %||% max(1L, as.integer(round(
      x$private_subject_count
    )))
  } else if (inherits(x, "pmx_structural_model")) {
    model <- x
    typical <- x$typical
    if (is.null(design)) {
      stop("`design` is required when generating from a structural model.",
           call. = FALSE)
    }
    n_subjects <- n_subjects %||% sum(design$cohort_sizes)
  } else {
    stop("`x` must be a `pmx_structural_model` or a `pmx_calibrated_model`.",
         call. = FALSE)
  }
  if (!inherits(design, "pmx_trial_design")) {
    stop("`design` must come from `pmx_trial_design()`.", call. = FALSE)
  }
  n_subjects <- as.integer(n_subjects)
  if (!is.finite(n_subjects) || n_subjects < 1L) {
    stop("`n_subjects` must be one positive integer.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  cohort <- .assign_cohorts(design, n_subjects)
  params <- .draw_subject_params(typical, model$iiv, n_subjects)
  dose_times <- .design_dose_times(design)
  endpoints <- model$endpoints

  pieces <- vector("list", n_subjects)
  for (i in seq_len(n_subjects)) {
    dose <- design$dose_levels[cohort[i]]
    p <- params[i, ]
    nominal <- .design_observation_times(design)
    clock <- .structural_clock(nominal, dose_times, design$visit_window)
    # Dropout truncates follow-up; it is a protocol assumption, not learned.
    if (dropout > 0 && stats::runif(1) < dropout) {
      keep <- clock$actual <= stats::runif(1) * max(clock$actual)
      if (sum(keep) >= 2L) {
        nominal <- nominal[keep]
        clock <- lapply(clock, `[`, keep)
      }
    }
    actual <- clock$actual
    doses <- rep(dose, length(dose_times))

    dose_rows <- data.frame(
      ID = i, TIME = dose_times, NTIME = dose_times, TAD = 0,
      OCC = seq_along(dose_times), DV = NA_real_, AMT = doses,
      RATE = if (design$duration > 0) doses / design$duration else 0,
      EVID = 1L, CMT = 1L, DVID = NA_character_, MDV = 1L,
      CENS = 0L, DOSE = dose, stringsAsFactors = FALSE
    )

    obs_rows <- list()
    for (ep in endpoints) {
      value <- if (ep == "cp") {
        .pk_profile(model, actual, doses, dose_times, p, design$duration)
      } else {
        .pd_profile(model, actual, doses, dose_times, p, design$duration)
      }
      if (model$residual_cv > 0) {
        value <- value * stats::rlnorm(
          length(value), 0, sqrt(log(1 + model$residual_cv^2))
        )
      }
      obs_rows[[ep]] <- data.frame(
        ID = i, TIME = actual, NTIME = nominal,
        TAD = clock$tad, OCC = clock$occasion,
        DV = value, AMT = 0, RATE = 0, EVID = 0L,
        CMT = if (ep == "cp") 2L else 3L, DVID = ep, MDV = 0L,
        CENS = 0L, DOSE = dose, stringsAsFactors = FALSE
      )
    }
    piece <- rbind(dose_rows, do.call(rbind, obs_rows))
    piece <- piece[order(piece$TIME, piece$EVID == 0L), , drop = FALSE]
    pieces[[i]] <- piece
  }

  out <- do.call(rbind, pieces)
  if (!is.null(lloq) && is.finite(lloq)) {
    below <- out$EVID == 0L & out$DVID == "cp" &
      is.finite(out$DV) & out$DV < lloq
    out$CENS[below] <- 1L
    out$DV[below] <- lloq
  }
  rownames(out) <- NULL
  attr(out, "pmx_source") <- if (inherits(x, "pmx_calibrated_model")) {
    "calibrated"
  } else "prior"
  out
}

#' Roles for tables produced by [pmx_generate()]
#'
#' @return A `pmx_roles` object matching the generated schema.
#' @export
pmx_generated_roles <- function() {
  pmx_roles(
    id = "ID", time = "TIME", nominal_time = "NTIME", tad = "TAD",
    occasion = "OCC", dv = "DV", amt = "AMT", rate = "RATE", evid = "EVID",
    cmt = "CMT", dvid = "DVID", mdv = "MDV", cens = "CENS",
    assigned_dose = "DOSE"
  )
}
