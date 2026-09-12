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
# prior-algorithm.Rmd, Step 5: Draw One Parameter Vector per Subject
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
# prior-algorithm.Rmd, Step 6: Place Dose and Observation Events on a Clock
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

# prior-algorithm.Rmd, Step 4: Assign Subjects to Cohorts
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

# A mixed-route model needs to know which dose went by which route, and the
# only place that can come from is the data: `pmx_trial_design()` has no route
# argument, and `.generate_structural()` calls `.pk_profile()` without one, so
# every dose would default to intravenous. For a `_mixed` form that also sets
# `f` to 1, silently discarding a bioavailability the constructor required --
# the output was byte-identical to the plain intravenous form at any `f`
# whatever (`REV-053`). `synpmx_model()` is the path that reads routes from a
# study and can use these forms.
.reject_mixed_route_model <- function(model, fn) {
  if (!inherits(model, "pmx_structural_model")) return(invisible(NULL))
  if (!model$pk %in% c("1cmt_mixed", "2cmt_mixed")) return(invisible(NULL))
  stop(.condition_text(
    "`", fn, "()` cannot generate from the `", model$pk, "` model.",
    why = paste(
      "A mixed-route model needs to know which dose went by which route, and",
      "a public trial design has no way to say. Declare the route the study",
      "actually used --", paste(sub("_mixed", "_iv", model$pk), "or",
                                sub("_mixed", "_oral", model$pk)),
      "-- or use `synpmx_model()`, which reads the route from the data."
    )), call. = FALSE)
}

# Taking a `pmx_roles()` means honouring it. Every generator that reads a study
# returns a table satisfying the roles it was handed, and these two must too --
# so a declaration naming a column they cannot produce is refused here rather
# than answered with a table that fails `validate_pmx()` under the caller's own
# declaration. The default roles name only columns the generator fills, so this
# is reachable only from a study's own declaration.
# prior-algorithm.Rmd, Step 3: Check the Output Declaration Can Be Filled
.reject_unfillable_roles <- function(roles, covariates) {
  # Column-naming roles with no counterpart in a structurally generated table.
  # `exclude` names columns to drop, and `routes`, `dose_endpoints` and
  # `endpoint_types` are mappings rather than columns, so none belong here.
  unfillable <- c("limit", "addl", "ii", "adm", "dose_covariate", "strata",
                  "keep")
  named <- unlist(roles[intersect(unfillable, names(roles))],
                  use.names = FALSE)
  wanted <- roles$covariates
  produced <- if (is.null(covariates)) character() else names(covariates)
  absent <- setdiff(wanted, produced)
  if (!length(named) && !length(absent)) return(invisible(NULL))

  why <- character()
  if (length(absent)) {
    why <- c(why, paste0(
      "The covariate(s) ", paste(absent, collapse = ", "),
      " are declared but not generated. These modes draw a covariate only ",
      "from a public distribution, so pass `covariates = pmx_covariates(",
      absent[[1L]], " = pmx_covariate(...))` as well."
    ))
  }
  if (length(named)) {
    why <- c(why, paste0(
      "The column(s) ", paste(named, collapse = ", "),
      " are declared and cannot be produced from a structural model and a ",
      "protocol. Leave those roles out of the declaration passed here."
    ))
  }
  stop(.condition_text(
    "`roles` declares columns this generator cannot fill.",
    why = paste(why, collapse = " ")), call. = FALSE)
}

# The generated table is built under the names in `pmx_generated_roles()`, and
# this renames it to whatever the caller declared. It is what makes
# `synpmx_prior()` and `synpmx_calibrated()` hand back the same thing
# `synpmx_avatar()`, `synpmx_pca()` and `synpmx_model()` do: a table wearing
# the study's own column names, so one `pmx_roles()` serves the generator,
# `compare_pmx_distributions()` and `synpmx_scorecard()` alike.
#
# Two consequences of adopting that convention rather than inventing one. A
# role the caller did not declare gets no column, because a study with no TAD
# column should not be handed one back -- and with the default roles every slot
# is named, so nothing is dropped. And `dvid` may name several columns, of
# which the first is used: the generated table holds one endpoint key and
# nothing here can split it across more.
# prior-algorithm.Rmd, Step 11: Name the Columns from the Declaration
.rename_to_roles <- function(out, roles) {
  if (!inherits(roles, "pmx_roles")) {
    stop("`roles` must come from `pmx_roles()`.", call. = FALSE)
  }
  generated <- pmx_generated_roles()
  # Walk the table's own columns rather than the role list, so the output keeps
  # the column order the generator produced whatever `roles` is.
  slot_of <- character()
  for (slot in names(generated)) {
    name <- generated[[slot]]
    if (length(name) == 1L) slot_of[[name]] <- slot
  }
  # `pmx_roles()` permits exactly one collision: `cmt` and `dvid` naming the
  # same column, because NONMEM's CMT routinely does both jobs -- the dosing
  # compartment on event rows, the endpoint key on observation rows. The
  # generated table keeps the two apart, and its CMT already tells the
  # endpoints apart (1 dose, 2 cp, 3 pd), so a shared declaration takes CMT
  # and the separate endpoint-key column is dropped rather than written twice.
  skip <- if (length(intersect(roles$cmt, roles$dvid))) generated$dvid
  targets <- character()
  values <- list()
  for (name in names(out)) {
    if (name %in% skip) next
    slot <- if (name %in% names(slot_of)) slot_of[[name]] else NA_character_
    if (is.na(slot)) {
      # Covariate columns are named by `pmx_covariates()` and pass through as
      # they are; they fill no role slot this table knows about.
      to <- name
    } else {
      to <- roles[[slot]]
      if (!length(to)) next
      to <- to[[1L]]
    }
    targets <- c(targets, to)
    values <- c(values, list(out[[name]]))
  }
  collisions <- unique(targets[duplicated(targets)])
  if (length(collisions)) {
    stop("`roles` puts more than one generated column in ",
         paste(collisions, collapse = ", "),
         ". Give each role its own column name.", call. = FALSE)
  }
  names(values) <- targets
  as.data.frame(values, stringsAsFactors = FALSE, check.names = FALSE)
}

#' Generate a synthetic PMX event table from a structural model
#'
#' Works in two modes. Supplied a [pmx_structural_model()] it generates purely
#' from public inputs, reads no confidential data, and makes no privacy claim.
#' Supplied a [.fit_calibrated()] result it uses the privately corrected
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
#' @param covariates Optional [pmx_covariates()] for prior-mode generation.
#'   Ignored for a calibrated model, which carries its own released covariate
#'   summaries.
#' @param roles A [pmx_roles()] naming the columns of the output. Defaults to
#'   the roles a calibrated model was fitted under, and otherwise to
#'   [pmx_generated_roles()].
#'
#' @return A data frame in PMX event-table form, under the names in `roles`.
#' @keywords internal
.generate_structural <- function(x, design = NULL, n_subjects = NULL, seed = NULL,
                         dropout = 0, lloq = NULL, covariates = NULL,
                         roles = NULL) {
  if (inherits(x, "pmx_calibrated_model")) {
    model <- x$model
    design <- design %||% x$design
    typical <- x$corrected_typical
    n_subjects <- n_subjects %||% max(1L, as.integer(round(
      x$private_subject_count
    )))
    covariates <- x$covariates
    covariate_summaries <- x$covariate_summaries
    # The release carries the declaration it was fitted under. The fallback is
    # the schema generation used before it carried one.
    roles <- roles %||% x$roles %||% pmx_generated_roles()
  } else if (inherits(x, "pmx_structural_model")) {
    model <- x
    typical <- x$typical
    if (is.null(design)) {
      stop("`design` is required when generating from a structural model.",
           call. = FALSE)
    }
    n_subjects <- n_subjects %||% sum(design$cohort_sizes)
    if (!is.null(covariates) && !inherits(covariates, "pmx_covariates")) {
      stop("`covariates` must come from `pmx_covariates()`.", call. = FALSE)
    }
    covariate_summaries <- NULL
    roles <- roles %||% pmx_generated_roles()
  } else {
    stop("`x` must be a `pmx_structural_model` or a `pmx_calibrated_model`.",
         call. = FALSE)
  }
  if (!inherits(design, "pmx_trial_design")) {
    stop("`design` must come from `pmx_trial_design()`.", call. = FALSE)
  }
  if (!inherits(roles, "pmx_roles")) {
    stop("`roles` must come from `pmx_roles()`.", call. = FALSE)
  }
  .reject_mixed_route_model(model, "synpmx_prior")
  .reject_unfillable_roles(roles, covariates)
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
    p <- params[i, ]
    # Per-occasion dose amounts. Equal for a parallel cohort, increasing for a
    # within-subject escalation.
    doses <- .design_dose_amounts(design, cohort[i])
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
    # The assigned dose for an observation is the amount of its qualifying dose,
    # so it tracks the escalation occasion by occasion.
    obs_dose <- doses[clock$occasion]

    dose_rows <- data.frame(
      ID = i, TIME = dose_times, NTIME = dose_times, TAD = 0,
      OCC = seq_along(dose_times), DV = NA_real_, AMT = doses,
      RATE = if (design$duration > 0) doses / design$duration else 0,
      EVID = 1L, CMT = 1L, DVID = NA_character_, MDV = 1L,
      CENS = 0L, DOSE = doses, stringsAsFactors = FALSE
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
        CENS = 0L, DOSE = obs_dose, stringsAsFactors = FALSE
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
  # Baseline covariates: one value per subject, constant across their rows.
  if (!is.null(covariates)) {
    cov_table <- .draw_covariate_table(covariates, covariate_summaries,
                                      n_subjects)
    cov_table$ID <- seq_len(n_subjects)
    out <- merge(out, cov_table, by = "ID", sort = FALSE)
    out <- out[order(out$ID, out$TIME, out$EVID == 0L), , drop = FALSE]
  }
  out <- .rename_to_roles(out, roles)
  rownames(out) <- NULL
  attr(out, "pmx_source") <- if (inherits(x, "pmx_calibrated_model")) {
    "calibrated"
  } else "prior"
  out
}

#' Roles for tables produced by [.generate_structural()]
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
