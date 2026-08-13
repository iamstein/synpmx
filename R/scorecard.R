# The scorecard ---------------------------------------------------------------
#
# `vignette("synthetic-data-checks")` states every check and the answer that
# counts as passing for it. This runs the ones that can be run and collects the
# answers into one table, so that deciding whether to ship a dataset does not
# require rereading the document or rewriting sixty lines of reporting code per
# study.
#
# Three rows deliberately have no pass mark. Doses per patient is the clearest:
# on a study with individualised dosing it can halve while every guarantee above
# it still reads 0, and any threshold picked for it would be wrong on some
# dataset. `"review"` is the honest verdict there, not a soft `"pass"`.

.scorecard_row <- function(check, question, reads, result, ok = NA) {
  data.frame(
    check = check, question = question, reads = reads,
    result = as.character(result),
    verdict = if (isTRUE(ok)) "pass" else if (isFALSE(ok)) "FAIL" else "review",
    stringsAsFactors = FALSE
  )
}

.scorecard_subjects <- function(data, roles) {
  length(unique(as.character(data[[roles$id]])))
}

.scorecard_endpoints <- function(data, roles) {
  sort(unique(.endpoint(data, roles)[.observation_rows(data, roles)]))
}

.scorecard_rows_per_patient <- function(data, roles, which) {
  rows <- if (which == "dose") .dose_rows(data, roles) else
    .observation_rows(data, roles)
  round(sum(rows) / .scorecard_subjects(data, roles), 1)
}

# Each subject's sorted observation times, or values, as one string. Two
# subjects share a string only if that whole vector is identical. Sorting is the
# stricter reading: a reordering of the same numbers still counts as a copy,
# since the reordering is not what protects the patient.
.scorecard_vectors <- function(data, roles, column) {
  observed <- .observation_rows(data, roles)
  tapply(suppressWarnings(as.numeric(data[[column]][observed])),
         as.character(data[[roles$id]][observed]),
         function(values) paste(sort(values), collapse = ","))
}

# What counts as a copy is a vector too *few real patients share*, which is the
# same rule `min_pattern_share` applies to visit sets. Reproducing a vector that
# a whole arm shares identifies nobody: on a fixed protocol grid every patient
# is sampled at the same times, and a check that simply asked "does any
# synthetic vector equal some real one" would fail every well-designed study
# while catching nothing.
.scorecard_copies <- function(source, synthetic, roles, column, floor = 2L) {
  source_vectors <- .scorecard_vectors(source, roles, column)
  holders <- table(source_vectors)
  exposed <- names(holders)[holders < floor]
  length(intersect(exposed, .scorecard_vectors(synthetic, roles, column)))
}

# One baseline value per subject, as character so a factor cannot collapse to
# its integer codes.
.scorecard_holders <- function(data, roles, column) {
  as.character(tapply(as.character(data[[column]]),
                      as.character(data[[roles$id]]),
                      function(x) x[[1L]]))
}

# `strata` are protocol assignments and are categorical whatever their storage
# type, so a numerically coded arm is not excluded here.
.scorecard_categorical <- function(data, roles) {
  numeric_covariate <- vapply(roles$covariates,
                              function(v) is.numeric(data[[v]]), logical(1))
  c(roles$strata, roles$covariates[!numeric_covariate])
}

#' The scorecard for one synthetic dataset
#'
#' Runs the checks from `vignette("synthetic-data-checks")` that can be run
#' automatically, and returns them as one table with the answer and whether that
#' answer passes. The vignette is the reference for what each check asks and why
#' its pass criterion is what it is.
#'
#' The `reads` column decides where the table may go. Rows marked `"source"` or
#' `"both"` were computed from real patient data, so the filled-in scorecard is
#' itself restricted output and belongs in the environment the source lives in.
#' Only the `"synthetic"` and `"run report"` rows can travel with the data.
#'
#' A `"review"` verdict is not a soft `"pass"`. It marks a row where no
#' threshold would be honest, and it has to be read.
#'
#' Three checks in the vignette are absent here because no function can produce
#' them: C2 (dose and observation ordering), the source-side rare-level census
#' behind B5, and the one that matters most -- whether the pipeline that will
#' consume the real study runs unchanged against this output.
#'
#' @param source Source PMX data.
#' @param synthetic Generated synthetic PMX data from [synpmx_avatar()],
#'   carrying its `"pmx_settings"` attribute.
#' @param roles Explicit roles from [pmx_roles()].
#' @param proximity An already-computed [compare_pmx_proximity()] result, to
#'   save recomputing it. Left `NULL` it is computed here, which is the slowest
#'   part of the scorecard.
#'
#' @return A `pmx_scorecard` data frame with columns `check`, `question`,
#'   `reads`, `result` and `verdict`, marked `"restricted_not_releasable"`.
#' @seealso [compare_pmx()], [pmx_masking_report()],
#'   `vignette("synthetic-data-checks")`.
#' @export
#' @examples
#' data <- pmx_simulated_fixture(30)
#' roles <- pmx_roles(
#'   id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
#'   cmt = "CMT", dvid = "DVID", covariates = "WT"
#' )
#' synthetic <- suppressWarnings(synpmx_avatar(data, roles, seed = 1))
#' pmx_scorecard(data, synthetic, roles)
pmx_scorecard <- function(source, synthetic, roles, proximity = NULL) {
  settings <- attr(synthetic, "pmx_settings")
  if (is.null(settings)) {
    stop("`synthetic` carries no `pmx_settings` attribute; it did not come ",
         "from `synpmx_avatar()`.", call. = FALSE)
  }
  if (is.null(proximity)) {
    proximity <- compare_pmx_proximity(source, synthetic, roles)
  }

  source_endpoints <- .scorecard_endpoints(source, roles)
  synthetic_endpoints <- .scorecard_endpoints(synthetic, roles)
  source_subjects <- .scorecard_subjects(source, roles)
  synthetic_subjects <- .scorecard_subjects(synthetic, roles)
  floor <- as.integer(settings$min_pattern_share %||% 2L)
  time_copies <- .scorecard_copies(source, synthetic, roles, roles$time, floor)
  dv_copies <- .scorecard_copies(source, synthetic, roles, roles$dv, floor)
  inside_null <- proximity$adversarial_accuracy >= proximity$null_lower &&
    proximity$adversarial_accuracy <= proximity$null_upper
  flagged <- flag_identifiable_subjects(synthetic, roles)

  rows <- list(
    .scorecard_row(
      "A1", "Synthetic table is a legal PMX dataset", "synthetic",
      validate_pmx(synthetic, roles)$valid,
      isTRUE(validate_pmx(synthetic, roles)$valid)
    ),
    .scorecard_row(
      "A2", "Source is legal under the declared roles", "source",
      validate_pmx(source, roles, strict = FALSE)$valid,
      isTRUE(validate_pmx(source, roles, strict = FALSE)$valid)
    ),
    .scorecard_row(
      "A3", "Every endpoint survived", "both",
      paste(length(synthetic_endpoints), "of", length(source_endpoints)),
      setequal(source_endpoints, synthetic_endpoints)
    ),
    .scorecard_row(
      "A4", "Cohort size survived", "both",
      paste(source_subjects, "->", synthetic_subjects),
      source_subjects == synthetic_subjects
    ),
    .scorecard_row(
      "A5", "Observations per patient", "both",
      paste(.scorecard_rows_per_patient(source, roles, "obs"), "->",
            .scorecard_rows_per_patient(synthetic, roles, "obs"))
    ),
    .scorecard_row(
      "A5", "Doses per patient", "both",
      paste(.scorecard_rows_per_patient(source, roles, "dose"), "->",
            .scorecard_rows_per_patient(synthetic, roles, "dose"))
    ),
    .scorecard_row(
      "B1a", "Avatars wearing one real patient's visit set", "run report",
      settings$identifying_visit_sets, settings$identifying_visit_sets == 0L
    ),
    .scorecard_row(
      "B1b", "Avatars wearing one real patient's dose schedule", "run report",
      settings$identifying_dose_schedules,
      settings$identifying_dose_schedules == 0L
    ),
    .scorecard_row(
      "B2", "Synthetic patients unusual within their stratum", "synthetic",
      paste(sum(flagged$flagged), "of", nrow(flagged))
    ),
    .scorecard_row(
      "B3", "Adversarial accuracy inside its null interval", "both",
      sprintf("%.3f in [%.3f, %.3f]", proximity$adversarial_accuracy,
              proximity$null_lower, proximity$null_upper),
      inside_null
    ),
    .scorecard_row(
      "B4a", "Generated time vectors copying an exposed real one", "both",
      time_copies, time_copies == 0L
    ),
    .scorecard_row(
      "B4b", "Generated DV vectors copying an exposed real one", "both",
      dv_copies, dv_copies == 0L
    )
  )

  categorical <- .scorecard_categorical(synthetic, roles)
  if (length(categorical)) {
    rarest <- min(vapply(categorical, function(column) {
      min(as.integer(table(.scorecard_holders(synthetic, roles, column))))
    }, integer(1)))
    rows <- c(rows, list(.scorecard_row(
      "B5", "Synthetic patients holding the least-held level", "synthetic",
      rarest, rarest > 1L
    )))
  }

  if (length(roles$strata)) {
    arm_size <- function(data) {
      first <- !duplicated(as.character(data[[roles$id]]))
      table(as.character(data[[roles$strata[[1L]]]])[first])
    }
    source_arms <- arm_size(source)
    synthetic_arms <- arm_size(synthetic)[names(source_arms)]
    synthetic_arms[is.na(synthetic_arms)] <- 0L
    matched <- sum(source_arms == synthetic_arms)
    rows <- c(rows, list(.scorecard_row(
      "C3", "Strata keeping their source size", "both",
      paste(matched, "of", length(source_arms)),
      matched == length(source_arms)
    )))
  }

  rows <- c(rows, list(.scorecard_row(
    "C4", "Dose regimens represented", "both",
    paste(settings$dose_regimens_represented, "of",
          settings$dose_regimens_source)
  )))

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out <- structure(out, class = c("pmx_scorecard", "data.frame"))
  .mark_release(out, "restricted_not_releasable")
}

#' @export
print.pmx_scorecard <- function(x, ...) {
  plain <- as.data.frame(x)
  cat("Scorecard: see vignette(\"synthetic-data-checks\") for what each asks\n\n")
  for (i in seq_len(nrow(plain))) {
    cat(sprintf("  %-5s %-48s %-11s %-24s %s\n",
                plain$check[[i]], plain$question[[i]], plain$reads[[i]],
                plain$result[[i]], plain$verdict[[i]]))
  }
  failed <- sum(plain$verdict == "FAIL")
  cat("\n", if (failed) paste0(failed, " FAIL, ") else "no failures, ",
      sum(plain$verdict == "review"), " to review.\n",
      "Rows reading `source` or `both` are restricted output.\n", sep = "")
  invisible(x)
}

#' @exportS3Method knitr::knit_print
knit_print.pmx_scorecard <- function(x, ...) {
  # No caption: the object is routinely subset to a row or two in prose, and a
  # caption written for the whole card is wrong on every such excerpt.
  knitr::knit_print(knitr::kable(as.data.frame(x), row.names = FALSE))
}
