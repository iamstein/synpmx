# The scorecard ---------------------------------------------------------------
#
# `vignette("scorecard-synthetic-data-checks")` states every check and the
# answer that counts as passing for it. This runs the ones that can be run and
# collects the answers into one table, so that deciding whether to ship a
# dataset does not require rereading the document or rewriting sixty lines of
# reporting code per study.
#
# Three rows deliberately have no pass mark. Doses per patient is the clearest:
# on a study with individualised dosing it can halve while every guarantee above
# it still reads 0, and any threshold picked for it would be wrong on some
# dataset. `"review"` is the honest verdict there, not a soft `"pass"`.
#
# Every row names the call that explains it, because a number without the tool
# that produced it is a dead end: "4 of 9" is not something a reader can act on
# until they can see which four. The calls are written against this function's
# own argument names -- `source`, `synthetic`, `roles` -- so they are
# copy-pasteable in a session that uses those names, and readable as a pointer
# in one that does not.

.scorecard_row <- function(check, question, reads, result, explore, ok = NA) {
  data.frame(
    check = check, question = question, reads = reads,
    result = as.character(result),
    verdict = if (isTRUE(ok)) "pass" else if (isFALSE(ok)) "FAIL" else "review",
    explore = explore,
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
#' Runs the checks from `vignette("scorecard-synthetic-data-checks")` that can
#' be run automatically, and returns them as one table with the answer, whether
#' that answer passes, and the call that explains it. The vignette is the
#' reference for what each check asks and why its pass criterion is what it is.
#'
#' The `reads` column decides where the table may go. Rows marked `"source"` or
#' `"both"` were computed from real patient data, so the filled-in scorecard is
#' itself restricted output and belongs in the environment the source lives in.
#' Only the `"synthetic"` and `"run settings"` rows can travel with the data.
#' `"run settings"` means the value is the generation run's own record of what
#' it did -- `attr(synthetic, "pmx_settings")` -- rather than a measurement
#' taken from either table.
#'
#' The `explore` column names the call to run when a row needs explaining. The
#' calls are written against this function's argument names (`source`,
#' `synthetic`, `roles`), so rename them to whatever the session calls those
#' objects. Every row carries one; printing lists them under the table, for the
#' rows that did not pass, rather than as a sixth column.
#'
#' Printing and knitting differ on purpose. `print()` is a console layout: the
#' verdict table, then the calls to run, then the B5b levels. Knitting a chunk
#' that returns this object emits `knitr::kable()` tables instead, so a `.Rmd`
#' or `.qmd` gets the whole card including the `explore` column. Running a chunk
#' interactively in an IDE shows the console form, since nothing is knitting.
#'
#' A `"review"` verdict is not a soft `"pass"`. It marks a row where no
#' threshold would be honest, and it has to be read.
#'
#' `"FAIL"` is reserved for the rows where the answer is always a defect: the
#' output is not a legal dataset (A1), it is not the study that went in (A3,
#' A6), or it reproduces one real patient's structure verbatim (B1a, B1b, B4a,
#' B4b). Everything else is `"review"`, including every row whose answer can
#' move for a legitimate reason -- a subject dropped for want of donors, a
#' cohort statistic at a small sample size, a source a validator objects to.
#'
#' Two checks in the vignette are absent here because no function can produce
#' them: C2 (dose and observation ordering), and the one that matters most --
#' whether the pipeline that will consume the real study runs unchanged against
#' this output.
#'
#' @param source Source PMX data.
#' @param synthetic Generated synthetic PMX data from [synpmx_avatar()],
#'   carrying its `"pmx_settings"` attribute.
#' @param roles Explicit roles from [pmx_roles()].
#' @param proximity An already-computed [compare_pmx_proximity()] result, to
#'   save recomputing it. Left `NULL` it is computed here, which is the slowest
#'   part of the scorecard.
#'
#' @return A `synpmx_scorecard` data frame with columns `check`, `question`,
#'   `reads`, `result`, `verdict` and `explore`, marked
#'   `"restricted_not_releasable"`.
#' @seealso [compare_pmx()], [pmx_masking_report()], [pmx_endpoint_types()],
#'   `vignette("scorecard-synthetic-data-checks")`.
#' @export
#' @examples
#' data <- pmx_simulated_fixture(30)
#' roles <- pmx_roles(
#'   id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
#'   cmt = "CMT", dvid = "DVID", covariates = "WT"
#' )
#' synthetic <- suppressWarnings(synpmx_avatar(data, roles, seed = 1))
#' synpmx_scorecard(data, synthetic, roles)
synpmx_scorecard <- function(source, synthetic, roles, proximity = NULL) {
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
  flagged <- flag_identifiable_subjects(synthetic, roles)

  rows <- list(
    .scorecard_row(
      "A1", "Synthetic table is a legal PMX dataset", "synthetic",
      validate_pmx(synthetic, roles)$valid,
      "validate_pmx(synthetic, roles)",
      isTRUE(validate_pmx(synthetic, roles)$valid)
    ),
    # Review rather than FAIL: this row is about the source and the roles
    # declared for it, not about the output. A real study that a validator
    # objects to is a normal thing to have -- the answer is often "yes, and here
    # is why" -- and the objection has to be read rather than scored.
    .scorecard_row(
      "A2", "Source is legal under the declared roles", "source",
      validate_pmx(source, roles, strict = FALSE)$valid,
      "validate_pmx(source, roles, strict = FALSE)",
      if (isTRUE(validate_pmx(source, roles, strict = FALSE)$valid)) TRUE
      else NA
    ),
    .scorecard_row(
      "A3", "Every endpoint survived", "both",
      paste(length(synthetic_endpoints), "of", length(source_endpoints)),
      "compare_pmx_distributions(source, synthetic, roles)",
      setequal(source_endpoints, synthetic_endpoints)
    ),
    # Equal is a pass; anything else is review, for the same reason C3 is.
    # `on_donor_shortfall = "drop"` removes a subject that could not be built
    # from enough donors, which is the correct answer and not a defect, and it
    # is the small cohorts that lose one -- a stratum holding a single patient
    # can legitimately come back empty. Scoring that a FAIL would mark the
    # generator wrong for declining to build on a patient it could not protect.
    .scorecard_row(
      "A4", "Cohort size survived", "both",
      paste(source_subjects, "->", synthetic_subjects),
      "pmx_masking_report(synthetic, source, roles)",
      if (source_subjects == synthetic_subjects) TRUE else NA
    ),
    .scorecard_row(
      "A5", "Observations per patient", "both",
      paste(.scorecard_rows_per_patient(source, roles, "obs"), "->",
            .scorecard_rows_per_patient(synthetic, roles, "obs")),
      "compare_pmx_distributions(source, synthetic, roles)"
    ),
    .scorecard_row(
      "A5", "Doses per patient", "both",
      paste(.scorecard_rows_per_patient(source, roles, "dose"), "->",
            .scorecard_rows_per_patient(synthetic, roles, "dose")),
      "pmx_masking_report(synthetic, source, roles)"
    ),
    .scorecard_row(
      "B1a", "Avatars wearing one real patient's visit set", "run settings",
      settings$identifying_visit_sets,
      "plot_pmx_schedule(source, roles)",
      settings$identifying_visit_sets == 0L
    ),
    .scorecard_row(
      "B1b", "Avatars wearing one real patient's dose schedule", "run settings",
      settings$identifying_dose_schedules,
      "skeleton_uniqueness(source, roles, coarsen_time = TRUE)",
      settings$identifying_dose_schedules == 0L
    ),
    .scorecard_row(
      "B2", "Synthetic patients unusual within their stratum", "synthetic",
      paste(sum(flagged$flagged), "of", nrow(flagged)),
      "flag_identifiable_subjects(synthetic, roles)"
    ),
    # Always review. Outside the interval means something, but not one thing:
    # below it is memorisation, above it is the two sets having separated, which
    # costs utility and discloses nothing. Either can also be a small-sample
    # artefact at pharmacometric cohort sizes, where the null interval is wide
    # and the statistic moves with the seed. The number and its interval are
    # printed so the direction can be read; no verdict is put on them.
    .scorecard_row(
      "B3", "Adversarial accuracy inside its null interval", "both",
      sprintf("%.3f in [%.3f, %.3f]", proximity$adversarial_accuracy,
              proximity$null_lower, proximity$null_upper),
      "compare_pmx_proximity(source, synthetic, roles)"
    ),
    .scorecard_row(
      "B4a", "Generated time vectors copying an exposed real one", "both",
      time_copies,
      "skeleton_uniqueness(source, roles, coarsen_time = TRUE)",
      time_copies == 0L
    ),
    # No function finds *which* vector was copied, so the pointer is the one
    # that measures the thing a copy is the extreme case of: whether generated
    # values sit closer to real ones than real ones sit to each other.
    .scorecard_row(
      "B4b", "Generated DV vectors copying an exposed real one", "both",
      dv_copies,
      "compare_pmx_proximity(source, synthetic, roles)",
      dv_copies == 0L
    )
  )

  # Only where there is something to ask. A study with no discrete endpoint
  # cannot pass or fail this, and a row reading "0 of 0" on every continuous
  # study is a row nobody reads on the one study that needs it.
  value_types <- .endpoint_value_types(source, roles)
  discrete <- names(value_types)[
    vapply(value_types, function(s) !identical(s$type, "continuous"),
           logical(1))
  ]
  if (length(discrete)) {
    violations <- .endpoint_type_violations(synthetic, roles, value_types)
    clean <- sum(violations[discrete] == 0L)
    rows <- c(rows, list(.scorecard_row(
      "A6", "Discrete endpoints keeping their source scale", "both",
      paste(clean, "of", length(discrete)),
      "pmx_endpoint_types(source, roles)",
      clean == length(discrete)
    )))
  }

  # NULL when the roles declare no categorical axis, which is also when there is
  # no B5 row at all.
  rare_levels <- NULL
  categorical <- .scorecard_categorical(synthetic, roles)
  if (length(categorical)) {
    # Which column and which level, not just the count. "1" on its own says a
    # patient is alone in holding something without saying what, which is not
    # something a reader can act on -- and the column it names is the one to
    # look at, since a level held by one patient is what singles that avatar
    # out.
    rarest <- do.call(rbind, lapply(categorical, function(column) {
      holders <- table(.scorecard_holders(synthetic, roles, column))
      data.frame(column = column,
                 level = names(holders)[[which.min(holders)]],
                 n = min(as.integer(holders)), stringsAsFactors = FALSE)
    }))
    rarest <- rarest[which.min(rarest$n), , drop = FALSE]
    rows <- c(rows, list(.scorecard_row(
      "B5a", "Patients holding the least-held categorical level", "synthetic",
      sprintf("%s = %s: %d", rarest$column, rarest$level, rarest$n),
      sprintf("table(synthetic$%s)", rarest$column),
      # Review rather than FAIL, because this is the weak form of the check and
      # it is wrong in both directions: a level held by thirty source patients
      # can land on one avatar and mean nothing, and a level held by two source
      # patients can be copied onto ten avatars and pass while the disclosure is
      # real. What the risk is made of is the source-side count, which B5b below
      # is. This flags a shape worth looking at.
      if (rarest$n > 1L) TRUE else NA
    )))

    # B5b is the strong form: not "is anyone alone in the output" but "did a
    # level too few REAL patients held get copied out at all". It reads the
    # source, so it is `both` rather than `synthetic` -- on a study where B5a
    # was the only categorical row, adding it makes the card restricted output.
    # Review rather than FAIL for now: on a small oncology cohort `RACE` will
    # light this up constantly, and the answer there is usually "do not carry
    # `RACE`" rather than "the generator is broken". It earns a verdict once it
    # has been read on real studies.
    census <- compare_pmx_rare_levels(source, synthetic, roles, floor = floor)
    exposed <- sum(census$exposed)
    leaked <- census[census$exposed & census$reached, , drop = FALSE]
    rare_levels <- as.data.frame(leaked)[, c("column", "level",
                                             "source_patients",
                                             "synthetic_patients")]
    rows <- c(rows, list(.scorecard_row(
      "B5b", "Rare source levels copied into the output", "both",
      # The count only. Which levels they were is a list, and it is printed
      # under the table and knitted as its own table, because a study that
      # trips this trips it several times over and one name in a cell is not
      # the answer.
      sprintf("%d of %d exposed", nrow(leaked), exposed),
      "compare_pmx_rare_levels(source, synthetic, roles)",
      if (!nrow(leaked)) TRUE else NA
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
    # Not a FAIL when the sizes move. `on_donor_shortfall = "drop"` removes a
    # subject that could not be built from enough donors, which is the correct
    # answer and lands on whichever arm was thinnest; an arm changing size is
    # therefore a reading about this cohort, not a defect. Only every arm
    # matching is a pass, so the row still has to be looked at.
    rows <- c(rows, list(.scorecard_row(
      "C3", "Strata keeping their source size", "both",
      paste(matched, "of", length(source_arms)),
      sprintf("table(synthetic$%s[!duplicated(synthetic$%s)])",
              roles$strata[[1L]], roles$id),
      if (matched == length(source_arms)) TRUE else NA
    )))
  }

  # C4 in the vignette is the open question "what was lost?", which has no pass
  # mark. This row measures one countable part of it, and that part does: every
  # distinct source regimen still appearing in the output is nothing lost, so it
  # is a pass. Fewer is `review` rather than `FAIL` -- declining to build on a
  # patient whose regimen nobody shares is the correct answer, and the row is
  # there to say what it cost.
  #
  # A "regimen" here is the set of dose TIMES a patient received, and nothing
  # else: `.dose_schedule_keys()` keys on times alone. Amounts are deliberately
  # not in it. On a weight-based study every patient's milligrams are their own,
  # so a key including amounts would report one regimen per patient on every run
  # and measure nothing. The row is named for what it counts.
  rows <- c(rows, list(.scorecard_row(
    "C4", "Distinct dose-time schedules represented", "both",
    paste(settings$dose_regimens_represented, "of",
          settings$dose_regimens_source),
    "pmx_masking_report(synthetic, source, roles)",
    if (settings$dose_regimens_represented ==
        settings$dose_regimens_source) TRUE else NA
  )))

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out <- structure(out, class = c("synpmx_scorecard", "data.frame"))
  # Carried beside the table rather than in it: B5b is one row whose answer is a
  # list, and printing that list is what makes the row actionable.
  attr(out, "rare_levels") <- rare_levels
  .mark_release(out, "restricted_not_releasable")
}

#' @export
print.synpmx_scorecard <- function(x, ...) {
  plain <- as.data.frame(x)
  # The table is the verdicts, and nothing else. `explore` was a sixth column of
  # long calls that pushed every line past 150 characters and wrapped in an
  # ordinary console, which made the whole card unreadable to save a second
  # glance. It goes below the table instead, one line per row that needs it, and
  # only those rows. The returned object still carries the column on every row.
  #
  # Widths come from the contents rather than from constants. A `result` holds a
  # study's own label -- "RACE = NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER" --
  # and clipping that to fit a fixed column turned the one row that needed
  # reading into the one row nobody could read. The table is as wide as the
  # study makes it, which on most studies is narrower than the old constants.
  columns <- c("check", "question", "reads", "result", "verdict")
  width <- vapply(columns, function(column) {
    max(nchar(c(column, plain[[column]])))
  }, integer(1))
  # `result` is the only column that can hold a study's own label, so it is the
  # only one capped. Padding is capped, but nothing is cut: a longer value runs
  # past its column and pushes that row's verdict right, which costs one ragged
  # line instead of a whole table sized for its worst label. The other columns
  # are package-authored and bounded, and capping those just made every long
  # question ragged for nothing.
  width[["result"]] <- min(width[["result"]], 45L)
  # Every column padded but the last, which nothing follows.
  format <- paste0("  ", paste(sprintf("%%-%ds", width[-length(width)]),
                               collapse = " "), " %s\n")
  line <- function(values) cat(do.call(sprintf, c(list(format), as.list(values))))
  cat("Scorecard: see vignette(\"scorecard-synthetic-data-checks\")",
      "for what each asks\n\n")
  line(columns)
  for (i in seq_len(nrow(plain))) {
    line(unlist(plain[i, columns], use.names = FALSE))
  }

  unresolved <- plain$verdict != "pass"
  if (any(unresolved)) {
    cat("\nTo explore, with `source`, `synthetic` and `roles` named as you",
        "have them:\n")
    for (i in which(unresolved)) {
      cat(sprintf("  %-5s %s\n", plain$check[[i]], plain$explore[[i]]))
    }
  }

  # B5b can only name one level in a table cell, and a study that trips it
  # usually trips it several times over. The whole list is carried on the card
  # so that reading it does not require running the census again.
  rare <- attr(x, "rare_levels")
  if (!is.null(rare) && nrow(rare)) {
    cat("\nRare source levels that reached the output (B5b):\n")
    for (i in seq_len(nrow(rare))) {
      cat(sprintf("  %s = %s -- %d source patient(s), %d avatar(s)\n",
                  rare$column[[i]], rare$level[[i]],
                  rare$source_patients[[i]], rare$synthetic_patients[[i]]))
    }
  }

  failed <- sum(plain$verdict == "FAIL")
  cat("\n", if (failed) paste0(failed, " FAIL, ") else "no failures, ",
      sum(plain$verdict == "review"), " to review.\n",
      "`run settings` rows come from the run's own record, ",
      "`attr(synthetic, \"pmx_settings\")`.\n",
      "Rows reading `source` or `both` are restricted output.\n", sep = "")
  invisible(x)
}

#' @exportS3Method knitr::knit_print
knit_print.synpmx_scorecard <- function(x, ...) {
  # No caption on the card itself: the object is routinely subset to a row or
  # two in prose, and a caption written for the whole card is wrong on every
  # such excerpt. The B5b detail is a second table rather than a footnote, so
  # that a knitted report says which levels without the reader running the
  # census again.
  out <- paste(knitr::kable(as.data.frame(x), row.names = FALSE),
               collapse = "\n")
  rare <- attr(x, "rare_levels")
  if (!is.null(rare) && nrow(rare)) {
    out <- c(out, paste(knitr::kable(
      rare, row.names = FALSE,
      caption = paste("RESTRICTED -- B5b: rare source levels that reached",
                      "the output. Each is one real patient's attribute,",
                      "copied.")
    ), collapse = "\n"))
  }
  knitr::asis_output(paste(out, collapse = "\n\n"))
}
