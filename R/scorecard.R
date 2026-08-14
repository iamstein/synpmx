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

.scorecard_row <- function(check, question, reads, result, explore, ok = NA,
                           verdict = NULL) {
  data.frame(
    check = check, question = question, reads = reads,
    result = as.character(result),
    verdict = verdict %||%
      if (isTRUE(ok)) "pass" else if (isFALSE(ok)) "FAIL" else "review",
    explore = explore,
    stringsAsFactors = FALSE
  )
}

# A row only the generating run can answer.
#
# Three of them -- B1a, B1b and C2 -- read `attr(synthetic, "pmx_settings")`
# rather than either table, and that is not a shortcut. Generated times are the
# coarsened grid plus resampled deviations, applied to dose rows as well, so an
# avatar's schedule cannot be matched back to a source patient's by exact key:
# snapping the finished table back onto the grid was tried and reported an
# avatar on `wbcSim` as holding a schedule it was never given. The run measured
# it before the deviations were applied, and that measurement is the better one.
#
# So a table with no such record -- another method's output, or this one's read
# back from a CSV -- gets the row marked `unavailable` rather than dropped. A
# card holding sixteen rows on one table and nineteen on another cannot be
# compared, and "not measured" is a different statement from "passed". `result`
# and `ok` are evaluated only when the record exists.
.scorecard_recorded <- function(check, question, settings, result, explore,
                                ok) {
  if (is.null(settings)) {
    return(.scorecard_row(check, question, "run settings",
                          "no run record", explore,
                          verdict = "unavailable"))
  }
  .scorecard_row(check, question, "run settings", result, explore, ok)
}

# A B1 result: how many avatars leaked, and out of which arms.
#
# Two arms at most in the cell. A study with five leaking arms has one answer --
# `unmaskable_strata()`, which the row already points at -- and five arm labels
# in a table cell is not it.
.scorecard_leak <- function(count, strata) {
  if (!count || !length(strata)) return(as.character(count))
  strata <- sort(strata, decreasing = TRUE)
  shown <- utils::head(strata, 2L)
  paste0(count, " in ",
         paste(sprintf("%s (%d)", names(shown), shown), collapse = ", "),
         if (length(strata) > 2L) sprintf(", +%d more", length(strata) - 2L)
         else "")
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

# D1 out of `compare_pmx_distributions()`, which returns a table per variable
# rather than anything a card can hold. The quantity worth putting in one cell
# is the spread: an avatar is a weighted mean of several donors, so
# between-subject variability shrinks on every run, and how far it shrank is the
# reading a modeller makes before deciding whether the output is usable.
#
# The variable reported is the one whose standard deviation moved FURTHEST in
# either direction, not the one that shrank most. A spread that grew is as much
# a reason to look -- it usually means an endpoint's noise was reinflated, or a
# covariate drew donors from across a bimodal cohort.
.scorecard_spread <- function(distributions) {
  # Three columns only: the endpoint table carries an `n_subjects` the covariate
  # table does not, and binding them whole is what an earlier version did.
  columns <- c("variable", "dataset", "sd")
  tables <- list(distributions$endpoints, distributions$covariates_numeric)
  tables <- lapply(tables[!vapply(tables, is.null, logical(1))],
                   function(df) as.data.frame(df)[, columns, drop = FALSE])
  rows <- do.call(rbind, tables)
  if (is.null(rows) || !nrow(rows)) return(NULL)
  from <- rows[rows$dataset == "source", , drop = FALSE]
  to <- rows[rows$dataset == "synthetic", , drop = FALSE]
  ratio <- to$sd[match(from$variable, to$variable)] / from$sd
  names(ratio) <- from$variable
  # A variable observed once has no `sd`, and one whose source `sd` is zero has
  # no ratio. Both are dropped rather than reported as Inf or NaN.
  ratio <- ratio[is.finite(ratio) & ratio > 0]
  if (!length(ratio)) return(NULL)
  furthest <- which.max(abs(log(ratio)))
  list(variable = names(ratio)[[furthest]], ratio = unname(ratio[[furthest]]),
       n = length(ratio))
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
#' # Scoring a table this package did not generate
#'
#' Everything measured from the two tables is measurable on any synthetic
#' dataset, whatever produced it, so a table carrying no `"pmx_settings"`
#' attribute is scored rather than refused: the three `"run settings"` rows
#' (B1a, B1b, C2) come back with the verdict `"unavailable"` and the rest of
#' the card is computed as usual. That covers another method's output, and this
#' package's own output read back from a file.
#'
#' Those three rows cannot be recomputed from the finished table, and that is a
#' property of the generator rather than an omission here. Generated times are
#' the coarsened visit grid plus resampled deviations -- applied to dose rows
#' too -- so an avatar's schedule no longer matches any source patient's key
#' exactly, and matching it back by snapping to the grid reports schedules that
#' were never given. The run measures both guarantees before the deviations are
#' applied. [unmaskable_strata()] is the part of that question answerable
#' without any run record: it reads the source alone and names the arms whose
#' patients no method could mask.
#'
#' The `explore` column names the call to run when a row needs explaining. The
#' calls are written against this function's argument names (`source`,
#' `synthetic`, `roles`), so rename them to whatever the session calls those
#' objects. Every row carries one; printing lists them under the table, for the
#' rows that did not pass, rather than as a sixth column.
#'
#' Printing and knitting differ on purpose. `print()` is a console layout: the
#' verdict table, then the calls to run, then the B5 levels. Knitting a chunk
#' that returns this object emits `knitr::kable()` tables instead, so a `.Rmd`
#' or `.qmd` gets the whole card including the `explore` column. Running a chunk
#' interactively in an IDE shows the console form, since nothing is knitting.
#'
#' A `"review"` verdict is not a soft `"pass"`. It marks a row where no
#' threshold would be honest, and it has to be read. Nor is `"unavailable"`:
#' it marks a row nothing was measured for.
#'
#' # Every card holds every row
#'
#' The same checks come back whatever the study declares. Where a study gives a
#' check nothing to ask -- no discrete endpoint for A6, no `strata` for C1 and
#' C3, no categorical axis for B5 -- the `result` says so and the
#' verdict is `"pass"`, rather than the row going missing. Two cards can then be
#' compared
#' row for row, and an absent row cannot be mistaken for one that passed.
#'
#' # Plot the data as well
#'
#' D1 reports the standard deviation that moved furthest between the two
#' tables, and a standard deviation cannot see a shape: one bell and two humps
#' with the same mean and spread give the same cell. Plot `DV` against time and
#' each covariate's distribution, source and synthetic on the same axes, before
#' deciding the output is usable. No function is offered for it -- every group
#' has plotting code it already trusts for its own study, and a generic one
#' would be a worse version of that.
#'
#' `"FAIL"` is reserved for the rows where the answer is always a defect: the
#' output is not a legal dataset (A1), it is not the study that went in (A3,
#' A6), or it reproduces one real patient's structure verbatim (B1a, B1b, B4a,
#' B4b). No other row can `"FAIL"`: the rest answer `"pass"` when there is
#' nothing to read and `"review"` when there is something whose meaning depends
#' on the study -- a subject dropped for want of donors, a cohort statistic at a
#' small sample size, a source a validator objects to. Four rows are `"review"`
#' whatever they land on, because no threshold on them would be honest: A5a,
#' A5b, B3 and D1.
#'
#' The check that matters most is absent here because no function can produce
#' it: whether the pipeline that will consume the real study runs unchanged
#' against this output.
#'
#' @param source Source PMX data.
#' @param synthetic Generated synthetic PMX data. From [synpmx_avatar()] it
#'   carries a `"pmx_settings"` attribute and the whole card can be filled in;
#'   without one, the three rows that need it read `"unavailable"`.
#' @param roles Explicit roles from [pmx_roles()].
#' @param proximity An already-computed [compare_pmx_proximity()] result, to
#'   save recomputing it. Left `NULL` it is computed here, which is the slowest
#'   part of the scorecard.
#'
#' @return A `synpmx_scorecard` data frame with columns `check`, `question`,
#'   `reads`, `result`, `verdict` and `explore`, marked
#'   `"restricted_not_releasable"`.
#' @seealso [synpmx_scorecard_datatable()] to colour the verdicts in an HTML
#'   report, [compare_pmx()], [pmx_masking_report()], [pmx_endpoint_types()],
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
  # No `pmx_settings` is not an error. Everything the card measures from the two
  # tables is measurable on any synthetic dataset, whatever produced it, and
  # refusing the whole card for the three rows that need the run's own record
  # made the one function meant to score an output unable to score anything but
  # this package's own -- including this package's own, once it had been through
  # `write.csv()`. Those three rows say so; the rest are computed as usual.
  settings <- attr(synthetic, "pmx_settings")
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

  # Always present, including on a study with nothing discrete to check. The row
  # was conditional, and a card whose rows depend on the study cannot be
  # compared across studies -- worse, an absent row reads as a row that passed,
  # when it means the question was never asked. A continuous study says so in
  # the result and passes, which is the true answer: no endpoint left its
  # source scale, because none had one to leave. Built here rather than
  # appended, because a row appended after the list prints after B4b, and A6
  # belongs with the other A rows.
  value_types <- .endpoint_value_types(source, roles)
  discrete <- names(value_types)[
    vapply(value_types, function(s) !identical(s$type, "continuous"),
           logical(1))
  ]
  a6_rows <- if (!length(discrete)) {
    list(.scorecard_row(
      "A6", "Discrete endpoints keeping their source scale", "both",
      "no discrete endpoint", "pmx_endpoint_types(source, roles)", TRUE
    ))
  } else {
    violations <- .endpoint_type_violations(synthetic, roles, value_types)
    clean <- sum(violations[discrete] == 0L)
    list(.scorecard_row(
      "A6", "Discrete endpoints keeping their source scale", "both",
      paste(clean, "of", length(discrete)),
      "pmx_endpoint_types(source, roles)",
      clean == length(discrete)
    ))
  }

  rows <- c(list(
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
    # Equal is a pass; anything else is review, for the same reason C1 is.
    # `on_donor_shortfall = "drop"` removes a subject that could not be built
    # from enough donors, which is the correct answer and not a defect, and it
    # is the small cohorts that lose one -- a stratum holding a single patient
    # can legitimately come back empty. Scoring that a FAIL would mark the
    # generator wrong for declining to build on a patient it could not protect.
    .scorecard_row(
      "A4", "Cohort size survived", "both",
      paste(source_subjects, "->", synthetic_subjects),
      'pmx_masking_report(synthetic, source, roles, section = "anchors")',
      if (source_subjects == synthetic_subjects) TRUE else NA
    ),
    .scorecard_row(
      "A5a", "Observations per patient", "both",
      paste(.scorecard_rows_per_patient(source, roles, "obs"), "->",
            .scorecard_rows_per_patient(synthetic, roles, "obs")),
      "compare_pmx_distributions(source, synthetic, roles)"
    ),
    .scorecard_row(
      "A5b", "Doses per patient", "both",
      paste(.scorecard_rows_per_patient(source, roles, "dose"), "->",
            .scorecard_rows_per_patient(synthetic, roles, "dose")),
      'pmx_masking_report(synthetic, source, roles, section = "dose_schedules")'
    )
  ), a6_rows, list(
    # Both rows name the arm rather than only the count. An avatar is built in
    # the arm it was allocated to and is never moved out of it -- that would
    # rewrite the arm sizes -- so one of these failing means one arm holds
    # nobody who can be masked, and which arm is the whole of what to do next.
    # `unmaskable_strata()` is the pointer for the same reason: it answers the
    # question from the source, before a run is made.
    .scorecard_recorded(
      "B1a", "Avatars with a visit set nobody else shares", settings,
      .scorecard_leak(settings$identifying_visit_sets,
                      settings$identifying_visit_set_strata),
      # Runnable either way, and the only one of the three pointers that is:
      # it reads the source alone, so it says which arms a table of unknown
      # provenance could be carrying a real schedule out of.
      "unmaskable_strata(source, roles)",
      settings$identifying_visit_sets == 0L
    ),
    .scorecard_recorded(
      "B1b", "Avatars with a dose schedule nobody else shares", settings,
      .scorecard_leak(settings$identifying_dose_schedules,
                      settings$identifying_dose_schedule_strata),
      "unmaskable_strata(source, roles)",
      settings$identifying_dose_schedules == 0L
    ),
    # 0 is a pass rather than a review: the row is a list of records to read,
    # and an empty list has nothing in it to read. Above 0 is `review` and never
    # `FAIL` -- the right count is not necessarily 0, since a study can contain
    # a patient who is genuinely unusual and whose avatar therefore is too.
    .scorecard_row(
      "B2", "Synthetic patients unusual within their stratum", "synthetic",
      paste(sum(flagged$flagged), "of", nrow(flagged)),
      "flag_identifiable_subjects(synthetic, roles)",
      verdict = if (sum(flagged$flagged) == 0L) "pass" else "review"
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
  ))

  # NULL when the roles declare no categorical axis. The B5 row is still
  # emitted in that case, saying so.
  rare_levels <- NULL
  categorical <- .scorecard_categorical(synthetic, roles)
  if (!length(categorical)) {
    rows <- c(rows, list(
      .scorecard_row(
        "B5", "Rare source levels copied into the output", "both",
        "no categorical covariate or stratum",
        "pmx_roles(strata = , covariates = )", TRUE
      )
    ))
  } else {
    # One B5 row, and it is the source-side one. A synthetic-side count -- "is
    # any level held by exactly one avatar" -- was a second row here and is
    # gone: it is wrong in both directions, since a level held by thirty source
    # patients can land on one avatar and mean nothing while a level held by two
    # can be copied onto ten avatars and pass while the disclosure is real. Its
    # only merit was being answerable without the source, and this function is
    # handed the source. `table(synthetic$COLUMN)` is that question for whoever
    # holds the released table and nothing else.
    #
    # Not "is anyone alone in the output" but "did a level too few REAL patients
    # held get copied out at all". Review rather than FAIL for now: on a small
    # oncology cohort `RACE` will light this up constantly, and the answer there
    # is usually "do not carry `RACE`" rather than "the generator is broken". It
    # earns a verdict once it has been read on real studies.
    census <- compare_pmx_rare_levels(source, synthetic, roles, floor = floor)
    exposed <- sum(census$exposed)
    leaked <- census[census$exposed & census$reached, , drop = FALSE]
    rare_levels <- as.data.frame(leaked)[, c("column", "level",
                                             "source_patients",
                                             "synthetic_patients")]
    rows <- c(rows, list(.scorecard_row(
      "B5", "Rare source levels copied into the output", "both",
      # The count only. Which levels they were is a list, and it is printed
      # under the table and knitted as its own table, because a study that
      # trips this trips it several times over and one name in a cell is not
      # the answer.
      sprintf("%d of %d exposed", nrow(leaked), exposed),
      "compare_pmx_rare_levels(source, synthetic, roles)",
      if (!nrow(leaked)) TRUE else NA
    )))
  }

  if (!length(roles$strata)) {
    # A study with no declared arm has no balance to preserve, so there is
    # nothing here to get wrong. The row says which of those two it is.
    rows <- c(rows, list(.scorecard_row(
      "C1", "Strata keeping their source size", "both",
      "no strata declared", "pmx_roles(strata = )", TRUE
    )))
  } else {
    arm_size <- function(data) {
      first <- !duplicated(as.character(data[[roles$id]]))
      table(as.character(data[[roles$strata[[1L]]]])[first])
    }
    source_arms <- arm_size(source)
    # Lined up by position: a blank arm label cannot be looked up by name, and
    # `[""]` gives NA, which then reads as an arm that lost every patient.
    at <- match(names(source_arms), names(arm_size(synthetic)))
    synthetic_arms <- as.integer(arm_size(synthetic))[at]
    synthetic_arms[is.na(synthetic_arms)] <- 0L
    matched <- sum(as.integer(source_arms) == synthetic_arms)
    # Not a FAIL when the sizes move. `on_donor_shortfall = "drop"` removes a
    # subject that could not be built from enough donors, which is the correct
    # answer and lands on whichever arm was thinnest; an arm changing size is
    # therefore a reading about this cohort, not a defect. Only every arm
    # matching is a pass, so the row still has to be looked at.
    rows <- c(rows, list(.scorecard_row(
      "C1", "Strata keeping their source size", "both",
      paste(matched, "of", length(source_arms)),
      "compare_pmx_strata_sizes(source, synthetic, roles)",
      if (matched == length(source_arms)) TRUE else NA
    )))
  }

  # C2 in the vignette is the open question "what was lost?", which has no pass
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
  #
  # `reads` is `run settings`, not `both`. Both counts come from the run's
  # record -- which source regimens it anchored on -- and the row said `both`
  # while reading neither table, which is a claim about where the filled-in card
  # may travel, not a cosmetic label.
  rows <- c(rows, list(.scorecard_recorded(
    "C2", "Distinct dose-time schedules represented", settings,
    paste(settings$dose_regimens_represented, "of",
          settings$dose_regimens_source),
    'pmx_masking_report(synthetic, source, roles, section = "dose_schedules")',
    if (settings$dose_regimens_represented ==
        settings$dose_regimens_source) TRUE else NA
  )))

  # C3 is the other half of the arm-level question C1 asks, and it is separate
  # because A3 does not cover it: A3 compares endpoint sets across the whole
  # cohort, which a placebo arm passes while holding no PK concentration of its
  # own. An avatar never leaves the arm it was anchored in, so an arm coming
  # back without an endpoint its source arm held is a real loss and nothing else
  # on the card would see it.
  #
  # `review` rather than `FAIL`, on the same reasoning as A4 and C1: an endpoint
  # held by one patient in an arm leaves with that patient when
  # `on_donor_shortfall = "drop"` removes them, which is the correct answer.
  # Only every arm matching is a pass.
  if (!length(roles$strata)) {
    rows <- c(rows, list(.scorecard_row(
      "C3", "Arms keeping their source endpoints", "both",
      "no strata declared", "pmx_roles(strata = )", TRUE
    )))
  } else {
    arm_endpoints <- function(data) {
      lapply(.strata_endpoint_holders(data, roles, roles$strata[[1L]]),
             function(endpoints) sort(unique(endpoints)))
    }
    source_endpoints_by_arm <- arm_endpoints(source)
    synthetic_endpoints_by_arm <- arm_endpoints(synthetic)
    # By name rather than by position: an arm that lost every patient has no
    # entry on the synthetic side at all, and `[[` on a missing name errors.
    matched_arms <- vapply(names(source_endpoints_by_arm), function(level) {
      at <- match(level, names(synthetic_endpoints_by_arm))
      setequal(source_endpoints_by_arm[[level]],
               if (is.na(at)) character() else synthetic_endpoints_by_arm[[at]])
    }, logical(1))
    rows <- c(rows, list(.scorecard_row(
      "C3", "Arms keeping their source endpoints", "both",
      paste(sum(matched_arms), "of", length(matched_arms)),
      "compare_pmx_strata_endpoints(source, synthetic, roles)",
      if (all(matched_arms)) TRUE else NA
    )))
  }

  # D1, and it is `review` on every study there will ever be. Spread shrinking
  # is what the algorithm does rather than something it got wrong, and spread
  # growing has no threshold either, so a pass mark here would be an invention.
  # The row exists to put the number in front of the reader and send them to the
  # table it came from.
  # `output = "tables"` is load-bearing rather than tidy: the default returns a
  # figure, a ggplot is a list, and `$endpoints` on one is NULL rather than an
  # error -- so the bare call would leave D1 reporting "no numeric variable to
  # compare" on every study, quietly.
  spread <- .scorecard_spread(
    compare_pmx_distributions(source, synthetic, roles, output = "tables")
  )
  rows <- c(rows, list(.scorecard_row(
    "D1", "Values landing in the same range", "both",
    if (is.null(spread)) "no numeric variable to compare" else
      sprintf("sd x%s on %s (furthest of %d)",
              format(signif(spread$ratio, 2)), spread$variable, spread$n),
    # Pointed at the tables rather than the default figure: this row quotes a
    # standard deviation, and a density curve does not carry one.
    'compare_pmx_distributions(source, synthetic, roles, output = "tables")',
    verdict = "review"
  )))

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out <- structure(out, class = c("synpmx_scorecard", "data.frame"))
  # Carried beside the table rather than in it: B5 is one row whose answer is a
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

  # B5 can only name one level in a table cell, and a study that trips it
  # usually trips it several times over. The whole list is carried on the card
  # so that reading it does not require running the census again.
  rare <- attr(x, "rare_levels")
  if (!is.null(rare) && nrow(rare)) {
    cat("\nRare source levels that reached the output (B5):\n")
    for (i in seq_len(nrow(rare))) {
      cat(sprintf("  %s = %s -- %d source patient(s), %d avatar(s)\n",
                  rare$column[[i]], .level_label(rare$level[[i]]),
                  rare$source_patients[[i]], rare$synthetic_patients[[i]]))
    }
  }

  failed <- sum(plain$verdict == "FAIL")
  unavailable <- sum(plain$verdict == "unavailable")
  cat("\n", if (failed) paste0(failed, " FAIL, ") else "no failures, ",
      sum(plain$verdict == "review"), " to review",
      # Said in the count line, not only in the rows: three rows quietly not
      # answered is the kind of thing a reader takes for three rows passed.
      if (unavailable) paste0(", ", unavailable, " unanswered") else "", ".\n",
      "`run settings` rows come from the run's own record, ",
      "`attr(synthetic, \"pmx_settings\")`",
      if (unavailable) ", which this table does not carry" else "", ".\n",
      "Rows reading `source` or `both` are restricted output.\n", sep = "")
  invisible(x)
}

#' @exportS3Method knitr::knit_print
knit_print.synpmx_scorecard <- function(x, ...) {
  # No caption on the card itself: the object is routinely subset to a row or
  # two in prose, and a caption written for the whole card is wrong on every
  # such excerpt. The B5 detail is a second table rather than a footnote, so
  # that a knitted report says which levels without the reader running the
  # census again.
  out <- paste(knitr::kable(as.data.frame(x), row.names = FALSE),
               collapse = "\n")
  rare <- attr(x, "rare_levels")
  if (!is.null(rare) && nrow(rare)) {
    out <- c(out, paste(knitr::kable(
      rare, row.names = FALSE,
      caption = paste("B5: rare source levels that reached",
                      "the output. Each is one real patient's attribute,",
                      "copied.")
    ), collapse = "\n"))
  }
  knitr::asis_output(paste(out, collapse = "\n\n"))
}

# One place for the palette, so the function and its documentation cannot
# disagree about which colour means what. The two dark shades are chosen to
# stay legible as text on white rather than to be the brightest orange and red
# available: plain `red` and `orange` are respectively glaring and washed out
# at this size, and the card is read, not glanced at. Grey says "nothing was
# measured here", which is a quieter statement than either.
.scorecard_verdict_colours <- c(FAIL = "#B00020", review = "#B45309",
                                unavailable = "#6C757D")

# The tints behind those two. Pale enough that the bold text on top stays the
# thing being read -- both clear 7:1 against their own foreground -- and pale
# enough that five tinted cells on a thirty-row card look like marks rather
# than like a warning banner. `"unavailable"` gets none: it is deliberately the
# quiet verdict, and a tint would put it back alongside the loud two.
.scorecard_verdict_fills <- c(FAIL = "#FDECEA", review = "#FFF4E5")

#' A scorecard as a coloured HTML table
#'
#' Displays a [synpmx_scorecard()] as an interactive table with the verdicts
#' coloured: `"FAIL"` in bold red on a light red, `"review"` in bold orange on
#' a light orange, `"unavailable"` in muted grey, and `"pass"` left as ordinary
#' text. A card is five verdicts among thirty-odd rows of prose, and the rows
#' that need reading are the ones that have to be findable without reading all
#' of it.
#'
#' `synpmx_scorecard()` computes the card and this displays it, so the object
#' is unchanged and can be subset, saved or printed as usual. The colouring is
#' the only thing added.
#'
#' What is emitted is what knitting the card itself emits, with the colouring
#' added: the card, then the B5 rare-level detail where a study has any.
#'
#' `DT` is a suggested package rather than a required one -- `synpmx` has no
#' hard dependencies -- so without it installed this says so and prints the
#' card in the console form instead. The verdicts are all still there; only the
#' colour is missing.
#'
#' The same restriction applies as to the card itself. Rows reading `"source"`
#' or `"both"` were computed from real patient data, so an HTML file holding
#' the whole table belongs in the environment the source lives in.
#'
#' @param x A [synpmx_scorecard()].
#' @param ... Passed to `DT::datatable()`. Paging is off and row numbers are
#'   suppressed by default, since the whole card is meant to be read at once
#'   and `check` already names each row.
#'
#' @return An `htmltools::tagList` holding the coloured card and the notes that
#'   knitting one carries. Without `DT` installed, `x` invisibly, having
#'   printed it.
#' @seealso [synpmx_scorecard()].
#' @export
#' @examples
#' data <- pmx_simulated_fixture(30)
#' roles <- pmx_roles(
#'   id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
#'   cmt = "CMT", dvid = "DVID", covariates = "WT"
#' )
#' synthetic <- suppressWarnings(synpmx_avatar(data, roles, seed = 1))
#' synpmx_scorecard_datatable(synpmx_scorecard(data, synthetic, roles))
synpmx_scorecard_datatable <- function(x, ...) {
  # `htmltools` is what `DT` itself is built on, so the second test only fails
  # on a broken installation; it is here so that a missing one is a message
  # rather than an error from inside the assembly below.
  if (!requireNamespace("DT", quietly = TRUE) ||
      !requireNamespace("htmltools", quietly = TRUE)) {
    message("DT is not installed, so the scorecard is printed uncoloured. ",
            "Install DT for the coloured table.")
    print(x)
    return(invisible(x))
  }
  verdicts <- names(.scorecard_verdict_colours)
  # Built through `do.call` so that the two defaults are defaults rather than
  # fixtures: passing `options` to a call that already names it is an error,
  # and a caller who wants paging or a caption should get it.
  arguments <- list(...)
  arguments$data <- as.data.frame(x)
  arguments$rownames <- arguments$rownames %||% FALSE
  arguments$options <- arguments$options %||% list(paging = FALSE)
  table <- do.call(DT::datatable, arguments)
  coloured <- DT::formatStyle(
    table, "verdict",
    color = DT::styleEqual(verdicts, unname(.scorecard_verdict_colours),
                           default = "inherit"),
    # `"unavailable"` is not bold. It marks a row nothing was measured for,
    # which is worth seeing but is not a finding to act on, and bolding it
    # would put it alongside the rows that are.
    fontWeight = DT::styleEqual(c("FAIL", "review"), c("bold", "bold"),
                                default = "normal"),
    # `"transparent"` rather than a white, so the table's own row striping and
    # hover still show through on the verdicts that are not tinted.
    backgroundColor = DT::styleEqual(names(.scorecard_verdict_fills),
                                     unname(.scorecard_verdict_fills),
                                     default = "transparent")
  )

  # The B5 levels come too. A reader who swaps a knitted card for a coloured
  # one is choosing a colour, not agreeing to lose the one row whose answer is
  # a list.
  parts <- list(coloured)
  rare <- attr(x, "rare_levels")
  if (!is.null(rare) && nrow(rare)) {
    parts <- c(parts, list(DT::datatable(
      rare, rownames = FALSE, options = list(paging = FALSE),
      caption = paste("B5: rare source levels that reached the output.",
                      "Each is one real patient's attribute, copied.")
    )))
  }
  # Spliced, not passed as one list: `tagList(list(...))` nests the list inside
  # a single element and the pieces stop rendering as siblings.
  do.call(htmltools::tagList, parts)
}
