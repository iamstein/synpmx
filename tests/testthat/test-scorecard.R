# `synpmx_scorecard()` is reporting code, so the load-bearing tests are the ones
# that prove it can say "FAIL". A scorecard that has only ever been seen to pass
# is an untested branch, not evidence.
#
# The rows it emits do NOT vary with the roles. A study that declares no
# `strata` and no categorical covariate -- the ordinary case at pharmacometric
# cohort sizes -- still gets C1 and B5, saying why they could not be
# asked, so two cards can be compared row for row.

sc_roles <- function(...) {
  pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
            cmt = "CMT", dvid = "DVID", covariates = "WT", ...)
}

sc_synthetic <- function(source, roles, seed = 1) {
  suppressWarnings(suppressMessages(synpmx_avatar(source, roles, seed = seed)))
}

sc_verdict <- function(card, check) card$verdict[card$check == check]

test_that("a clean run passes the guarantees and marks its judgement calls", {
  source <- pmx_simulated_fixture(40)
  roles <- sc_roles()
  synthetic <- sc_synthetic(source, roles)

  card <- synpmx_scorecard(source, synthetic, roles)

  expect_s3_class(card, "synpmx_scorecard")
  expect_identical(names(card),
                   c("check", "question", "reads", "result", "verdict",
                     "explore"))
  # Every row has to name something runnable. A blank `explore` is a row that
  # hands a reader a number and then abandons them.
  expect_true(all(nzchar(card$explore)))
  # B3's live value is deliberately not asserted here: on a 40-patient fixture
  # the null interval is wide and the statistic lands outside it in the
  # *utility* direction, which is a reading about this fixture rather than a
  # defect in the scorecard. The test below drives all three directions.
  # The four structural guarantees must be exact, not "review".
  expect_identical(sc_verdict(card, "A1"), "pass")
  expect_identical(sc_verdict(card, "A4"), "pass")
  expect_identical(sc_verdict(card, "B1a"), "pass")
  expect_identical(sc_verdict(card, "B4a"), "pass")
  # A study that came through with its event counts intact has nothing for a
  # reader to decide, so A5a and A5b pass rather than asking to be read.
  expect_identical(sc_verdict(card, "A5a"), "pass")
  expect_identical(sc_verdict(card, "A5b"), "pass")
})

test_that("A5a and A5b pass within 5% and review beyond it, never FAIL", {
  source <- pmx_simulated_fixture(40)
  roles <- sc_roles()
  synthetic <- sc_synthetic(source, roles)

  # Every subject holds the same number of observations here, so removing all
  # of one subject's moves the cohort mean by a known fraction of itself: one
  # subject in 40 is 2.5%, inside the tolerance, and three is 7.5%, outside it.
  # Avatars are given new identifiers, so the subjects to thin are the
  # synthetic table's own.
  subjects <- unique(synthetic$ID)
  thin <- function(ids) {
    synthetic[!(synthetic$EVID == 0 & synthetic$ID %in% ids), , drop = FALSE]
  }

  inside <- synpmx_scorecard(source, thin(subjects[1L]), roles)
  expect_identical(sc_verdict(inside, "A5a"), "pass")

  outside <- synpmx_scorecard(source, thin(subjects[1:3]), roles)
  expect_identical(sc_verdict(outside, "A5a"), "review")

  # However far it moves. An avatar whose dose course was cut back to reach the
  # B1b guarantee is the generator doing its job, and what is left is the
  # reader's judgement, so neither row can reach `FAIL`.
  far <- synpmx_scorecard(source, thin(subjects[1:20]), roles)
  expect_identical(sc_verdict(far, "A5a"), "review")
  expect_false(any(sc_verdict(far, "A5a") == "FAIL"))
})

test_that("B3 reviews only below its null interval, and says which side", {
  source <- pmx_simulated_fixture(40)
  roles <- sc_roles()
  synthetic <- sc_synthetic(source, roles)
  proximity <- compare_pmx_proximity(source, synthetic, roles)

  # The statistic is forced rather than provoked: a synthetic cohort that
  # actually memorises would have to be built to land below a null this wide,
  # and what is under test is the reading of the number, not the generator.
  at <- function(value) {
    forced <- proximity
    forced$adversarial_accuracy <- value
    as.data.frame(
      synpmx_scorecard(source, synthetic, roles, proximity = forced)
    )
  }
  b3 <- function(card) card$result[card$check == "B3"]

  inside <- at(mean(c(proximity$null_lower, proximity$null_upper)))
  expect_identical(sc_verdict(inside, "B3"), "pass")
  expect_match(b3(inside), " in \\[")

  # Memorisation: the one direction that answers what section B asks.
  below <- at(proximity$null_lower - 0.05)
  expect_identical(sc_verdict(below, "B3"), "review")
  expect_match(b3(below), " below \\[")

  # Separation costs utility and discloses nothing, so it passes -- and the
  # result still says which way it went, which is why the side is printed.
  above <- at(proximity$null_upper + 0.05)
  expect_identical(sc_verdict(above, "B3"), "pass")
  expect_match(b3(above), " above \\[")

  # A cohort too small to compare carries the reason, not "NA in [NA, NA]".
  expect_identical(
    .scorecard_proximity_result(
      data.frame(adversarial_accuracy = NA_real_,
                 verdict = "too few subjects to compare")
    ),
    "too few subjects to compare"
  )
})

test_that("B2 passes on an empty list and reviews a non-empty one", {
  source <- pmx_simulated_fixture(40)
  roles <- sc_roles()
  synthetic <- sc_synthetic(source, roles)

  card <- synpmx_scorecard(source, synthetic, roles)
  expect_identical(card$result[card$check == "B2"], "0 of 40")
  expect_identical(sc_verdict(card, "B2"), "pass")

  # One avatar followed three times as long as anybody else. The row becomes a
  # list of one record to read, and never a `FAIL`.
  odd <- unique(synthetic[[roles$id]])[[1L]]
  late <- synthetic[synthetic[[roles$id]] == odd, ][1L, ]
  late[[roles$time]] <- max(synthetic[[roles$time]]) * 3
  late[[roles$evid]] <- 0
  stretched <- rbind(synthetic, late)
  attr(stretched, "pmx_settings") <- attr(synthetic, "pmx_settings")

  card <- synpmx_scorecard(source, stretched, roles)
  expect_identical(sc_verdict(card, "B2"), "review")
  expect_match(card$result[card$check == "B2"], "^[1-9][0-9]* of 40$")
})

test_that("a verbatim copy of the source fails the exact-copy row", {
  source <- pmx_simulated_fixture(40)
  roles <- sc_roles()
  synthetic <- sc_synthetic(source, roles)
  # Keep the settings attribute -- the scorecard needs it -- but hand it back
  # the source's own event table under new identifiers. This is the positive
  # control: B4b exists to catch exactly this.
  copied <- source
  copied[[roles$id]] <- paste0("copy_", copied[[roles$id]])
  attr(copied, "pmx_settings") <- attr(synthetic, "pmx_settings")

  card <- synpmx_scorecard(source, copied, roles)

  expect_identical(sc_verdict(card, "B4b"), "FAIL")
  # B4a stays quiet, and that is correct rather than a miss. This fixture puts
  # every patient on one protocol grid, so a time vector belongs to the design
  # and not to a patient. What B4a watches is a vector too few real patients
  # share, which is the same rule `min_pattern_share` applies to visit sets.
  expect_identical(sc_verdict(card, "B4a"), "pass")
})

test_that("a lost endpoint fails and a lost patient is review", {
  source <- pmx_simulated_fixture(40)
  roles <- sc_roles()
  synthetic <- sc_synthetic(source, roles)
  settings <- attr(synthetic, "pmx_settings")

  dropped <- synthetic[synthetic$DVID == synthetic$DVID[[1L]], ]
  keep <- utils::head(unique(as.character(dropped[[roles$id]])), -1L)
  dropped <- dropped[as.character(dropped[[roles$id]]) %in% keep, ]
  attr(dropped, "pmx_settings") <- settings

  card <- synpmx_scorecard(source, dropped, roles)

  # An endpoint cannot go missing for a legitimate reason, so A3 is a failure.
  # A patient can: `on_donor_shortfall = "drop"` removes one that could not be
  # built, and a one-patient stratum can legitimately come back empty, so A4
  # says review however far the count moved.
  expect_identical(sc_verdict(card, "A3"), "FAIL")
  expect_identical(sc_verdict(card, "A4"), "review")
})

test_that("only the rows that are always a defect can say FAIL", {
  # The whole list, pinned. Adding a FAIL to any other row is a decision about
  # what the package calls broken, and it should have to change this test.
  source <- pmx_simulated_fixture(40)
  source$ARM <- ifelse(as.integer(factor(source$ID)) %% 2L == 0L, "A", "B")
  roles <- sc_roles(strata = "ARM")
  card <- synpmx_scorecard(source, sc_synthetic(source, roles), roles)

  can_fail <- c("A1", "A3", "A6", "B1a", "B1b", "B4a", "B4b")
  expect_true(all(card$check[card$verdict == "FAIL"] %in% can_fail))
  # And the softened rows are present, so this is not passing by their absence.
  expect_true(all(c("A2", "A4", "B3", "B5", "C1", "C2") %in% card$check))
})

test_that("an objection to the source is review, not FAIL", {
  source <- pmx_simulated_fixture(20)
  roles <- sc_roles()
  synthetic <- sc_synthetic(source, roles)
  # Break the source after generating from it: time must not run backwards
  # within a subject. Generation would refuse this, which is why the row exists
  # -- it reports on the source you handed in, whatever the validator makes of
  # it, and being handed a real study a validator objects to is ordinary.
  broken <- source
  broken[[roles$time]][[2L]] <- broken[[roles$time]][[1L]] - 1

  card <- synpmx_scorecard(broken, synthetic, roles)

  expect_identical(sc_verdict(card, "A2"), "review")
})

test_that("B5 marks a level that reaches the output on one patient only", {
  source <- pmx_simulated_fixture(40)
  ids <- unique(as.character(source$ID))
  # One patient carries a level nobody else has. `strata` are copied from the
  # anchor verbatim, so it reaches the output on whoever anchors to them.
  source$ARM <- ifelse(as.character(source$ID) == ids[[1L]], "rare", "common")
  roles <- sc_roles(strata = "ARM")
  synthetic <- sc_synthetic(source, roles)
  skip_if_not(any(as.character(synthetic$ARM) == "rare"),
              "the rare level did not reach this run's output")

  card <- synpmx_scorecard(source, synthetic, roles)

  # Review rather than FAIL: on a small cohort a rare `RACE` lights this up
  # constantly and the answer is usually to stop carrying the covariate.
  expect_identical(sc_verdict(card, "B5"), "review")
  # The count alone is not actionable, so the levels ride along on the card.
  expect_identical(card$result[card$check == "B5"], "1 of 1 exposed")
  rare <- attr(card, "rare_levels")
  expect_identical(rare$column, "ARM")
  expect_identical(rare$level, "rare")
})

test_that("a stratum that changed size is review, not FAIL", {
  source <- pmx_simulated_fixture(40)
  source$ARM <- ifelse(as.integer(factor(source$ID)) %% 2L == 0L, "A", "B")
  roles <- sc_roles(strata = "ARM")
  # Anchors are drawn with replacement, so switching the balance guarantee off
  # is what makes an arm change size. Dropping a subject for want of donors does
  # the same thing on a real study, and neither is a defect -- which is why this
  # row cannot say FAIL.
  synthetic <- suppressWarnings(suppressMessages(
    synpmx_avatar(source, roles, seed = 3, preserve_strata_balance = FALSE)
  ))
  arms <- table(as.character(synthetic$ARM)[!duplicated(synthetic$ID)])
  skip_if(identical(as.integer(arms), c(20L, 20L)),
          "this seed happened to reproduce both arm sizes exactly")

  card <- synpmx_scorecard(source, synthetic, roles)

  expect_identical(sc_verdict(card, "C1"), "review")
})

test_that("an arm that lost an endpoint is C3, and A3 does not see it", {
  # A3 compares endpoint sets across the whole cohort, so an arm losing an
  # endpoint the rest of the study still holds passes it. C3 is the only row
  # that looks arm by arm, and the vignette described this check for weeks
  # before anything computed it.
  source <- pmx_simulated_fixture(40)
  source$ARM <- ifelse(as.integer(factor(source$ID)) %% 2L == 0L, "A", "B")
  roles <- sc_roles(strata = "ARM")
  synthetic <- sc_synthetic(source, roles, seed = 4)

  endpoint <- sort(unique(as.character(synthetic$DVID)))[[1L]]
  lost <- synthetic[!(synthetic$ARM == "A" & synthetic$DVID == endpoint &
                        synthetic$EVID == 0), ]
  attr(lost, "pmx_settings") <- attr(synthetic, "pmx_settings")

  card <- suppressWarnings(synpmx_scorecard(source, lost, roles))

  expect_identical(sc_verdict(card, "C3"), "review")
  expect_identical(card$result[card$check == "C3"], "1 of 2")
  expect_identical(sc_verdict(card, "A3"), "pass")
  expect_identical(sc_verdict(synpmx_scorecard(source, synthetic, roles), "C3"),
                   "pass")
})

test_that("a table with no run record is scored, not refused", {
  # The card used to stop() without `pmx_settings`, which made the one function
  # meant to score an output unable to score anything but this package's own --
  # including this package's own, once it had been through `write.csv()`. The
  # three rows that need the run's record say so; the rest are measured from the
  # two tables as usual.
  source <- pmx_simulated_fixture(24)
  roles <- sc_roles()
  synthetic <- sc_synthetic(source, roles)
  full <- synpmx_scorecard(source, synthetic, roles)

  bare <- synthetic
  attributes(bare) <- attributes(bare)[c("names", "class", "row.names")]
  card <- synpmx_scorecard(source, bare, roles)

  recorded <- c("B1a", "B1b", "C2")
  expect_identical(card$verdict[card$check %in% recorded],
                   rep("unavailable", length(recorded)))
  expect_identical(card$result[card$check %in% recorded],
                   rep("no run record", length(recorded)))
  # Same rows in the same order, and every other verdict unchanged: a card that
  # holds different rows on different tables cannot be compared across them.
  expect_identical(card$check, full$check)
  expect_identical(card$verdict[!card$check %in% recorded],
                   full$verdict[!full$check %in% recorded])
  # `unavailable` is not `pass`: the count line has to say so.
  expect_output(print(card), "3 unanswered")
})

test_that("C2 reads the run's record, and says that it does", {
  # Both of its counts come from `pmx_settings`, so labelling the row `both`
  # claimed the filled-in card was restricted output on the strength of a row
  # that reads neither table.
  source <- pmx_simulated_fixture(20)
  roles <- sc_roles()
  synthetic <- sc_synthetic(source, roles)

  card <- synpmx_scorecard(source, synthetic, roles)

  expect_identical(card$reads[card$check == "C2"], "run settings")
})

test_that("every card holds the same checks, whatever the roles declare", {
  # A card whose rows depend on the study cannot be compared with another card,
  # and an absent row reads as a check that passed when it means the question
  # was never asked. A6, B5 and C1 were each conditional on the study
  # having something for them to ask.
  source <- pmx_simulated_fixture(40)
  full_roles <- sc_roles()
  bare_roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                          evid = "EVID", cmt = "CMT", dvid = "DVID")

  full <- synpmx_scorecard(source, sc_synthetic(source, full_roles), full_roles)
  bare <- synpmx_scorecard(source, sc_synthetic(source, bare_roles), bare_roles)

  expect_identical(bare$check, full$check)
  # And each says why it could not be asked, rather than reporting a number.
  not_asked <- c(A6 = "no discrete endpoint",
                 B5 = "no categorical covariate or stratum",
                 C1 = "no strata declared",
                 C3 = "no strata declared")
  for (check in names(not_asked)) {
    expect_identical(bare$result[bare$check == check], not_asked[[check]])
    expect_identical(bare$verdict[bare$check == check], "pass")
  }
})

test_that("D1 reports the spread that moved furthest, and never passes", {
  # The utility tier was on the card in the vignette and absent from the
  # function, so a reader working from the printed card never saw the one
  # quantity the algorithm moves on purpose.
  source <- pmx_simulated_fixture(40)
  roles <- sc_roles()
  synthetic <- sc_synthetic(source, roles)

  card <- synpmx_scorecard(source, synthetic, roles)

  expect_identical(card$verdict[card$check == "D1"], "review")
  expect_identical(card$reads[card$check == "D1"], "both")
  expect_match(card$result[card$check == "D1"],
               "^sd x[0-9.]+ on .+ \\(furthest of [0-9]+\\)$")
  expect_identical(
    card$explore[card$check == "D1"],
    'compare_pmx_distributions(source, synthetic, roles, output = "tables")'
  )
})

test_that("D1 says so rather than erroring when no spread can be compared", {
  # A single observation per endpoint has no `sd`, so the ratio is NaN on every
  # variable. The row still has to appear -- a card that drops a row on some
  # studies cannot be compared across them.
  source <- pmx_simulated_fixture(20)
  roles <- sc_roles()
  synthetic <- sc_synthetic(source, roles)
  distributions <- compare_pmx_distributions(source, synthetic, roles,
                                             output = "tables")
  distributions$endpoints$sd <- NA_real_
  distributions$covariates_numeric$sd <- NA_real_

  expect_null(.scorecard_spread(distributions))
})

test_that("the scorecard is restricted output and reuses a given proximity", {
  source <- pmx_simulated_fixture(40)
  roles <- sc_roles()
  synthetic <- sc_synthetic(source, roles)
  proximity <- compare_pmx_proximity(source, synthetic, roles)

  card <- synpmx_scorecard(source, synthetic, roles, proximity = proximity)

  expect_identical(attr(card, "release_status"), "restricted_not_releasable")
  expect_true(any(card$reads == "source"))
  expect_match(card$result[card$check == "B3"],
               sprintf("^%.3f", proximity$adversarial_accuracy))
})

test_that("the datatable colours every verdict the card can carry", {
  # The colouring is a lookup keyed on the verdict strings, so the test worth
  # having is that the keys and the verdicts cannot drift apart: a renamed or
  # added verdict would silently stop being coloured rather than error.
  skip_if_not_installed("DT")
  source <- pmx_simulated_fixture(30)
  roles <- sc_roles()
  # Settings stripped, so the card carries "unavailable" rows as well.
  synthetic <- sc_synthetic(source, roles)
  attr(synthetic, "pmx_settings") <- NULL
  card <- synpmx_scorecard(source, synthetic, roles)

  expect_true(any(card$verdict == "unavailable"))
  expect_true(all(setdiff(unique(card$verdict), "pass") %in%
                    names(.scorecard_verdict_colours)))

  shown <- synpmx_scorecard_datatable(card)

  expect_s3_class(shown, "shiny.tag.list")
  expect_s3_class(shown[[1]], "datatables")
  rendered <- paste(unlist(shown[[1]]$x), collapse = " ")
  expect_true(all(vapply(c(.scorecard_verdict_colours, .scorecard_verdict_fills),
                         function(colour) grepl(colour, rendered, fixed = TRUE),
                         logical(1))))
  # Every tint is a verdict the text palette also names, so the two cannot
  # drift into tinting something that is left uncoloured.
  expect_true(all(names(.scorecard_verdict_fills) %in%
                    names(.scorecard_verdict_colours)))
})

test_that("the datatable keeps the B5 detail that knitting the card shows", {
  # `knit_print()` emits the rare-level table as well as the card, because the
  # card can only name one level in a cell. Colouring must not be a way to
  # quietly drop the rest of them.
  skip_if_not_installed("DT")
  source <- pmx_simulated_fixture(30)
  roles <- sc_roles()
  card <- synpmx_scorecard(source, sc_synthetic(source, roles), roles)
  attr(card, "rare_levels") <- data.frame(
    column = "RACE", level = "OTHER",
    source_patients = 1L, synthetic_patients = 2L,
    stringsAsFactors = FALSE
  )

  shown <- synpmx_scorecard_datatable(card)

  expect_length(Filter(function(part) inherits(part, "datatables"), shown), 2L)
  expect_match(paste(unlist(shown), collapse = " "), "OTHER")
})

test_that("the datatable says so and prints the card when DT is missing", {
  source <- pmx_simulated_fixture(20)
  roles <- sc_roles()
  card <- synpmx_scorecard(source, sc_synthetic(source, roles), roles)
  local_mocked_bindings(requireNamespace = function(...) FALSE,
                        .package = "base")

  printed <- capture.output(
    expect_message(result <- synpmx_scorecard_datatable(card),
                   "DT is not installed")
  )

  expect_s3_class(result, "synpmx_scorecard")
  expect_true(any(grepl("verdict", printed)))
})

test_that("B4a is not applicable where attendance is drawn per visit", {
  data <- pmx_simulated_fixture(30)
  roles <- pmx_roles(id = "ID", time = "TIME", nominal_time = "NTIME",
                     dv = "DV", amt = "AMT", evid = "EVID", cmt = "CMT",
                     dvid = "DVID", mdv = "MDV")
  synthetic <- synpmx_pca(data, roles, seed = 3)
  card <- as.data.frame(synpmx_scorecard(data, synthetic, roles))
  b4a <- card[card$check == "B4a", ]
  expect_identical(b4a$verdict, "not applicable")
  expect_match(b4a$result, "drawn per visit")
  # B4b is the row that answers the disclosure question, and it is untouched.
  expect_true(card$verdict[card$check == "B4b"] %in% c("pass", "FAIL"))
})

test_that("B4a is still computed where attendance is copied from donors", {
  data <- pmx_simulated_fixture(30)
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", cmt = "CMT", dvid = "DVID")
  synthetic <- suppressWarnings(synpmx_avatar(data, roles, seed = 1))
  card <- as.data.frame(synpmx_scorecard(data, synthetic, roles))
  expect_false(card$verdict[card$check == "B4a"] == "not applicable")
})

test_that("a table with no source attribute is measured the AVATAR way", {
  # Through `write.csv()` the attribute is gone, and the row falls back rather
  # than silently reading `not applicable` on a table nothing is known about.
  data <- pmx_simulated_fixture(30)
  roles <- pmx_roles(id = "ID", time = "TIME", nominal_time = "NTIME",
                     dv = "DV", amt = "AMT", evid = "EVID", cmt = "CMT",
                     dvid = "DVID", mdv = "MDV")
  synthetic <- synpmx_pca(data, roles, seed = 3)
  attr(synthetic, "pmx_source") <- NULL
  card <- as.data.frame(synpmx_scorecard(data, synthetic, roles))
  expect_false(card$verdict[card$check == "B4a"] == "not applicable")
})
