# `synpmx_scorecard()` is reporting code, so the load-bearing tests are the ones
# that prove it can say "FAIL". A scorecard that has only ever been seen to pass
# is an untested branch, not evidence.
#
# The rows it emits do NOT vary with the roles. A study that declares no
# `strata` and no categorical covariate -- the ordinary case at pharmacometric
# cohort sizes -- still gets C1, B5a and B5b, saying why they could not be
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
  # B3 is deliberately not asserted here: on a 40-patient fixture the null
  # interval is wide and the statistic can land outside it in the *utility*
  # direction, which is a reading about this fixture rather than a defect in
  # the scorecard.
  # The four structural guarantees must be exact, not "review".
  expect_identical(sc_verdict(card, "A1"), "pass")
  expect_identical(sc_verdict(card, "A4"), "pass")
  expect_identical(sc_verdict(card, "B1a"), "pass")
  expect_identical(sc_verdict(card, "B4a"), "pass")
  # And the rows where no threshold would be honest must not claim one.
  expect_identical(sc_verdict(card, "A5a"), "review")
  expect_identical(sc_verdict(card, "A5b"), "review")
  expect_identical(sc_verdict(card, "B2"), "review")
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
  expect_true(all(c("A2", "A4", "B3", "B5a", "B5b", "C1", "C2") %in%
                    card$check))
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

test_that("B5a marks a level that reaches the output on one patient only", {
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

  # Review rather than FAIL: this is the weak, synthetic-side form of the
  # check, and it is wrong in both directions. What the risk is made of is how
  # many SOURCE patients held the level, which nothing computes yet.
  expect_identical(sc_verdict(card, "B5a"), "review")
  # The count alone is not actionable: the row has to say which column and
  # which level, and point at the call that tabulates that column.
  expect_identical(card$result[card$check == "B5a"], "ARM = rare: 1")
  expect_identical(card$explore[card$check == "B5a"], "table(synthetic$ARM)")
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
  # was never asked. A6, B5a, B5b and C1 were each conditional on the study
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
                 B5a = "no categorical covariate or stratum",
                 B5b = "no categorical covariate or stratum",
                 C1 = "no strata declared")
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
  expect_identical(card$explore[card$check == "D1"],
                   "compare_pmx_distributions(source, synthetic, roles)")
  # Reading the card is not enough: two distributions with the same mean and
  # spread can be different shapes, and nothing here draws them.
  expect_output(print(card), "Plot source and synthetic")
})

test_that("D1 says so rather than erroring when no spread can be compared", {
  # A single observation per endpoint has no `sd`, so the ratio is NaN on every
  # variable. The row still has to appear -- a card that drops a row on some
  # studies cannot be compared across them.
  source <- pmx_simulated_fixture(20)
  roles <- sc_roles()
  synthetic <- sc_synthetic(source, roles)
  distributions <- compare_pmx_distributions(source, synthetic, roles)
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
