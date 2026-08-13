# `pmx_scorecard()` is reporting code, so the load-bearing tests are the ones
# that prove it can say "FAIL". A scorecard that has only ever been seen to pass
# is an untested branch, not evidence.
#
# The rows it emits vary with the roles: no `strata` means no C3, and no
# categorical covariate means no B5. Both are tested, because a study with
# neither is the ordinary case at pharmacometric cohort sizes.

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

  card <- pmx_scorecard(source, synthetic, roles)

  expect_s3_class(card, "pmx_scorecard")
  expect_identical(names(card),
                   c("check", "question", "reads", "result", "verdict"))
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
  expect_true(all(sc_verdict(card, "A5") == "review"))
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

  card <- pmx_scorecard(source, copied, roles)

  expect_identical(sc_verdict(card, "B4b"), "FAIL")
  # B4a stays quiet, and that is correct rather than a miss. This fixture puts
  # every patient on one protocol grid, so a time vector belongs to the design
  # and not to a patient. What B4a watches is a vector too few real patients
  # share, which is the same rule `min_pattern_share` applies to visit sets.
  expect_identical(sc_verdict(card, "B4a"), "pass")
})

test_that("a lost endpoint and a lost patient are reported as failures", {
  source <- pmx_simulated_fixture(40)
  roles <- sc_roles()
  synthetic <- sc_synthetic(source, roles)
  settings <- attr(synthetic, "pmx_settings")

  dropped <- synthetic[synthetic$DVID == synthetic$DVID[[1L]], ]
  keep <- utils::head(unique(as.character(dropped[[roles$id]])), -1L)
  dropped <- dropped[as.character(dropped[[roles$id]]) %in% keep, ]
  attr(dropped, "pmx_settings") <- settings

  card <- pmx_scorecard(source, dropped, roles)

  expect_identical(sc_verdict(card, "A3"), "FAIL")
  expect_identical(sc_verdict(card, "A4"), "FAIL")
})

test_that("the optional rows appear only when the roles declare them", {
  source <- pmx_simulated_fixture(40)
  source$ARM <- ifelse(as.integer(factor(source$ID)) %% 2L == 0L, "A", "B")

  with_strata <- sc_roles(strata = "ARM")
  card <- pmx_scorecard(source, sc_synthetic(source, with_strata), with_strata)
  expect_true("C3" %in% card$check)
  expect_true("B5" %in% card$check)

  without <- sc_roles()
  bare <- pmx_scorecard(source, sc_synthetic(source, without), without)
  # `WT` is numeric, so with no `strata` there is no categorical axis at all.
  expect_false("C3" %in% bare$check)
  expect_false("B5" %in% bare$check)
})

test_that("B5 fails when a level reaches the output on one patient only", {
  source <- pmx_simulated_fixture(40)
  ids <- unique(as.character(source$ID))
  # One patient carries a level nobody else has. `strata` are copied from the
  # anchor verbatim, so it reaches the output on whoever anchors to them.
  source$ARM <- ifelse(as.character(source$ID) == ids[[1L]], "rare", "common")
  roles <- sc_roles(strata = "ARM")
  synthetic <- sc_synthetic(source, roles)
  skip_if_not(any(as.character(synthetic$ARM) == "rare"),
              "the rare level did not reach this run's output")

  card <- pmx_scorecard(source, synthetic, roles)

  expect_identical(sc_verdict(card, "B5"), "FAIL")
})

test_that("the scorecard refuses data that did not come from the generator", {
  source <- pmx_simulated_fixture(20)
  roles <- sc_roles()
  expect_error(pmx_scorecard(source, source, roles), "pmx_settings")
})

test_that("the scorecard is restricted output and reuses a given proximity", {
  source <- pmx_simulated_fixture(40)
  roles <- sc_roles()
  synthetic <- sc_synthetic(source, roles)
  proximity <- compare_pmx_proximity(source, synthetic, roles)

  card <- pmx_scorecard(source, synthetic, roles, proximity = proximity)

  expect_identical(attr(card, "release_status"), "restricted_not_releasable")
  expect_true(any(card$reads == "source"))
  expect_match(card$result[card$check == "B3"],
               sprintf("^%.3f", proximity$adversarial_accuracy))
})
