# The census exists because the synthetic-side check (B5a) is wrong in both
# directions, so the tests that matter are the two cases where the two answers
# disagree: a level one real patient held reaching the output, and a level many
# real patients held landing on a single avatar.

rare_source <- function(n = 40L, solo = 1L, pair = 2:3) {
  source <- pmx_simulated_fixture(n)
  ids <- unique(as.character(source$ID))
  source$ARM <- "common"
  source$ARM[as.character(source$ID) %in% ids[solo]] <- "solo"
  source$ARM[as.character(source$ID) %in% ids[pair]] <- "pair"
  source
}

rare_roles <- function(covariates = "WT", ...) {
  pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
            cmt = "CMT", dvid = "DVID", covariates = covariates, ...)
}

rare_synthetic <- function(source, roles, seed = 2) {
  suppressWarnings(suppressMessages(synpmx_avatar(source, roles, seed = seed)))
}

test_that("the census counts holders on both sides and marks the exposed", {
  source <- rare_source()
  roles <- rare_roles(strata = "ARM")
  synthetic <- rare_synthetic(source, roles)

  census <- compare_pmx_rare_levels(source, synthetic, roles)

  expect_s3_class(census, "pmx_rare_levels")
  expect_identical(names(census),
                   c("column", "level", "source_patients",
                     "synthetic_patients", "exposed", "reached"))
  # `solo` is held by one real patient and `pair` by two, so at the default
  # floor of 2 exactly one level is exposed.
  expect_identical(census$source_patients[census$level == "solo"], 1L)
  expect_identical(census$source_patients[census$level == "pair"], 2L)
  expect_true(census$exposed[census$level == "solo"])
  expect_false(census$exposed[census$level == "pair"])
  expect_false(census$exposed[census$level == "common"])
  # It reads the source on both sides, so it cannot leave the environment.
  expect_identical(attr(census, "release_status"), "restricted_not_releasable")
})

test_that("the floor comes from the run, and can be raised by hand", {
  source <- rare_source()
  roles <- rare_roles(strata = "ARM")
  synthetic <- rare_synthetic(source, roles)

  expect_identical(attr(compare_pmx_rare_levels(source, synthetic, roles),
                        "floor"), 2L)
  # At a floor of 3 the two-patient level is exposed as well.
  raised <- compare_pmx_rare_levels(source, synthetic, roles, floor = 3)
  expect_true(raised$exposed[raised$level == "pair"])
  expect_error(compare_pmx_rare_levels(source, synthetic, roles, floor = 0),
               "one integer of 1 or more")
})

test_that("no categorical axis means an empty census rather than an error", {
  source <- pmx_simulated_fixture(20)
  roles <- rare_roles()   # `WT` is numeric, and no strata are declared

  census <- compare_pmx_rare_levels(source, rare_synthetic(source, roles),
                                    roles)

  expect_identical(nrow(census), 0L)
  expect_output(print(census), "nothing to census")
})

test_that("B5b reports a rare source level that was copied into the output", {
  source <- rare_source()
  roles <- rare_roles(strata = "ARM")
  synthetic <- rare_synthetic(source, roles)
  skip_if_not(any(as.character(synthetic$ARM) == "solo"),
              "the one-patient level did not reach this run's output")

  card <- synpmx_scorecard(source, synthetic, roles)

  # Review, not FAIL: on a small cohort a rare `RACE` lights this up constantly
  # and the answer is usually to stop carrying the covariate.
  expect_identical(card$verdict[card$check == "B5b"], "review")
  # The cell counts; the levels themselves ride along on the card, because
  # which level it was is what decides what to do about it.
  expect_identical(card$result[card$check == "B5b"], "1 of 1 exposed")
  rare <- attr(card, "rare_levels")
  expect_identical(rare$column, "ARM")
  expect_identical(rare$level, "solo")
  expect_identical(rare$source_patients, 1L)
  expect_identical(card$explore[card$check == "B5b"],
                   "compare_pmx_rare_levels(source, synthetic, roles)")
  # And it reads the source, which is what makes the card restricted here.
  expect_identical(card$reads[card$check == "B5b"], "both")
})

test_that("B5b passes when every level in the output is widely held", {
  source <- rare_source(solo = integer(), pair = integer())
  roles <- rare_roles(strata = "ARM")
  synthetic <- rare_synthetic(source, roles)

  card <- synpmx_scorecard(source, synthetic, roles)

  expect_identical(card$verdict[card$check == "B5b"], "pass")
  expect_identical(card$result[card$check == "B5b"], "0 of 0 exposed")
})

test_that("the census counts a blank level like any other", {
  # A blank cell is a level a real dataset carries, and `[[` on a name never
  # matches it, which used to take the census out of bounds rather than count it.
  # A blank *stratum* is rejected outright, so this is a covariate: the census
  # has to survive one either way.
  source <- rare_source()
  source$SEX <- "F"
  source$SEX[as.character(source$ID) == unique(as.character(source$ID))[[4L]]] <- ""
  roles <- rare_roles(covariates = c("WT", "SEX"))
  synthetic <- rare_synthetic(source, roles)

  census <- compare_pmx_rare_levels(source, synthetic, roles)

  expect_true("" %in% census$level)
  expect_identical(census$source_patients[census$level == ""], 1L)
  expect_true(census$exposed[census$level == ""])
  # Printed as something rather than as nothing.
  expect_output(print(census), "<blank>")
})
