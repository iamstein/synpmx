# B5 is the source-side census. A synthetic-side count -- is any level held by
# exactly one avatar -- was a second row on the card and is gone, because it is
# wrong in both directions. The tests that matter are the two cases where the
# two answers disagreed: a level one real patient held reaching the output, and
# a level many real patients held landing on a single avatar.

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

test_that("B5 reports a rare source level that was copied into the output", {
  source <- rare_source()
  roles <- rare_roles(strata = "ARM")
  synthetic <- rare_synthetic(source, roles)
  skip_if_not(any(as.character(synthetic$ARM) == "solo"),
              "the one-patient level did not reach this run's output")

  card <- synpmx_scorecard(source, synthetic, roles)

  # Review, not FAIL: on a small cohort a rare `RACE` lights this up constantly
  # and the answer is usually to stop carrying the covariate.
  expect_identical(card$verdict[card$check == "B5"], "review")
  # The cell counts; the levels themselves ride along on the card, because
  # which level it was is what decides what to do about it.
  expect_identical(card$result[card$check == "B5"], "1 of 1 exposed")
  rare <- attr(card, "rare_levels")
  expect_identical(rare$column, "ARM")
  expect_identical(rare$level, "solo")
  expect_identical(rare$source_patients, 1L)
  expect_identical(card$explore[card$check == "B5"],
                   "compare_pmx_rare_levels(source, synthetic, roles)")
  # And it reads the source, which is what makes the card restricted here.
  expect_identical(card$reads[card$check == "B5"], "both")
})

test_that("B5 passes when every level in the output is widely held", {
  source <- rare_source(solo = integer(), pair = integer())
  roles <- rare_roles(strata = "ARM")
  synthetic <- rare_synthetic(source, roles)

  card <- synpmx_scorecard(source, synthetic, roles)

  expect_identical(card$verdict[card$check == "B5"], "pass")
  expect_identical(card$result[card$check == "B5"], "0 of 0 exposed")
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

# The claim `avatar-algorithm.Rmd` makes about categorical covariates, pinned
# here because the vignette states it in prose rather than computing it.
#
# A categorical covariate has nothing to average, so a synthetic patient's
# category is always some real patient's category, copied. What decides whether
# a level escapes is the number of patients holding it, and the mechanism is
# geometric rather than protective: `.build_profiles()` one-hot encodes every
# level, so a sole holder sits alone on that axis, is nobody's nearest
# neighbour, and is almost never a donor. A second holder makes the two each
# other's nearest neighbour and the level travels between them.
propagation_fixture <- function(n_holders, n = 40L) {
  set.seed(1)
  do.call(rbind, lapply(seq_len(n), function(i) {
    times <- c(0, 1, 2, 4, 8, 24)
    data.frame(
      ID = i, TIME = c(0, times), NTIME = c(0, times),
      DV = c(NA, round(10 * exp(-0.1 * times) +
                         stats::rnorm(length(times), 0, 0.3), 3)),
      AMT = c(100, rep(0, length(times))),
      EVID = c(1L, rep(0L, length(times))),
      CMT = c(1L, rep(2L, length(times))),
      WT = round(stats::runif(1, 60, 80), 1),
      RARE = if (i <= n_holders) "rare-level" else "common",
      stringsAsFactors = FALSE
    )
  }))
}

avatars_carrying <- function(n_holders, n_subjects = 200L) {
  source <- propagation_fixture(n_holders)
  roles <- pmx_roles(
    id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
    cmt = "CMT", nominal_time = "NTIME", covariates = c("WT", "RARE")
  )
  synthetic <- suppressWarnings(suppressMessages(
    synpmx_avatar(source, roles, n_subjects = n_subjects, seed = 9)
  ))
  length(unique(synthetic$ID[synthetic$RARE == "rare-level"]))
}

test_that("a level held by one patient stays in, and two holders leak it", {
  expect_identical(avatars_carrying(1L), 0L)
  expect_gt(avatars_carrying(2L), 0L)
  # And it is the holder count that does it, not rarity as such: by ten holders
  # the level is ordinary and tracks its source frequency.
  expect_gt(avatars_carrying(10L), avatars_carrying(2L))
})
