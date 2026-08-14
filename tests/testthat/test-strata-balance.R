# Strata balance preservation, and the identifiability screen scored within
# stratum. Both exist because the cohort is the wrong comparison group as soon
# as a study assigns patients to anything.

strata_fixture <- function(sizes = c(A = 20, B = 12, C = 2, D = 1), seed = 4) {
  set.seed(seed)
  arms <- rep(names(sizes), times = sizes)
  subject <- function(i, arm) {
    times <- c(0, 1, 2, 4, 8, 24)
    data.frame(
      ID = i, TIME = c(0, times), NTIME = c(0, times),
      DV = c(NA, round(10 * exp(-0.1 * times) + rnorm(6, 0, 0.3), 3)),
      AMT = c(100, rep(0, 6)), EVID = c(1L, rep(0L, 6)),
      CMT = c(1L, rep(2L, 6)), WT = round(runif(1, 60, 80), 1),
      ARM = arm, stringsAsFactors = FALSE
    )
  }
  do.call(rbind, Map(subject, seq_along(arms), arms))
}

strata_fixture_roles <- function() {
  pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
            cmt = "CMT", nominal_time = "NTIME", covariates = "WT",
            strata = "ARM")
}

arm_counts <- function(data) {
  table(factor(data$ARM[!duplicated(data$ID)], levels = c("A", "B", "C", "D")))
}

test_that("strata at or above the floor keep their exact share on every seed", {
  source <- strata_fixture()
  roles <- strata_fixture_roles()
  for (seed in 1:4) {
    synthetic <- suppressWarnings(suppressMessages(
      synpmx_avatar(source, roles, seed = seed)
    ))
    counts <- arm_counts(synthetic)
    expect_equal(as.integer(counts[["A"]]), 20L)
    expect_equal(as.integer(counts[["B"]]), 12L)
    # Cohort size is preserved whatever the small strata did.
    expect_equal(length(unique(synthetic$ID)), 35L)
    settings <- attr(synthetic, "pmx_settings")
    expect_equal(settings$strata_balanced, 2L)
    expect_equal(settings$strata_stochastic, 2L)
  }
})

test_that("strata below the floor are left stochastic, so their size is not disclosed", {
  source <- strata_fixture()
  roles <- strata_fixture_roles()
  small <- vapply(1:8, function(seed) {
    synthetic <- suppressWarnings(suppressMessages(
      synpmx_avatar(source, roles, seed = seed)
    ))
    as.integer(arm_counts(synthetic)[["D"]])
  }, integer(1))
  # D holds one real patient. Reproducing that exactly on every seed would say
  # so; the count has to move.
  expect_gt(length(unique(small)), 1L)
})

test_that("balance holds when the cohort size differs from the source", {
  source <- strata_fixture(c(A = 20, B = 20))
  roles <- strata_fixture_roles()
  synthetic <- suppressWarnings(suppressMessages(
    synpmx_avatar(source, roles, n_subjects = 20, seed = 1)
  ))
  counts <- arm_counts(synthetic)
  expect_equal(as.integer(counts[["A"]]), 10L)
  expect_equal(as.integer(counts[["B"]]), 10L)
  expect_equal(length(unique(synthetic$ID)), 20L)
})

test_that("preserve_strata_balance = FALSE restores the unstratified draw", {
  source <- strata_fixture(c(A = 20, B = 20))
  roles <- strata_fixture_roles()
  counts <- vapply(1:6, function(seed) {
    synthetic <- suppressWarnings(suppressMessages(
      synpmx_avatar(source, roles, seed = seed,
                    preserve_strata_balance = FALSE)
    ))
    as.integer(arm_counts(synthetic)[["A"]])
  }, integer(1))
  expect_gt(length(unique(counts)), 1L)
  expect_error(
    synpmx_avatar(source, roles, preserve_strata_balance = NA),
    "must be TRUE or FALSE"
  )
})

test_that("an avatar never leaves its anchor's stratum, balanced or not", {
  source <- strata_fixture(c(A = 20, B = 12))
  roles <- strata_fixture_roles()
  for (preserve in c(TRUE, FALSE)) {
    synthetic <- suppressWarnings(suppressMessages(
      synpmx_avatar(source, roles, seed = 2,
                    preserve_strata_balance = preserve)
    ))
    expect_setequal(unique(synthetic$ARM), unique(source$ARM))
    # One stratum per subject, still.
    per_subject <- tapply(synthetic$ARM, synthetic$ID,
                          function(x) length(unique(x)))
    expect_true(all(per_subject == 1L))
  }
})

test_that("the identifiability screen scores within stratum, not across arms", {
  skip_if_not_installed("xgxr")
  source <- as.data.frame(
    get(utils::data(list = "case1_pkpd", package = "xgxr"))
  )
  roles <- pmx_roles(
    id = "ID", time = "TIME", dv = "LIDV", amt = "AMT", evid = "EVID",
    cmt = "CMT", dvid = "NAME", nominal_time = "NOMTIME",
    strata = c("TRTACT", "DOSE"), covariates = "WEIGHTB"
  )
  synthetic <- suppressWarnings(suppressMessages(
    synpmx_avatar(source, roles, seed = 808)
  ))

  no_strata <- pmx_roles(
    id = "ID", time = "TIME", dv = "LIDV", amt = "AMT", evid = "EVID",
    cmt = "CMT", dvid = "NAME", nominal_time = "NOMTIME",
    covariates = "WEIGHTB"
  )
  # Six arms from 3 mg to 300 mg. Give one patient in a low arm the dose a
  # HIGHER arm was assigned: cohort-wide that dose is thirty other patients'
  # dose, so the gap to the nearest subject is 0 and nothing is reported. It is
  # only unusual next to the arm it was given in, which is the comparison group
  # the screen has to use.
  odd <- unique(synthetic$ID[synthetic$TRTACT == "10 mg"])[[1L]]
  high <- max(synthetic$AMT[synthetic$TRTACT == "300 mg"], na.rm = TRUE)
  synthetic$AMT[synthetic$ID == odd & synthetic$EVID != 0] <- high

  cohort_wide <- flag_identifiable_subjects(synthetic, no_strata)
  expect_false(cohort_wide$flagged[cohort_wide$subject_id == as.character(odd)])

  stratified <- flag_identifiable_subjects(synthetic, roles)
  expect_true(stratified$flagged[stratified$subject_id == as.character(odd)])
  expect_match(flagged <- stratified$outlier_axes[
    stratified$subject_id == as.character(odd)
  ], "dose magnitude")
})

test_that("a genuine within-stratum outlier is still flagged", {
  skip_if_not_installed("xgxr")
  source <- as.data.frame(
    get(utils::data(list = "case1_pkpd", package = "xgxr"))
  )
  roles <- pmx_roles(
    id = "ID", time = "TIME", dv = "LIDV", amt = "AMT", evid = "EVID",
    cmt = "CMT", dvid = "NAME", nominal_time = "NOMTIME",
    strata = c("TRTACT", "DOSE"), covariates = "WEIGHTB"
  )
  synthetic <- suppressWarnings(suppressMessages(
    synpmx_avatar(source, roles, seed = 808)
  ))
  odd <- unique(synthetic$ID[synthetic$TRTACT == "300 mg"])[[1L]]
  synthetic$AMT[synthetic$ID == odd & synthetic$EVID != 0] <- 600

  flags <- flag_identifiable_subjects(synthetic, roles)
  expect_true(flags$flagged[flags$subject_id == as.character(odd)])
  expect_match(flags$outlier_axes[flags$subject_id == as.character(odd)],
               "dose magnitude")
})

test_that("a stratum too small to score falls back to the cohort", {
  source <- strata_fixture(c(A = 20, B = 12, C = 2))
  roles <- strata_fixture_roles()
  synthetic <- suppressWarnings(suppressMessages(
    synpmx_avatar(source, roles, seed = 1)
  ))
  # C holds two patients, below the scoring floor, so it must not produce a
  # zero-scale stratum that flags everything or nothing by accident.
  expect_no_error(flag_identifiable_subjects(synthetic, roles))
})
