# compare_pmx_proximity() is the measurement for donor blending -- the one
# masking mechanism that acts on the values rather than the structure, and the
# only one that had no measurement attached.
#
# The load-bearing test here is the positive control. A privacy metric that has
# never been seen to fire is an untested branch, not evidence, so the first test
# hands it data it must reject.

px_roles <- function() {
  pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
            cmt = "CMT", covariates = "WT")
}

px_source <- function(n = 40L, seed = 1) {
  set.seed(seed)
  times <- c(0.5, 1, 2, 4, 8)
  weight <- round(stats::runif(n, 50, 100), 1)
  clearance <- stats::runif(n, 0.1, 0.4)
  do.call(rbind, lapply(seq_len(n), function(i) {
    data.frame(
      ID = as.character(i), TIME = c(0, times),
      DV = c(0, 10 * exp(-clearance[i] * times)),
      AMT = c(100, rep(0, length(times))),
      EVID = c(1L, rep(0L, length(times))),
      CMT = c(1L, rep(2L, length(times))),
      WT = weight[i], stringsAsFactors = FALSE
    )
  }))
}

test_that("a verbatim copy is caught (positive control)", {
  source <- px_source()
  roles <- px_roles()
  # The "synthetic" data is half the source, relabelled. Every synthetic subject
  # therefore sits exactly on a real one, so each one's nearest neighbour is its
  # own twin in the other set -- which is what the statistic is built to see.
  copied <- source[as.integer(source$ID) <= 20L, , drop = FALSE]
  copied$ID <- paste0("copy_", copied$ID)

  report <- compare_pmx_proximity(source, copied, roles, replicates = 30,
                                  seed = 2)
  expect_lt(report$adversarial_accuracy, report$null_lower)
  expect_match(report$verdict, "below the null")
  # And the distance quantile agrees: a copy is at zero distance from its twin.
  expect_lt(report$synthetic_to_source_q05, report$source_to_source_q05)
})

test_that("an ordinary AVATAR run is not flagged", {
  source <- px_source()
  roles <- px_roles()
  synthetic <- suppressWarnings(synpmx_avatar(source, roles, seed = 3))
  report <- compare_pmx_proximity(source, synthetic, roles, replicates = 30,
                                  seed = 2)
  expect_gte(report$adversarial_accuracy, report$null_lower)
  expect_false(grepl("below the null", report$verdict))
})

test_that("the statistic is scale-free where a raw distance is not", {
  # Doubling the cohort halves nothing about blending's safety, but it does pull
  # every nearest-neighbour distance down, which is why a raw distance to the
  # closest record cannot be the headline. The ratio-based statistic should not
  # move nearly as much.
  roles <- px_roles()
  small <- px_source(n = 30L, seed = 5)
  large <- px_source(n = 90L, seed = 5)
  report <- function(source) {
    synthetic <- suppressWarnings(synpmx_avatar(source, roles, seed = 4))
    compare_pmx_proximity(source, synthetic, roles, replicates = 20, seed = 2)
  }
  a <- report(small)
  b <- report(large)
  expect_lt(b$source_to_source_q05, a$source_to_source_q05)   # distances shrink
  expect_lt(abs(b$adversarial_accuracy - a$adversarial_accuracy), 0.35)
})

test_that("the null is built by the same code path as the observed value", {
  # Both arms of the comparison and both arms of the null run at one size, so
  # every small-sample artefact is present in each and cancels.
  source <- px_source(n = 30L)
  roles <- px_roles()
  synthetic <- suppressWarnings(synpmx_avatar(source, roles, seed = 6))
  report <- compare_pmx_proximity(source, synthetic, roles, replicates = 20,
                                  seed = 2)
  expect_equal(report$n_compared, 15L)          # min(nrow(real) %/% 2, n_synth)
  expect_lt(report$null_lower, report$null_upper)
  expect_true(report$null_lower > 0 && report$null_upper < 1)
})

test_that("too small a cohort reports that rather than a number", {
  source <- px_source(n = 4L)
  roles <- px_roles()
  synthetic <- suppressWarnings(synpmx_avatar(source, roles, seed = 7))
  report <- compare_pmx_proximity(source, synthetic, roles, replicates = 5)
  expect_true(is.na(report$adversarial_accuracy))
  expect_match(report$verdict, "too few subjects")
})

test_that("the report is marked restricted and leaves the RNG alone", {
  source <- px_source(n = 20L)
  roles <- px_roles()
  synthetic <- suppressWarnings(synpmx_avatar(source, roles, seed = 8))
  set.seed(99)
  before <- .Random.seed
  report <- compare_pmx_proximity(source, synthetic, roles, replicates = 10)
  expect_identical(.Random.seed, before)
  expect_identical(attr(report, "release_status"), "restricted_not_releasable")
  expect_error(
    compare_pmx_proximity(source, synthetic, roles, replicates = 0),
    "`replicates` must be one positive integer"
  )
})
