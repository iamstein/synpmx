# An anchor that cannot be given a masking visit set is swapped for another
# one. The swap used to draw from the whole eligible pool, so the avatar took
# the replacement's `strata` values and quietly changed arm -- while the run
# still reported the balance `.strata_targets()` had planned.

test_that("a re-anchored avatar stays in its own stratum", {
  strata_key <- c("A", "A", "A", "B", "B", "B")
  allowed <- seq_along(strata_key)

  move <- synpmx:::.reanchor_candidates(allowed, anchor = 1L, strata_key, "A")

  expect_identical(move$candidates, c(2L, 3L))
  expect_false(move$crossed)
})

test_that("the whole pool is the fallback, and the crossing is reported", {
  # A stratum whose only member is this anchor has nobody else to offer, and
  # emitting one real patient's visit set is worse than moving arm -- but the
  # move has to be counted, because it is why C3 will not add up.
  strata_key <- c("SOLO", "B", "B", "B")
  allowed <- seq_along(strata_key)

  move <- synpmx:::.reanchor_candidates(allowed, anchor = 1L, strata_key,
                                        "SOLO")

  expect_identical(move$candidates, c(2L, 3L, 4L))
  expect_true(move$crossed)
})

test_that("with no strata declared every other anchor is a candidate", {
  allowed <- 1:4

  move <- synpmx:::.reanchor_candidates(allowed, anchor = 2L, NULL, NULL)

  expect_identical(move$candidates, c(1L, 3L, 4L))
  expect_false(move$crossed)
})

test_that("an anchor excluded from the pool is not offered back to itself", {
  strata_key <- c("A", "A")
  move <- synpmx:::.reanchor_candidates(1:2, anchor = 1L, strata_key, "A")
  expect_false(1L %in% move$candidates)
})

test_that("the run records how many avatars changed arm", {
  source <- pmx_simulated_fixture(40)
  ids <- unique(as.character(source$ID))
  source$TRT <- rep(c("A", "B"), each = 20L)[match(as.character(source$ID), ids)]
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", cmt = "CMT", dvid = "DVID",
                     covariates = "WT", strata = "TRT")

  synthetic <- suppressWarnings(suppressMessages(
    synpmx_avatar(source, roles, seed = 1)
  ))

  settings <- attr(synthetic, "pmx_settings")
  expect_identical(settings$strata_crossed, 0L)
  # And with nothing crossing, the arms are exactly their source sizes.
  sizes <- compare_pmx_strata_sizes(source, synthetic, roles)
  expect_identical(sizes$synthetic_patients, sizes$source_patients)
})
