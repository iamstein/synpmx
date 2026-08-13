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

test_that("a stratum with no safe member is left rather than drawn from again", {
  # SIM-044. Keeping the move inside the arm was right, but "the stratum has
  # nobody else to offer" was read as "nobody at all" rather than "nobody who
  # can be masked". An arm in which every dose schedule is unique has plenty of
  # members, so the move stayed inside it and picked another unmaskable anchor.
  strata_key <- c("A", "A", "B", "B", "B")
  safe <- c(TRUE, TRUE, FALSE, FALSE, FALSE)

  move <- synpmx:::.reanchor_candidates(seq_along(strata_key), anchor = 3L,
                                        strata_key, "B", safe)

  expect_identical(move$candidates, c(1L, 2L))
  expect_true(move$crossed)
})

test_that("a safe member of the avatar's own arm outranks a safe one outside", {
  strata_key <- c("A", "A", "B", "B")
  safe <- c(TRUE, TRUE, FALSE, TRUE)

  move <- synpmx:::.reanchor_candidates(seq_along(strata_key), anchor = 3L,
                                        strata_key, "B", safe)

  expect_identical(move$candidates, 4L)
  expect_false(move$crossed)
})

test_that("with nobody safe anywhere the arm is kept rather than spent", {
  # Crossing cannot fix a cohort in which no anchor can be masked, so paying
  # the arm sizes for it buys nothing.
  strata_key <- c("A", "A", "B", "B")
  safe <- rep(FALSE, 4L)

  move <- synpmx:::.reanchor_candidates(seq_along(strata_key), anchor = 3L,
                                        strata_key, "B", safe)

  expect_identical(move$candidates, 4L)
  expect_false(move$crossed)
})

test_that("an arm whose dose schedules are all unique leaks nothing", {
  # The end-to-end form of SIM-044, recomputed from the run rather than from
  # `.reanchor_candidates()`: arm B doses every patient on their own days, so
  # no member of B is a safe anchor and all twelve of its avatars have to be
  # built on arm A instead.
  pieces <- lapply(1:24, function(s) {
    arm <- if (s <= 12L) "A" else "B"
    dose_times <- if (arm == "A") c(0, 24, 48) else c(s - 12, s + 12, s + 36)
    obs_times <- c(1, 4, 8, 25, 28, 49, 52)
    data.frame(
      ID = s, TIME = c(dose_times, obs_times),
      DV = c(rep(0, 3L), 10 + s * 0.1 + seq_along(obs_times)),
      AMT = c(rep(100, 3L), rep(0, length(obs_times))),
      EVID = c(rep(1L, 3L), rep(0L, length(obs_times))),
      CMT = 1L, DVID = "cp", WT = 70 + s, TRT = arm,
      stringsAsFactors = FALSE
    )
  })
  source <- do.call(rbind, pieces)
  source <- source[order(source$ID, source$TIME), ]
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", cmt = "CMT", dvid = "DVID",
                     covariates = "WT", strata = "TRT")

  # `coarsen_time = FALSE` because the grid would snap those dose days back
  # together and there would be nothing to mask.
  for (seed in 1:3) {
    synthetic <- suppressWarnings(suppressMessages(
      synpmx_avatar(source, roles, seed = seed, coarsen_time = FALSE)
    ))
    settings <- attr(synthetic, "pmx_settings")
    expect_identical(settings$identifying_dose_schedules, 0L)
    expect_identical(settings$identifying_visit_sets, 0L)
    # Recomputed from the finished table: no avatar holds a set of dose times
    # that exactly one real patient has.
    rare <- names(which(table(synpmx:::.dose_schedule_keys(source, roles)) < 2L))
    expect_length(
      intersect(synpmx:::.dose_schedule_keys(synthetic, roles), rare), 0L
    )
    # And the arm it could not hold is reported as crossed, not held silently.
    expect_gt(settings$strata_crossed, 0L)
  }
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
