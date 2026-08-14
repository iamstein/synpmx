# An anchor that cannot be given a masking visit set is swapped for another
# one. The swap used to draw from the whole eligible pool, so the avatar took
# the replacement's `strata` values and quietly changed arm -- while the run
# still reported the balance `.strata_targets()` had planned. The move now stays
# in the arm always: where the arm cannot mask its own avatars, the run says so
# by name and B1a/B1b read FAIL, rather than the arm sizes being spent on it.

test_that("a re-anchored avatar stays in its own stratum", {
  strata_key <- c("A", "A", "A", "B", "B", "B")
  allowed <- seq_along(strata_key)

  move <- synpmx:::.reanchor_candidates(allowed, anchor = 1L, strata_key, "A")

  expect_identical(move$candidates, c(2L, 3L))
})

test_that("a stratum with nobody else to offer offers nobody", {
  # A one-member arm has no other anchor, and the answer is to keep the one it
  # has -- not to reach into another arm, which would put a patient in an arm
  # they were never in and leave the run reporting a balance it did not build.
  strata_key <- c("SOLO", "B", "B", "B")
  allowed <- seq_along(strata_key)

  move <- synpmx:::.reanchor_candidates(allowed, anchor = 1L, strata_key,
                                        "SOLO")

  expect_length(move$candidates, 0L)
})

test_that("with no strata declared every other anchor is a candidate", {
  allowed <- 1:4

  move <- synpmx:::.reanchor_candidates(allowed, anchor = 2L, NULL, NULL)

  expect_identical(move$candidates, c(1L, 3L, 4L))
})

test_that("an anchor excluded from the pool is not offered back to itself", {
  strata_key <- c("A", "A")
  move <- synpmx:::.reanchor_candidates(1:2, anchor = 1L, strata_key, "A")
  expect_false(1L %in% move$candidates)
})

test_that("a stratum with no safe member offers nobody rather than anybody", {
  # SIM-044. Keeping the move inside the arm was right, but "the stratum has
  # nobody else to offer" was read as "nobody at all" rather than "nobody who
  # can be masked". An arm in which every dose schedule is unique has plenty of
  # members, so the move stayed inside it and picked another unmaskable anchor
  # -- twelve times, and then emitted one anyway. Re-drawing among patients who
  # are equally unmaskable only changes which real patient is exposed.
  strata_key <- c("A", "A", "B", "B", "B")
  safe <- c(TRUE, TRUE, FALSE, FALSE, FALSE)

  move <- synpmx:::.reanchor_candidates(seq_along(strata_key), anchor = 3L,
                                        strata_key, "B", safe)

  expect_length(move$candidates, 0L)
})

test_that("only the safe members of the avatar's own arm are offered", {
  strata_key <- c("A", "A", "B", "B")
  safe <- c(TRUE, TRUE, FALSE, TRUE)

  move <- synpmx:::.reanchor_candidates(seq_along(strata_key), anchor = 3L,
                                        strata_key, "B", safe)

  expect_identical(move$candidates, 4L)
})

test_that("a tier nobody qualifies for falls through to the next one", {
  # Where every visit set in the study is unique, "safe on both counts" is
  # nobody. Reading that as "nothing to prefer" dropped the dose condition too
  # and drew from the unmaskable members again, so the tiers are tried in turn.
  strata_key <- c("A", "A", "B", "B", "B")
  both <- rep(FALSE, 5L)
  dose_safe <- c(TRUE, TRUE, FALSE, FALSE, TRUE)

  move <- synpmx:::.reanchor_candidates(seq_along(strata_key), anchor = 3L,
                                        strata_key, "B",
                                        list(both, dose_safe))

  expect_identical(move$candidates, 5L)
})

test_that("an arm whose dose schedules are all unique fails by name", {
  # SIM-044 end to end, recomputed from the run rather than from
  # `.reanchor_candidates()`: arm B doses every patient on their own days, so no
  # member of B is a safe anchor and no avatar allocated to B can be masked.
  # Arm A is not borrowed from -- the arm sizes are what a reader takes the
  # study to be -- so the run reports the leak against arm B and B1b fails.
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

  # Answerable before generating: arm B has no safe anchor, arm A has twelve.
  arms <- unmaskable_strata(source, roles, coarsen_time = FALSE)
  expect_identical(arms$safe_anchors[arms$stratum == "B"], 0L)
  expect_gt(arms$safe_anchors[arms$stratum == "A"], 0L)

  # `coarsen_time = FALSE` because the grid would snap those dose days back
  # together and there would be nothing to mask.
  for (seed in 1:3) {
    synthetic <- suppressWarnings(suppressMessages(
      synpmx_avatar(source, roles, seed = seed, coarsen_time = FALSE)
    ))
    settings <- attr(synthetic, "pmx_settings")
    # The leak is arm B's alone, and it is named rather than counted.
    expect_identical(names(settings$identifying_dose_schedule_strata), "B")
    expect_identical(settings$identifying_visit_sets, 0L)
    # Every avatar holding a rare dose schedule is in arm B, and the arms are
    # still their source sizes -- which is the trade this makes.
    rare <- names(which(table(synpmx:::.dose_schedule_keys(source, roles)) < 2L))
    leaked <- names(which(
      synpmx:::.dose_schedule_keys(synthetic, roles) %in% rare
    ))
    first <- !duplicated(as.character(synthetic$ID))
    arm_of <- stats::setNames(as.character(synthetic$TRT)[first],
                              as.character(synthetic$ID)[first])
    expect_true(all(arm_of[leaked] == "B"))
    sizes <- compare_pmx_strata_sizes(source, synthetic, roles)
    expect_identical(sizes$synthetic_patients, sizes$source_patients)

    card <- synpmx_scorecard(source, synthetic, roles)
    expect_identical(card$verdict[card$check == "B1b"], "FAIL")
    expect_match(card$result[card$check == "B1b"], "in B \\(")
  }
})

test_that("a study of unique visit sets still masks what it can", {
  # The tier fall-through end to end. Every patient here misses a different
  # visit, so no anchor is visit-safe and the combined flag is empty for the
  # whole cohort -- which used to take the dose condition down with it, leaving
  # avatars in arms A-C on dose days those arms could have masked.
  visits <- seq(0, 12L) * 24
  pieces <- lapply(1:16, function(s) {
    arm <- c("A", "B", "C", "D")[[(s - 1L) %/% 4L + 1L]]
    dose_times <- if (arm == "D") {
      (seq_len(10L) * 24)[-(s %% 4L + 2L)]   # each patient's own dose days
    } else {
      seq_len(10L) * 24
    }
    obs_times <- visits[-(s %% 11L + 2L)]    # each patient's own missed visit
    data.frame(
      ID = s, TIME = c(dose_times, obs_times),
      DV = c(rep(0, length(dose_times)), 10 + s * 0.1 + seq_along(obs_times)),
      AMT = c(rep(100, length(dose_times)), rep(0, length(obs_times))),
      EVID = c(rep(1L, length(dose_times)), rep(0L, length(obs_times))),
      CMT = 1L, DVID = "cp", WT = 70 + s, TRT = arm,
      stringsAsFactors = FALSE
    )
  })
  source <- do.call(rbind, pieces)
  source <- source[order(source$ID, source$TIME), ]
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", cmt = "CMT", dvid = "DVID",
                     covariates = "WT", strata = "TRT")

  for (seed in 1:3) {
    synthetic <- suppressWarnings(suppressMessages(
      synpmx_avatar(source, roles, seed = seed, coarsen_time = FALSE)
    ))
    settings <- attr(synthetic, "pmx_settings")
    # Arm D cannot be masked and says so; nothing leaks out of A, B or C.
    expect_identical(names(settings$identifying_dose_schedule_strata), "D")
    sizes <- compare_pmx_strata_sizes(source, synthetic, roles)
    expect_identical(sizes$synthetic_patients, sizes$source_patients)
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

  # The arms are exactly their source sizes, which is now unconditional: no
  # mechanism moves an avatar between arms.
  sizes <- compare_pmx_strata_sizes(source, synthetic, roles)
  expect_identical(sizes$synthetic_patients, sizes$source_patients)
})
