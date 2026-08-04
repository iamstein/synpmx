# Two mechanisms that both rest on the same decomposition as `coarsen_time`:
# separate the protocol part, which a whole stratum shares, from the individual
# part, which is already blended or can be resampled.
#
#   dose        the mg/kg multiplier is protocol; the patient's weight is
#               individual, and AVATAR already blends it. So keep the multiplier
#               and recompute the amount from the avatar's own weight.
#   attendance  the visit grid is protocol; which visits were attended is
#               individual. So sample the pattern from ones several subjects
#               share instead of copying the anchor's.
#
# `strata` is the stratum for both. It is deliberately NOT a blending
# barrier -- only route is -- so these tests also pin that donors still cross it.

da_roles <- function(..., properties = NULL) {
  pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
            cmt = "CMT", covariates = "WT", strata = properties,
            ...)
}

# Exactly `mg_per_kg` mg/kg, so any departure in the output is the generator's.
da_source <- function(n = 20L, mg_per_kg = 5, arm = NULL, seed = 1) {
  set.seed(seed)
  weight <- round(stats::runif(n, 55, 95), 1)
  times <- c(0.5, 1, 2, 4, 8)
  do.call(rbind, lapply(seq_len(n), function(i) {
    row <- data.frame(
      ID = as.character(i), TIME = c(0, times),
      DV = c(0, 5 * exp(-0.2 * times)),
      AMT = c(mg_per_kg * weight[i], rep(0, length(times))),
      EVID = c(1L, rep(0L, length(times))),
      CMT = c(1L, rep(2L, length(times))),
      WT = weight[i], stringsAsFactors = FALSE
    )
    if (!is.null(arm)) row$ARM <- arm
    row
  }))
}

test_that("weight-based dosing is detected and the avatar's own weight is used", {
  source <- da_source()
  roles <- da_roles()
  expect_true(all(abs(
    source$AMT[source$EVID != 0] / source$WT[source$EVID != 0] - 5
  ) < 1e-9))

  synthetic <- suppressWarnings(
    synpmx_avatar(source, roles, n_subjects = 20, seed = 7)
  )
  dosed <- synthetic[synthetic$EVID != 0, , drop = FALSE]

  expect_identical(attr(synthetic, "pmx_settings")$dose_basis, "WT")
  # The coherence half: every avatar obeys the protocol it claims to follow.
  # Before this, AMT came from the anchor while WT came from the donor blend, so
  # avatars landed anywhere from 4.4 to 5.3 mg/kg.
  expect_equal(dosed$AMT / dosed$WT, rep(5, nrow(dosed)), tolerance = 1e-8)
  # The privacy half: no avatar carries a real patient's real dose.
  expect_equal(sum(round(dosed$AMT, 6) %in%
                     round(source$AMT[source$EVID != 0], 6)), 0L)
})

test_that("dosing that is not covariate-proportional is left alone", {
  source <- da_source()
  # Break proportionality: amounts now vary independently of weight.
  set.seed(2)
  dose_rows <- which(source$EVID != 0)
  source$AMT[dose_rows] <- source$AMT[dose_rows] *
    stats::runif(length(dose_rows), 0.6, 1.6)
  synthetic <- suppressWarnings(
    synpmx_avatar(source, da_roles(), n_subjects = 10, seed = 3)
  )
  # Detection must fail closed: recomputing amounts on a study that is not
  # dose-proportional would invent a relationship the protocol never had.
  expect_true(is.na(attr(synthetic, "pmx_settings")$dose_basis))
})

test_that("flat dosing needs no basis and gets none", {
  source <- da_source()
  source$AMT[source$EVID != 0] <- 100
  synthetic <- suppressWarnings(
    synpmx_avatar(source, da_roles(), n_subjects = 10, seed = 3)
  )
  expect_true(is.na(attr(synthetic, "pmx_settings")$dose_basis))
  expect_true(all(synthetic$AMT[synthetic$EVID != 0] == 100))
})

test_that("several dose levels are recognised with or without a declared arm", {
  # Ratios are clustered, not averaged within a declared stratum, so the levels
  # are found from the numbers themselves. Declaring the arm is then a
  # convenience rather than a requirement.
  source <- rbind(da_source(n = 10L, mg_per_kg = 1, arm = "a", seed = 1),
                  da_source(n = 10L, mg_per_kg = 2, arm = "b", seed = 2),
                  da_source(n = 10L, mg_per_kg = 3, arm = "c", seed = 3))
  source$ID <- paste0(source$ARM, "_", source$ID)

  for (roles in list(da_roles(keep = "ARM"), da_roles(properties = "ARM"))) {
    synthetic <- suppressWarnings(
      synpmx_avatar(source, roles, n_subjects = 30, seed = 5)
    )
    settings <- attr(synthetic, "pmx_settings")
    expect_identical(settings$dose_basis, "WT")
    expect_equal(sort(settings$dose_levels), c(1, 2, 3), tolerance = 1e-8)
    dosed <- synthetic[synthetic$EVID != 0, , drop = FALSE]
    expect_equal(dosed$AMT / dosed$WT,
                 c(a = 1, b = 2, c = 3)[dosed$ARM],
                 tolerance = 1e-8, ignore_attr = TRUE)
  }
})

test_that("intra-patient escalation is recognised", {
  # The case a stratum cannot reach: the level changes *within* subject, so any
  # subject-constant grouping sees a non-constant ratio and fails closed. A
  # ratio does not care which subject or occasion it came from.
  set.seed(4)
  weight <- round(stats::runif(12, 55, 95), 1)
  ladder <- c(1, 2, 4)
  source <- do.call(rbind, lapply(seq_along(weight), function(i) {
    do.call(rbind, lapply(seq_along(ladder), function(step) {
      start <- (step - 1) * 24
      times <- start + c(1, 2, 4, 8)
      data.frame(
        ID = as.character(i), TIME = c(start, times),
        DV = c(0, 5 * exp(-0.2 * (times - start))),
        AMT = c(ladder[step] * weight[i], rep(0, length(times))),
        EVID = c(1L, rep(0L, length(times))),
        CMT = c(1L, rep(2L, length(times))),
        WT = weight[i], stringsAsFactors = FALSE
      )
    }))
  }))
  roles <- da_roles()
  synthetic <- suppressWarnings(
    synpmx_avatar(source, roles, n_subjects = 12, seed = 8)
  )
  settings <- attr(synthetic, "pmx_settings")
  expect_identical(settings$dose_basis, "WT")
  expect_equal(sort(settings$dose_levels), ladder, tolerance = 1e-8)

  dosed <- synthetic[synthetic$EVID != 0, , drop = FALSE]
  # Every avatar climbs the same ladder against its own blended weight, and each
  # subject keeps three ascending doses.
  expect_true(all(round(dosed$AMT / dosed$WT, 6) %in% ladder))
  by_subject <- split(dosed$AMT / dosed$WT, as.character(dosed$ID))
  expect_true(all(vapply(by_subject, function(v) {
    identical(round(sort(v), 6), ladder)
  }, logical(1))))
  expect_equal(sum(round(dosed$AMT, 6) %in%
                     round(source$AMT[source$EVID != 0], 6)), 0L)
})

test_that("strata is a stratum, never a blending barrier", {
  # Route is the only absolute barrier. A stratum with too few subjects to reach
  # the donor floor must still borrow across strata rather than fall back to a
  # near-copy of one real patient.
  source <- rbind(da_source(n = 12L, arm = "a", seed = 1),
                  da_source(n = 2L, arm = "b", seed = 2))
  source$ID <- paste0(source$ARM, "_", source$ID)
  synthetic <- suppressWarnings(synpmx_avatar(
    source, da_roles(properties = "ARM"), n_subjects = 14, seed = 6
  ))
  expect_true(validate_pmx(synthetic, da_roles(properties = "ARM"))$valid)
  # Every avatar is still blended to the floor, arm "b"'s two subjects included.
  expect_gte(attr(synthetic, "pmx_settings")$min_effective_donors, 1)
  expect_true(all(c("a", "b") %in% synthetic$ARM))
})

# Attendance ------------------------------------------------------------------

# Ten subjects on one schedule, plus two who each miss a different visit. The
# full pattern is held by ten, so it qualifies at any floor up to ten; the two
# partial patterns are held by one subject each and must never be reproduced at
# the default floor of 2.
at_source <- function() {
  times <- c(0.5, 1, 2, 4, 8)
  subject <- function(id, drop = integer()) {
    kept <- setdiff(times, times[drop])
    data.frame(
      ID = as.character(id), TIME = c(0, kept),
      DV = c(0, 5 * exp(-0.2 * kept)),
      AMT = c(100, rep(0, length(kept))),
      EVID = c(1L, rep(0L, length(kept))),
      CMT = c(1L, rep(2L, length(kept))),
      WT = 70 + id, stringsAsFactors = FALSE
    )
  }
  rbind(do.call(rbind, lapply(1:10, subject)),
        subject(11, drop = 3), subject(12, drop = c(4, 5)))
}

test_that("a pattern held by one subject is never reproduced", {
  source <- at_source()
  roles <- da_roles()
  lone <- vapply(c("11", "12"), function(id) {
    rows <- source$ID == id & source$EVID == 0
    paste(sort(source$TIME[rows]), collapse = ",")
  }, character(1))

  synthetic <- suppressWarnings(
    synpmx_avatar(source, roles, n_subjects = 30, seed = 11)
  )
  produced <- tapply(synthetic$TIME[synthetic$EVID == 0],
                     as.character(synthetic$ID[synthetic$EVID == 0]),
                     function(v) length(v))
  # Only the ten-subject pattern qualifies, so every avatar has its visit count.
  expect_true(all(produced == 5L))
  expect_false(any(vapply(lone, function(k) {
    length(strsplit(k, ",")[[1]]) == 4L && all(produced == 4L)
  }, logical(1))))
  expect_true(validate_pmx(synthetic, roles)$valid)
})

test_that("min_pattern_share = 1 restores anchor copying exactly", {
  source <- at_source()
  roles <- da_roles()
  copied <- suppressWarnings(
    synpmx_avatar(source, roles, n_subjects = 30, seed = 11,
                  min_pattern_share = 1)
  )
  settings <- attr(copied, "pmx_settings")
  expect_equal(settings$pattern_sampled_fraction, 0)
  # With copying restored, the rare four-visit patterns come through again --
  # which is the behaviour the default exists to prevent.
  counts <- tapply(copied$TIME[copied$EVID == 0],
                   as.character(copied$ID[copied$EVID == 0]), length)
  expect_true(any(counts == 4L))
})

test_that("dose rows are untouched by attendance sampling", {
  source <- at_source()
  synthetic <- suppressWarnings(
    synpmx_avatar(source, da_roles(), n_subjects = 20, seed = 12)
  )
  # Sampling replaces observation rows only. A regimen no protocol permits is
  # the objection that ruled out randomized dose dropping, so the dose schedule
  # stays exactly as the anchor had it.
  doses <- tapply(synthetic$AMT[synthetic$EVID != 0],
                  as.character(synthetic$ID[synthetic$EVID != 0]), length)
  expect_true(all(doses == 1L))
  expect_true(all(synthetic$AMT[synthetic$EVID != 0] == 100))
})

# Visit sets that no patient shares -----------------------------------------
#
# Three outcomes, and the fixtures differ only in how much ROOM there is on the
# grid. `.place_attendance()` may not hand an avatar a visit set that a real
# patient holds, so what matters is whether any legal arrangement is left over.

visit_sets <- function(d) {
  observed <- d$EVID == 0
  vapply(split(d$TIME[observed], as.character(d$ID[observed])),
         function(x) paste(sort(x), collapse = ","), character(1))
}

# n patients, each discontinuing after a different number of visits on a grid of
# `visits` points. Every exact pattern has one holder and so does every
# (kind, count) shape.
staggered_source <- function(n, visits) {
  grid <- c(0, seq_len(visits - 1L) * 24)
  do.call(rbind, lapply(seq_len(n), function(i) {
    kept <- grid[seq_len(min(i + 1L, visits))]
    data.frame(
      ID = as.character(i), TIME = c(0, kept),
      DV = c(0, 5 * exp(-0.02 * kept)),
      AMT = c(100, rep(0, length(kept))),
      EVID = c(1L, rep(0L, length(kept))),
      CMT = c(1L, rep(2L, length(kept))),
      WT = 70 + i, stringsAsFactors = FALSE
    )
  }))
}

test_that("staggered discontinuation is rescued by the kind-only fallback", {
  # SIM-039. No exact pattern is shared AND no (kind, count) shape is either,
  # because every dropout depth has one holder. The group used to get no pool at
  # all and every avatar kept its anchor's own visit set -- the one thing the
  # mechanism exists to prevent. "These patients dropped out" is shared by all
  # of them, and with room on the grid that is enough to build from.
  source <- staggered_source(21L, 30L)
  synthetic <- suppressWarnings(suppressMessages(
    synpmx_avatar(source, da_roles(), n_subjects = 21, seed = 13)
  ))
  settings <- attr(synthetic, "pmx_settings")

  expect_equal(settings$pattern_sampled_fraction, 1)
  expect_equal(sum(visit_sets(synthetic) %in% visit_sets(source)), 0L)
  # Reaching a legal placement required moving the number of missed visits: a
  # discontinuation at a given depth has exactly one possible placement, so
  # without that the draw fails however many times it is retried.
  expect_gt(settings$pattern_shifted_fraction, 0)
})

test_that("a saturated grid alerts instead of quietly copying the anchor", {
  # The union here is exactly six cells and the six patients occupy every depth
  # over it -- one complete attender and trailing 1 through 5 -- so no legal
  # arrangement is left at any depth AND no set is shared widely enough to
  # substitute. Emitting the anchor's own set is then the only thing left, and
  # since that set is unique to one real patient it must be loud.
  grid <- c(0.5, 1, 2, 4, 8, 12, 24)
  source <- do.call(rbind, lapply(1:6, function(i) {
    kept <- utils::head(grid, length(grid) - i)
    data.frame(
      ID = as.character(i), TIME = c(0, kept),
      DV = c(0, 5 * exp(-0.2 * kept)),
      AMT = c(100, rep(0, length(kept))),
      EVID = c(1L, rep(0L, length(kept))),
      CMT = c(1L, rep(2L, length(kept))),
      WT = 70 + i, stringsAsFactors = FALSE
    )
  }))
  raised <- raised_alerts(
    synpmx_avatar(source, da_roles(), n_subjects = 6, seed = 13)
  )
  expect_true(any(grepl("carrying a visit set nobody else shares", raised,
                        fixed = TRUE)))
})

test_that("a rare anchor set is swapped for a shared one rather than emitted", {
  # Ten patients on the complete schedule and one who missed a middle visit.
  # That one patient's set is theirs alone, so an avatar anchored on them must
  # not come out wearing it -- and here there IS a widely held set to use
  # instead, so nothing needs to alert and nothing identifying is emitted.
  grid <- c(0.5, 1, 2, 4, 8, 12, 24)
  source <- do.call(rbind, lapply(1:11, function(i) {
    kept <- if (i == 11L) grid[-4L] else grid
    data.frame(
      ID = as.character(i), TIME = c(0, kept),
      DV = c(0, 5 * exp(-0.2 * kept)),
      AMT = c(100, rep(0, length(kept))),
      EVID = c(1L, rep(0L, length(kept))),
      CMT = c(1L, rep(2L, length(kept))),
      WT = 70 + i, stringsAsFactors = FALSE
    )
  }))
  synthetic <- suppressWarnings(suppressMessages(
    synpmx_avatar(source, da_roles(), n_subjects = 30, seed = 4)
  ))
  settings <- attr(synthetic, "pmx_settings")

  # The disclosure count is the one that must be zero.
  expect_equal(settings$pattern_identifying_fraction, 0)
  # And no avatar holds the lone patient's set.
  rare_set <- paste(sort(c(0.5, 1, 2, 8, 12, 24)), collapse = ",")
  observed <- synthetic$EVID == 0
  produced <- vapply(split(synthetic$TIME[observed],
                           as.character(synthetic$ID[observed])),
                     function(x) paste(sort(x), collapse = ","), character(1))
  expect_false(any(produced == rare_set))
})

test_that("a stratum with nothing shared at all still alerts rather than copying", {
  # The residue the kind-only fallback cannot reach: two patients, one complete
  # and one who missed a middle visit. The kinds are "complete" and "block",
  # one holder each, so not even the coarsest abstraction clears a floor of 2
  # and the group gets no pool.
  grid <- c(0.5, 1, 2, 4, 8, 12, 24)
  source <- do.call(rbind, lapply(1:2, function(i) {
    kept <- if (i == 1L) grid else grid[-3L]
    data.frame(
      ID = as.character(i), TIME = c(0, kept),
      DV = c(0, 5 * exp(-0.2 * kept)),
      AMT = c(100, rep(0, length(kept))),
      EVID = c(1L, rep(0L, length(kept))),
      CMT = c(1L, rep(2L, length(kept))),
      WT = 70 + i, stringsAsFactors = FALSE
    )
  }))
  raised <- raised_alerts(
    synpmx_avatar(source, da_roles(), n_subjects = 4, seed = 13)
  )
  expect_true(any(grepl("no visit set is shared by 2 or more patients",
                        raised, fixed = TRUE)))
})

# Six subjects, each missing exactly one *interior* visit, so every exact
# pattern has a single holder while all six share the shape "one visit missed,
# not at the end".
#
# The grid is ten visits and not six on purpose. The shape fallback can only
# rescue a pattern if there is somewhere to put the miss that no real subject
# already occupies: `.place_attendance()` rejects any placement reproducing a
# pattern too rare to reuse. With six visits and six subjects every placement is
# somebody's, the draw fails, and each avatar keeps its anchor's own visit set
# -- correct behaviour, but it demonstrates the opposite of what these two tests
# are about. Ten visits leaves free placements.
one_miss_grid <- function() c(0.5, 1, 2, 4, 6, 8, 10, 12, 16, 24)

one_miss_source <- function() {
  grid <- one_miss_grid()
  do.call(rbind, lapply(1:6, function(i) {
    kept <- grid[-(i + 1L)]
    data.frame(
      ID = as.character(i), TIME = c(0, kept),
      DV = c(0, 5 * exp(-0.2 * kept)),
      AMT = c(100, rep(0, length(kept))),
      EVID = c(1L, rep(0L, length(kept))),
      CMT = c(1L, rep(2L, length(kept))),
      WT = 70 + i, stringsAsFactors = FALSE
    )
  }))
}

test_that("a shape rescues patterns no single subject shares", {
  # Six subjects each miss exactly one visit, but a different one, so every
  # *exact* pattern is held by one subject and the strict rule discards all six.
  # As shapes they are one group of six: "one visit missed, mid-study".
  grid <- one_miss_grid()
  source <- one_miss_source()
  synthetic <- suppressWarnings(
    synpmx_avatar(source, da_roles(), n_subjects = 12, seed = 21)
  )
  settings <- attr(synthetic, "pmx_settings")
  # Nothing is discarded, where the exact-cell rule alone would have lost all six.
  expect_equal(settings$patterns_dropped, 0L)
  expect_equal(settings$pattern_sampled_fraction, 1)
  # Each avatar still misses exactly one visit: the shape is preserved even
  # though which visit was missed is not.
  counts <- tapply(synthetic$TIME[synthetic$EVID == 0],
                   as.character(synthetic$ID[synthetic$EVID == 0]), length)
  expect_true(all(counts == length(grid) - 1L))
  expect_true(validate_pmx(synthetic, da_roles())$valid)
})

test_that("a generated placement never reproduces a pattern too rare to reuse", {
  # The guarantee the rejection check defends. Every exact pattern here is held
  # by one subject, so every one of them is off limits; the shape they share is
  # not. No avatar may come out holding any source subject's exact visit set.
  grid <- one_miss_grid()
  source <- one_miss_source()
  roles <- da_roles()
  coarsened <- synpmx:::.coarsen_source_time(source, roles)$source
  rare <- vapply(unique(as.character(coarsened$ID)), function(id) {
    paste(sort(coarsened$TIME[coarsened$ID == id & coarsened$EVID == 0]),
          collapse = ",")
  }, character(1))

  synthetic <- suppressWarnings(
    synpmx_avatar(source, roles, n_subjects = 40, seed = 22)
  )
  produced <- vapply(unique(as.character(synthetic$ID)), function(id) {
    paste(sort(synthetic$TIME[synthetic$ID == id & synthetic$EVID == 0]),
          collapse = ",")
  }, character(1))
  # Times are re-noised after placement, so an exact string match is not the
  # test; the visit *set* is. Compare on the grid the placement chose.
  expect_gt(attr(synthetic, "pmx_settings")$pattern_generated_fraction, 0)
  expect_true(validate_pmx(synthetic, roles)$valid)
  expect_equal(length(produced), 40L)
})

test_that("the shape is drawn before the arrangement, so the mix is preserved", {
  # Ten complete attenders and two who dropped out. The cohort is mostly
  # complete, and sampling must keep it that way rather than flattening the mix.
  grid <- c(0.5, 1, 2, 4, 8)
  subject <- function(id, keep_n) {
    kept <- utils::head(grid, keep_n)
    data.frame(
      ID = as.character(id), TIME = c(0, kept),
      DV = c(0, 5 * exp(-0.2 * kept)),
      AMT = c(100, rep(0, length(kept))),
      EVID = c(1L, rep(0L, length(kept))),
      CMT = c(1L, rep(2L, length(kept))),
      WT = 70 + id, stringsAsFactors = FALSE
    )
  }
  source <- rbind(do.call(rbind, lapply(1:10, subject, keep_n = 5)),
                  subject(11, keep_n = 3), subject(12, keep_n = 3))
  synthetic <- suppressWarnings(
    synpmx_avatar(source, da_roles(), n_subjects = 60, seed = 23)
  )
  counts <- tapply(synthetic$TIME[synthetic$EVID == 0],
                   as.character(synthetic$ID[synthetic$EVID == 0]), length)
  complete <- mean(counts == 5L)
  # Ten of twelve source subjects are complete attenders; the synthetic cohort
  # should be near that, not uniform across shapes.
  expect_gt(complete, 0.6)
  expect_lt(complete, 0.95)
})

test_that("the cost of the floor is counted and reported", {
  source <- at_source()
  synthetic <- suppressWarnings(
    synpmx_avatar(source, da_roles(), n_subjects = 20, seed = 14)
  )
  settings <- attr(synthetic, "pmx_settings")
  # Subjects 11 and 12 each hold a pattern nobody else does, so both are
  # excluded and the run must say so rather than losing them silently.
  expect_equal(settings$patterns_dropped, 2L)
  expect_equal(settings$subjects_with_dropped_pattern, 2L)

  raised <- raised_alerts(
    synpmx_avatar(source, da_roles(), n_subjects = 20, seed = 14)
  )
  # The alert has to name the count and say plainly what was removed. It is
  # written for a reader who has never heard the phrase "attendance pattern",
  # so the assertion pins the meaning rather than any one turn of phrase.
  expect_true(any(grepl(
    "2 of 3 distinct visit sets, held by 2 patients", raised, fixed = TRUE)))
  expect_true(any(grepl("are given to no avatar", raised, fixed = TRUE)))

  # Nothing is dropped when the mechanism is off.
  off <- suppressWarnings(
    synpmx_avatar(source, da_roles(), n_subjects = 20, seed = 14,
                  min_pattern_share = 1)
  )
  expect_equal(attr(off, "pmx_settings")$patterns_dropped, 0L)
})

test_that("min_pattern_share is validated", {
  expect_error(
    synpmx_avatar(at_source(), da_roles(), min_pattern_share = 0),
    "`min_pattern_share` must be one integer of 1 or more"
  )
  expect_error(
    synpmx_avatar(at_source(), da_roles(), min_pattern_share = 2.5),
    "`min_pattern_share` must be one integer of 1 or more"
  )
})

# SIM-035 -- attendance sampling cloned one template row per endpoint and
# rewrote only TIME, so every sampled visit inherited the *first* visit's
# nominal time, occasion, and any other row-varying column. The visit metadata
# has to travel with the visit.
test_that("sampled attendance carries each visit's own row metadata", {
  source <- pmx_simulated_fixture(24)
  roles <- pmx_roles(
    id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
    cmt = "CMT", dvid = "DVID", mdv = "MDV", nominal_time = "NTIME",
    tad = "TAD", occasion = "OCC", covariates = "WT"
  )

  synthetic <- suppressWarnings(suppressMessages(
    synpmx_avatar(source, roles, seed = 7)
  ))

  # The fixture is observed exactly on its nominal grid, so the two columns
  # agree row for row in the source and must still agree in the output.
  expect_true(all(abs(source$TIME - source$NTIME) < 1e-9))
  expect_true(all(abs(synthetic$TIME - synthetic$NTIME) < 1e-9))

  # Occasion is not recomputed anywhere, so it is the column that shows a
  # stamped template most plainly: the second dose interval must still be
  # occasion 2.
  second <- synthetic$TIME >= 12
  expect_true(any(second))
  expect_true(all(synthetic$OCC[second] == 2L))
  expect_true(all(synthetic$OCC[!second] == 1L))

  # The defect collapsed the nominal column onto one or two values.
  expect_gt(length(unique(synthetic$NTIME)), 10L)
})

# SIM-036 -- `cmt` says where a dose goes, `dvid` says which endpoint a row is,
# and nothing infers the second from the first. A dataset whose observations sat
# in two compartments with no `dvid` had every row labelled one endpoint, so two
# endpoints read at one visit collapsed to a single attendance cell and one of
# them was rebuilt away. 108 rows in, 60 out, no warning.
ek_source <- function(n = 12L) {
  do.call(rbind, lapply(seq_len(n), function(subject) {
    evid <- c(1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L)
    data.frame(
      ID = subject,
      TIME = c(0, 0, 0, 1, 1, 4, 4, 8, 8),
      DV = c(0, 0, 50, 8, 40, 4, 30, 1, 45) * (1 + 0.1 * sin(subject)),
      AMT = ifelse(evid != 0, 100, 0),
      EVID = evid,
      CMT = c(1L, 2L, 3L, 2L, 3L, 2L, 3L, 2L, 3L),
      WT = 70 + subject,
      stringsAsFactors = FALSE
    )
  }))
}

test_that("multi-compartment observations without dvid are refused", {
  source <- ek_source()
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", cmt = "CMT", covariates = "WT")

  expect_error(synpmx_avatar(source, roles, seed = 1),
               "Observations occupy 2 compartments")

  # One compartment is unambiguous and must still be allowed.
  single <- source[source$CMT != 3, ]
  expect_no_error(
    suppressWarnings(suppressMessages(synpmx_avatar(single, roles, seed = 1)))
  )
})

test_that("one column may be both cmt and dvid, and keeps both endpoints", {
  source <- ek_source()
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", cmt = "CMT", dvid = "CMT",
                     covariates = "WT")

  synthetic <- suppressWarnings(suppressMessages(
    synpmx_avatar(source, roles, seed = 1)
  ))

  # No row is lost and no endpoint disappears.
  expect_equal(nrow(synthetic), nrow(source))
  expect_equal(sort(unique(synthetic$CMT)), sort(unique(source$CMT)))
  expect_equal(as.integer(table(synthetic$CMT)),
               as.integer(table(source$CMT)))

  # Separated, not pooled: each endpoint keeps its own transform and its own
  # scale rather than sharing one fitted across both.
  expect_length(attr(synthetic, "pmx_settings")$endpoint_transforms, 2L)
  observed <- synthetic$EVID == 0
  expect_lt(max(synthetic$DV[observed & synthetic$CMT == 2]),
            min(synthetic$DV[observed & synthetic$CMT == 3]))
})

test_that("a repeated endpoint and time is rebuilt as two rows, not one", {
  # Genuine duplicate records, with dvid declared: the attendance cell for
  # (cp, t = 1) is held twice and must come back twice.
  source <- do.call(rbind, lapply(1:12, function(subject) {
    evid <- c(1L, 0L, 0L, 0L, 0L)
    data.frame(
      ID = subject, TIME = c(0, 1, 1, 4, 8),
      DV = c(0, 8, 8.2, 4, 1) * (1 + 0.1 * sin(subject)),
      AMT = ifelse(evid != 0, 100, 0), EVID = evid,
      CMT = c(1L, 2L, 2L, 2L, 2L),
      DVID = factor("cp", levels = "cp"), WT = 70 + subject,
      stringsAsFactors = FALSE
    )
  }))
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", cmt = "CMT", dvid = "DVID",
                     covariates = "WT")

  synthetic <- suppressWarnings(suppressMessages(
    synpmx_avatar(source, roles, seed = 1)
  ))

  expect_equal(nrow(synthetic), nrow(source))
  per_subject <- table(synthetic$ID[synthetic$TIME == 1 & synthetic$EVID == 0])
  expect_true(all(per_subject == 2L))
})

# An absent `dvid` is legitimate -- a single-endpoint study needs no endpoint
# key -- but it is also how a multi-endpoint dataset gets silently pooled into
# one endpoint, and where `cmt` is not declared nothing can detect that. The
# note is the only defense, so it must fire whenever `dvid` is absent.
test_that("an absent dvid is announced, and named when dvid is declared", {
  source <- ek_source()
  observed <- source$EVID == 0

  # Single endpoint, no dvid: the note fires and says what is assumed.
  single <- source[source$CMT != 3, ]
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", cmt = "CMT", covariates = "WT")
  expect_message(
    suppressWarnings(synpmx_avatar(single, roles, seed = 1)),
    "no `dvid` declared"
  )

  # Repeated subject-and-time observations are reported when present, since
  # they are the visible hint that two endpoints may be in play. They are not
  # refused: replicate measurements are ordinary.
  no_roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                        evid = "EVID", covariates = "WT")
  expect_message(
    suppressWarnings(synpmx_avatar(source, no_roles, seed = 1)),
    "share a subject and time"
  )
  expect_no_error(
    suppressWarnings(suppressMessages(synpmx_avatar(source, no_roles, seed = 1)))
  )

  # Declaring the endpoint key silences it.
  with_key <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                        evid = "EVID", cmt = "CMT", dvid = "CMT",
                        covariates = "WT")
  expect_no_message(
    suppressWarnings(synpmx_avatar(source, with_key, seed = 1))
  )
})

# Declared dose basis --------------------------------------------------------
#
# `.detect_dose_basis()` has to guess that a study is dose-proportional, and it
# guesses conservatively. A study that escalates within a patient and rounds the
# milligrams -- the ordinary oncology shape -- fails that test, correctly, and
# then every avatar's amount is its anchor's milligrams over its own blended
# weight, which is nobody's protocol. `dose_covariate` says it outright.

escalating_source <- function(n = 21L) {
  # Three doses per patient, escalating "mostly 2x" -- not every patient by the
  # same factor -- and dispensed in 25 mg vials. Both of those are ordinary, and
  # together they are what stops the ratios collapsing onto a few levels.
  steps <- list(c(1, 2, 4), c(1, 2, 3), c(1, 1.5, 3), c(1, 2, 4))
  do.call(rbind, lapply(seq_len(n), function(i) {
    wt <- 52 + i * 2.3
    amt <- round(steps[[(i %% 4L) + 1L]] * wt / 25) * 25
    dose_t <- c(0, 336, 1008)
    obs_t <- c(24, 360, 1032)
    out <- rbind(
      data.frame(ID = i, TIME = dose_t, DV = NA_real_, AMT = amt, EVID = 1L,
                 CMT = 1L, WT = wt),
      data.frame(ID = i, TIME = obs_t, DV = round(5 + i / 10, 2), AMT = 0,
                 EVID = 0L, CMT = 2L, WT = wt))
    out[order(out$TIME, -out$EVID), ]
  }))
}

escalating_roles <- function(declare = FALSE) {
  args <- list(id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
               cmt = "CMT", covariates = "WT")
  if (declare) args$dose_covariate <- "WT"
  do.call(pmx_roles, args)
}

test_that("inference declines rounded intra-patient escalation, and says why", {
  synthetic <- suppressWarnings(suppressMessages(
    synpmx_avatar(escalating_source(), escalating_roles(), seed = 2)
  ))
  settings <- attr(synthetic, "pmx_settings")
  expect_true(is.na(settings$dose_basis))
  # The reason has to name the covariate that was tried, or "not detected" and
  # "never attempted" are the same blank in a report.
  expect_match(settings$dose_basis_note, "WT")
})

test_that("a declared dose_covariate holds each dose row's own ratio", {
  source <- escalating_source()
  synthetic <- suppressWarnings(suppressMessages(
    synpmx_avatar(source, escalating_roles(declare = TRUE), seed = 2)
  ))
  settings <- attr(synthetic, "pmx_settings")
  expect_equal(settings$dose_basis, "WT")
  expect_true(settings$dose_basis_declared)

  # Each avatar's mg/kg must match its own anchor's, dose row for dose row.
  # Nothing is snapped to a shared level, which is what lets a study where
  # patients escalate by different factors come through intact.
  ratio <- function(d) {
    dosed <- d[d$EVID == 1L, ]
    lapply(split(round(dosed$AMT / dosed$WT, 2), as.character(dosed$ID)),
           unname)
  }
  source_ratios <- unique(ratio(source))
  for (r in ratio(synthetic)) {
    expect_true(any(vapply(source_ratios,
                           function(s) isTRUE(all.equal(s, r, tolerance = 0.02)),
                           logical(1))))
    expect_true(r[[2L]] > r[[1L]] && r[[3L]] > r[[2L]])   # still escalating
  }
})

test_that("dose_covariate is validated against the roles and the data", {
  # It must be blended, so it must be a covariate.
  expect_error(
    pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
              covariates = "WT", dose_covariate = "BSA"),
    "must also be named in `covariates`"
  )
  # Naming it twice is not a role collision.
  expect_no_error(
    pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
              covariates = "WT", dose_covariate = "WT")
  )
  # A covariate that is missing where a dose is given cannot be divided by.
  source <- escalating_source()
  source$WT[source$EVID == 1L][3L] <- NA_real_
  expect_error(
    suppressMessages(synpmx_avatar(source, escalating_roles(declare = TRUE),
                                   seed = 2)),
    "missing or non-positive on 1 dose row"
  )
})

# SIM-042. The guarantee, asserted on the finished table rather than inferred
# from the mechanisms that built it. Every leak in this area so far leaked
# because each mechanism was correct on its own terms and nobody asked the whole
# question afterwards, so the whole question is asked here.
test_that("no avatar is emitted with a visit set only one real patient has", {
  # One patient measuring a different set of endpoints from everybody else is a
  # schedule group of one, which can never hold a shared visit set. Nothing can
  # be substituted for it, so an avatar anchored there used to wear that
  # patient's exact pattern.
  visits <- c(0, 1, 2, 7, 14, 28, 56, 84)
  mk <- function(i) {
    pk <- if (i == 12L) NULL else
      data.frame(NTIME = visits[c(1:4, 4 + which(i %% c(2, 3, 4) == 0))],
                 YTYPE = 1L, CMT = 2L)
    obs <- rbind(data.frame(NTIME = visits, YTYPE = 2L, CMT = 3L), pk)
    obs <- obs[order(obs$NTIME, obs$YTYPE), ]
    rbind(
      data.frame(ID = i, TIME = 0, NTIME = 0, DV = NA_real_, AMT = 100,
                 EVID = 1L, CMT = 1L, YTYPE = 0L, WT = 60 + i),
      data.frame(ID = i, TIME = obs$NTIME, NTIME = obs$NTIME,
                 DV = round(5 + i / 10, 2), AMT = 0, EVID = 0L,
                 CMT = obs$CMT, YTYPE = obs$YTYPE, WT = 60 + i))
  }
  source <- do.call(rbind, lapply(1:12, mk))
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", cmt = "CMT", dvid = "YTYPE",
                     nominal_time = "NTIME", covariates = "WT")

  synthetic <- suppressWarnings(suppressMessages(
    synpmx_avatar(source, roles, n_subjects = 24, seed = 5)
  ))
  settings <- attr(synthetic, "pmx_settings")

  # The number the guarantee is stated in, measured on the output.
  expect_equal(settings$identifying_visit_sets, 0L)
  # Independently of the settings, recomputed here from the tables themselves.
  cells <- function(d) {
    observed <- d$EVID == 0
    vapply(split(paste0(d$YTYPE[observed], "@", d$TIME[observed]),
                 as.character(d$ID[observed])),
           function(x) paste(sort(unique(x)), collapse = ";"), character(1))
  }
  source_sets <- cells(source)
  rare <- names(which(table(source_sets) < 2L))
  expect_equal(sum(cells(synthetic) %in% source_sets[source_sets %in% rare]), 0L)
  # The lone-endpoint-set patient is still a donor and still in the cohort; it
  # is the avatar that moved, not the patient that was removed.
  expect_equal(settings$anchors_available, 12L)
})

# SIM-043. The dose side of the same guarantee. Dose events are copied from the
# anchor verbatim -- resampling them would emit regimens the protocol never
# permitted -- so an avatar's dose schedule simply IS its anchor's, and the only
# lever is which patients may anchor.
test_that("no avatar inherits a dose schedule only one real patient has", {
  # Ten patients on the protocol schedule, two with an extra interim dose that
  # nobody else received.
  mk <- function(i) {
    doses <- if (i > 10L) c(0, 168, 252, 336) else c(0, 168, 336)
    obs <- c(1, 24, 169, 192, 337, 360)
    out <- rbind(
      data.frame(ID = i, TIME = doses, DV = NA_real_, AMT = 100, EVID = 1L,
                 CMT = 1L, WT = 60 + i),
      data.frame(ID = i, TIME = obs, DV = round(5 + i / 10, 2), AMT = 0,
                 EVID = 0L, CMT = 2L, WT = 60 + i))
    out[order(out$TIME, -out$EVID), ]
  }
  # Patient 11 and 12 differ from each other too, so each is alone.
  source <- do.call(rbind, lapply(1:12, mk))
  source$TIME[source$ID == 12L & source$TIME == 252] <- 264
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", cmt = "CMT", covariates = "WT")

  synthetic <- suppressWarnings(suppressMessages(
    synpmx_avatar(source, roles, n_subjects = 40, seed = 3)
  ))
  settings <- attr(synthetic, "pmx_settings")
  expect_equal(settings$identifying_dose_schedules, 0L)

  # Recomputed from the tables, independently of the recorded setting: no
  # avatar holds the dose-time set of a patient who was alone in having it.
  dose_sets <- function(d) {
    dosed <- d$EVID == 1L
    vapply(split(d$TIME[dosed], as.character(d$ID[dosed])),
           function(t) paste(sort(unique(t)), collapse = ","), character(1))
  }
  source_sets <- dose_sets(source)
  alone <- names(which(table(source_sets) < 2L))
  expect_gt(length(alone), 0L)                     # the fixture really is loaded
  expect_equal(sum(dose_sets(synthetic) %in% alone), 0L)
  # Those two patients are not removed from the cohort -- they still anchor
  # nothing, but they still act as donors.
  expect_equal(settings$anchors_available, 12L)
})

# SIM-044. Dosing that stopped early, as a shape. A schedule that is a PREFIX of
# the full one can be re-truncated to a different depth, which is the one dose
# edit that is protocol-valid: the result is a regimen the study permitted.
onc_source <- function() {
  mk <- function(i, nd) {
    dt <- seq(0, by = 336, length.out = 10)[seq_len(nd)]
    out <- rbind(
      data.frame(ID = i, TIME = dt, DV = NA_real_, AMT = 100, EVID = 1L,
                 CMT = 1L, WT = 60 + i),
      data.frame(ID = i, TIME = dt + 24, DV = 5 + i / 10, AMT = 0, EVID = 0L,
                 CMT = 2L, WT = 60 + i))
    out[order(out$TIME, -out$EVID), ]
  }
  # 28 patients complete all ten doses; four stop at depths nobody else used.
  do.call(rbind, c(lapply(1:28, function(i) mk(i, 10)),
                   list(mk(29, 4), mk(30, 6), mk(31, 7), mk(32, 9))))
}

onc_depths <- function(d) as.integer(table(d$ID[d$EVID == 1L]))

test_that("early discontinuation survives, at a depth nobody used", {
  source <- onc_source()
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", cmt = "CMT", covariates = "WT")
  synthetic <- suppressWarnings(suppressMessages(
    synpmx_avatar(source, roles, n_subjects = 32, seed = 5)
  ))
  settings <- attr(synthetic, "pmx_settings")

  # The guarantee holds.
  expect_equal(settings$identifying_dose_schedules, 0L)
  held_once <- as.integer(names(which(table(onc_depths(source)) == 1L)))
  expect_gt(length(held_once), 0L)
  expect_equal(sum(onc_depths(synthetic) %in% held_once), 0L)

  # And early stopping is still represented, rather than every avatar being
  # pushed onto the full schedule -- which is what re-anchoring alone did.
  expect_true(any(onc_depths(synthetic) < 10L))
  expect_gt(settings$dose_truncated_fraction, 0)
})

test_that("a truncation target is never deeper than the schedule it replaces", {
  # Truncation drops dose rows and cannot invent one, so a deeper target would
  # silently leave the schedule where it was -- and put the avatar back on the
  # depth exactly one patient had used.
  plan <- synpmx:::.dose_truncation_plan(
    onc_source(),
    pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
              cmt = "CMT", covariates = "WT"),
    2L
  )
  moved <- !is.na(plan$target)
  expect_true(any(moved))
  expect_true(all(plan$target[moved] < plan$depth[moved]))
})
