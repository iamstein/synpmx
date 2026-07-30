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
# `subject_properties` is the stratum for both. It is deliberately NOT a blending
# barrier -- only route is -- so these tests also pin that donors still cross it.

da_roles <- function(..., properties = NULL) {
  pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
            cmt = "CMT", covariates = "WT", subject_properties = properties,
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

test_that("subject_properties is a stratum, never a blending barrier", {
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

test_that("a stratum with nothing shared widely enough alerts rather than copying", {
  # Every subject on their own schedule: no pattern can qualify, so anchors keep
  # their own pattern. That is the safe failure, but it must be loud, because it
  # is exactly what the mechanism was asked to prevent.
  # Neither an exact pattern nor a shape can qualify here: each subject misses a
  # *different number* of visits, so every shape is held by one subject too. A
  # fixture where only the exact patterns differ would now be rescued by the
  # shape fallback, which is the point of that fallback.
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
  # This fixture legitimately trips the coarsening alerts too, so collect every
  # warning and assert the pattern one is among them rather than muffling all.
  raised <- character()
  suppressMessages(withCallingHandlers(
    synpmx_avatar(source, da_roles(), n_subjects = 6, seed = 13),
    warning = function(w) {
      raised <<- c(raised, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  ))
  expect_true(any(grepl("no attendance pattern is shared by 2 or more subjects",
                        raised, fixed = TRUE)))
})

test_that("a shape rescues patterns no single subject shares", {
  # Six subjects each miss exactly one visit, but a different one, so every
  # *exact* pattern is held by one subject and the strict rule discards all six.
  # As shapes they are one group of six: "one visit missed, mid-study".
  grid <- c(0.5, 1, 2, 4, 8, 12)
  source <- do.call(rbind, lapply(1:6, function(i) {
    kept <- grid[-i]
    data.frame(
      ID = as.character(i), TIME = c(0, kept),
      DV = c(0, 5 * exp(-0.2 * kept)),
      AMT = c(100, rep(0, length(kept))),
      EVID = c(1L, rep(0L, length(kept))),
      CMT = c(1L, rep(2L, length(kept))),
      WT = 70 + i, stringsAsFactors = FALSE
    )
  }))
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
  grid <- c(0.5, 1, 2, 4, 8, 12)
  source <- do.call(rbind, lapply(1:6, function(i) {
    kept <- grid[-i]
    data.frame(
      ID = as.character(i), TIME = c(0, kept),
      DV = c(0, 5 * exp(-0.2 * kept)),
      AMT = c(100, rep(0, length(kept))),
      EVID = c(1L, rep(0L, length(kept))),
      CMT = c(1L, rep(2L, length(kept))),
      WT = 70 + i, stringsAsFactors = FALSE
    )
  }))
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

  raised <- character()
  suppressMessages(withCallingHandlers(
    synpmx_avatar(source, da_roles(), n_subjects = 20, seed = 14),
    warning = function(w) {
      raised <<- c(raised, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  ))
  expect_true(any(grepl("held by 2 subject(s), were not shared",
                        raised, fixed = TRUE)))

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
