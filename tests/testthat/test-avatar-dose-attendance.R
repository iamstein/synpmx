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

test_that("a declared stratum lets several dose levels be recognised", {
  # Two arms at different mg/kg. Pooled, the ratio is not constant and detection
  # correctly refuses; declared as a stratum, each arm is proportional and both
  # are recomputed against their own multiplier. This is what `subject_properties`
  # buys, and why a multi-level study has to declare its dose group.
  source <- rbind(da_source(n = 12L, mg_per_kg = 5, arm = "low", seed = 1),
                  da_source(n = 12L, mg_per_kg = 9, arm = "high", seed = 2))
  source$ID <- paste0(source$ARM, "_", source$ID)

  pooled <- suppressWarnings(synpmx_avatar(
    source, da_roles(keep = "ARM"), n_subjects = 12, seed = 5
  ))
  expect_true(is.na(attr(pooled, "pmx_settings")$dose_basis))

  stratified <- suppressWarnings(synpmx_avatar(
    source, da_roles(properties = "ARM"), n_subjects = 24, seed = 5
  ))
  expect_identical(attr(stratified, "pmx_settings")$dose_basis, "WT")
  dosed <- stratified[stratified$EVID != 0, , drop = FALSE]
  expected <- ifelse(dosed$ARM == "low", 5, 9)
  expect_equal(dosed$AMT / dosed$WT, expected, tolerance = 1e-8)
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
# partial patterns are held by one subject each and must never be reproduced.
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
  # The times themselves are shared -- coarsening puts everyone on one grid --
  # but each subject attended a different subset of it, so every *pattern* is
  # held by exactly one subject and none can qualify. This is the residual
  # coarsening cannot reach, in its purest form.
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
  expect_true(any(grepl("no attendance pattern is shared by 3 or more subjects",
                        raised, fixed = TRUE)))
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
