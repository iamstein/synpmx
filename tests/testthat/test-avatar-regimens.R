# Structural invariants for the regimen shapes a real study actually has:
# more than one drug, more than one dosing compartment, and dosing patterns that
# differ between arms. None of these are exercised by the public datasets, which
# are all single-drug and single-regimen, so they are pinned here.
#
# The rule they all test is the same one: M3 redraws which *DV observations* an
# avatar has and never touches dose events, so every avatar's regimen must be one
# some real patient actually received. A synthetic regimen nobody was on would be
# a protocol violation invented by the generator.

rg_roles <- function(...) {
  pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
            cmt = "CMT", covariates = "WT", ...)
}

# Dose-gap signature per subject, in days, as a single string.
rg_gaps <- function(data, id) {
  times <- sort(data$TIME[data$ID == id & data$EVID != 0])
  paste(round(diff(times) / 24), collapse = ",")
}

rg_subject <- function(id, dose_times, cmt = 1L, amt = 100, rate = 0) {
  data.frame(TIME = dose_times, AMT = amt, EVID = 1L, CMT = cmt, RATE = rate,
             stringsAsFactors = FALSE)
}

test_that("two drugs into two compartments keep their dose structure", {
  source <- do.call(rbind, lapply(1:24, function(i) {
    doses <- rbind(
      rg_subject(i, c(0, 24, 48), cmt = 1L, amt = 100),          # oral drug
      rg_subject(i, 0, cmt = 3L, amt = 50, rate = 25)            # IV infusion
    )
    obs <- data.frame(
      TIME = c(1, 4, 12, 25, 28, 49, 52, 1, 6, 24, 48), AMT = 0, EVID = 0L,
      CMT = c(rep(2L, 7), rep(4L, 4)), RATE = 0, stringsAsFactors = FALSE)
    out <- rbind(doses, obs)
    out$ID <- i
    out$DV <- ifelse(out$EVID != 0, 0,
                     ifelse(out$CMT == 2L,
                            8 * exp(-0.1 * (out$TIME %% 24)) * (1 + 0.1 * sin(i)),
                            40 * exp(-0.02 * out$TIME) * (1 + 0.1 * cos(i))))
    out$YTYPE <- ifelse(out$EVID != 0, 0L, ifelse(out$CMT == 2L, 1L, 2L))
    out$WT <- 70 + i
    out[order(out$TIME, -out$EVID), ]
  }))
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     rate = "RATE", evid = "EVID", cmt = "CMT", dvid = "YTYPE",
                     covariates = "WT")

  synthetic <- suppressWarnings(suppressMessages(
    synpmx_avatar(source, roles, seed = 11)
  ))

  expect_true(validate_pmx(synthetic, roles)$valid)
  # Both dosing compartments survive, in the same proportion.
  expect_equal(as.integer(table(synthetic$CMT[synthetic$EVID != 0])),
               as.integer(table(source$CMT[source$EVID != 0])))
})

test_that("a route arm below the donor floor is dropped, loudly", {
  # 20 patients on one drug, 4 on both. The four-patient combination is its own
  # route, cannot reach k + 1 donors, and by default is dropped -- so the second
  # drug disappears from the output entirely. Pinned because it is a surprising
  # consequence of a deliberate rule, not because it is desirable.
  mk <- function(i, both) {
    doses <- rg_subject(i, c(0, 24), cmt = 1L)
    if (both) doses <- rbind(doses, rg_subject(i, 0, cmt = 3L, amt = 50, rate = 25))
    obs <- data.frame(TIME = c(1, 4, 12, 25, 30), AMT = 0, EVID = 0L, CMT = 2L,
                      RATE = 0, stringsAsFactors = FALSE)
    out <- rbind(doses, obs)
    out$ID <- i
    out$DV <- ifelse(out$EVID != 0, 0, 8 * exp(-0.08 * out$TIME) * (1 + 0.1 * sin(i)))
    out$WT <- 70 + i
    out[order(out$TIME, -out$EVID), ]
  }
  source <- do.call(rbind, c(lapply(1:20, mk, both = FALSE),
                             lapply(21:24, mk, both = TRUE)))
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     rate = "RATE", evid = "EVID", cmt = "CMT",
                     covariates = "WT")

  raised <- character()
  synthetic <- withCallingHandlers(
    suppressMessages(synpmx_avatar(source, roles, seed = 13)),
    warning = function(w) {
      raised <<- c(raised, conditionMessage(w))
      invokeRestart("muffleWarning")
    })

  expect_equal(attr(synthetic, "pmx_settings")$routes, 2L)
  expect_equal(attr(synthetic, "pmx_settings")$anchors_route_excluded, 4L)
  expect_equal(sum(synthetic$CMT == 3L), 0L)          # the second drug is gone
  expect_true(any(grepl("below the donor floor", raised)))
})

test_that("mixed regimens stay coherent and no avatar gets an invented one", {
  # Continuous QD against QD with three weeks on and one week off -- the shape a
  # real oncology study has, and the one where an invented regimen would be a
  # protocol violation rather than a cosmetic flaw.
  continuous <- seq(0, 24 * 55, by = 24)
  holiday <- c(seq(0, 24 * 20, by = 24), seq(24 * 28, 24 * 48, by = 24))
  mk <- function(i, dose_times) {
    obs <- sort(unique(c(dose_times + 2, seq(0, 24 * 55, by = 24 * 7) + 4)))
    out <- rbind(
      data.frame(TIME = dose_times, AMT = 100, EVID = 1L, CMT = 1L),
      data.frame(TIME = obs, AMT = 0, EVID = 0L, CMT = 2L))
    out$ID <- i
    dosing <- vapply(out$TIME, function(t) {
      any(t - dose_times >= 0 & t - dose_times < 24)
    }, TRUE)
    out$DV <- ifelse(out$EVID != 0, 0,
                     ifelse(dosing, 6, 1.2) * (1 + 0.1 * sin(i)))
    out$WT <- 70 + i
    out[order(out$TIME, -out$EVID), ]
  }
  source <- do.call(rbind, c(lapply(1:15, mk, dose_times = continuous),
                             lapply(16:30, mk, dose_times = holiday)))
  roles <- rg_roles()

  synthetic <- suppressWarnings(suppressMessages(
    synpmx_avatar(source, roles, seed = 21)
  ))

  expect_true(validate_pmx(synthetic, roles)$valid)

  # The heart of it: every avatar's dose-gap sequence is one a real patient had.
  # A "hybrid" -- three weeks on, then continuous -- would be a regimen the
  # protocol never permitted, assembled by the generator.
  source_gaps <- unique(vapply(unique(source$ID), rg_gaps, character(1),
                               data = source))
  synthetic_gaps <- vapply(unique(synthetic$ID), rg_gaps, character(1),
                           data = synthetic)
  expect_true(all(synthetic_gaps %in% source_gaps))
  expect_setequal(unique(synthetic_gaps), source_gaps)

  # Dose counts are the source's two values, never something in between.
  counts <- as.integer(table(synthetic$ID[synthetic$EVID != 0]))
  expect_true(all(counts %in% c(42L, 56L)))

  # The drug holiday is visible in the output: during the off week, avatars on
  # the interrupted regimen sit far below the continuously dosed ones. If donor
  # borrowing across regimens washed this out, the mode would be useless for the
  # studies it matters in.
  interrupted <- vapply(unique(synthetic$ID), function(id) {
    grepl("8", rg_gaps(synthetic, id), fixed = TRUE)
  }, TRUE)
  off_week <- synthetic$EVID == 0 &
    synthetic$TIME >= 24 * 21 & synthetic$TIME <= 24 * 27
  on_holiday <- synthetic$ID %in% unique(synthetic$ID)[interrupted]
  expect_lt(stats::median(synthetic$DV[off_week & on_holiday]),
            0.5 * stats::median(synthetic$DV[off_week & !on_holiday]))
})

test_that("two drugs with different PK profiles keep their own shapes", {
  # Drug A oral into CMT 1, fast, with an absorption peak around 2 h.
  # Drug B IV into CMT 3, slow, monotone decline over three days.
  # Declaring `dvid` is what earns each endpoint its own transform; without it
  # both would be pooled into one and the shapes would blur together.
  source <- do.call(rbind, lapply(1:24, function(i) {
    doses <- rbind(
      data.frame(TIME = c(0, 24, 48), AMT = 100, EVID = 1L, CMT = 1L, RATE = 0),
      data.frame(TIME = 0, AMT = 50, EVID = 1L, CMT = 3L, RATE = 25))
    obs <- rbind(
      data.frame(TIME = c(0.5, 2, 6, 12, 24.5, 26, 30, 48.5, 50, 54, 72),
                 AMT = 0, EVID = 0L, CMT = 2L, RATE = 0),
      data.frame(TIME = c(1, 6, 24, 48, 72), AMT = 0, EVID = 0L, CMT = 4L,
                 RATE = 0))
    out <- rbind(doses, obs)
    out$ID <- i
    scale <- 1 + 0.12 * sin(i)
    tad <- out$TIME %% 24
    out$DV <- ifelse(
      out$EVID != 0, 0,
      ifelse(out$CMT == 2L,
             10 * (exp(-0.10 * tad) - exp(-1.2 * tad)) * scale + 0.2,
             40 * exp(-0.02 * out$TIME) * scale))
    out$YTYPE <- ifelse(out$EVID != 0, 0L, ifelse(out$CMT == 2L, 1L, 2L))
    out$WT <- 70 + i
    out[order(out$TIME, -out$EVID), ]
  }))
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     rate = "RATE", evid = "EVID", cmt = "CMT", dvid = "YTYPE",
                     covariates = "WT")

  synthetic <- suppressWarnings(suppressMessages(
    synpmx_avatar(source, roles, seed = 31)
  ))

  expect_true(validate_pmx(synthetic, roles)$valid)
  # One transform per endpoint, never one shared across both scales.
  expect_length(attr(synthetic, "pmx_settings")$endpoint_transforms, 2L)

  decline <- function(data, ytype) {
    keep <- data$EVID == 0 & data$YTYPE == ytype
    early <- stats::median(data$DV[keep & data$TIME <= 3])
    late <- stats::median(data$DV[keep & data$TIME >= 48])
    late / early
  }
  # The two drugs decline at very different rates, and each avatar must keep the
  # rate belonging to its own endpoint. Pooling them would drag both toward a
  # common middle.
  expect_equal(decline(synthetic, 1L), decline(source, 1L), tolerance = 0.15)
  expect_equal(decline(synthetic, 2L), decline(source, 2L), tolerance = 0.15)
  expect_gt(decline(synthetic, 1L), 2 * decline(synthetic, 2L))

  # Drug A's absorption peak survives at roughly its source position.
  peak <- function(data) {
    keep <- data$EVID == 0 & data$YTYPE == 1L
    stats::median((data$TIME[keep][order(-data$DV[keep])][1:5]) %% 24)
  }
  expect_lt(abs(peak(synthetic) - peak(source)), 2)
})
