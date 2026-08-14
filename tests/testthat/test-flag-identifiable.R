# flag_identifiable_subjects(): a per-subject screen for structurally unusual
# subjects on four axes -- follow-up time, number of doses, dose magnitude, and
# DV value. Each test isolates one axis by holding the other three uniform.

mk <- function(id, obs_times = c(0.5, 1, 2, 4), n_dose = 1L, dose = 100,
               peak = 5) {
  dose_times <- if (n_dose >= 1L) seq(0, by = 12, length.out = n_dose) else
    numeric(0)
  time <- c(dose_times, obs_times)
  data.frame(
    ID = as.character(id),
    TIME = time,
    DV = c(rep(0, length(dose_times)), peak * exp(-0.3 * obs_times)),
    AMT = c(rep(dose, length(dose_times)), rep(0, length(obs_times))),
    EVID = c(rep(1L, length(dose_times)), rep(0L, length(obs_times))),
    CMT = c(rep(1L, length(dose_times)), rep(2L, length(obs_times))),
    stringsAsFactors = FALSE
  )
}

id_roles <- function() {
  pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
            evid = "EVID", cmt = "CMT")
}

flagged_axes <- function(report, id) {
  report$outlier_axes[report$subject_id == id]
}

test_that("an unusually long follow-up is flagged on the time axis", {
  data <- rbind(
    do.call(rbind, lapply(1:5, mk)),
    mk(6, obs_times = c(0.5, 1, 2, 40))     # long follow-up
  )
  report <- flag_identifiable_subjects(data, id_roles())

  expect_true(report$flagged[report$subject_id == "6"])
  expect_match(flagged_axes(report, "6"), "follow-up time")
  expect_equal(sum(report$flagged), 1L)
})

test_that("an unusual number of doses is flagged on the dose-count axis", {
  data <- rbind(
    do.call(rbind, lapply(1:5, mk)),
    mk(6, n_dose = 5L)                       # five doses vs one
  )
  report <- flag_identifiable_subjects(data, id_roles())

  expect_true(report$flagged[report$subject_id == "6"])
  expect_match(flagged_axes(report, "6"), "number of doses")
  expect_equal(sum(report$flagged), 1L)
})

test_that("a rare dose magnitude is flagged on the dose axis", {
  data <- rbind(
    do.call(rbind, lapply(1:5, mk)),
    mk(6, dose = 1000)                       # 10x dose
  )
  report <- flag_identifiable_subjects(data, id_roles())

  expect_true(report$flagged[report$subject_id == "6"])
  expect_match(flagged_axes(report, "6"), "dose magnitude")
  expect_equal(sum(report$flagged), 1L)
})

test_that("an extreme DV peak is flagged on the DV axis", {
  data <- rbind(
    do.call(rbind, lapply(1:5, mk)),
    mk(6, peak = 50)                         # 10x peak
  )
  report <- flag_identifiable_subjects(data, id_roles())

  expect_true(report$flagged[report$subject_id == "6"])
  expect_match(flagged_axes(report, "6"), "DV value")
  expect_equal(sum(report$flagged), 1L)
})

test_that("a homogeneous cohort flags no one", {
  data <- do.call(rbind, lapply(1:6, mk))
  report <- flag_identifiable_subjects(data, id_roles())

  expect_equal(sum(report$flagged), 0L)
  expect_true(all(report$outlier_axes == ""))
})

test_that("the report is sorted most-unusual first and marked restricted", {
  data <- rbind(do.call(rbind, lapply(1:5, mk)), mk(6, peak = 50))
  report <- flag_identifiable_subjects(data, id_roles())

  expect_equal(report$subject_id[1], "6")            # riskiest first
  expect_equal(attr(report, "release_status"), "restricted_not_releasable")
  expect_equal(attr(report, "n_flagged"), sum(report$flagged))
})

test_that("threshold is validated and printing returns invisibly", {
  data <- do.call(rbind, lapply(1:6, mk))
  roles <- id_roles()
  expect_error(flag_identifiable_subjects(data, roles, threshold = 0),
               "threshold")
  report <- flag_identifiable_subjects(data, roles)
  expect_output(out <- withVisible(print(report)), "identifiability check")
  expect_false(out$visible)
})

test_that("a real AVATAR run over a public fixture is screenable end to end", {
  data <- pmx_simulated_fixture(30)
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", cmt = "CMT", dvid = "DVID",
                     covariates = "WT")
  synthetic <- suppressWarnings(synpmx_avatar(data, roles, seed = 7))

  report <- flag_identifiable_subjects(synthetic, roles)
  expect_equal(nrow(report), length(unique(synthetic$ID)))
  expect_true(all(c("subject_id", "follow_up_time", "n_doses", "max_dose",
                    "max_dv", "outlier_axes", "flagged") %in% names(report)))
})

# SIM-041. The median absolute deviation is zero whenever more than half a
# cohort shares one exact value, which on trial data is the ordinary case, not a
# degenerate one: most patients complete the protocol and end at the same visit.
# Treating every departure from the median as infinitely extreme then flagged
# every patient who stopped even slightly early.
test_that("a cohort clustered on one follow-up time does not flag everyone", {
  z <- synpmx:::.modified_z(c(rep(145, 11), 8, 8, 10, 11, 37, 82, 93, 109, 109))
  expect_true(all(is.finite(z)))
  expect_equal(sum(abs(z) > 3.5), 0L)
  # Two patients agreeing to three decimal places must not both be outliers,
  # which is what gave the defect away on a real study.
  expect_equal(z[[length(z)]], z[[length(z) - 1L]])
})

test_that("the zero-MAD fallback still fires on a genuine extreme", {
  # wbcSim's shape: a cohort at ~650 hours with one patient at 4580.
  z <- synpmx:::.modified_z(c(rep(650, 18), 1730, 4580))
  expect_gt(abs(z[[length(z)]]), 3.5)
  # A count axis behaves the same way: three doses for nearly everyone, and a
  # single-dose patient is genuinely distinctive.
  doses <- synpmx:::.modified_z(c(rep(3, 19), 1, 2))
  expect_gt(abs(doses[[20L]]), 3.5)
})

test_that("no variation at all yields no outliers rather than infinities", {
  z <- synpmx:::.modified_z(rep(7, 10))
  expect_true(all(z == 0))
})

# SIM-042. A modified z says a subject is unusual for its cohort; it cannot say
# the difference is big enough to single anybody out. On a study whose subjects
# all completed the protocol, follow-up times sit within an hour or two of each
# other, the median absolute deviation of that is minutes, and a subject forty
# minutes from the median scores past 3.5 while being identifiable to nobody.
# `min_relative_gap` is the second condition: the distance to the NEAREST other
# subject, in units of the group's median.
test_that("an outlier separated by a trivial amount is not flagged", {
  # Everybody finishes at about 216 hours, one subject an hour and a half later.
  finish <- c(215.1, 215.4, 215.6, 215.7, 215.8, 215.9, 216.0, 216.1, 217.5)
  data <- do.call(rbind, Map(function(id, last) {
    mk(id, obs_times = c(0.5, 1, 2, last))
  }, seq_along(finish), finish))

  expect_equal(sum(flag_identifiable_subjects(data, id_roles())$flagged), 0L)
  # The z alone does call it an outlier, which is the behaviour being fixed.
  loose <- flag_identifiable_subjects(data, id_roles(), min_relative_gap = 0)
  expect_true(loose$flagged[loose$subject_id == "9"])
})

test_that("two subjects sharing an extreme value single out neither", {
  alone <- rbind(do.call(rbind, lapply(1:5, mk)),
                 mk(6, obs_times = c(0.5, 1, 2, 40)))
  paired <- rbind(alone, mk(7, obs_times = c(0.5, 1, 2, 40)))

  expect_true(flag_identifiable_subjects(alone, id_roles())$flagged[
    flag_identifiable_subjects(alone, id_roles())$subject_id == "6"
  ])
  # Same value, now held by two: the gap to the nearest other subject is 0.
  expect_equal(sum(flag_identifiable_subjects(paired, id_roles())$flagged), 0L)
})

test_that("min_relative_gap is validated and recorded", {
  data <- do.call(rbind, lapply(1:6, mk))
  expect_error(flag_identifiable_subjects(data, id_roles(),
                                          min_relative_gap = -1),
               "min_relative_gap")
  report <- flag_identifiable_subjects(data, id_roles(),
                                       min_relative_gap = 0.25)
  expect_equal(attr(report, "min_relative_gap"), 0.25)
})
