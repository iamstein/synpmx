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
