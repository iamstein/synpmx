# remediate_identifiable_subjects(): truncate a lone long follow-up, drop the
# rest. Reuses the single-axis fixtures from the flag tests.

rm_mk <- function(id, obs_times = c(0.5, 1, 2, 4), n_dose = 1L, dose = 100,
                  peak = 5) {
  dose_times <- seq(0, by = 12, length.out = n_dose)
  time <- c(dose_times, obs_times)
  data.frame(
    ID = as.character(id), TIME = time,
    DV = c(rep(0, length(dose_times)), peak * exp(-0.3 * obs_times)),
    AMT = c(rep(dose, length(dose_times)), rep(0, length(obs_times))),
    EVID = c(rep(1L, length(dose_times)), rep(0L, length(obs_times))),
    CMT = c(rep(1L, length(dose_times)), rep(2L, length(obs_times))),
    stringsAsFactors = FALSE
  )
}

rm_roles <- function() {
  pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
            evid = "EVID", cmt = "CMT")
}

test_that("a time-only outlier is truncated, not dropped", {
  data <- rbind(do.call(rbind, lapply(1:6, rm_mk)),
                rm_mk(7, obs_times = c(0.5, 1, 2, 40)))
  out <- suppressMessages(remediate_identifiable_subjects(data, rm_roles()))

  expect_equal(attr(out, "truncated"), "7")
  expect_length(attr(out, "dropped"), 0L)
  expect_true("7" %in% out$ID)                       # kept
  late <- out$TIME[out$ID == "7" & out$EVID == 0]
  expect_true(max(late) <= attr(out, "horizon"))     # tail removed
  expect_true(validate_pmx(out, rm_roles())$valid)
})

test_that("a DV outlier is dropped, not truncated", {
  data <- rbind(do.call(rbind, lapply(1:6, rm_mk)), rm_mk(7, peak = 50))
  out <- suppressMessages(remediate_identifiable_subjects(data, rm_roles()))

  expect_equal(attr(out, "dropped"), "7")
  expect_length(attr(out, "truncated"), 0L)
  expect_false("7" %in% out$ID)
})

test_that("a subject flagged for both time and DV is dropped, not truncated", {
  data <- rbind(do.call(rbind, lapply(1:6, rm_mk)),
                rm_mk(7, obs_times = c(0.5, 1, 2, 40), peak = 50))
  out <- suppressMessages(remediate_identifiable_subjects(data, rm_roles()))

  expect_equal(attr(out, "dropped"), "7")
  expect_length(attr(out, "truncated"), 0L)
})

test_that("the time and other options are honored", {
  data <- rbind(do.call(rbind, lapply(1:6, rm_mk)),
                rm_mk(7, obs_times = c(0.5, 1, 2, 40)),
                rm_mk(8, peak = 50))

  dropped_all <- suppressMessages(
    remediate_identifiable_subjects(data, rm_roles(), time = "drop")
  )
  expect_setequal(attr(dropped_all, "dropped"), c("7", "8"))
  expect_length(attr(dropped_all, "truncated"), 0L)

  kept_other <- suppressMessages(
    remediate_identifiable_subjects(data, rm_roles(), other = "keep")
  )
  expect_true("8" %in% kept_other$ID)                # DV outlier kept
  expect_equal(attr(kept_other, "truncated"), "7")   # time outlier truncated
})

test_that("a homogeneous cohort is returned unchanged", {
  data <- do.call(rbind, lapply(1:6, rm_mk))
  out <- suppressMessages(remediate_identifiable_subjects(data, rm_roles()))

  expect_length(attr(out, "dropped"), 0L)
  expect_length(attr(out, "truncated"), 0L)
  expect_equal(nrow(out), nrow(data))
})

test_that("dropped subjects are replaced when a source is supplied", {
  source_data <- do.call(rbind, lapply(1:8, rm_mk))
  roles <- rm_roles()
  synth <- suppressWarnings(
    synpmx_avatar(source_data, roles, n_subjects = 8L, seed = 4)
  )
  n0 <- length(unique(synth$ID))

  # Force one synthetic subject to be an extreme-DV outlier.
  victim <- unique(synth$ID)[1]
  synth$DV[synth$ID == victim & synth$EVID == 0] <- 500

  out <- suppressMessages(remediate_identifiable_subjects(
    synth, roles, source = source_data, seed = 99
  ))

  expect_gte(attr(out, "replaced"), 1L)
  expect_false(victim %in% out$ID)                    # outlier gone
  expect_equal(length(unique(out$ID)), n0)            # cohort size restored
  expect_true(validate_pmx(out, roles)$valid)
  # the refilled cohort has no remaining flagged subjects
  expect_equal(sum(flag_identifiable_subjects(out, roles)$flagged), 0L)
})

test_that("without a source, dropped subjects are not replaced", {
  data <- rbind(do.call(rbind, lapply(1:6, rm_mk)), rm_mk(7, peak = 50))
  out <- suppressMessages(remediate_identifiable_subjects(data, rm_roles()))
  expect_equal(attr(out, "replaced"), 0L)
  expect_equal(length(unique(out$ID)), 6L)            # 7 dropped, not refilled
})
