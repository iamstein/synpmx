# Steps 1 and 2 of the PMX model generator: which endpoint is the drug
# concentration, and what design produced it. Fixtures are built to land on each
# branch, including the ones where the source does not settle the question and
# both answers survive into the candidate set.

.design_roles <- function(...) {
  pmx_roles(id = "ID", time = "TIME", nominal_time = "NTIME", dv = "DV",
            amt = "AMT", evid = "EVID", cmt = "CMT", dvid = "DVID",
            mdv = "MDV", ...)
}

# One arm, `n` subjects, sampled at `times` after a dose at 0, with `dv_at` a
# function of time. Everything below is this study with one property changed.
.pk_fixture <- function(n = 24, times = c(0.25, 0.5, 1, 2, 4, 8, 12, 24),
                        dv_at = function(t, i) 10 * (exp(-0.1 * t) -
                                                       exp(-1.2 * t)) * i,
                        amount = 100, second_endpoint = NULL, cmt_obs = 2L) {
  pieces <- lapply(seq_len(n), function(subject) {
    scale <- 1 + 0.2 * sin(subject)
    dose <- data.frame(TIME = 0, NTIME = 0, DV = 0,
                       AMT = amount[(subject - 1L) %% length(amount) + 1L],
                       EVID = 1L, CMT = 1L, DVID = "cp", MDV = 1L)
    obs <- data.frame(
      TIME = times, NTIME = times,
      DV = dv_at(times, scale) *
        (amount[(subject - 1L) %% length(amount) + 1L] / amount[1L]),
      AMT = 0, EVID = 0L, CMT = cmt_obs, DVID = "cp", MDV = 0L
    )
    rows <- rbind(dose, obs)
    if (!is.null(second_endpoint)) {
      rows <- rbind(rows, data.frame(
        TIME = times, NTIME = times, DV = second_endpoint(times, scale),
        AMT = 0, EVID = 0L, CMT = 3L, DVID = "biomarker", MDV = 0L
      ))
    }
    rows$ID <- subject
    rows
  })
  out <- do.call(rbind, pieces)
  out$DVID <- factor(out$DVID)
  rownames(out) <- NULL
  out
}

test_that("a one-endpoint oral study is classified and its route detected", {
  data <- .pk_fixture()
  roles <- .design_roles()
  obs <- .model_observations(data, roles)
  classified <- .model_classify_endpoints(data, roles, obs)
  expect_identical(classified$pk, "cp")
  expect_identical(classified$decided_by, "inferred")

  design <- .model_detect_design(data, roles, obs, classified$pk)
  expect_identical(design$route, "oral")
  expect_true(design$richness$rich)
  # One compartment is the whole default set even where the sampling would
  # support a distribution phase; `pk = "2cmt_oral"` is how you ask for one.
  expect_identical(design$candidates, "1cmt_oral")
})

test_that("a sample at the moment of the dose makes it intravenous", {
  data <- .pk_fixture(times = c(0, 0.5, 1, 2, 4, 8, 12, 24),
                      dv_at = function(t, i) 10 * exp(-0.15 * t) * i)
  roles <- .design_roles()
  obs <- .model_observations(data, roles)
  design <- .model_detect_design(data, roles, obs, "cp")
  expect_identical(design$route, "iv")
  expect_identical(design$candidates, "1cmt_iv")
})

test_that("a declared rate makes it an infusion, with no search over routes", {
  data <- .pk_fixture(dv_at = function(t, i) 10 * exp(-0.15 * t) * i)
  data$RATE <- ifelse(data$EVID != 0L, 50, 0)
  roles <- .design_roles(rate = "RATE")
  obs <- .model_observations(data, roles)
  design <- .model_detect_design(data, roles, obs, "cp")
  expect_identical(design$route, "infusion")
  expect_identical(design$candidates, "1cmt_infusion")
  expect_match(design$reason, "rate")
})

test_that("a study sampled only at troughs would not support two compartments", {
  data <- .pk_fixture(times = c(12, 24))
  roles <- .design_roles()
  obs <- .model_observations(data, roles)
  design <- .model_detect_design(data, roles, obs, "cp")
  expect_false(design$richness$rich)
})

test_that("no design offers a two-compartment candidate on its own", {
  # The candidate set is one-compartment whatever the sampling shows. A
  # distribution phase is a refinement of a shape the one-compartment model
  # already has, and fitting for it costs five times as long.
  for (times in list(c(0.25, 0.5, 1, 2, 4, 8, 12, 24), c(12, 24),
                     c(0, 0.5, 1, 2, 4, 8, 12, 24))) {
    data <- .pk_fixture(times = times)
    roles <- .design_roles()
    obs <- .model_observations(data, roles)
    design <- .model_detect_design(data, roles, obs, "cp")
    expect_false(any(grepl("^2cmt", design$candidates)))
  }
})

test_that("a study nobody samples early enough offers both routes", {
  # Every profile begins after the peak, so every one of them declines from its
  # first point whichever route it came from. Nothing here is evidence, and the
  # candidate set says so rather than committing to the reading the schedule
  # produced. This is the `warfarin` case, where 22 of 32 patients are first
  # sampled at 24 hours.
  data <- .pk_fixture(times = c(24, 36, 48, 72),
                      dv_at = function(t, i) 10 * exp(-0.05 * t) * i)
  roles <- .design_roles()
  obs <- .model_observations(data, roles)
  design <- .model_detect_design(data, roles, obs, "cp")
  expect_identical(design$route, "both")
  expect_setequal(design$candidates, c("1cmt_iv", "1cmt_oral"))
  expect_match(design$reason, "declines from its first sample")
})

test_that("a biomarker beside a concentration is not mistaken for it", {
  # The biomarker rises from a baseline and stays up: present before the dose,
  # and no rise-and-fall. The concentration is neither.
  data <- .pk_fixture(second_endpoint = function(t, i) 80 + 2 * t * i)
  data$NTIME[data$DVID == "biomarker"] <- data$TIME[data$DVID == "biomarker"]
  # Give the biomarker a genuine pre-dose observation, which is what a baseline
  # has and a drug concentration does not.
  baseline <- data[data$DVID == "biomarker" &
                     abs(data$TIME - 0.25) < 1e-8, , drop = FALSE]
  baseline$TIME <- -1
  baseline$NTIME <- -1
  data <- rbind(data, baseline)
  roles <- .design_roles()
  obs <- .model_observations(data, roles)
  signals <- .model_endpoint_signals(data, roles, obs)
  expect_false(signals$post_dose[signals$endpoint == "biomarker"])
  expect_true(signals$post_dose[signals$endpoint == "cp"])
  classified <- .model_classify_endpoints(data, roles, obs)
  expect_identical(classified$pk, "cp")
  expect_identical(classified$pd, "biomarker")
})

test_that("dose proportionality separates a concentration from a biomarker", {
  # Two arms at 100 and 400. The concentration scales with the dose; the
  # biomarker does not.
  data <- .pk_fixture(amount = c(100, 400),
                      second_endpoint = function(t, i) 80 - 0.5 * t * i)
  roles <- .design_roles()
  obs <- .model_observations(data, roles)
  signals <- .model_endpoint_signals(data, roles, obs)
  expect_true(signals$proportional[signals$endpoint == "cp"])
  expect_false(signals$proportional[signals$endpoint == "biomarker"])
  expect_identical(.model_classify_endpoints(data, roles, obs)$pk, "cp")
})

test_that("a dose that varies per patient is not read as a set of dose levels", {
  # Dosed by body weight: every patient gets a different amount, so reading each
  # one as its own level would compare a median of one against a median of one
  # over a dose ratio near 1, which almost anything passes.
  data <- .pk_fixture(amount = 100 + seq_len(24))
  roles <- .design_roles()
  obs <- .model_observations(data, roles)
  expect_length(.model_dose_levels(obs), 0L)
  expect_true(is.na(.model_endpoint_signals(data, roles, obs)$proportional))
})

test_that("the compartment signal breaks a tie the required signals leave", {
  # Both endpoints pass the required signals; only one sits in the compartment
  # one above the dosing compartment nobody observes.
  data <- .pk_fixture(second_endpoint = function(t, i) 9 * (exp(-0.1 * t) -
                                                              exp(-1.2 * t)) * i)
  roles <- .design_roles()
  obs <- .model_observations(data, roles)
  signals <- .model_endpoint_signals(data, roles, obs)
  expect_true(signals$compartment[signals$endpoint == "cp"])
  expect_false(signals$compartment[signals$endpoint == "biomarker"])
  expect_identical(.model_classify_endpoints(data, roles, obs)$pk, "cp")
})

test_that("no endpoint looking like a concentration is refused, not guessed", {
  # A baseline biomarker on its own: present before the dose, so the one
  # required signal that can be computed here fails.
  data <- .pk_fixture(dv_at = function(t, i) 80 + 2 * t * i)
  pre <- data[data$EVID == 0L & abs(data$TIME - 0.25) < 1e-8, , drop = FALSE]
  pre$TIME <- -1
  pre$NTIME <- -1
  data <- rbind(data, pre)
  roles <- .design_roles()
  obs <- .model_observations(data, roles)
  expect_error(.model_classify_endpoints(data, roles, obs),
               "No endpoint looks like a drug concentration")
  expect_error(.model_classify_endpoints(data, roles, obs), "endpoint_roles")
})

test_that("endpoint_roles overrides the signals rather than hinting at them", {
  data <- .pk_fixture(second_endpoint = function(t, i) 80 + 2 * t * i)
  roles <- .design_roles()
  obs <- .model_observations(data, roles)
  classified <- .model_classify_endpoints(data, roles, obs,
                                          endpoint_roles = c(pk = "biomarker"))
  expect_identical(classified$pk, "biomarker")
  expect_identical(classified$decided_by, "declared")
  expect_identical(classified$pd, "cp")
})

test_that("endpoint_roles naming something unfittable says which", {
  data <- .pk_fixture()
  roles <- .design_roles()
  obs <- .model_observations(data, roles)
  expect_error(
    .model_classify_endpoints(data, roles, obs, endpoint_roles = c(pk = "nope")),
    "not an endpoint in this data"
  )
})

test_that("a binary endpoint is never a candidate and never a PD time course", {
  data <- .pk_fixture(second_endpoint = function(t, i) as.numeric(t > 4))
  roles <- .design_roles()
  obs <- .model_observations(data, roles)
  classified <- .model_classify_endpoints(data, roles, obs)
  expect_identical(classified$pk, "cp")
  expect_false("biomarker" %in% classified$pd)
  expect_true("biomarker" %in% classified$discrete)
})

test_that("an endpoint recorded as whole numbers is still a time course", {
  # `warfarin`'s prothrombin activity is typed `integer` and is a continuous
  # quantity that was rounded, which is what `.snap_endpoint_values()` puts back.
  data <- .pk_fixture(second_endpoint = function(t, i) round(80 + 2 * t * i))
  roles <- .design_roles()
  obs <- .model_observations(data, roles)
  classified <- .model_classify_endpoints(data, roles, obs)
  expect_identical(classified$pd, "biomarker")
})

test_that("the signals are read on the nominal grid, not the recorded clock", {
  # Same study, recorded times jittered per patient. On recorded times no two
  # patients share a time, the median profile is one point per time, and the
  # shape signal reads noise.
  data <- .pk_fixture()
  set.seed(11)
  jitter <- stats::runif(nrow(data), -0.05, 0.05)
  data$TIME <- ifelse(data$EVID == 0L, data$TIME + jitter, data$TIME)
  roles <- .design_roles()
  obs <- .model_observations(data, roles)
  expect_true(.model_signal_shape(obs, "cp", 1L))
  expect_identical(.model_detect_design(data, roles, obs, "cp")$route, "oral")

  # And the recorded axis is carried too, because estimation needs it.
  expect_false(isTRUE(all.equal(obs$tad, obs$actual_tad)))
})

test_that("the richest dose interval is the one the signals are read in", {
  data <- .pk_fixture()
  second <- data[data$EVID == 0L, , drop = FALSE]
  second$TIME <- second$TIME + 48
  second$NTIME <- second$NTIME + 48
  dose <- data[data$EVID != 0L, , drop = FALSE]
  dose$TIME <- 48
  dose$NTIME <- 48
  # A second interval with one sample per patient cannot be the richest.
  data <- rbind(data, dose, second[seq_len(24), , drop = FALSE])
  roles <- .design_roles()
  obs <- .model_observations(data, roles)
  expect_identical(.model_richest_interval(obs), 1L)
})

test_that("a PD time course is fitted on the axis it will be generated on", {
  # Two bugs met here, both invisible on a single-dose study. A PD endpoint is
  # evaluated at the visit grid's times, which are study time from the first
  # dose, and it used to be fitted against time after dose. On a daily regimen
  # almost every sample sits at the same time after its own dose, so the curve
  # was fitted against one point and then drawn over months.
  times <- c(0, 24, 48, 168, 336, 672)
  pieces <- lapply(seq_len(24), function(subject) {
    doses <- data.frame(TIME = seq(0, 672, by = 24), NTIME = seq(0, 672, by = 24),
                        DV = 0, AMT = 100, EVID = 1L, CMT = 1L, DVID = "pd",
                        MDV = 1L)
    obs <- data.frame(TIME = times, NTIME = times,
                      DV = 100 - 0.08 * times + 0.5 * sin(subject),
                      AMT = 0, EVID = 0L, CMT = 2L, DVID = "pd", MDV = 0L)
    rows <- rbind(doses, obs)
    rows$ID <- subject
    rows
  })
  data <- do.call(rbind, pieces)
  data$DVID <- factor(data$DVID)
  roles <- .design_roles()
  obs <- .model_observations(data, roles)

  # Time after dose spans one dosing interval; study time spans the study.
  expect_lt(diff(range(obs$tad, na.rm = TRUE)), 24)
  expect_gt(diff(range(obs$aligned)), 600)

  shape <- .model_fit_pd(obs, "pd")
  expect_identical(shape$pd, "linear")
  expect_equal(unname(shape$typical[["slope"]]), -0.08, tolerance = 0.05)
})

test_that("an arm with no dose records does not take the PD fit down", {
  # A placebo arm records its administrations with `AMT = 0`, which
  # `.dose_rows()` drops whenever any positive amount exists in the study, so
  # every one of its rows sits outside any dose interval. Asking for each
  # subject's earliest observation by time after dose then returns nothing.
  times <- c(0, 24, 48, 168)
  make <- function(subject, amount) {
    doses <- data.frame(TIME = 0, NTIME = 0, DV = 0, AMT = amount, EVID = 1L,
                        CMT = 1L, DVID = "pd", MDV = 1L)
    obs <- data.frame(TIME = times, NTIME = times,
                      DV = 90 - 0.05 * times * (amount > 0) + 0.3 * sin(subject),
                      AMT = 0, EVID = 0L, CMT = 2L, DVID = "pd", MDV = 0L)
    rows <- rbind(doses, obs)
    rows$ID <- subject
    rows
  }
  data <- do.call(rbind, c(lapply(1:12, make, amount = 100),
                           lapply(13:24, make, amount = 0)))
  data$DVID <- factor(data$DVID)
  roles <- .design_roles()
  obs <- .model_observations(data, roles)
  placebo <- obs[obs$subject %in% as.character(13:24), ]
  expect_true(all(is.na(placebo$tad)))
  expect_true(all(is.finite(placebo$aligned)))
  expect_s3_class(.model_fit_pd(obs, "pd")$candidates, "data.frame")
})
