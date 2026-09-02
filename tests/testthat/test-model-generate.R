# Generating from a fitted model. Everything here runs on base R against a
# hand-constructed `pmx_fitted_model`, because the generator has no fitter in
# it: `.pk_profile()` evaluates the profile and the dosing and visit models are
# summaries of the source. Nothing in this file needs nlmixr2.

.generate_roles <- function(...) {
  pmx_roles(id = "ID", time = "TIME", nominal_time = "NTIME", dv = "DV",
            amt = "AMT", evid = "EVID", cmt = "CMT", dvid = "DVID",
            mdv = "MDV", covariates = "WT", ...)
}

# A study on a weekly cycle, so that a dose ladder and the three rates have
# somewhere to live. `reduce_from` puts every patient onto half dose from that
# cycle on, which is what gives the dosing model a ladder to read.
.cycle_fixture <- function(n = 24, cycles = 6, reduce_from = NULL,
                           stop_after = NULL) {
  times <- (seq_len(cycles) - 1L) * 168
  samples <- c(0, 1, 4, 12, 24, 72, 167)
  pieces <- lapply(seq_len(n), function(subject) {
    last <- if (is.null(stop_after)) cycles else
      if (subject %% 3L == 0L) stop_after else cycles
    amounts <- rep(100, cycles)
    if (!is.null(reduce_from) && subject %% 2L == 0L) {
      amounts[reduce_from:cycles] <- 50
    }
    doses <- data.frame(TIME = times[seq_len(last)], NTIME = times[seq_len(last)],
                        DV = 0, AMT = amounts[seq_len(last)], EVID = 1L,
                        CMT = 1L, DVID = "cp", MDV = 1L)
    observation_times <- as.numeric(outer(samples, times[seq_len(last)], "+"))
    obs <- data.frame(
      TIME = observation_times, NTIME = observation_times,
      DV = 5 * exp(-0.02 * (observation_times %% 168)) + 0.1 * subject,
      AMT = 0, EVID = 0L, CMT = 2L, DVID = "cp", MDV = 0L
    )
    rows <- rbind(doses, obs)
    rows$ID <- subject
    rows$WT <- 70 + 10 * sin(subject)
    rows
  })
  out <- do.call(rbind, pieces)
  out$DVID <- factor(out$DVID)
  out <- out[order(out$ID, out$TIME, out$EVID == 0L), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# A fitted model with the apparatus built from a real source and the estimated
# half supplied by hand. This is exactly what the build order calls a complete
# generator with no fitter in it.
.hand_built_fit <- function(data, roles, fixed = c(cl = 5, v = 50, ka = 1.2),
                            structural = "1cmt_oral", cv = 0.15) {
  subject_group <- .model_subject_arms(data, roles)
  observations <- .model_observations(data, roles)
  classified <- .model_classify_endpoints(data, roles, observations,
                                          endpoint_roles = c(pk = "cp"))
  cells <- .model_cells(data, roles, c(classified$pk, classified$pd,
                                       classified$discrete), 3L)
  planned <- data
  planned[[roles$time]] <- as.numeric(data[[roles$nominal_time]])
  arm_models <- .arm_models(planned, roles, cells, subject_group, 3L)
  omega <- diag(rep(0.09, length(fixed)))
  dimnames(omega) <- list(names(fixed), names(fixed))

  .pmx_fitted_model(
    structural = structural,
    candidates = data.frame(model = structural, converged = TRUE, aic = 1,
                            note = "", stringsAsFactors = FALSE),
    parameters = list(fixed = fixed, omega = omega,
                      residual = list(kind = "proportional", cv = cv)),
    endpoints = list(pk = classified$pk, pd = character(),
                     discrete = classified$discrete, decided_by = "declared"),
    arms = list(arms = arm_models$arms, sizes = arm_models$sizes),
    dosing = arm_models$dosing, visits = arm_models$visits,
    schema = .source_schema(data, roles, classified$pk, subject_group),
    roles = roles,
    settings = list(min_arm_patients = 3L), n_source = length(subject_group),
    cells = cells, covariates = .covariate_model(data, roles, subject_group),
    discrete = .discrete_model(data, roles, cells, subject_group)
  )
}

test_that("a hand-built model generates a legal dataset in the source's shape", {
  data <- .cycle_fixture()
  roles <- .generate_roles()
  fit <- .hand_built_fit(data, roles)
  synthetic <- synpmx_model_generate(fit, n_subjects = 20, seed = 3)

  expect_true(validate_pmx(synthetic, roles)$valid)
  expect_length(unique(synthetic$ID), 20L)
  expect_setequal(names(synthetic), names(data))
  expect_true(all(synthetic$DV[synthetic$EVID == 0L] >= 0))
})

test_that("generation is reproducible from its seed and reads no patient data", {
  data <- .cycle_fixture()
  roles <- .generate_roles()
  fit <- .hand_built_fit(data, roles)
  expect_identical(synpmx_model_generate(fit, n_subjects = 12, seed = 5),
                   synpmx_model_generate(fit, n_subjects = 12, seed = 5))
  # The second stage's arguments are the model and a count. Nothing it can
  # reach holds a row of the source.
  expect_false(any(vapply(fit, is.data.frame, logical(1))[
    c("dosing", "visits", "covariates", "discrete")
  ]))
})

test_that("no generated subject reproduces a source subject's values", {
  data <- .cycle_fixture()
  roles <- .generate_roles()
  synthetic <- synpmx_model_generate(.hand_built_fit(data, roles),
                                     n_subjects = 24, seed = 9)
  source_values <- vapply(split(data$DV[data$EVID == 0L],
                                data$ID[data$EVID == 0L]),
                          function(x) paste(sprintf("%.6g", x), collapse = "|"),
                          character(1))
  generated <- vapply(split(synthetic$DV[synthetic$EVID == 0L],
                            synthetic$ID[synthetic$EVID == 0L]),
                      function(x) paste(sprintf("%.6g", x), collapse = "|"),
                      character(1))
  expect_length(intersect(source_values, generated), 0L)
})

test_that("a drawn dose reduction reaches the concentrations", {
  # The largest fidelity gain over `synpmx_pca()`, and the check the design
  # asks for by name: the schedule is drawn first and the profile computed from
  # it, so a subject who steps down a level has a lower exposure from that cycle
  # on. Comparing dosing records alone would pass whether or not that is true.
  data <- .cycle_fixture(reduce_from = 4L)
  roles <- .generate_roles()
  fit <- .hand_built_fit(data, roles)
  expect_gt(length(fit$dosing[[1L]]$levels), 1L)
  expect_gt(fit$dosing[[1L]]$reduction, 0)

  synthetic <- synpmx_model_generate(fit, n_subjects = 60, seed = 21)
  by_subject <- split(synthetic, synthetic$ID)
  ratios <- vapply(by_subject, function(part) {
    doses <- part[part$EVID != 0L, , drop = FALSE]
    if (nrow(doses) < 2L) return(NA_real_)
    dropped <- which(diff(doses$AMT) < -1e-8)
    if (!length(dropped)) return(NA_real_)
    at <- doses$TIME[dropped[1L] + 1L]
    observations <- part[part$EVID == 0L, , drop = FALSE]
    before <- observations$DV[observations$TIME < at]
    after <- observations$DV[observations$TIME >= at]
    if (!length(before) || !length(after)) return(NA_real_)
    stats::median(after) / stats::median(before)
  }, numeric(1))
  ratios <- ratios[is.finite(ratios)]
  expect_gt(length(ratios), 5L)
  # Exposure after the reduction is materially lower, not merely different.
  expect_lt(stats::median(ratios), 0.9)
})

test_that("a study nobody departs from reproduces its planned schedule", {
  data <- .cycle_fixture()
  roles <- .generate_roles()
  fit <- .hand_built_fit(data, roles)
  expect_identical(fit$dosing[[1L]]$levels, 1)
  expect_identical(fit$dosing[[1L]]$reduction, 0)
  synthetic <- synpmx_model_generate(fit, n_subjects = 15, seed = 4)
  amounts <- unique(synthetic$AMT[synthetic$EVID != 0L])
  expect_length(amounts, 1L)
})

test_that("a skipped cycle leaves the trough it implies", {
  data <- .cycle_fixture(stop_after = 3L)
  roles <- .generate_roles()
  fit <- .hand_built_fit(data, roles)
  expect_gt(fit$dosing[[1L]]$discontinuation, 0)
  synthetic <- synpmx_model_generate(fit, n_subjects = 40, seed = 6)
  spans <- vapply(split(synthetic, synthetic$ID), function(part) {
    sum(part$EVID != 0L)
  }, integer(1))
  expect_gt(length(unique(spans)), 1L)
})

test_that("the visit model refuses a grid no arm shares", {
  data <- .cycle_fixture(n = 24)
  # Give every patient their own observation times: no cell is held by three.
  data$NTIME[data$EVID == 0L] <- data$NTIME[data$EVID == 0L] +
    data$ID[data$EVID == 0L] / 1000
  roles <- .generate_roles()
  expect_error(.model_cells(data, roles, "cp", 3L),
               "held by at least 3 patients")
})

test_that("covariates are drawn from the arm's model, not copied", {
  data <- .cycle_fixture()
  roles <- .generate_roles()
  synthetic <- synpmx_model_generate(.hand_built_fit(data, roles),
                                     n_subjects = 24, seed = 12)
  drawn <- vapply(split(synthetic$WT, synthetic$ID), function(x) x[1L],
                  numeric(1))
  source_weights <- vapply(split(data$WT, data$ID), function(x) x[1L],
                           numeric(1))
  expect_length(intersect(sprintf("%.8g", drawn),
                          sprintf("%.8g", source_weights)), 0L)
  expect_lt(abs(mean(drawn) - mean(source_weights)), 5)
})

test_that("a covariate effect moves the profile it is declared on", {
  data <- .cycle_fixture()
  roles <- .generate_roles()
  fit <- .hand_built_fit(data, roles)
  fit$covariate_effects <- list(
    cl = list(covariate = "WT", reference = 70, exponent = 0.75)
  )
  light <- .subject_parameters(fit$parameters$fixed, fit$covariate_effects,
                               c(cl = 0, v = 0, ka = 0), list(WT = 50))
  heavy <- .subject_parameters(fit$parameters$fixed, fit$covariate_effects,
                               c(cl = 0, v = 0, ka = 0), list(WT = 100))
  expect_lt(light[["cl"]], heavy[["cl"]])
  expect_identical(light[["v"]], heavy[["v"]])
})

test_that("random effects are drawn from the covariance, not read off anybody", {
  omega <- matrix(c(0.16, 0.06, 0.06, 0.09), 2, 2,
                  dimnames = list(c("cl", "v"), c("cl", "v")))
  set.seed(2)
  draws <- .draw_random_effects(omega, 20000L)
  expect_equal(diag(stats::cov(draws)), diag(omega), tolerance = 0.05)
  expect_equal(stats::cov(draws)[1L, 2L], omega[1L, 2L], tolerance = 0.05)
  expect_identical(colnames(draws), c("cl", "v"))
})

test_that("generation refuses anything that is not a fitted model", {
  expect_error(synpmx_model_generate(list()), "must come from")
  data <- .cycle_fixture()
  roles <- .generate_roles()
  fit <- .hand_built_fit(data, roles)
  expect_error(synpmx_model_generate(fit, n_subjects = 0), "positive integer")
})

test_that("the generated dataset carries its fitted model", {
  data <- .cycle_fixture()
  roles <- .generate_roles()
  fit <- .hand_built_fit(data, roles)
  synthetic <- synpmx_model_generate(fit, n_subjects = 10, seed = 1)
  expect_s3_class(attr(synthetic, "pmx_fitted_model"), "pmx_fitted_model")
  expect_identical(attr(synthetic, "pmx_source"), "model")
})

test_that("a discrete endpoint holds up where a visit recorded one level", {
  # `sample(x, 1)` reads a length-one numeric `x` as `seq_len(x)`, so a visit
  # every patient recorded the same level at asked for a draw from that
  # level's own value and failed with "incorrect number of probabilities".
  # `xgxr::mad` is the study that found it; the fixture below is the smallest
  # reproduction. SIM-059.
  data <- .cycle_fixture(n = 12)
  observations <- data[data$EVID == 0L, , drop = FALSE]
  severity <- observations
  severity$DVID <- "severity"
  severity$CMT <- 3L
  # A three-level scale over the study, but every patient records 3 at the
  # first visit -- so that one cell's marginal holds a single level whose
  # value is larger than one, which is a draw from `seq_len(3)` under the old
  # code.
  severity$DV <- 1 + (severity$ID %% 3L)
  severity$DV[severity$NTIME == min(severity$NTIME)] <- 3
  data <- rbind(data, severity)
  data$DVID <- factor(data$DVID)
  data <- data[order(data$ID, data$TIME, data$EVID == 0L), , drop = FALSE]
  rownames(data) <- NULL

  roles <- .generate_roles()
  fit <- .hand_built_fit(data, roles)
  expect_true("severity" %in% fit$endpoints$discrete)

  synthetic <- synpmx_model_generate(fit, n_subjects = 12, seed = 7)
  drawn <- synthetic$DV[synthetic$DVID == "severity" & synthetic$EVID == 0L]
  expect_true(length(drawn) > 0L)
  expect_true(all(drawn %in% c(1, 2, 3)))
})
