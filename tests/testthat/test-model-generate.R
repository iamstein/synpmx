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
    cells = cells, covariates = .covariate_model(data, roles),
    discrete = .discrete_model(data, roles, cells, subject_group,
                               setdiff(unique(cells$endpoint), classified$pk))
  )
}

# SIM-073. An infusion that generates as a bolus is a different study. The
# duration is what makes it an infusion, and it comes from the rate the arm was
# given the dose at rather than from a parameter nobody set.
test_that("an infusion generates as an infusion, with its rate written back", {
  data <- .cycle_fixture()
  data$RATE <- ifelse(data$EVID != 0L, data$AMT / 2, 0)
  roles <- .generate_roles(rate = "RATE")
  fit <- .hand_built_fit(data, roles, structural = "1cmt_infusion",
                         fixed = c(cl = 5, v = 50))
  synthetic <- as.data.frame(synpmx_model_generate(fit, n_subjects = 8, seed = 3))
  doses <- synthetic[synthetic$EVID != 0L, , drop = FALSE]

  expect_true(all(doses$RATE > 0))
  # Two hours, which is the duration the source's rate implies.
  expect_equal(unique(round(doses$AMT / doses$RATE, 6)), 2)
  expect_true(all(doses$AMT > 0))

  # And the concentration is the infusion's, not the bolus it used to be: an
  # infusion has a lower peak than the same dose given instantly.
  bolus <- .hand_built_fit(data, roles, structural = "1cmt_iv",
                           fixed = c(cl = 5, v = 50))
  peak <- function(fit) {
    d <- as.data.frame(synpmx_model_generate(fit, n_subjects = 8, seed = 3))
    max(d$DV[d$EVID == 0L], na.rm = TRUE)
  }
  expect_lt(peak(fit), peak(bolus))
})

# SIM-072. NONMEM writes the end of an infusion as a mirror record, `AMT` and
# `RATE` both negated. Reading every event as an administration planned a
# schedule of alternating doses and anti-doses, and generated them.
test_that("an infusion stop record is not planned as a dose", {
  data <- .cycle_fixture()
  data$RATE <- ifelse(data$EVID != 0L, data$AMT, 0)
  stops <- data[data$EVID != 0L, , drop = FALSE]
  stops$TIME <- stops$TIME + 1
  stops$NTIME <- stops$NTIME + 1
  stops$AMT <- -stops$AMT
  stops$RATE <- -stops$RATE
  data <- rbind(data, stops)
  data <- data[order(data$ID, data$TIME, data$EVID == 0L), , drop = FALSE]

  roles <- .generate_roles(rate = "RATE")
  fit <- .hand_built_fit(data, roles, structural = "1cmt_infusion",
                         fixed = c(cl = 5, v = 50))
  expect_true(all(fit$dosing[[1L]]$planned$amt >= 0))
  synthetic <- as.data.frame(synpmx_model_generate(fit, n_subjects = 6, seed = 4))
  expect_true(all(synthetic$AMT >= 0))
})

# The subject's arm has to reach their PD value, or fitting per arm buys
# nothing. Built by hand: two arms, two shapes, and the arm each subject was
# assigned decides which one is evaluated.
test_that("generation evaluates the arm's own PD shape where there is one", {
  data <- .cycle_fixture()
  data$ARM <- ifelse(as.integer(data$ID) %% 2L == 0L, "high", "low")
  roles <- .generate_roles(strata = "ARM")
  fit <- .hand_built_fit(data, roles)

  flat <- list(pd = "constant", typical = c(baseline = 10), baseline_cv = 0,
               residual = list(kind = "additive", sd = 0))
  steep <- list(pd = "constant", typical = c(baseline = 100), baseline_cv = 0,
                residual = list(kind = "additive", sd = 0))
  fit$endpoints$pd <- "effect"
  fit$pd <- list(effect = c(flat, list(arms = list(low = flat, high = steep))))
  fit$cells <- rbind(fit$cells, transform(fit$cells[1, , drop = FALSE],
                                          endpoint = "effect"))
  for (arm in names(fit$visits)) {
    fit$visits[[arm]]$probability <- c(fit$visits[[arm]]$probability, 1)
  }

  synthetic <- synpmx_model_generate(fit, n_subjects = 40, seed = 9)
  effect <- synthetic[synthetic$DVID == "effect", , drop = FALSE]
  by_arm <- tapply(effect$DV, effect$ARM, mean)
  expect_equal(unname(by_arm[["low"]]), 10)
  expect_equal(unname(by_arm[["high"]]), 100)
})

# Everything on the object is an input to generation, so printing it has to
# report all of them: a reader who has to call a second function to see half the
# simulation will read half of it.
test_that("printing a fit reports every input the generator simulates from", {
  data <- .cycle_fixture(reduce_from = 4L, stop_after = 3L)
  roles <- .generate_roles()
  fit <- .hand_built_fit(data, roles)
  out <- paste(utils::capture.output(print(fit)), collapse = "\n")
  for (heading in c("structural model", "fixed effects", "between-subject",
                    "residual error", "cohort",
                    "dose changes", "visit grid", "visit attendance",
                    "covariates",
                    "columns emitted")) {
    expect_match(out, heading, fixed = TRUE)
  }
  # The arm key is joined with a control character internally and never shown
  # with one.
  expect_false(grepl("\r", out, fixed = TRUE))
})

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

test_that("covariates are drawn from the study's model, not copied", {
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

# SIM-061. The proportional multiplier is lognormal, so a residual estimated
# from a misfitted structural model widens the band instead of dropping values
# onto zero.
test_that("how often a value reaches zero does not depend on the residual CV", {
  data <- .cycle_fixture()
  roles <- .generate_roles()
  zero_rate <- function(cv) {
    fit <- .hand_built_fit(data, roles, cv = cv)
    fit$quantification_floor <- NULL
    synthetic <- synpmx_model_generate(fit, n_subjects = 60, seed = 4)
    dv <- synthetic$DV[synthetic$EVID == 0L & !is.na(synthetic$DV)]
    mean(dv == 0)
  }
  # Whatever reaches zero here is the profile underflowing late in a dose
  # interval, which the floor is what handles. The residual must not add to it:
  # under `1 + N(0, cv)` the rate would climb by `pnorm(-1 / cv)`, which is
  # 7.7 points between these two, and each of those was a value clamped to zero.
  expect_equal(zero_rate(0.7), zero_rate(0.05), tolerance = 0.005)
})

test_that("the multiplier keeps the coefficient of variation it was given", {
  residual <- list(kind = "proportional", cv = 0.6)
  set.seed(11)
  drawn <- .add_residual_error(rep(100, 2e5), residual)
  expect_true(all(drawn > 0))
  expect_equal(stats::sd(drawn) / mean(drawn), 0.6, tolerance = 0.02)
  expect_equal(stats::median(drawn), 100, tolerance = 0.02)
})

# SIM-061. The floor: a study that declared no censoring column still had an
# assay, and nothing is emitted below half the smallest value it reported.
test_that("nothing is emitted below the quantification floor", {
  data <- .cycle_fixture()
  roles <- .generate_roles()
  floor_value <- .model_quantification_floor(data, roles, "cp")
  fit <- .hand_built_fit(data, roles, cv = 0.5)
  fit$quantification_floor <- floor_value

  synthetic <- synpmx_model_generate(fit, n_subjects = 60, seed = 4)
  dv <- synthetic$DV[synthetic$EVID == 0L & !is.na(synthetic$DV)]
  expect_true(all(dv >= floor_value$cp))
  expect_false(any(dv == 0))
})

test_that("a per-visit marginal is built only where generation reads one", {
  # SIM-083. `.discrete_model()` built a marginal for every cell of every
  # endpoint, and a marginal is the source's own values with their
  # frequencies. `.model_generate()` reads one only for an endpoint that is
  # neither a concentration nor a fitted shape, so on a study whose endpoints
  # are all time courses every one of those marginals was unreadable -- and
  # the fit object, which is written to `inst/extdata/` and shipped, carried
  # thousands of real measurements nothing could ever draw from.
  data <- .cycle_fixture()
  roles <- .generate_roles()
  cells <- .model_cells(data, roles, unique(data$DVID), 1L)
  group <- .model_subject_arms(data, roles)

  none <- .discrete_model(data, roles, cells, group, drawn = character(0))
  expect_true(all(vapply(unlist(none, recursive = FALSE), is.null,
                         logical(1))))

  one <- unique(cells$endpoint)[[1L]]
  some <- .discrete_model(data, roles, cells, group, drawn = one)
  kept <- cells$endpoint[!vapply(some[[1L]], is.null, logical(1))]
  expect_true(length(kept) > 0L)
  expect_equal(unique(kept), one)
})

test_that("a fit carrying no floor still generates", {
  data <- .cycle_fixture()
  roles <- .generate_roles()
  fit <- .hand_built_fit(data, roles)
  fit$quantification_floor <- NULL
  expect_no_error(synpmx_model_generate(fit, n_subjects = 12, seed = 2))
})

# SIM-061. A floor doing this to much of a dataset is hiding a bad fit, so it
# says what it caught rather than quietly raising the values.
test_that("a floor catching a large share of the output warns", {
  data <- .cycle_fixture()
  roles <- .generate_roles()
  fit <- .hand_built_fit(data, roles, cv = 0.5)
  fit$quantification_floor <- .model_quantification_floor(data, roles, "cp")
  expect_warning(synthetic <- synpmx_model_generate(fit, n_subjects = 40,
                                                    seed = 4),
                 "fell below the smallest value")
  caught <- attr(synthetic, "pmx_floored")
  expect_true(caught[["raised"]] > 0 && caught[["raised"]] <= caught[["seen"]])
})

test_that("a fit with no floor records nothing and warns about nothing", {
  data <- .cycle_fixture()
  roles <- .generate_roles()
  fit <- .hand_built_fit(data, roles)
  fit$quantification_floor <- NULL
  synthetic <- synpmx_model_generate(fit, n_subjects = 12, seed = 2)
  expect_null(attr(synthetic, "pmx_floored"))
})

# Two drugs in one study ------------------------------------------------------
#
# A combination is the case every PK endpoint sharing every dose record is wrong
# for. Undeclared, drug B's doses enter drug A's fit and drive drug A's
# generated profile; `dose_endpoints` is what tells the two apart, because the
# data cannot -- a metabolite has no dose records of its own.

.combination_fixture <- function(n = 12, cycles = 4) {
  cycle_times <- (seq_len(cycles) - 1L) * 168
  samples <- c(0, 1, 4, 24, 167)
  pieces <- lapply(seq_len(n), function(subject) {
    # Drug A weekly at 100 mg (ADM 1), drug B on the first day of every other
    # week at 900 mg (ADM 2). Different amounts and different schedules, which
    # is what pooling them destroys.
    a_times <- cycle_times
    b_times <- cycle_times[seq(1L, cycles, by = 2L)]
    doses <- rbind(
      data.frame(TIME = a_times, NTIME = a_times, DV = 0, AMT = 100,
                 EVID = 1L, CMT = 1L, ADM = 1L, DVID = "A conc", MDV = 1L),
      data.frame(TIME = b_times, NTIME = b_times, DV = 0, AMT = 900,
                 EVID = 1L, CMT = 3L, ADM = 2L, DVID = "B conc", MDV = 1L)
    )
    times <- as.numeric(outer(samples, cycle_times, "+"))
    obs <- rbind(
      data.frame(TIME = times, NTIME = times,
                 DV = 5 * exp(-0.02 * (times %% 168)) + 0.1 * subject,
                 AMT = 0, EVID = 0L, CMT = 2L, ADM = NA_integer_,
                 DVID = "A conc", MDV = 0L),
      data.frame(TIME = times, NTIME = times,
                 DV = 40 * exp(-0.01 * (times %% 336)) + 0.5 * subject,
                 AMT = 0, EVID = 0L, CMT = 4L, ADM = NA_integer_,
                 DVID = "B conc", MDV = 0L)
    )
    rows <- rbind(doses, obs)
    rows$ID <- subject
    rows$WT <- 70 + 10 * sin(subject)
    rows
  })
  out <- do.call(rbind, pieces)
  out <- out[order(out$ID, out$TIME, out$EVID == 0L), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.combination_roles <- function(...) {
  pmx_roles(id = "ID", time = "TIME", nominal_time = "NTIME", dv = "DV",
            amt = "AMT", evid = "EVID", cmt = "CMT", dvid = "DVID",
            mdv = "MDV", covariates = "WT", adm = "ADM",
            routes = c("1" = "extravascular", "2" = "extravascular"), ...)
}

# Two structural models, two parameter sets, and the apparatus built from the
# source exactly as `.hand_built_fit()` does it for one endpoint.
.combination_fit <- function(data, roles) {
  subject_group <- .model_subject_arms(data, roles)
  observations <- .model_observations(data, roles)
  classified <- .model_classify_endpoints(
    data, roles, observations, endpoint_roles = c(pk = c("A conc", "B conc")))
  cells <- .model_cells(data, roles, classified$pk, 3L)
  planned <- data
  planned[[roles$time]] <- as.numeric(data[[roles$nominal_time]])
  groups <- if (is.null(roles$dose_endpoints)) NULL else
    .dose_endpoint(planned, roles)
  arm_models <- .arm_models(planned, roles, cells, subject_group, 3L, groups)
  model <- function(fixed) {
    omega <- diag(rep(0.09, length(fixed)))
    dimnames(omega) <- list(names(fixed), names(fixed))
    list(structural = "1cmt_oral",
         parameters = list(fixed = fixed, omega = omega,
                           residual = list(kind = "proportional", cv = 0.1)),
         effects = list())
  }
  pk_models <- list(`A conc` = model(c(cl = 5, v = 50, ka = 1.2)),
                    `B conc` = model(c(cl = 20, v = 200, ka = 0.5)))
  dose_rows <- .dose_rows(data, roles)
  drives <- .dose_endpoint(data, roles)
  .pmx_fitted_model(
    structural = "1cmt_oral",
    candidates = data.frame(model = "1cmt_oral", converged = TRUE, aic = 1,
                            note = "", stringsAsFactors = FALSE),
    parameters = pk_models[[1L]]$parameters,
    pk_models = pk_models,
    endpoints = list(pk = classified$pk, pd = character(),
                     discrete = character(), decided_by = "declared"),
    arms = list(arms = arm_models$arms, sizes = arm_models$sizes),
    dosing = arm_models$dosing, visits = arm_models$visits,
    schema = .source_schema(data, roles, classified$pk, subject_group),
    roles = roles, settings = list(min_arm_patients = 3L),
    n_source = length(subject_group), cells = cells,
    dose_records = vapply(stats::setNames(classified$pk, classified$pk),
                          function(endpoint) {
                            sum(if (is.null(roles$dose_endpoints)) dose_rows else
                              dose_rows & !is.na(drives) & drives == endpoint)
                          }, integer(1)),
    covariates = .covariate_model(data, roles),
    discrete = .discrete_model(data, roles, cells, subject_group, character())
  )
}

test_that("each endpoint is fitted against its own drug's doses", {
  data <- .combination_fixture()
  roles <- .combination_roles(
    dose_endpoints = c("1" = "A conc", "2" = "B conc"))
  a <- .model_estimation_data(data, roles, "A conc")
  b <- .model_estimation_data(data, roles, "B conc")
  expect_setequal(a$AMT[a$EVID != 0L], 100)
  expect_setequal(b$AMT[b$EVID != 0L], 900)
  # Undeclared, both fits see both drugs, which is the parent-and-metabolite
  # reading and the reason the declaration exists.
  shared <- .combination_roles()
  pooled <- .model_estimation_data(data, shared, "A conc")
  expect_setequal(pooled$AMT[pooled$EVID != 0L], c(100, 900))
})

test_that("a declared combination generates both drugs' dose records", {
  data <- .combination_fixture()
  roles <- .combination_roles(
    dose_endpoints = c("1" = "A conc", "2" = "B conc"))
  fit <- .combination_fit(data, roles)
  synthetic <- synpmx_model_generate(fit, n_subjects = 8, seed = 5)
  doses <- synthetic[synthetic$EVID != 0L, , drop = FALSE]
  expect_setequal(as.integer(doses$ADM), c(1L, 2L))
  # Each drug keeps its own amount, its own schedule and its own compartment.
  expect_setequal(doses$AMT[doses$ADM == 1L], 100)
  expect_setequal(doses$AMT[doses$ADM == 2L], 900)
  expect_setequal(doses$CMT[doses$ADM == 1L], 1L)
  expect_setequal(doses$CMT[doses$ADM == 2L], 3L)
  expect_gt(sum(doses$ADM == 1L), sum(doses$ADM == 2L))
})

# The whole point of the split: drug B's dose can move without moving drug A's
# concentration. Pooled, it cannot -- both endpoints read one schedule.
test_that("a drug's profile answers to its own doses only", {
  data <- .combination_fixture()
  roles <- .combination_roles(
    dose_endpoints = c("1" = "A conc", "2" = "B conc"))
  fit <- .combination_fit(data, roles)
  base <- synpmx_model_generate(fit, n_subjects = 8, seed = 5)

  heavier <- fit
  for (arm in names(heavier$dosing)) {
    heavier$dosing[[arm]][["B conc"]]$planned$amt <-
      heavier$dosing[[arm]][["B conc"]]$planned$amt * 10
  }
  moved <- synpmx_model_generate(heavier, n_subjects = 8, seed = 5)

  a_of <- function(d) d$DV[d$EVID == 0L & d$DVID == "A conc"]
  b_of <- function(d) d$DV[d$EVID == 0L & d$DVID == "B conc"]
  expect_equal(a_of(moved), a_of(base))
  expect_gt(mean(b_of(moved)), 5 * mean(b_of(base)))
})

# A NONMEM dataset has no administration column: the compartment a dose enters
# is the administration id, so `adm` and `cmt` may name one column.
test_that("a NONMEM study separates its drugs by compartment", {
  data <- .combination_fixture()
  roles <- pmx_roles(
    id = "ID", time = "TIME", nominal_time = "NTIME", dv = "DV", amt = "AMT",
    evid = "EVID", cmt = "CMT", dvid = "DVID", mdv = "MDV", covariates = "WT",
    adm = "CMT", routes = c("1" = "extravascular", "3" = "extravascular"),
    dose_endpoints = c("1" = "A conc", "3" = "B conc")
  )
  a <- .model_estimation_data(data, roles, "A conc")
  b <- .model_estimation_data(data, roles, "B conc")
  expect_setequal(a$AMT[a$EVID != 0L], 100)
  expect_setequal(b$AMT[b$EVID != 0L], 900)
})

test_that("a dose_endpoints mapping is checked against the study", {
  data <- .combination_fixture()
  expect_error(
    .model_check_dose_endpoints(
      data, .combination_roles(dose_endpoints = c("1" = "A conc")),
      c("A conc", "B conc")),
    "does not say which endpoint")
  expect_error(
    .model_check_dose_endpoints(
      data,
      .combination_roles(dose_endpoints = c("1" = "A conc", "2" = "A conc")),
      c("A conc", "B conc")),
    "no dose record")
  expect_error(
    .model_check_dose_endpoints(
      data,
      .combination_roles(dose_endpoints = c("1" = "A conc", "2" = "nobody")),
      c("A conc", "B conc")),
    "no structural model")
  expect_error(pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                         evid = "EVID", dose_endpoints = c("1" = "A conc")),
               "needs `adm`")
})

# The dosing history each observation is read against. Time after dose, the
# dose interval and the first dose amount feed the starting values, the route
# reading and the dose-proportionality signal, so a B sample measured from A's
# last dose puts every one of them on a mixture of two drugs.
test_that("an observation's dosing history is its own drug's", {
  data <- .combination_fixture()
  roles <- .combination_roles(
    dose_endpoints = c("1" = "A conc", "2" = "B conc"))
  obs <- .model_observations(data, roles)
  b <- obs[obs$endpoint == "B conc", , drop = FALSE]
  expect_setequal(b$first_dose_amt, 900)
  # Drug B is dosed on the first day of every other week, so a sample at 169 h
  # is 169 hours after its own drug's dose and 1 hour after the other's.
  expect_equal(b$tad[abs(b$time - 169) < 1e-8][[1L]], 169)
  a <- obs[obs$endpoint == "A conc", , drop = FALSE]
  expect_setequal(a$first_dose_amt, 100)
  expect_equal(a$tad[abs(a$time - 169) < 1e-8][[1L]], 1)

  # Undeclared, both endpoints read every dose, which is the reading a parent
  # and its metabolite need.
  pooled <- .model_observations(data, .combination_roles())
  pooled_b <- pooled[pooled$endpoint == "B conc", , drop = FALSE]
  expect_equal(pooled_b$tad[abs(pooled_b$time - 169) < 1e-8][[1L]], 1)
})

# The same trap `endpoint_types` fell into: an entry keyed by administration id
# with an ENDPOINT NAME for a value must not reach the places that unlist the
# roles object as a list of columns, or every reader of the study errors with
# "role columns not found" on the endpoint's own name.
test_that("`dose_endpoints` is not treated as a column name", {
  data <- .combination_fixture()
  roles <- .combination_roles(
    dose_endpoints = c("1" = "A conc", "2" = "B conc"))
  expect_true(.assert_roles(data, roles))
  expect_false(any(c("A conc", "B conc") %in% .retained_role_columns(roles)))
  expect_silent(validate_pmx(data, roles))
})
