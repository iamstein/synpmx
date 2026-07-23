# Version 3: public structural model, correction release, generation.

.v3_model <- function(pd = "none") {
  typical <- c(cl = 10, v = 70, ka = 1)
  if (pd == "exponential") {
    typical <- c(typical, baseline = 100, plateau = 40, rate = 0.02)
  } else if (pd == "linear") {
    typical <- c(typical, baseline = 100, slope = 0.3)
  } else if (pd == "constant") {
    typical <- c(typical, baseline = 100)
  }
  pmx_structural_model("1cmt_oral", typical, pd = pd, source = "unit test")
}

.v3_design <- function(...) {
  pmx_trial_design(c(10, 30, 100), c(6, 6, 6), c(0, .5, 1, 2, 4, 8, 12, 24),
                   source = "unit test protocol", ...)
}

test_that("public inputs require provenance", {
  expect_error(
    pmx_structural_model("1cmt_iv", c(cl = 1, v = 1)),
    "source"
  )
  expect_error(pmx_prior(c(1 / 4, 4)), "source")
  expect_error(pmx_prior(c(4, 1 / 4), "x"), "increasing")
  expect_error(
    pmx_trial_design(10, 6, c(0, 1)),
    "source"
  )
  # A model missing a parameter its structure needs must fail loudly.
  expect_error(
    pmx_structural_model("1cmt_oral", c(cl = 1, v = 1), source = "x"),
    "missing required parameters"
  )
})

test_that("an unused rxode2 model warns instead of silently doing nothing", {
  skip_if_not_installed("rxode2")
  # REV-020: `rx` is accepted but never consumed by profile evaluation, so the
  # curve is the analytic one either way. The warning must say so.
  expect_warning(
    model <- pmx_structural_model(
      "1cmt_oral", c(cl = 10, v = 70, ka = 1), source = "unit test",
      rx = "d/dt(central) = -ke * central"
    ),
    "not yet used"
  )
  expect_equal(
    .pk_profile(model, c(0, 1, 4), 100, 0),
    .pk_profile(.v3_model(), c(0, 1, 4), 100, 0)
  )
})

test_that("analytic PK is dose-proportional and accumulates", {
  model <- .v3_model()
  t <- c(0, 1, 2, 4, 8, 24)
  one <- .pk_profile(model, t, 10, 0)
  ten <- .pk_profile(model, t, 100, 0)
  # Linear PK: ten times the dose is exactly ten times the concentration.
  expect_equal(ten, one * 10, tolerance = 1e-10)

  # Accumulation: one hour after the second dose must exceed one hour after the
  # first, because residual drug from dose one is still present. Comparing at
  # exactly the dose time would show no difference, since an oral dose
  # contributes nothing at zero elapsed time.
  first_peak <- .pk_profile(model, 1, 100, 0)
  second_peak <- .pk_profile(model, 25, rep(100, 2), c(0, 24))
  expect_gt(second_peak, first_peak)

  # Concentration must be zero before the dose and nonnegative throughout.
  expect_equal(.pk_profile(model, 0, 100, 5), 0)
  expect_true(all(.pk_profile(model, seq(0, 48, 2), 100, 0) >= 0))
})

test_that("the ka == ke singularity is handled", {
  model <- pmx_structural_model("1cmt_oral", c(cl = 10, v = 10, ka = 1),
                                source = "x")   # ke = 1 = ka
  values <- .pk_profile(model, c(0, 1, 4, 12), 100, 0)
  expect_true(all(is.finite(values)))
  expect_gt(values[2L], 0)
})

test_that("infusion rises to a plateau then declines", {
  model <- pmx_structural_model("1cmt_infusion", c(cl = 10, v = 70),
                                source = "x")
  during <- .pk_profile(model, c(1, 2, 3), 100, 0, duration = 4)
  after <- .pk_profile(model, c(5, 8, 12), 100, 0, duration = 4)
  expect_true(all(diff(during) > 0))
  expect_true(all(diff(after) < 0))
})

test_that("prior-mode generation makes a valid PMX table with no budget", {
  table <- .generate_structural(.v3_model(), .v3_design(), seed = 1)
  expect_identical(attr(table, "pmx_source"), "prior")
  expect_true(validate_pmx(table, pmx_generated_roles())$valid)
  expect_equal(length(unique(table$ID)), 18L)

  observations <- table[table$EVID == 0 & table$DVID == "cp", ]
  cmax <- tapply(observations$DV, observations$DOSE, max)
  # Higher dose must give higher exposure. This is structural, not learned.
  expect_true(all(diff(cmax) > 0))

  # Timing invariants: TAD is never negative and never reaches the next dose.
  expect_true(all(table$TAD >= 0))
  expect_true(all(table$TIME >= 0))
})

test_that("repeated dosing keeps every sample inside its own occasion", {
  design <- pmx_trial_design(100, 8, c(0, 1, 4, 12), n_doses = 4,
                             dose_interval = 24, visit_window = 0.2,
                             source = "x")
  table <- .generate_structural(.v3_model(), design, seed = 3)
  expect_true(validate_pmx(table, pmx_generated_roles())$valid)
  expect_true(all(table$TAD >= 0))
  # SIM-005: a dose-relative time must stay strictly before the next dose.
  expect_true(all(table$TAD < 24))
  expect_setequal(unique(table$OCC), 1:4)
})

test_that("the PD endpoint produces a coherent time course", {
  model <- .v3_model("exponential")
  design <- pmx_trial_design(100, 10, c(0, 2, 8, 24, 48, 96),
                             source = "x")
  table <- .generate_structural(model, design, seed = 5)
  pd <- table[table$EVID == 0 & table$DVID == "pd", ]
  expect_true(nrow(pd) > 0)
  expect_true(all(is.finite(pd$DV)))
  expect_true(all(pd$DV > 0))
})

test_that("censoring uses the Monolix convention", {
  table <- .generate_structural(.v3_model(), .v3_design(), seed = 1, lloq = 0.05)
  censored <- table[table$CENS == 1L, ]
  expect_true(nrow(censored) > 0)
  expect_true(all(censored$DV == 0.05))
  expect_true(all(table$DV[table$EVID == 0 & table$CENS == 0 &
                             table$DVID == "cp"] >= 0.05 - 1e-12))
})

test_that("pre-flight reports f and refuses a worthless configuration", {
  priors <- pmx_priors(pk = pmx_prior(c(1 / 4, 4), "x"))
  good <- pmx_preflight(priors, epsilon = 0.5, n_subjects = 20)
  expect_equal(good$d, 2L)
  expect_equal(good$f, 2 / (0.5 * 20))
  expect_equal(good$verdict, "worthwhile")

  # Tiny cohort, small epsilon: the noise is as wide as the prior.
  bad <- pmx_preflight(priors, epsilon = 0.1, n_subjects = 6)
  expect_gte(bad$f, 1)
  expect_equal(bad$verdict, "worthless")

  # f is a fraction of the prior width, so a wider prior gives the same f but
  # a worse fold-error.
  wide <- pmx_preflight(pmx_priors(pk = pmx_prior(c(1 / 100, 100), "x")),
                        epsilon = 0.5, n_subjects = 20)
  expect_equal(wide$f, good$f)
  expect_gt(wide$table$expected_fold_error, good$table$expected_fold_error)
})

test_that("calibration recovers a known correction", {
  skip_if_not(dp_backend_status()$available, "OpenDP unavailable")
  truth <- pmx_structural_model("1cmt_oral", c(cl = 25, v = 70, ka = 1),
                                source = "simulated truth")
  design <- .v3_design()
  confidential <- .generate_structural(truth, design, n_subjects = 60, seed = 7)

  predicted <- .v3_model()                       # says cl = 10, 2.5x wrong
  priors <- pmx_priors(pk = pmx_prior(c(1 / 4, 4), "scaling literature"))
  fit <- .fit_calibrated(confidential, pmx_generated_roles(), predicted,
                            design, priors, epsilon = 2, backend = "opendp")

  expect_s3_class(fit, "pmx_calibrated_model")
  expect_true(fit$privacy$formal_dp)
  # The correction should land near 2.5x. The tolerance is wide because this is
  # a privacy mechanism, not an estimator.
  expect_gt(fit$corrections$pk$factor, 1.5)
  expect_lt(fit$corrections$pk$factor, 4)
  expect_gt(fit$corrected_typical[["cl"]], 15)
  expect_lt(fit$corrected_typical[["cl"]], 40)
})

test_that("accounting never exceeds the requested budget", {
  skip_if_not(dp_backend_status()$available, "OpenDP unavailable")
  design <- .v3_design()
  confidential <- .generate_structural(.v3_model(), design, n_subjects = 30, seed = 9)
  priors <- pmx_priors(pk = pmx_prior(c(1 / 4, 4), "x"))
  fit <- .fit_calibrated(confidential, pmx_generated_roles(), .v3_model(),
                            design, priors, epsilon = 1, backend = "opendp")
  entries <- fit$privacy$accounting$entries
  expect_equal(nrow(entries), 2L)                # count + one correction
  expect_lte(fit$privacy$accounting$realized_epsilon, 1)
  expect_true(all(entries$epsilon > 0))
  expect_true(all(entries$sensitivity == 1))

  # The fitted object must carry no raw records.
  expect_false(any(grepl("raw|donor|template|profile",
                         tolower(.recursive_names(fit)))))
})

test_that("calibrated diagnostics do not retain exact usable-subject counts", {
  design <- .v3_design()
  data <- .generate_structural(.v3_model(), design, n_subjects = 20, seed = 9)
  priors <- pmx_priors(pk = pmx_prior(c(1 / 4, 4), "x"))
  fit <- .fit_calibrated(
    data, pmx_generated_roles(), .v3_model(), design, priors,
    epsilon = 1, backend = "public", public_source = TRUE
  )

  expect_false("n_used" %in% names(fit$corrections$pk))
  expect_equal(fit$preflight$d, 2)
  expect_false(any(grepl("n_used", .recursive_names(fit))))
})

test_that("generation from a calibrated model is post-processing", {
  skip_if_not(dp_backend_status()$available, "OpenDP unavailable")
  design <- .v3_design()
  confidential <- .generate_structural(.v3_model(), design, n_subjects = 30, seed = 9)
  priors <- pmx_priors(pk = pmx_prior(c(1 / 4, 4), "x"))
  fit <- .fit_calibrated(confidential, pmx_generated_roles(), .v3_model(),
                            design, priors, epsilon = 1, backend = "opendp")

  before <- fit$privacy$accounting$realized_epsilon
  a <- .generate_structural(fit, seed = 1)
  b <- .generate_structural(fit, seed = 2)
  expect_equal(fit$privacy$accounting$realized_epsilon, before)
  expect_identical(attr(a, "pmx_source"), "calibrated")
  expect_true(validate_pmx(a, pmx_generated_roles())$valid)
  expect_false(isTRUE(all.equal(a$DV, b$DV)))

  # Generated times must never reproduce a source subject's exact vector.
  source_times <- split(confidential$TIME, confidential$ID)
  generated_times <- split(a$TIME, a$ID)
  expect_false(any(vapply(generated_times, function(g) {
    any(vapply(source_times, function(s) isTRUE(all.equal(g, s)), logical(1)))
  }, logical(1))))
})

test_that("the public backend is refused for undeclared sources", {
  design <- .v3_design()
  table <- .generate_structural(.v3_model(), design, n_subjects = 10, seed = 1)
  priors <- pmx_priors(pk = pmx_prior(c(1 / 4, 4), "x"))
  expect_error(
    .fit_calibrated(table, pmx_generated_roles(), .v3_model(), design,
                       priors, epsilon = 1, backend = "public"),
    "public_source"
  )
  fit <- .fit_calibrated(table, pmx_generated_roles(), .v3_model(), design,
                            priors, epsilon = 1, backend = "public",
                            public_source = TRUE)
  expect_false(fit$privacy$formal_dp)
})

test_that("a pd prior without a pd model is rejected", {
  design <- .v3_design()
  table <- .generate_structural(.v3_model(), design, n_subjects = 10, seed = 1)
  expect_error(
    .fit_calibrated(
      table, pmx_generated_roles(), .v3_model(), design,
      pmx_priors(pk = pmx_prior(c(1 / 4, 4), "x"),
                 pd = pmx_prior(c(1 / 4, 4), "x")),
      epsilon = 1, backend = "public", public_source = TRUE
    ),
    "no PD component"
  )
})

test_that("fitting and generating agree on what `typical` means", {
  # `typical` is the median of a lognormal parameter, and the released
  # correction is a mean on the log scale, so it targets the same quantity.
  # Centering the arithmetic mean when generating would leave a systematic
  # exp(sigma^2 / 2) gap that no amount of N or epsilon removes.
  design <- .v3_design()
  truth <- pmx_structural_model("1cmt_oral", c(cl = 12, v = 70, ka = 1),
                                source = "round-trip truth")
  data <- .generate_structural(truth, design, n_subjects = 400, seed = 21)

  priors <- pmx_priors(pk = pmx_prior(c(1 / 4, 4), "x"))
  # Noiseless backend isolates estimator bias from privacy noise.
  fit <- .fit_calibrated(data, pmx_generated_roles(), .v3_model(), design,
                            priors, epsilon = 1, backend = "public",
                            public_source = TRUE)
  recovered <- fit$corrected_typical[["cl"]]
  expect_gt(recovered, 12 / 1.15)
  expect_lt(recovered, 12 * 1.15)
})

test_that("median-centred variability keeps the population median on target", {
  set.seed(4)
  drawn <- .draw_subject_params(c(cl = 10), c(cl = 0.4), 20000)
  expect_equal(stats::median(drawn[, "cl"]), 10, tolerance = 0.03)
  # Arithmetic mean is deliberately above the median for a lognormal.
  expect_gt(mean(drawn[, "cl"]), 10)
})

test_that("pre-flight caps the fold-error at the prior half-width", {
  priors <- pmx_priors(pk = pmx_prior(c(1 / 4, 4), "x"))
  # f = 1.33 would give exp(1.33 * 2.08) = 16x uncapped, but clipping means a
  # release can never land outside the prior.
  hopeless <- pmx_preflight(priors, epsilon = 0.25, n_subjects = 6)
  expect_gt(hopeless$f, 1)
  expect_lte(hopeless$table$expected_fold_error,
             exp(priors$pk$span / 2) + 1e-9)
})

test_that("two-compartment PK is biphasic and conserves Dose/CL", {
  model <- pmx_structural_model("2cmt_iv", c(cl = 10, v = 30, q = 15, v2 = 100),
                                source = "x")
  grid <- seq(0, 400, by = 0.05)
  profile <- .pk_profile(model, grid, 100, 0)
  expect_equal(.trapezoid(grid, profile), 100 / 10, tolerance = 1e-3)
  expect_true(all(diff(profile) < 0))

  # Distribution phase: early decline is steeper than terminal decline.
  early <- .pk_profile(model, c(0, 1), 100, 0)
  late <- .pk_profile(model, c(24, 25), 100, 0)
  expect_gt(log(early[1] / early[2]), log(late[1] / late[2]))

  oral <- pmx_structural_model("2cmt_oral",
                               c(cl = 10, v = 30, q = 15, v2 = 100, ka = 1),
                               source = "x")
  values <- .pk_profile(oral, c(0, 0.5, 1, 2, 8, 24), 100, 0)
  expect_equal(values[1L], 0)
  expect_true(all(values >= 0))
  expect_gt(max(values), 0)
})

test_that("simple PD shapes are closed-form and go the right way", {
  spec <- function(pd, extra) pmx_structural_model(
    "1cmt_oral", c(cl = 10, v = 70, ka = 1, baseline = 100, extra),
    pd = pd, source = "x"
  )
  grid <- c(0, 6, 24, 48, 96, 168)

  flat <- .pd_profile(spec("constant", NULL), grid, 100, 0)
  expect_true(all(flat == 100))

  up <- .pd_profile(spec("linear", c(slope = 0.3)), grid, 100, 0)
  expect_true(all(diff(up) > 0))
  expect_equal(up[1L], 100)

  decay <- .pd_profile(spec("exponential", c(plateau = 40, rate = 0.02)),
                       grid, 100, 0)
  expect_true(all(diff(decay) < 0))
  expect_equal(decay[1L], 100)
  expect_gt(decay[length(decay)], 40)          # approaches, never crosses

  # The same shape covers growth when the plateau is above the baseline.
  growth <- .pd_profile(spec("exponential", c(plateau = 300, rate = 0.02)),
                        grid, 100, 0)
  expect_true(all(diff(growth) > 0))

  # A simple shape must not depend on dose.
  low <- .pd_profile(spec("exponential", c(plateau = 40, rate = 0.02)),
                     grid, 10, 0)
  expect_equal(low, decay)
})

test_that("the simple PD level correction survives residual error", {
  skip_if_not(dp_backend_status()$available, "OpenDP unavailable")
  design <- pmx_trial_design(c(30, 100), c(1, 1), c(0, 6, 24, 48, 96, 168),
                             source = "x")
  spec <- function(baseline, cv) pmx_structural_model(
    "1cmt_oral",
    c(cl = 10, v = 70, ka = 1, baseline = baseline,
      plateau = baseline * 0.4, rate = 0.02),
    pd = "exponential", source = "x", residual_cv = cv,
    iiv = c(cl = 0.3, v = 0.2, baseline = 0.25)
  )
  priors <- pmx_priors(pk = pmx_prior(c(1 / 4, 4), "x"),
                       pd = pmx_prior(c(1 / 10, 10), "x"))
  # A ratio of means, unlike the exposure-driven deviation statistic, is not
  # meaningfully biased by residual error.
  for (cv in c(0, 0.15, 0.3)) {
    data <- .generate_structural(spec(100, cv), design, n_subjects = 200, seed = 6)
    fit <- suppressWarnings(.fit_calibrated(
      data, pmx_generated_roles(), spec(40, cv), design, priors,
      epsilon = 50, backend = "opendp"
    ))
    expect_gt(fit$corrections$pd$factor, 2.2)
    expect_lt(fit$corrections$pd$factor, 2.9)
  }
})

test_that("the simple PD correction scales the whole curve", {
  skip_if_not(dp_backend_status()$available, "OpenDP unavailable")
  design <- pmx_trial_design(100, 10, c(0, 6, 24, 48, 96), source = "x")
  spec <- function(baseline) pmx_structural_model(
    "1cmt_oral",
    c(cl = 10, v = 70, ka = 1, baseline = baseline,
      plateau = baseline * 0.4, rate = 0.02),
    pd = "exponential", source = "x"
  )
  data <- .generate_structural(spec(100), design, n_subjects = 150, seed = 8)
  fit <- suppressWarnings(.fit_calibrated(
    data, pmx_generated_roles(), spec(40), design,
    pmx_priors(pk = pmx_prior(c(1 / 4, 4), "x"),
               pd = pmx_prior(c(1 / 10, 10), "x")),
    epsilon = 50, backend = "opendp"
  ))
  # Baseline and plateau move together, so the shape is preserved and only the
  # level is calibrated.
  ratio <- fit$corrected_typical[["plateau"]] /
    fit$corrected_typical[["baseline"]]
  expect_equal(ratio, 0.4, tolerance = 1e-8)
  expect_gt(fit$corrected_typical[["baseline"]], 80)
  expect_lt(fit$corrected_typical[["baseline"]], 125)
})

test_that("within-subject dose escalation increases exposure by occasion", {
  model <- pmx_structural_model("1cmt_oral", c(cl = 10, v = 70, ka = 1),
                                source = "x")
  design <- pmx_trial_design(
    dose_escalation = c(10, 30, 100), dose_times = c(0, 7, 14),
    sampling = c(0, 1, 4, 24, 72, 167), source = "x"
  )
  expect_equal(.design_dose_times(design), c(0, 7, 14))
  expect_equal(.design_dose_amounts(design), c(10, 30, 100))

  table <- .generate_structural(model, design, n_subjects = 6, seed = 1)
  expect_true(validate_pmx(table, pmx_generated_roles())$valid)

  # AMT and the assigned-dose column both follow the escalation by occasion.
  doses <- table[table$EVID == 1, ]
  expect_equal(as.numeric(tapply(doses$AMT, doses$OCC, unique)),
               c(10, 30, 100))
  obs <- table[table$EVID == 0, ]
  expect_equal(as.numeric(tapply(obs$DOSE, obs$OCC, unique)), c(10, 30, 100))

  # Exposure rises with dose across occasions.
  cmax <- tapply(obs$DV, obs$OCC, max)
  expect_true(all(diff(cmax) > 0))
})

test_that("escalation and parallel design are mutually exclusive", {
  expect_error(
    pmx_trial_design(dose_levels = 10, dose_escalation = c(10, 30),
                     sampling = c(0, 1), source = "x"),
    "not both"
  )
  expect_error(
    pmx_trial_design(dose_escalation = c(10, 30), dose_times = c(0, 7, 14),
                     sampling = c(0, 1), source = "x"),
    "must give 2"
  )
})

test_that("calibration works on an escalating regimen", {
  skip_if_not(dp_backend_status()$available, "OpenDP unavailable")
  design <- pmx_trial_design(dose_escalation = c(10, 30, 100),
                             dose_times = c(0, 7, 14),
                             sampling = c(0, 1, 4, 24, 72, 167), source = "x")
  truth <- pmx_structural_model("1cmt_oral", c(cl = 22, v = 70, ka = 1),
                                source = "truth")
  pred <- pmx_structural_model("1cmt_oral", c(cl = 10, v = 70, ka = 1),
                               source = "pred")
  data <- .generate_structural(truth, design, n_subjects = 60, seed = 3)
  fit <- suppressWarnings(.fit_calibrated(
    data, pmx_generated_roles(), pred, design,
    pmx_priors(pk = pmx_prior(c(1 / 4, 4), "x")), epsilon = 2, backend = "opendp"
  ))
  expect_gt(fit$corrected_typical[["cl"]], 14)
  expect_lt(fit$corrected_typical[["cl"]], 32)
})

test_that("per-cohort escalation gives each cohort its own sequence", {
  model <- pmx_structural_model("1cmt_oral", c(cl = 10, v = 70, ka = 1),
                                source = "x")
  design <- pmx_trial_design(
    dose_escalation = list(c(1, 2, 4), c(2, 4, 8), c(4, 8, 16)),
    cohort_sizes = c(6, 6, 6), dose_times = c(0, 7, 14),
    sampling = c(0, 1, 4, 24, 72, 167), source = "x"
  )
  expect_length(design$escalation, 3L)

  table <- .generate_structural(model, design, seed = 1)
  expect_true(validate_pmx(table, pmx_generated_roles())$valid)
  expect_equal(length(unique(table$ID)), 18L)

  # Exactly the three protocol sequences appear, six subjects each.
  doses <- table[table$EVID == 1, ]
  seqs <- tapply(doses$AMT, doses$ID, function(a) paste(a, collapse = "-"))
  expect_setequal(unique(seqs), c("1-2-4", "2-4-8", "4-8-16"))
  expect_true(all(table(seqs) == 6))
})

test_that("escalation sequences must share a length", {
  expect_error(
    pmx_trial_design(dose_escalation = list(c(1, 2, 4), c(2, 4)),
                     dose_times = c(0, 7, 14), sampling = c(0, 1), source = "x"),
    "same number of doses"
  )
})

test_that("covariate declarations validate their inputs", {
  expect_error(pmx_covariate(range = c(40, 120)), "source")
  expect_error(pmx_covariate(source = "x"), "exactly one")
  expect_error(pmx_covariate(range = c(1, 2), levels = c("a"), source = "x"),
               "exactly one")
  expect_error(pmx_covariate(range = c(120, 40), source = "x"), "increasing")
  expect_error(pmx_covariate(levels = c("a", "a"), source = "x"), "unique")
  expect_error(pmx_covariates(pmx_covariate(range = c(1, 2), source = "x")),
               "named")
})

test_that("prior-mode covariates appear, constant within subject, in range", {
  covariates <- pmx_covariates(
    WT  = pmx_covariate(range = c(40, 120), source = "x"),
    SEX = pmx_covariate(levels = c("M", "F"), source = "x")
  )
  table <- .generate_structural(.v3_model(), .v3_design(), seed = 1,
                        covariates = covariates)
  expect_true(all(c("WT", "SEX") %in% names(table)))
  expect_true(validate_pmx(table, pmx_generated_roles())$valid)

  # One value per subject, held across their rows.
  expect_true(all(tapply(table$WT, table$ID,
                         function(v) length(unique(v))) == 1L))
  expect_true(all(table$WT >= 40 & table$WT <= 120))
  expect_true(all(table$SEX %in% c("M", "F")))
})

test_that("each covariate adds exactly one budget slice with sensitivity one", {
  skip_if_not(dp_backend_status()$available, "OpenDP unavailable")
  design <- .v3_design()
  data <- .generate_structural(.v3_model(), design, n_subjects = 60, seed = 3)
  data$WT <- ave(data$ID, data$ID, FUN = function(i) 70 + i[1L] %% 30)
  data$SEX <- ifelse(data$ID %% 2 == 0, "M", "F")

  covariates <- pmx_covariates(
    WT  = pmx_covariate(range = c(40, 120), source = "x"),
    SEX = pmx_covariate(levels = c("M", "F"), source = "x")
  )
  fit <- .fit_calibrated(
    data, pmx_generated_roles(), .v3_model(), design,
    pmx_priors(pk = pmx_prior(c(1 / 4, 4), "x")),
    epsilon = 2, covariates = covariates, backend = "opendp"
  )
  entries <- fit$privacy$accounting$entries
  expect_equal(nrow(entries), 4L)               # count + pk + WT + SEX
  expect_true(all(entries$sensitivity == 1))
  expect_lte(fit$privacy$accounting$realized_epsilon, 2)

  # Generation from the fit carries the columns without further budget.
  synthetic <- .generate_structural(fit, seed = 11)
  expect_true(all(c("WT", "SEX") %in% names(synthetic)))
  expect_true(validate_pmx(synthetic, pmx_generated_roles())$valid)
})

test_that("pmx_covariates_auto validates and defaults to 1-99 clipping", {
  cov <- pmx_covariates_auto(c("WT", "SEX"))
  expect_s3_class(cov, "pmx_covariates")
  expect_equal(cov$WT$type, "bootstrap")
  expect_equal(cov$WT$clip, c(0.01, 0.99))
  expect_true(.covariates_have_bootstrap(cov))

  expect_null(pmx_covariates_auto("WT", clip = NULL)$WT$clip)
  expect_error(pmx_covariates_auto(character()), "non-empty")
  expect_error(pmx_covariates_auto(c("a", "a")), "unique")
  expect_error(pmx_covariates_auto("WT", clip = c(0.9, 0.1)), "increasing")
})

test_that("bootstrap covariates resample from data without spending budget", {
  skip_if_not(dp_backend_status()$available, "OpenDP unavailable")
  design <- .v3_design()
  data <- .generate_structural(.v3_model(), design, n_subjects = 80, seed = 3)
  # A continuous column with a clear range, and a skewed categorical.
  per_id <- function(f) f(length(unique(data$ID)))[match(data$ID,
                                                         unique(data$ID))]
  set.seed(2)
  data$EGFR <- per_id(function(k) round(rnorm(k, 90, 20)))
  data$RACE <- per_id(function(k) sample(c("A", "B", "C"), k, replace = TRUE,
                                         prob = c(0.7, 0.2, 0.1)))

  cov <- pmx_covariates_auto(c("EGFR", "RACE"))
  # Bootstrap covariates warn that they are not DP; that is expected here.
  expect_warning(
    fit <- .fit_calibrated(
      data, pmx_generated_roles(), .v3_model(), design,
      pmx_priors(pk = pmx_prior(c(1 / 4, 4), "x")),
      epsilon = 1, covariates = cov, backend = "opendp"
    ),
    "differentially private"
  )
  # Bootstrap covariates are NOT budgeted: only count + pk are released.
  expect_equal(nrow(fit$privacy$accounting$entries), 2L)
  expect_false(fit$privacy$covariates_private)

  synthetic <- .generate_structural(fit, seed = 5)
  expect_true(all(c("EGFR", "RACE") %in% names(synthetic)))
  expect_true(validate_pmx(synthetic, pmx_generated_roles())$valid)

  # Continuous values stay inside the clipped source range; the extremes are
  # trimmed by the default 1st/99th-percentile clip.
  src_egfr <- tapply(data$EGFR, data$ID, unique)
  expect_gte(min(synthetic$EGFR), min(src_egfr))
  expect_lte(max(synthetic$EGFR), max(src_egfr))

  # Categorical proportions are broadly preserved (the majority stays majority).
  expect_equal(names(which.max(table(synthetic$RACE))), "A")
})

test_that("clip = NULL exposes the raw min and max, synadam-style", {
  summary <- .bootstrap_summary(c(1, 2, 3, 4, 100), clip = NULL)
  expect_equal(summary$range, c(1, 100))
  clipped <- .bootstrap_summary(c(1, 2, 3, 4, 100), clip = c(0.01, 0.99))
  expect_lt(clipped$range[2L], 100)             # the 100 outlier is trimmed
})

test_that("a bootstrap covariate cannot be generated in prior mode", {
  expect_error(
    .generate_structural(.v3_model(), .v3_design(),
                 covariates = pmx_covariates_auto("WT")),
    "needs the data"
  )
})

# Per-occasion sampling ------------------------------------------------------
#
# A bare `sampling` vector is applied after every dose, which is right for a
# uniform protocol and wrong for most real ones: repeated-dose studies commonly
# sample richly on the first and last days and take troughs in between.

test_that("a bare sampling vector still repeats after every dose", {
  design <- pmx_trial_design(100, 6, c(0, 1, 4), n_doses = 3,
                             dose_interval = 24, source = "unit test")
  expect_length(design$sampling, 3L)
  expect_true(all(vapply(design$sampling, identical, logical(1), c(0, 1, 4))))
  # 3 doses x 3 times, none coinciding.
  expect_length(.design_observation_times(design), 9L)
})

test_that("a sampling list gives one schedule per dose", {
  rich <- c(0, 1, 2, 4, 12)
  design <- pmx_trial_design(
    100, 6, sampling = list(rich, 0, NULL, rich), n_doses = 4,
    dose_interval = 24, source = "unit test"
  )
  times <- .design_observation_times(design)
  # Dose 1 rich (0..12), dose 2 trough at 24, dose 3 nothing, dose 4 rich
  # (72..84). Nothing coincides, so every declared time survives.
  expect_length(times, length(rich) * 2L + 1L)
  expect_true(all(times[times > 12 & times < 72] == 24))
  # The gap after dose 3 is real: no observation between 24 and 72.
  expect_false(any(times > 24 & times < 72))
})

test_that("NULL and numeric(0) both mean no samples after that dose", {
  rich <- c(0, 1, 2)
  a <- pmx_trial_design(100, 6, list(rich, NULL, rich), n_doses = 3,
                        dose_interval = 24, source = "unit test")
  b <- pmx_trial_design(100, 6, list(rich, numeric(0), rich), n_doses = 3,
                        dose_interval = 24, source = "unit test")
  expect_identical(.design_observation_times(a), .design_observation_times(b))
})

test_that("a sampling list must match the dose count", {
  expect_error(
    pmx_trial_design(100, 6, list(c(0, 1), c(0, 1)), n_doses = 3,
                     source = "unit test"),
    "one element per dose"
  )
})

test_that("a design with fewer than two observation times is rejected", {
  expect_error(
    pmx_trial_design(100, 6, list(0, NULL, NULL), n_doses = 3,
                     source = "unit test"),
    "at least two observation times"
  )
})

test_that("per-occasion sampling reaches the generated data", {
  rich <- c(0, 0.5, 1, 2, 4, 9, 24)
  design <- pmx_trial_design(
    320, 12, sampling = list(rich, 0, NULL, NULL, NULL, NULL, rich),
    n_doses = 7, dose_interval = 24, source = "unit test"
  )
  model <- pmx_structural_model("1cmt_oral", c(cl = 6, v = 35, ka = 1.5),
                                source = "unit test")
  generated <- synpmx_prior(model, design, n_subjects = 6, seed = 1)
  observed <- generated[generated$EVID == 0, ]
  per_subject <- nrow(observed) / 6
  # A bare vector would give 7 profiles; this protocol declares two plus a
  # trough, so the study is sampled roughly three times less densely.
  expect_lt(per_subject, 20)
  # Doses 3 to 6 span (48, 144): the protocol says nothing is collected there.
  middle <- observed$NTIME > 50 & observed$NTIME < 140
  expect_false(any(middle))
})
