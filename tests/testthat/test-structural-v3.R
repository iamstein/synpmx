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
  table <- pmx_generate(.v3_model(), .v3_design(), seed = 1)
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
  table <- pmx_generate(.v3_model(), design, seed = 3)
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
  table <- pmx_generate(model, design, seed = 5)
  pd <- table[table$EVID == 0 & table$DVID == "pd", ]
  expect_true(nrow(pd) > 0)
  expect_true(all(is.finite(pd$DV)))
  expect_true(all(pd$DV > 0))
})

test_that("censoring uses the Monolix convention", {
  table <- pmx_generate(.v3_model(), .v3_design(), seed = 1, lloq = 0.05)
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
  confidential <- pmx_generate(truth, design, n_subjects = 60, seed = 7)

  predicted <- .v3_model()                       # says cl = 10, 2.5x wrong
  priors <- pmx_priors(pk = pmx_prior(c(1 / 4, 4), "scaling literature"))
  fit <- fit_calibrated_pmx(confidential, pmx_generated_roles(), predicted,
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
  confidential <- pmx_generate(.v3_model(), design, n_subjects = 30, seed = 9)
  priors <- pmx_priors(pk = pmx_prior(c(1 / 4, 4), "x"))
  fit <- fit_calibrated_pmx(confidential, pmx_generated_roles(), .v3_model(),
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

test_that("generation from a calibrated model is post-processing", {
  skip_if_not(dp_backend_status()$available, "OpenDP unavailable")
  design <- .v3_design()
  confidential <- pmx_generate(.v3_model(), design, n_subjects = 30, seed = 9)
  priors <- pmx_priors(pk = pmx_prior(c(1 / 4, 4), "x"))
  fit <- fit_calibrated_pmx(confidential, pmx_generated_roles(), .v3_model(),
                            design, priors, epsilon = 1, backend = "opendp")

  before <- fit$privacy$accounting$realized_epsilon
  a <- pmx_generate(fit, seed = 1)
  b <- pmx_generate(fit, seed = 2)
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
  table <- pmx_generate(.v3_model(), design, n_subjects = 10, seed = 1)
  priors <- pmx_priors(pk = pmx_prior(c(1 / 4, 4), "x"))
  expect_error(
    fit_calibrated_pmx(table, pmx_generated_roles(), .v3_model(), design,
                       priors, epsilon = 1, backend = "public"),
    "public_source"
  )
  fit <- fit_calibrated_pmx(table, pmx_generated_roles(), .v3_model(), design,
                            priors, epsilon = 1, backend = "public",
                            public_source = TRUE)
  expect_false(fit$privacy$formal_dp)
})

test_that("a pd prior without a pd model is rejected", {
  design <- .v3_design()
  table <- pmx_generate(.v3_model(), design, n_subjects = 10, seed = 1)
  expect_error(
    fit_calibrated_pmx(
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
  data <- pmx_generate(truth, design, n_subjects = 400, seed = 21)

  priors <- pmx_priors(pk = pmx_prior(c(1 / 4, 4), "x"))
  # Noiseless backend isolates estimator bias from privacy noise.
  fit <- fit_calibrated_pmx(data, pmx_generated_roles(), .v3_model(), design,
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
    data <- pmx_generate(spec(100, cv), design, n_subjects = 200, seed = 6)
    fit <- suppressWarnings(fit_calibrated_pmx(
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
  data <- pmx_generate(spec(100), design, n_subjects = 150, seed = 8)
  fit <- suppressWarnings(fit_calibrated_pmx(
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

  table <- pmx_generate(model, design, n_subjects = 6, seed = 1)
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
  data <- pmx_generate(truth, design, n_subjects = 60, seed = 3)
  fit <- suppressWarnings(fit_calibrated_pmx(
    data, pmx_generated_roles(), pred, design,
    pmx_priors(pk = pmx_prior(c(1 / 4, 4), "x")), epsilon = 2, backend = "opendp"
  ))
  expect_gt(fit$corrected_typical[["cl"]], 14)
  expect_lt(fit$corrected_typical[["cl"]], 32)
})
