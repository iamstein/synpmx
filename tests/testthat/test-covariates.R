# Public covariate declarations and the reference catalogue -------------------
#
# `pmx_covariate()` separates two things the first version conflated: `range`
# is the clipping bound that makes a private release possible, and `median`
# plus `cv` say where the population sits inside it. Before the split, the
# spread was `diff(range) / 6`, so widening a safety bound widened the
# distribution and every covariate was symmetric.

.cov_model <- function() {
  pmx_structural_model("1cmt_oral", c(cl = 6, v = 35, ka = 1.5),
                       source = "unit test")
}

.cov_design <- function(n) {
  pmx_trial_design(320, n, c(0, 1, 4), source = "unit test protocol")
}

.cov_roles <- function(names) {
  pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
            cmt = "CMT", covariates = names)
}

# One value per subject, which is what a baseline covariate is.
.cov_draw <- function(covariates, n = 2000, seed = 1) {
  data <- synpmx_prior(.cov_model(), .cov_design(n), .cov_roles(names(covariates)),
                       n_subjects = n, seed = seed, covariates = covariates)
  data[!duplicated(data$ID), names(covariates), drop = FALSE]
}

test_that("a stated median and CV are what comes out, not the range midpoint", {
  drawn <- .cov_draw(pmx_covariates(
    WT = pmx_covariate(range = c(35, 160), median = 75, cv = 0.18,
                       source = "unit test")
  ))
  expect_equal(stats::median(drawn$WT), 75, tolerance = 0.05)
  expect_equal(stats::sd(drawn$WT) / mean(drawn$WT), 0.18, tolerance = 0.1)
  # Lognormal by default, so the tail runs to the right rather than symmetric.
  skew <- mean((drawn$WT - mean(drawn$WT))^3) / stats::sd(drawn$WT)^3
  expect_gt(skew, 0.2)
})

test_that("the range no longer sets the spread", {
  # Same population, a bound twice as wide. Before the split this doubled the
  # standard deviation and moved the centre from 75 to 97.5.
  tight <- .cov_draw(pmx_covariates(
    WT = pmx_covariate(range = c(50, 110), median = 75, cv = 0.18,
                       source = "unit test")
  ))
  wide <- .cov_draw(pmx_covariates(
    WT = pmx_covariate(range = c(35, 160), median = 75, cv = 0.18,
                       source = "unit test")
  ))
  expect_equal(stats::median(tight$WT), stats::median(wide$WT),
               tolerance = 0.05)
})

test_that("a range with no stated distribution keeps the old behaviour", {
  drawn <- .cov_draw(pmx_covariates(
    WT = pmx_covariate(range = c(40, 130), source = "unit test")
  ))
  expect_equal(mean(drawn$WT), 85, tolerance = 0.05)
  expect_equal(stats::sd(drawn$WT), 15, tolerance = 0.1)
})

test_that("level probabilities are honoured, and default to uniform", {
  drawn <- .cov_draw(pmx_covariates(
    SEX = pmx_covariate(levels = c("M", "F"), prob = c(0.7, 0.3),
                        source = "unit test"),
    ARM = pmx_covariate(levels = c("a", "b", "c", "d"), source = "unit test")
  ))
  expect_equal(mean(drawn$SEX == "M"), 0.7, tolerance = 0.1)
  expect_equal(mean(drawn$ARM == "a"), 0.25, tolerance = 0.2)
})

test_that("an integer covariate is drawn whole", {
  drawn <- .cov_draw(pmx_covariates(
    AGE = pmx_covariate(range = c(18, 90), median = 45, cv = 0.3,
                        distribution = "normal", integer = TRUE,
                        source = "unit test")
  ), n = 200)
  expect_true(all(drawn$AGE == round(drawn$AGE)))
  expect_gte(min(drawn$AGE), 18)
  expect_lte(max(drawn$AGE), 90)
})

test_that("a covariate refuses a declaration it cannot draw from", {
  expect_error(pmx_covariate(range = c(40, 130), cv = 0.2, source = "x"),
               "needs `median`")
  expect_error(
    pmx_covariate(range = c(40, 130), median = 200, cv = 0.2, source = "x"),
    "outside `range`"
  )
  expect_error(pmx_covariate(range = c(40, 130), median = 75, cv = 0,
                             source = "x"), "positive")
  expect_error(pmx_covariate(range = c(40, 130), distribution = "normal",
                             source = "x"), "needs `median` and `cv`")
  expect_error(pmx_covariate(range = c(40, 130), prob = c(1, 2), source = "x"),
               "categorical")
  expect_error(pmx_covariate(levels = c("M", "F"), median = 75, source = "x"),
               "continuous")
  expect_error(
    pmx_covariate(levels = c("M", "F"), prob = c(1, 2, 3), source = "x"),
    "one non-negative number per level"
  )
})

test_that("the reference catalogue builds usable covariates", {
  covariates <- pmx_covariates_reference(c("WT", "AGE", "SEX", "RACE"))
  expect_s3_class(covariates, "pmx_covariates")
  expect_named(covariates, c("WT", "AGE", "SEX", "RACE"))
  # Provenance travels with each declaration and names its own basis, so a
  # value the protocol is supposed to decide says so rather than reading as a
  # survey measurement.
  expect_match(covariates$WT$source, "NHANES")
  expect_match(covariates$AGE$source, "placeholder")
  expect_match(covariates$SEX$source, "placeholder")
  drawn <- .cov_draw(covariates)
  expect_equal(stats::median(drawn$WT), 84, tolerance = 0.05)
  expect_true(all(drawn$AGE == round(drawn$AGE)))
  expect_equal(mean(drawn$RACE == "White"), 0.65, tolerance = 0.1)
})

test_that("a column name maps onto a reference quantity", {
  # A named vector renames, and the spellings that occur in the public
  # datasets resolve without one.
  renamed <- pmx_covariates_reference(c(BWT = "WT"))
  expect_named(renamed, "BWT")
  expect_equal(renamed$BWT$median, 84)
  expect_named(pmx_covariates_reference(c("WEIGHTB", "HGT")),
               c("WEIGHTB", "HGT"))
  expect_equal(pmx_covariates_reference("WEIGHTB")$WEIGHTB$median, 84)

  expect_error(pmx_covariates_reference("BODYWEIGHT"), "No reference")
  expect_error(pmx_covariates_reference(c("WT", "WT")), "each output column")
  expect_error(pmx_covariates_reference("WT", medians = c(AGE = 50)),
               "not being built")
  expect_error(pmx_covariates_reference("WT", cvs = c(AGE = 0.2)),
               "not being built")
})

test_that("a replaced median or CV is honoured and stops claiming the source", {
  cohort <- pmx_covariates_reference(c(WEIGHTB = "WT"),
                                     medians = c(WEIGHTB = 75),
                                     cvs = c(WEIGHTB = 0.16))
  expect_equal(cohort$WEIGHTB$median, 75)
  expect_equal(cohort$WEIGHTB$cv, 0.16)
  expect_match(cohort$WEIGHTB$source, "median and CV replaced by the caller")

  # A median far from the default widens the bound with it, or every draw
  # would clip to one edge.
  heavy <- pmx_covariates_reference(c(WEIGHTB = "WT"),
                                    medians = c(WEIGHTB = 117))
  expect_gt(heavy$WEIGHTB$range[2L], 117)
  expect_equal(heavy$WEIGHTB$cv, 0.25)
  expect_match(heavy$WEIGHTB$source, "median replaced by the caller")
  expect_false(grepl("CV replaced", heavy$WEIGHTB$source))

  drawn <- .cov_draw(cohort)
  expect_equal(stats::median(drawn$WEIGHTB), 75, tolerance = 0.05)
  expect_equal(stats::sd(drawn$WEIGHTB) / mean(drawn$WEIGHTB), 0.16,
               tolerance = 0.15)
})

test_that("the reference table lists every catalogue entry with its basis", {
  table <- pmx_covariate_reference_table()
  expect_setequal(table$covariate, names(.covariate_reference))
  expect_true(all(nzchar(table$quantity)))
  expect_true(all(nzchar(table$basis)))
  # The units are the ones this package writes, and the table is where a
  # reader checks them, so none may be blank for a continuous covariate.
  continuous <- table$distribution != "categorical"
  expect_false(any(is.na(table$unit[continuous])))
  # Anything the protocol is meant to decide is marked, not dressed up as a
  # measurement.
  expect_match(table$basis[table$covariate %in% c("AGE", "SEX", "RACE")],
               "placeholder", all = TRUE)
})

# The gap the documentation claims, held as a test. NHANES describes the
# general adult population; eligibility criteria cut its tails, so a trial
# cohort runs lighter and tighter. Medians measured from the studies in the
# public-data surveys, so neither the reference value nor the claim about it
# can drift without this failing.
test_that("the survey reference runs heavier and wider than a trial cohort", {
  reference <- .covariate_reference$WT
  cohorts <- c(theo_md = 70.5, mavoglurant = 82.1, mad = 78.9,
               mixroute_sim = 72.0)
  expect_true(all(reference$median > cohorts))
  # Heavier, but a starting point rather than a different population: within a
  # fifth of every ordinary adult cohort measured.
  expect_true(all((reference$median - cohorts) / cohorts < 0.20))
  # Wider too. Those cohorts run at CVs of 14 to 18 per cent.
  expect_gt(reference$cv, 0.18)
  # A cohort selected for obesity runs the other way, which is the other
  # reason `medians` exists.
  expect_lt(reference$median, 117.0)
})

test_that("a calibrated release centres the covariate on what it measured", {
  # What a release carries is an arithmetic mean; what a lognormal is centred
  # on is its median. Generating without converting between them put the
  # column's mean above the quantity that was actually released.
  model <- .cov_model()
  design <- .cov_design(400)
  roles <- .cov_roles("WT")
  covariates <- pmx_covariates_reference("WT")
  source <- synpmx_prior(model, design, roles, n_subjects = 400, seed = 1,
                         covariates = covariates)
  observed <- mean(source$WT[!duplicated(source$ID)])

  # A large epsilon so the released centre is the measured one and this tests
  # the conversion rather than the noise.
  synthetic <- synpmx_calibrated(
    data = source, roles = roles, model = model, design = design,
    priors = pmx_priors(pk = pmx_prior(c(1 / 4, 4), source = "unit test")),
    epsilon = 50, seed = 2, covariates = covariates,
    backend = "public", public_source = TRUE
  )
  drawn <- synthetic$WT[!duplicated(synthetic$ID)]
  expect_equal(mean(drawn), observed, tolerance = 0.03)

  # One budget slice per covariate, whatever its distribution says: the
  # subject count, the PK correction, and weight.
  expect_equal(attr(synthetic, "synpmx_release")$preflight$d, 3)
})
