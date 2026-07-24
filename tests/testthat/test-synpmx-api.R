# The four generation modes as a public surface.
#
# The engines themselves are covered by test-avatar.R, test-structural-v3.R,
# test-generation.R, and test-private-fitting.R. What is tested here is the
# contract the four public functions add on top: every mode returns a dataset,
# the confidential modes carry their release with them, and drawing further
# datasets from that release spends nothing.

api_model <- function() {
  pmx_structural_model("1cmt_oral", c(cl = 10, v = 70, ka = 1),
                       source = "unit test")
}

api_design <- function() {
  pmx_trial_design(c(10, 30, 100), c(6, 6, 6), c(0, .5, 1, 2, 4, 8, 12, 24),
                   source = "unit test protocol")
}

api_priors <- function() {
  pmx_priors(pk = pmx_prior(c(1 / 4, 4), "unit test"))
}

# A public stand-in for confidential data, so the calibrated mode can be
# exercised without a DP backend. `public_source = TRUE` makes the release
# noiseless and explicitly claimless; that is the point of the fixture.
api_source <- function() {
  .generate_structural(api_model(), api_design(), n_subjects = 30, seed = 7)
}

api_calibrated <- function(data = api_source(), ...) {
  suppressWarnings(synpmx_calibrated(
    data = data, roles = pmx_generated_roles(), model = api_model(),
    design = api_design(), priors = api_priors(), epsilon = 1,
    backend = "public", public_source = TRUE, ...
  ))
}

test_that("every mode returns a data frame, not a generator", {
  source <- pmx_simulated_fixture(20)
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", dvid = "DVID", cmt = "CMT", mdv = "MDV",
                     covariates = "WT")
  avatar <- suppressWarnings(synpmx_avatar(source, roles, seed = 101))
  expect_s3_class(avatar, "data.frame")
  expect_gt(nrow(avatar), 0)

  prior <- synpmx_prior(api_model(), api_design(), n_subjects = 6, seed = 202)
  expect_s3_class(prior, "data.frame")
  expect_gt(nrow(prior), 0)

  expect_s3_class(api_calibrated(), "data.frame")
})

test_that("synpmx_prior refuses a calibrated release", {
  # The mode names describe what the output is built from, so handing a release
  # to the public-inputs-only mode is a category error, not a shortcut.
  release <- attr(api_calibrated(), "synpmx_release", exact = TRUE)
  expect_error(synpmx_prior(release, api_design()), "pmx_structural_model")
})

test_that("a confidential dataset carries its own release", {
  syn <- api_calibrated()
  release <- attr(syn, "synpmx_release", exact = TRUE)
  expect_s3_class(release, "pmx_calibrated_model")
  expect_identical(.release_of(syn), release)
})

test_that("datasets from public-input modes carry no release", {
  prior <- synpmx_prior(api_model(), api_design(), n_subjects = 6, seed = 202)
  expect_null(attr(prior, "synpmx_release", exact = TRUE))
  expect_error(.release_of(prior), "carries no privacy release")
})

test_that("regenerating from a dataset spends no further budget", {
  syn <- api_calibrated()
  release <- attr(syn, "synpmx_release", exact = TRUE)

  again <- synpmx_generate(syn, seed = 999)
  expect_s3_class(again, "data.frame")
  # The same release object, so the same epsilon: nothing new was spent.
  expect_identical(attr(again, "synpmx_release", exact = TRUE), release)
  # A different seed is post-processing, so the data itself does differ.
  expect_false(identical(syn, again))
})

test_that("n_datasets draws several from one release", {
  many <- api_calibrated(n_datasets = 3L)
  expect_type(many, "list")
  expect_length(many, 3L)
  releases <- lapply(many, attr, "synpmx_release", exact = TRUE)
  # One fit, three datasets: every dataset points at the same release.
  expect_true(all(vapply(releases, identical, logical(1), releases[[1L]])))
  expect_false(identical(many[[1L]], many[[2L]]))
})

test_that("synpmx_generate accepts a release object as readily as a dataset", {
  syn <- api_calibrated()
  release <- attr(syn, "synpmx_release", exact = TRUE)
  expect_equal(
    synpmx_generate(release, seed = 555),
    synpmx_generate(syn, seed = 555),
    ignore_attr = TRUE
  )
})

test_that("a repeated budget-spending fit warns about the total spent", {
  # Nothing else in the package would notice a second fit, and the cost is
  # silent and irreversible, so the warning names the running total.
  data <- api_source()
  api_calibrated(data)
  expect_warning(
    synpmx_calibrated(
      data = data, roles = pmx_generated_roles(), model = api_model(),
      design = api_design(), priors = api_priors(), epsilon = 1,
      backend = "public", public_source = TRUE
    ),
    "budget-spending fit"
  )
})

test_that("the empirical mode returns a dataset carrying its release", {
  syn <- suppressWarnings(synpmx_empirical(
    data = private_fixture(), roles = private_roles(),
    endpoints = private_endpoints(), epsilon = 5, delta = 0,
    bounds = private_bounds(), public_design = private_design(private_fixture()),
    contribution_limits = private_limits(),
    budget_allocation = private_budget(), backend = "public",
    public_source = TRUE
  ))
  expect_s3_class(syn, "data.frame")
  release <- attr(syn, "synpmx_release", exact = TRUE)
  expect_s3_class(release, "private_pmx_model")
  # The accounting helpers read straight off the dataset.
  expect_s3_class(privacy_report(syn), "pmx_privacy_report")
})

# REV-023: the DP engines' unaudited status must be an explicit, session-level
# acknowledgment, not just documentation, so it cannot be reached by accident.

test_that("the DP engines refuse to run without synpmx_enable_dp_engines()", {
  synpmx_disable_dp_engines()
  on.exit(suppressMessages(synpmx_enable_dp_engines()), add = TRUE)

  data <- api_source()
  expect_error(
    synpmx_calibrated(
      data = data, roles = pmx_generated_roles(), model = api_model(),
      design = api_design(), priors = api_priors(), epsilon = 1
    ),
    "synpmx_enable_dp_engines"
  )
  expect_error(
    synpmx_empirical(
      data = private_fixture(), roles = private_roles(),
      endpoints = private_endpoints(), epsilon = 5, delta = 0,
      bounds = private_bounds(), public_design = private_design(private_fixture()),
      contribution_limits = private_limits(), budget_allocation = private_budget()
    ),
    "synpmx_enable_dp_engines"
  )
})

test_that("backend = \"public\" is exempt from the DP-engines gate", {
  synpmx_disable_dp_engines()
  on.exit(suppressMessages(synpmx_enable_dp_engines()), add = TRUE)
  expect_s3_class(api_calibrated(), "data.frame")
})

test_that("synpmx_enable_dp_engines() unlocks the gate for the session", {
  skip_if_not(dp_backend_status()$available, "OpenDP unavailable")
  synpmx_disable_dp_engines()
  expect_message(synpmx_enable_dp_engines(), "DP engines enabled")
  on.exit(suppressMessages(synpmx_enable_dp_engines()), add = TRUE)
  expect_s3_class(api_calibrated(backend = "opendp",
                                 public_source = FALSE), "data.frame")
})
