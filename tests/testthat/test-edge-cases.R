test_that("missing and incomplete public configuration fails clearly", {
  source <- private_fixture()
  expect_error(pmx_bounds(c(0, 1), list()), "named list")
  expect_error(pmx_contribution_limits(1, 1, 1, 1, 1), "at least two")
  expect_error(pmx_budget_allocation(.5, .5, .5, 0, 0, 0), "sum")
  expect_error(
    pmx_budget_allocation(.1, .15, .15, .1, .5, 1e-12), "sum"
  )
  tampered_budget <- private_budget()
  tampered_budget$endpoints <- 1.5
  expect_error(
    fit_private_pmx(
      source, private_roles(), private_endpoints(), 5, 0, private_bounds(),
      private_design(source), private_limits(), tampered_budget,
      backend = "public", public_source = TRUE
    ),
    "modified"
  )
  wrong <- private_design(source)
  wrong$schema$columns <- wrong$schema$columns[-1L]
  expect_error(
    fit_private_pmx(
      source, private_roles(), private_endpoints(), 5, 0, private_bounds(),
      wrong, private_limits(), private_budget(), backend = "public",
      public_source = TRUE
    ), "schema"
  )
  expect_error(
    pmx_public_design(
      pmx_schema(source),
      endpoint_occasion_grids = list(cp = list("1.5" = c(0, 1)))
    ),
    "positive occasion number"
  )
  expect_error(
    pmx_public_design(
      pmx_schema(source),
      endpoint_occasion_grids = list(cp = list("1" = c(1, 0)))
    ),
    "strictly increasing"
  )
})

test_that("small studies and excessive dimension warn without weakening privacy", {
  source <- private_fixture(2)
  expect_warning(
    model <- fit_private_pmx(
      source, private_roles(), private_endpoints(), 5, 0,
      private_bounds(), private_design(source), private_limits(),
      private_budget(), backend = "public", public_source = TRUE
    ),
    "below six|dimension"
  )
  expect_equal(model$privacy$epsilon, 5)
  expect_equal(model$privacy$delta, 0)
  generated <- generate_pmx(model, 2, 3)
  expect_true(validate_pmx(
    generated, private_roles(), private_endpoints()
  )$valid)
})

test_that("six-, twelve-, and larger simulated studies retain broad utility", {
  for (n in c(6L, 12L)) {
    source <- private_fixture(n)
    model <- fit_public_fixture(source)
    synthetic <- generate_pmx(model, n_subjects = n, seed = 700 + n)
    expect_equal(model$population$private_subject_count, n)
    expect_equal(length(unique(synthetic$ID)), n)
    expect_true(validate_pmx(synthetic, private_roles(), private_endpoints())$valid)
  }

  source <- pmx_simulated_fixture(60L)
  model <- fit_public_fixture(source)
  synthetic <- generate_pmx(model, n_subjects = 60L, seed = 760)
  source_observed <- source[source$EVID == 0, , drop = FALSE]
  synthetic_observed <- synthetic[synthetic$EVID == 0, , drop = FALSE]
  source_center <- tapply(source_observed$DV, source_observed$DVID, mean)
  synthetic_center <- tapply(synthetic_observed$DV, synthetic_observed$DVID, mean)
  expect_lt(abs(log1p(synthetic_center[["cp"]]) -
                  log1p(source_center[["cp"]])), 1)
  expect_lt(abs(synthetic_center[["pd"]] - source_center[["pd"]]), 45)
  expect_true(validate_pmx(synthetic, private_roles(), private_endpoints())$valid)
})

test_that("restricted comparisons label every source-derived component", {
  source <- private_fixture()
  synthetic <- generate_pmx(fit_public_fixture(source), 3, 90)
  comparison <- compare_pmx(source, synthetic, private_roles(), private_endpoints())
  expect_s3_class(comparison, "pmx_comparison")
  expect_true(all(comparison$release_status$release_status[
    comparison$release_status$component != "validation.synthetic"
  ] == "restricted_not_releasable"))
  expect_equal(attr(comparison$validation$synthetic, "release_status"),
               "releasable_post_processing")
})

test_that("an empirical audit flags a deliberately broken mechanism", {
  audit <- function(mechanism, left, right, epsilon, repetitions = 2000L) {
    x <- replicate(repetitions, mechanism(left))
    y <- replicate(repetitions, mechanism(right))
    event <- sort(unique(c(x, y)))[1L]
    px <- mean(x == event)
    py <- mean(y == event)
    px > exp(epsilon) * py + 0.02 || py > exp(epsilon) * px + 0.02
  }
  broken <- function(data) length(data)
  expect_true(audit(broken, 1:5, 1:6, epsilon = 1))
  # This is only a bug-finding test; passing an audit would not prove DP.
})
