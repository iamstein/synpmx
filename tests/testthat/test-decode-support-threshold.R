# Regression tests for SIM-020 / REV-001.
#
# Released presence fields are unnormalized subject counts, not probabilities.
# Decoding them against a bare constant lets Laplace noise around zero pass as
# "this grid cell had support", after which the separately noised value sum
# clamps to zero and the cell decodes to the bottom of the working domain.

test_that("the support threshold scales with the release perturbation", {
  # Noiseless releases keep the historical constant exactly, so public-fixture
  # behavior is unchanged.
  expect_equal(.support_threshold(count = 60, noise_scale = 0), 0.25)
  expect_equal(.support_threshold(count = 2000, noise_scale = 0), 0.25)

  # With noise, the gate must grow with the noise scale.
  expect_equal(.support_threshold(count = 2000, noise_scale = 20), 60)
  expect_gt(.support_threshold(count = 2000, noise_scale = 40),
            .support_threshold(count = 2000, noise_scale = 20))

  # ... but must never exceed half the cohort, or it would delete genuinely
  # occupied cells whenever noise is large relative to N.
  expect_equal(.support_threshold(count = 8, noise_scale = 20), 4)
  expect_lte(.support_threshold(count = 40, noise_scale = 1e6), 20)

  # Degenerate inputs fall back to the constant rather than propagating.
  expect_equal(.support_threshold(count = NA_real_, noise_scale = NA_real_),
               0.25)
  expect_equal(.support_threshold(count = -5, noise_scale = -5), 0.25)
})

test_that("noise in an unsupported cell does not decode to the domain floor", {
  fixture <- .threshold_fixture()
  endpoints <- fixture$endpoints
  cells <- length(endpoints$cp$grid)
  floor_value <- .endpoint_working_bounds(endpoints$cp)[1L]

  # Cells 1-3 hold real support from a 600-subject cohort. Cells 4-5 hold no
  # support at all; cell 4 carries a plausible Laplace excursion for a release
  # with sensitivity 40 at epsilon 0.4 (scale 100 -> an excursion of 40 is
  # entirely ordinary).
  released <- c(
    stats::setNames(c(600, 600, 600, 40, 0), fixture$map$cp$presence),
    stats::setNames(c(300, 420, 360, 0, 0), fixture$map$cp$value)
  )

  decoded <- .decode_trajectories(
    released, fixture$map, endpoints, count = 600, noise_scale = 100
  )$cp
  expect_length(decoded$mean_working, cells)
  expect_true(all(is.finite(decoded$mean_working)))
  expect_false(any(decoded$mean_working <= floor_value + 1e-8))

  # The supported cells must still decode to their released means, so the gate
  # suppresses noise without discarding signal.
  supported <- .from_unit(c(300, 420, 360) / 600,
                          .endpoint_working_bounds(endpoints$cp))
  expect_equal(decoded$mean_working[1:3], supported, tolerance = 1e-8)
})

test_that("the bare constant is what produced the domain-floor artifact", {
  fixture <- .threshold_fixture()
  floor_value <- .endpoint_working_bounds(fixture$endpoints$cp)[1L]
  released <- c(
    stats::setNames(c(600, 600, 600, 40, 0), fixture$map$cp$presence),
    stats::setNames(c(300, 420, 360, 0, 0), fixture$map$cp$value)
  )
  # noise_scale = 0 reproduces the pre-fix gate. Cell 4 passes it, divides
  # zero by forty, and pins to the bottom of the working domain.
  decoded <- .decode_trajectories(
    released, fixture$map, fixture$endpoints, count = 600, noise_scale = 0
  )$cp
  expect_true(any(decoded$mean_working <= floor_value + 1e-8))
})
