test_that("the public censoring fixture covers all supported conventions", {
  fixture <- pmx_censoring_fixture()
  roles <- pmx_roles(
    id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
    dvid = "DVID", mdv = "MDV", cens = "CENS", limit = "LIMIT"
  )
  expect_true(validate_pmx(fixture, roles)$valid)
  expect_setequal(unique(fixture$CENS), c(-1L, 0L, 1L))
  expect_true(any(fixture$CENS == 1L & is.finite(fixture$LIMIT)))

  broken <- fixture
  broken$CENS[2] <- 2L
  expect_false(validate_pmx(broken, roles)$valid)
  broken <- fixture
  broken$LIMIT[5] <- 5
  expect_false(validate_pmx(broken, roles)$valid)
})

test_that("latent trajectories produce coherent generated censoring", {
  source <- private_fixture()
  endpoints <- private_endpoints()
  endpoints$cp$censoring <- list(left = 3, right = 15)
  endpoints$pd$censoring <- list(interval = c(45, 55))
  model <- suppressWarnings(fit_private_pmx(
    source, private_roles(), endpoints, 5, 0, private_bounds(),
    private_design(source), private_limits(), private_budget(),
    backend = "public", public_source = TRUE
  ))
  synthetic <- generate_pmx(model, 10, 202)
  observed <- synthetic$EVID == 0
  expect_true(any(synthetic$CENS[observed] != 0))
  expect_true(all(synthetic$CENS %in% c(-1L, 0L, 1L)))
  interval <- observed & synthetic$CENS == 1L & is.finite(synthetic$LIMIT)
  expect_true(all(synthetic$LIMIT[interval] <= synthetic$DV[interval]))
  expect_true(validate_pmx(synthetic, private_roles(), endpoints)$valid)
})
