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
  mock <- generate_pmx(model, 10, 202)
  observed <- mock$EVID == 0
  expect_true(any(mock$CENS[observed] != 0))
  expect_true(all(mock$CENS %in% c(-1L, 0L, 1L)))
  interval <- observed & mock$CENS == 1L & is.finite(mock$LIMIT)
  expect_true(all(mock$LIMIT[interval] <= mock$DV[interval]))
  expect_true(validate_pmx(mock, private_roles(), endpoints)$valid)
})
