test_that("roles are explicit and validated", {
  roles <- fixture_roles()
  expect_s3_class(roles, "pmx_roles")
  expect_equal(roles$id, "ID")
  expect_error(
    pmx_roles("ID", "TIME", "DV", evid = "EVID", covariates = "ID"),
    "multiple roles"
  )
  expect_error(
    validate_pmx(pmx_fixture(), pmx_roles("unknown", "TIME", "DV",
                                         evid = "EVID"), strict = TRUE),
    "not found"
  )
})

test_that("validation returns useful reports and strict mode fails", {
  source <- pmx_fixture()
  roles <- fixture_roles()
  report <- validate_pmx(source, roles)
  expect_s3_class(report, "pmx_validation")
  expect_true(report$valid)
  expect_equal(report$summary$subjects, 7L)
  expect_true(all(c("check", "status", "message") %in% names(report$checks)))

  broken <- source
  broken$WT[2L] <- broken$WT[2L] + 1
  expect_false(validate_pmx(broken, roles)$valid)
  expect_error(validate_pmx(broken, roles, strict = TRUE), "varies within")
})

test_that("optional roles can be absent", {
  source <- pmx_fixture()[, c("ID", "TIME", "DV", "AMT", "EVID", "CMT", "WT")]
  roles <- pmx_roles(
    "ID", "TIME", "DV", "AMT", "EVID", "CMT", covariates = "WT"
  )
  mock <- mock_pmx(source, roles, n_subjects = 2, seed = 1)
  expect_true(validate_pmx(mock, roles)$valid)
  expect_null(attr(mock, "pmx_settings")$roles$dvid)
  expect_null(attr(mock, "pmx_settings")$roles$mdv)
})
