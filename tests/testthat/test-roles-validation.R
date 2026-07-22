test_that("roles and endpoint metadata are explicit", {
  roles <- private_roles()
  expect_s3_class(roles, "pmx_roles")
  expect_equal(roles$nominal_time, "NTIME")
  expect_equal(roles$cens, "CENS")
  expect_error(
    pmx_roles("ID", "TIME", "DV", evid = "EVID", covariates = "ID"),
    "multiple roles"
  )
  endpoint <- pmx_endpoint(
    dvid = "cp", alignment = "dose_relative", transform = "auto",
    shape = "occasion"
  )
  expect_s3_class(endpoint, "pmx_endpoint")
  expect_error(pmx_endpoint(alignment = "unknown"), "arg")
  expect_error(pmx_endpoint(shape = "global"), "alignment")
  expect_error(pmx_endpoint(alignment = "study_time"), "shape")
})

test_that("validation covers timing, endpoints, and baseline constancy", {
  source <- private_fixture()
  report <- validate_pmx(source, private_roles(), private_endpoints())
  expect_s3_class(report, "pmx_validation")
  expect_true(report$valid)
  expect_equal(report$summary$subjects, 8L)

  broken <- source
  broken$WT[2L] <- broken$WT[2L] + 1
  expect_false(validate_pmx(broken, private_roles(), private_endpoints())$valid)
  expect_error(
    validate_pmx(broken, private_roles(), private_endpoints(), strict = TRUE),
    "varies within"
  )
})

test_that("validation supports reset occasion clocks and coherent properties", {
  source <- private_fixture()
  reset <- source
  reset$TIME[reset$OCC == 2L] <- reset$TIME[reset$OCC == 2L] - 12
  expect_true(validate_pmx(reset, private_roles(), private_endpoints())$valid)

  source$ARM <- ifelse(source$ID %% 2L, "A", "B")
  role_args <- unclass(private_roles())
  role_args$subject_properties <- "ARM"
  roles <- do.call(pmx_roles, role_args)
  expect_true(validate_pmx(source, roles, private_endpoints())$valid)
  source$ARM[2L] <- "B"
  expect_false(validate_pmx(source, roles, private_endpoints())$valid)
})

test_that("Version 1 synthesis modes are not public", {
  exports <- getNamespaceExports("pmxSynthData")
  expect_false("mock_pmx" %in% exports)
  expect_false(any(grepl("avatar|blend|template", exports, ignore.case = TRUE)))
  expect_true(all(c(
    "fit_private_pmx", "generate_pmx", "privacy_report",
    "validate_private_model", "subject_property_summary"
  ) %in% exports))
})
