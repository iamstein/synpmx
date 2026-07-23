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

test_that("all engines are exported side by side", {
  # Version 4 restored the AVATAR engine as the primary method. The DP (v2) and
  # structural (v3) engines are kept as superseded alternatives, so all three
  # public entry points coexist.
  exports <- getNamespaceExports("synpmx")
  expect_true("synthesize_pmx" %in% exports)             # AVATAR (v4, primary)
  expect_true("fit_private_pmx" %in% exports)      # DP dense grid (v2)
  expect_true("fit_calibrated_pmx" %in% exports)   # structural correction (v3)
  expect_true(all(c(
    "generate_pmx", "privacy_report", "validate_private_model"
  ) %in% exports))
})
