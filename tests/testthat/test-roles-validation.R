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

test_that("the four generation modes are exported side by side", {
  # One function per mode, each returning a synthetic dataset. The fit and
  # generate primitives behind the two confidential modes stay internal; the
  # only supported way to spend budget is through these four.
  exports <- getNamespaceExports("synpmx")
  expect_true(all(c(
    "synpmx_avatar",       # real templates, blended trajectories; no DP claim
    "synpmx_prior",        # public model and protocol only; reads no data
    "synpmx_calibrated",   # public model, magnitude privately corrected
    "synpmx_empirical"     # dense noised population summaries
  ) %in% exports))
  expect_true(all(c(
    "synpmx_generate", "privacy_report", "validate_private_model"
  ) %in% exports))
  expect_false(any(c(
    ".fit_private", ".fit_calibrated", ".generate_private",
    ".generate_structural"
  ) %in% exports))
})

# Real assignment and schedule columns have gaps. Validation used to refuse both,
# which stopped a run over data that the generator handles perfectly well: a
# missing nominal time falls back to the inferred grid row by row, and a missing
# stratum is simply its own level. Found 2026-07-29 on a real study.

gap_roles <- function(...) {
  pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
            cmt = "CMT", dvid = "DVID", covariates = "WT", ...)
}

test_that("a nominal time with gaps is accepted and falls back per row", {
  data <- pmx_simulated_fixture(20)
  data$NTIME[c(3L, 40L, 77L)] <- NA_real_
  roles <- gap_roles(nominal_time = "NTIME")

  report <- validate_pmx(data, roles)
  expect_true(report$valid)
  expect_match(
    report$checks$message[report$checks$check == "nominal_time"],
    "3 row\\(s\\) have none"
  )
  # Rows that have a nominal time snap to it; the rest use the inferred grid.
  synthetic <- suppressWarnings(synpmx_avatar(data, roles, seed = 1))
  expect_identical(attr(synthetic, "pmx_settings")$time_grid, "mixed")
  expect_true(validate_pmx(synthetic, roles)$valid)
})

test_that("a nominal time that is entirely missing is refused", {
  data <- pmx_simulated_fixture(10)
  data$NTIME <- NA_real_
  report <- validate_pmx(data, gap_roles(nominal_time = "NTIME"))
  expect_false(report$valid)
  expect_match(
    report$checks$message[report$checks$check == "nominal_time"],
    "entirely missing"
  )
})

test_that("a stratum with gaps warns and becomes its own level", {
  data <- pmx_simulated_fixture(20)
  data$TRT <- rep(c("A", "B"), each = nrow(data) / 2)
  data$TRT[data$ID %in% c("5", "6")] <- NA_character_
  roles <- gap_roles(subject_properties = "TRT")

  report <- validate_pmx(data, roles)
  # A warning, not an error: visible without stopping a run over data that is
  # merely incomplete.
  expect_true(report$valid)
  expect_identical(
    report$checks$status[report$checks$check == "subject_property_TRT"],
    "warning"
  )
  synthetic <- suppressWarnings(synpmx_avatar(data, roles, seed = 1))
  expect_true(validate_pmx(synthetic, roles)$valid)
  expect_true(anyNA(synthetic$TRT))
})

test_that("a stratum that varies within subject is still an error", {
  # Missing is incomplete data; varying means the column is not a subject-level
  # assignment at all, and no stratum can be built from it.
  data <- pmx_simulated_fixture(10)
  data$TRT <- rep(c("A", "B"), length.out = nrow(data))
  report <- validate_pmx(data, gap_roles(subject_properties = "TRT"))
  expect_false(report$valid)
  expect_match(
    report$checks$message[report$checks$check == "subject_property_TRT"],
    "varies within"
  )
})

test_that("a declared column that does not exist is named, not ignored", {
  # Deliberately still fatal. Only role-named columns survive generation, so
  # silently skipping a role that points at nothing would drop data on a typo.
  data <- pmx_simulated_fixture(10)
  expect_error(
    synpmx_avatar(data, gap_roles(keep = c("STUDYID", "TRT")), seed = 1),
    "Role columns not found in `data`: STUDYID, TRT"
  )
})

test_that("cmt and dvid may name the same column, other collisions may not", {
  # NONMEM's CMT is the dosing compartment on event rows and the endpoint key on
  # observation rows. Making the user copy the column to itself under another
  # name declares nothing extra.
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", cmt = "CMT", dvid = "CMT")
  expect_equal(roles$cmt, "CMT")
  expect_equal(roles$dvid, "CMT")

  # The overlap is exactly cmt/dvid; nothing else is loosened.
  expect_error(
    pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
              cmt = "CMT", covariates = "CMT"),
    "A column cannot have multiple roles"
  )
  expect_error(
    pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
              dvid = "CMT", keep = "CMT"),
    "A column cannot have multiple roles"
  )
})
