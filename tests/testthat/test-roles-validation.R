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
  role_args$strata <- "ARM"
  roles <- do.call(pmx_roles, role_args)
  expect_true(validate_pmx(source, roles, private_endpoints())$valid)
  source$ARM[2L] <- "B"
  expect_false(validate_pmx(source, roles, private_endpoints())$valid)
})

test_that("the generation modes are exported side by side", {
  # One function per mode, each returning a synthetic dataset. The fit and
  # generate primitives behind the two confidential modes stay internal; the
  # only supported way to spend budget is through the four listed here.
  # `synpmx_pca()` is the exception by design: it makes no DP claim, so its two
  # stages are public and can be run separately.
  exports <- getNamespaceExports("synpmx")
  expect_true(all(c(
    "synpmx_avatar",       # real templates, blended trajectories; no DP claim
    "synpmx_pca",          # fitted basis, drawn scores; no DP claim
    "synpmx_prior",        # public model and protocol only; reads no data
    "synpmx_calibrated",   # public model, magnitude privately corrected
    "synpmx_empirical"     # dense noised population summaries
  ) %in% exports))
  expect_true(all(c(
    "synpmx_pca_summarize", "synpmx_pca_generate"
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

test_that("a stratum with gaps is an error", {
  # It used to be a warning, on the reasoning that a missing value could stand
  # as its own stratum level. It cannot: the unlabelled group is then sized,
  # balanced and reported as though it were an arm.
  data <- pmx_simulated_fixture(20)
  data$TRT <- rep(c("A", "B"), each = nrow(data) / 2)
  data$TRT[data$ID %in% c("5", "6")] <- NA_character_
  roles <- gap_roles(strata = "TRT")

  report <- validate_pmx(data, roles)
  expect_false(report$valid)
  expect_identical(
    report$checks$status[report$checks$check == "stratum_TRT"],
    "error"
  )
  expect_match(report$checks$message[report$checks$check == "stratum_TRT"],
               "missing for 2 subject")
  expect_error(synpmx_avatar(data, roles, seed = 1), "must be recorded")
})

test_that("a blank stratum is caught the same way, and named as blank", {
  # `""` is how the gap usually arrives, and it used to pass every check here:
  # neither missing nor varying.
  data <- pmx_simulated_fixture(20)
  data$TRT <- rep(c("A", "B"), each = nrow(data) / 2)
  data$TRT[data$ID %in% c("5", "6")] <- ""
  roles <- gap_roles(strata = "TRT")

  report <- validate_pmx(data, roles)
  expect_false(report$valid)
  expect_match(report$checks$message[report$checks$check == "stratum_TRT"],
               "blank in 2 of them")
})

test_that("a stratum that varies within subject is still an error", {
  # Missing is incomplete data; varying means the column is not a subject-level
  # assignment at all, and no stratum can be built from it.
  data <- pmx_simulated_fixture(10)
  data$TRT <- rep(c("A", "B"), length.out = nrow(data))
  report <- validate_pmx(data, gap_roles(strata = "TRT"))
  expect_false(report$valid)
  expect_match(
    report$checks$message[report$checks$check == "stratum_TRT"],
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

# The `tad` role is an output, not an input: `synpmx_avatar()` overwrites it.
# `validate_pmx()` is the only place the declared values are read, and it
# reports rather than refuses, because the source may be right and the
# derivation wrong for that study.
tad_source <- function() {
  do.call(rbind, lapply(1:6, function(i) {
    out <- rbind(
      data.frame(ID = i, TIME = c(0, 24), DV = NA_real_, AMT = 100, EVID = 1L,
                 CMT = 1L, TAD = 0, WT = 60 + i),
      data.frame(ID = i, TIME = c(1, 6, 25, 30), DV = c(9, 5, 8, 4),
                 AMT = 0, EVID = 0L, CMT = 2L, TAD = c(1, 6, 1, 6),
                 WT = 60 + i))
    out[order(out$TIME, -out$EVID), ]
  }))
}

tad_roles <- function(...) {
  pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
            cmt = "CMT", tad = "TAD", covariates = "WT", ...)
}

test_that("a TAD that agrees with the times passes and says so", {
  report <- validate_pmx(tad_source(), tad_roles())
  expect_true(report$valid)
  agreement <- report$checks[report$checks$check == "tad_agreement", ]
  expect_equal(nrow(agreement), 1L)
  expect_equal(agreement$status, "pass")
})

test_that("a TAD on another convention warns without invalidating", {
  # A study measuring TAD from a nominal dose time an hour earlier than the
  # recorded one. Nothing here is malformed -- it is a different convention --
  # so it must warn and stay valid.
  source <- tad_source()
  source$TAD[source$EVID == 0L] <- source$TAD[source$EVID == 0L] + 1
  report <- validate_pmx(source, tad_roles())
  expect_true(report$valid)
  agreement <- report$checks[report$checks$check == "tad_agreement", ]
  expect_equal(agreement$status, "warning")
  expect_match(agreement$message, "24 of 24 observation row")
  # And the warning must be visible, which it was not: print() returned on the
  # valid branch before reaching the warnings block.
  expect_output(print(report), "Warnings \\(not fatal\\)")
})

test_that("the agreement check is skipped, loudly, when addl or ii is declared", {
  source <- tad_source()
  source$ADDL <- 0L
  source$II <- 0
  report <- validate_pmx(source, tad_roles(addl = "ADDL", ii = "II"))
  agreement <- report$checks[report$checks$check == "tad_agreement", ]
  expect_equal(agreement$status, "warning")
  expect_match(agreement$message, "does not expand them")
})

test_that("generation overwrites TAD from the generated times", {
  source <- tad_source()
  source$TAD[source$EVID == 0L] <- 999           # nonsense the generator ignores
  synthetic <- suppressWarnings(suppressMessages(
    synpmx_avatar(source, tad_roles(), n_subjects = 6, seed = 2)
  ))
  expect_false(any(synthetic$TAD == 999))
  expect_true(all(synthetic$TAD >= 0))
  # It equals the derivation, by construction.
  expect_equal(synthetic$TAD, synpmx:::.derived_tad(synthetic, tad_roles()),
               tolerance = 1e-8)
})
