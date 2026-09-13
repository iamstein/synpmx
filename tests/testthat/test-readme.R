# README.md is hand-maintained: there is no README.Rmd knitting it, so nothing
# re-runs its example when behavior changes. This test is what keeps the example
# honest -- it declares exactly the roles shown under "Example with
# `synpmx_model()`" and checks that every argument the README passes still
# exists on the functions it passes them to.
#
# If this fails, the README is telling readers something the package no longer
# does. Fix the package or update the README -- do not relax the test.
#
# What it does not do is run the example. `synpmx_model()` fits a population
# model, which compiles, and `R CMD check` must never compile one -- the same
# reason `pmxmodel-demo.Rmd` knits against a stored fit. So the test pins the
# inputs and the interface, which is where a rename or a dropped role would
# show up, and the demo vignette and the public-data survey exercise the fit
# itself.

# The README's example is a real public dataset, so the test uses that dataset
# rather than a fixture standing in for it.
readme_study <- function() {
  study <- as.data.frame(get(utils::data(list = "case1_pkpd", package = "xgxr")))
  # CENS here flags the PK assay limit only, as the README's comment says.
  study$CENS[study$NAME == "PD - Continuous"] <- 0
  study
}

# Every role the README names, in the order it names them. A role dropped from
# `pmx_roles()` or renamed fails here rather than in a reader's session.
readme_roles <- function() {
  pmx_roles(
    id             = "ID",
    time           = "TIME",
    dv             = "LIDV",
    evid           = "EVID",
    amt            = "AMT",
    cmt            = "CMT",
    dvid           = "NAME",
    mdv            = NULL,
    rate           = NULL,
    nominal_time   = "NOMTIME",
    tad            = NULL,
    occasion       = NULL,
    cens           = "CENS",
    limit          = NULL,
    addl           = NULL,
    ii             = NULL,
    covariates     = "WEIGHTB",
    strata         = c("TRTACT", "DOSE"),
    dose_covariate = NULL,
    endpoint_types = NULL,
    keep           = "STUDY"
  )
}

test_that("the README declaration is accepted for the dataset it shows", {
  skip_if_not_installed("xgxr")
  study <- readme_study()
  roles <- readme_roles()

  #> [1] TRUE
  expect_true(validate_pmx(study, roles)$valid)

  # The README says every column that is not described is dropped, and names
  # `STUDY` as the one carried through verbatim.
  expect_true(all(c("WEIGHTB", "TRTACT", "DOSE", "STUDY") %in% names(study)))

  # It also says the cohort is large enough to fit, which is a gate rather
  # than a suggestion: `synpmx_model_estimate()` refuses under 20 subjects.
  expect_gte(length(unique(study$ID)), formals(synpmx_model_estimate)$min_subjects)
})

test_that("the README passes arguments these functions still have", {
  # The one-call form.
  expect_true(all(c("data", "roles", "n_subjects", "seed") %in%
                    names(formals(synpmx_model))))
  # The two-stage form the README shows underneath it.
  expect_true(all(c("data", "roles", "seed") %in%
                    names(formals(synpmx_model_estimate))))
  expect_true(all(c("n_subjects", "seed") %in%
                    names(formals(synpmx_model_generate))))
  # And the reader is told to call this on the fit.
  expect_true(is.function(model_report))
})
