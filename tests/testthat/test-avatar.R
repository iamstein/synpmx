# Version 4: the AVATAR-style synthesizer, synthesize_pmx().

test_that("synthesize_pmx preserves schema and produces fresh subjects", {
  source <- pmx_simulated_fixture(40)
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", dvid = "DVID", cmt = "CMT", mdv = "MDV",
                     covariates = c("WT", "AGE", "SEX"))
  synthetic <- suppressWarnings(synthesize_pmx(source, roles, n_subjects = 20,
                                               seed = 1))

  expect_true(validate_pmx(synthetic, roles)$valid)
  expect_equal(length(unique(synthetic$ID)), 20L)
  # New IDs: no generated subject reuses a source identifier.
  expect_length(intersect(synthetic$ID, source$ID), 0L)
  # Schema is restored: same columns and classes.
  expect_setequal(names(synthetic), names(source))
  expect_equal(vapply(synthetic[names(source)], class, character(1)),
               vapply(source, class, character(1)))
  expect_true(all(c("WT", "AGE", "SEX") %in% names(synthetic)))
})

test_that("generation is reproducible by seed and varies without it", {
  source <- pmx_simulated_fixture(30)
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", dvid = "DVID", cmt = "CMT", mdv = "MDV",
                     covariates = "WT")
  a <- suppressWarnings(synthesize_pmx(source, roles, n_subjects = 15, seed = 7))
  b <- suppressWarnings(synthesize_pmx(source, roles, n_subjects = 15, seed = 7))
  c <- suppressWarnings(synthesize_pmx(source, roles, n_subjects = 15, seed = 9))
  expect_equal(a, b)
  expect_false(isTRUE(all.equal(a$DV, c$DV)))
})

test_that("the default cohort size matches the source", {
  source <- pmx_simulated_fixture(24)
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", dvid = "DVID", cmt = "CMT", mdv = "MDV")
  synthetic <- suppressWarnings(synthesize_pmx(source, roles, seed = 1))
  expect_equal(length(unique(synthetic$ID)), 24L)
})

test_that("the caller's RNG state is left untouched", {
  source <- pmx_simulated_fixture(20)
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", dvid = "DVID", cmt = "CMT", mdv = "MDV")
  set.seed(123)
  before <- stats::runif(1)
  set.seed(123)
  invisible(suppressWarnings(synthesize_pmx(source, roles, seed = 5)))
  after <- stats::runif(1)
  expect_equal(before, after)
})

# The five nlmixr2data demonstrations, exercised through AVATAR. -------------

.avatar_datasets <- function() {
  list(
    theo_md = pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                        evid = "EVID", cmt = "CMT", covariates = "WT"),
    warfarin = pmx_roles(id = "id", time = "time", dv = "dv", amt = "amt",
                         evid = "evid", dvid = "dvid",
                         covariates = c("wt", "age", "sex")),
    wbcSim = pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                       evid = "EVID", cmt = "CMT"),
    nimoData = pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                         evid = "EVID", rate = "RATE", mdv = "MDV", tad = "TAD",
                         occasion = "OCC", covariates = c("BSA", "AGE", "HGT"),
                         subject_properties = "DOS", exclude = "WGT"),
    mavoglurant = pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                            evid = "EVID", cmt = "CMT", rate = "RATE",
                            mdv = "MDV", occasion = "OCC",
                            assigned_dose = "DOSE",
                            covariates = c("AGE", "SEX", "WT", "HT"))
  )
}

test_that("AVATAR synthesizes every nlmixr2data demonstration", {
  skip_if_not_installed("nlmixr2data")
  for (name in names(.avatar_datasets())) {
    roles <- .avatar_datasets()[[name]]
    env <- new.env()
    utils::data(list = name, package = "nlmixr2data", envir = env)
    source <- get(name, envir = env)

    synthetic <- suppressWarnings(synthesize_pmx(source, roles, seed = 1))
    expect_true(validate_pmx(synthetic, roles)$valid,
                info = paste(name, "should validate"))
    # Cohort size is preserved and identifiers are fresh.
    expect_equal(length(unique(synthetic[[roles$id]])),
                 length(unique(source[[roles$id]])), info = name)
    expect_length(intersect(synthetic[[roles$id]], source[[roles$id]]), 0L)
    # Declared endpoints all survive.
    if (!is.null(roles$dvid)) {
      expect_setequal(unique(synthetic[[roles$dvid]][!is.na(
        synthetic[[roles$dvid]])]),
        unique(source[[roles$dvid]][!is.na(source[[roles$dvid]])]))
    }
  }
})
