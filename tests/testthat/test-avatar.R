# Version 4: the AVATAR-style synthesizer, synpmx_avatar().

test_that("synpmx_avatar preserves declared schema and produces fresh subjects", {
  source <- pmx_simulated_fixture(40)
  # Declaring every fixture column is what keeps the whole schema: undeclared
  # columns are dropped by design.
  roles <- pmx_roles(id = "ID", time = "TIME", nominal_time = "NTIME",
                     tad = "TAD", occasion = "OCC", dv = "DV", amt = "AMT",
                     rate = "RATE", evid = "EVID", cmt = "CMT", dvid = "DVID",
                     mdv = "MDV", cens = "CENS", limit = "LIMIT",
                     covariates = c("WT", "AGE", "SEX"))
  synthetic <- suppressWarnings(synpmx_avatar(source, roles, n_subjects = 20,
                                               seed = 1))

  expect_true(validate_pmx(synthetic, roles)$valid)
  expect_equal(length(unique(synthetic$ID)), 20L)
  # New IDs: no generated subject reuses a source identifier.
  expect_length(intersect(synthetic$ID, source$ID), 0L)
  # Every declared column is restored, with its class.
  expect_setequal(names(synthetic), names(source))
  expect_equal(vapply(synthetic[names(source)], class, character(1)),
               vapply(source, class, character(1)))
  expect_true(all(c("WT", "AGE", "SEX") %in% names(synthetic)))
})

test_that("synpmx_avatar drops undeclared columns and keeps `keep` ones", {
  source <- pmx_simulated_fixture(20)
  source$USUBJID <- sprintf("SECRET-%04d", source$ID)   # undeclared identifier
  source$ARM <- ifelse(source$ID <= 10, "A", "B")        # kept verbatim
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", dvid = "DVID", cmt = "CMT", mdv = "MDV",
                     covariates = "WT", keep = "ARM")
  synthetic <- suppressWarnings(
    synpmx_avatar(source, roles, n_subjects = 15, seed = 2)
  )
  # The undeclared identifier does not survive; the kept column does.
  expect_false("USUBJID" %in% names(synthetic))
  expect_true("ARM" %in% names(synthetic))
  # A kept value is copied verbatim, so it stays constant within a subject that
  # came from one anchor, and uses only real source levels.
  per_subject <- tapply(synthetic$ARM, synthetic$ID,
                        function(x) length(unique(x)))
  expect_true(all(per_subject == 1L))
  expect_true(all(synthetic$ARM %in% c("A", "B")))
})

test_that("generation is reproducible by seed and varies without it", {
  source <- pmx_simulated_fixture(30)
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", dvid = "DVID", cmt = "CMT", mdv = "MDV",
                     covariates = "WT")
  a <- suppressWarnings(synpmx_avatar(source, roles, n_subjects = 15, seed = 7))
  b <- suppressWarnings(synpmx_avatar(source, roles, n_subjects = 15, seed = 7))
  c <- suppressWarnings(synpmx_avatar(source, roles, n_subjects = 15, seed = 9))
  expect_equal(a, b)
  expect_false(isTRUE(all.equal(a$DV, c$DV)))
})

test_that("the default cohort size matches the source", {
  source <- pmx_simulated_fixture(24)
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", dvid = "DVID", cmt = "CMT", mdv = "MDV")
  synthetic <- suppressWarnings(synpmx_avatar(source, roles, seed = 1))
  expect_equal(length(unique(synthetic$ID)), 24L)
})

test_that("factor IDs stay factors and get fresh levels", {
  source <- pmx_simulated_fixture(12)
  source$ID <- factor(source$ID)
  roles <- pmx_roles(
    id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
    dvid = "DVID", cmt = "CMT", mdv = "MDV"
  )
  synthetic <- suppressWarnings(synpmx_avatar(source, roles, n_subjects = 6,
                                                seed = 12))

  expect_true(is.factor(synthetic$ID))
  expect_false(anyNA(synthetic$ID))
  expect_true(all(as.character(synthetic$ID) %in% levels(synthetic$ID)))
  expect_length(intersect(as.character(synthetic$ID),
                          as.character(source$ID)), 0L)
})

test_that("AVATAR rejects differential-privacy-only roles", {
  source <- pmx_simulated_fixture(12)
  base <- list(id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
               dvid = "DVID", cmt = "CMT", mdv = "MDV")
  expect_error(
    synpmx_avatar(source, do.call(pmx_roles,
                                  c(base, list(subject_properties = "SEX")))),
    "does not use it"
  )
  expect_error(
    synpmx_avatar(source, do.call(pmx_roles,
                                  c(base, list(exclude = "AGE")))),
    "leave the column undeclared"
  )
})

test_that("dose magnitudes constrain AVATAR donors", {
  source <- do.call(rbind, lapply(1:2, function(id) data.frame(
    ID = as.integer(id), TIME = c(0, 1, 2), DV = c(0, id, id / 2),
    AMT = c(100 * id, 0, 0), EVID = c(1L, 0L, 0L), CMT = c(1L, 2L, 2L)
  )))
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", cmt = "CMT")
  synthetic <- suppressWarnings(synpmx_avatar(
    source, roles, n_subjects = 100, seed = 22,
    subject_noise_sd = 0, residual_noise_sd = 0
  ))
  per_subject <- lapply(split(synthetic, synthetic$ID), function(x) {
    data.frame(dose = max(x$AMT), mean_dv = mean(x$DV[x$EVID == 0]))
  })
  summary <- aggregate(mean_dv ~ dose, do.call(rbind, per_subject), mean)
  expect_lt(summary$mean_dv[summary$dose == 100],
            summary$mean_dv[summary$dose == 200])
})

test_that("the caller's RNG state is left untouched", {
  source <- pmx_simulated_fixture(20)
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", dvid = "DVID", cmt = "CMT", mdv = "MDV")
  set.seed(123)
  before <- stats::runif(1)
  set.seed(123)
  invisible(suppressWarnings(synpmx_avatar(source, roles, seed = 5)))
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
                         keep = "DOS"),
    mavoglurant = pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                            evid = "EVID", cmt = "CMT", rate = "RATE",
                            mdv = "MDV", occasion = "OCC", keep = "DOSE",
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

    synthetic <- suppressWarnings(synpmx_avatar(source, roles, seed = 1))
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

# Redundant endpoint labels: a numeric key beside a character one. -----------

.two_label_source <- function(consistent = TRUE) {
  mk <- function(id) {
    t <- c(0, 2, 8)
    name_obs <- if (consistent) "cp" else if (id %% 2 == 0) "cp" else "pca"
    rbind(
      data.frame(ID = id, TIME = 0, DV = NA, AMT = 100, EVID = 1L,
                 YTYPE = 0L, NAME = "dose"),
      data.frame(ID = id, TIME = t, DV = exp(-0.1 * t), AMT = 0, EVID = 0L,
                 YTYPE = 1L, NAME = name_obs))
  }
  do.call(rbind, lapply(1:6, mk))
}

test_that("two consistent endpoint-key columns validate and are both kept", {
  source <- .two_label_source(consistent = TRUE)
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", dvid = c("YTYPE", "NAME"))
  expect_true(validate_pmx(source, roles)$valid)

  synthetic <- suppressWarnings(synpmx_avatar(source, roles, n_subjects = 4,
                                               seed = 1))
  expect_true(all(c("YTYPE", "NAME") %in% names(synthetic)))
  # The two labels stay a clean 1:1 mapping in the output.
  pairs <- unique(synthetic[synthetic$EVID == 0, c("YTYPE", "NAME")])
  expect_equal(nrow(pairs), length(unique(synthetic$NAME[synthetic$EVID == 0])))
})

test_that("inconsistent endpoint-key columns fail validation", {
  source <- .two_label_source(consistent = FALSE)
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", dvid = c("YTYPE", "NAME"))
  report <- validate_pmx(source, roles)
  expect_false(report$valid)
  consistency <- report$checks$status[report$checks$check == "dvid_consistency"]
  expect_identical(consistency, "error")
})

test_that("the primary dvid column drives endpoint grouping", {
  source <- .two_label_source(consistent = TRUE)
  # YTYPE first: grouping uses the numeric key, NAME rides along.
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", dvid = c("YTYPE", "NAME"))
  synthetic <- suppressWarnings(synpmx_avatar(source, roles, n_subjects = 4,
                                               seed = 1))
  expect_true(validate_pmx(synthetic, roles)$valid)
})
