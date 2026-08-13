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
  # `strata` is deliberately NOT in this set. AVATAR reads it as the
  # assigned stratum -- grouping the dose basis and the attendance-pattern pool
  # -- so declaring it is valid rather than an error. It is still not a blending
  # barrier; only route is.
  expect_silent(suppressMessages(suppressWarnings(
    synpmx_avatar(source, do.call(pmx_roles,
                                  c(base, list(strata = "SEX"))),
                  n_subjects = 5, seed = 1)
  )))
  expect_error(
    synpmx_avatar(source, do.call(pmx_roles,
                                  c(base, list(exclude = "AGE")))),
    "leave the column undeclared"
  )
})

test_that("AVATAR keeps donors dose-matched when the dose group meets the floor", {
  # Six subjects per dose, so each dose is its own signature group of size >= k
  # and no cross-dose borrowing is needed: the REV-019 dose-matching guarantee
  # still holds for well-populated groups. (Small groups deliberately borrow
  # across doses; see test-avatar-pooling.R.)
  mk <- function(id, dose) data.frame(
    ID = as.integer(id), TIME = c(0, 1, 2),
    DV = c(0, dose / 100, dose / 200),
    AMT = c(dose, 0, 0), EVID = c(1L, 0L, 0L), CMT = c(1L, 2L, 2L)
  )
  source <- do.call(rbind, c(
    lapply(1:6, mk, dose = 100),
    lapply(7:12, mk, dose = 200)
  ))
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

# Donor interpolation (SIM-046) ----------------------------------------------

test_that("a donor is not extrapolated past the times it was observed at", {
  # SIM-046. This used to rescale the target's range onto the donor's, so a
  # donor first measured at 24 h answered a 0.5 h query with its 24 h value --
  # an elimination concentration standing in for an absorption one. On
  # `warfarin`, where 19 of 32 subjects have no `cp` sample before 24 h, that
  # replaced the absorption limb with a flat plateau.
  donor <- list(time = c(24, 36, 48, 72), value = c(8.5, 7.5, 6.3, 4.4))
  out <- synpmx:::.interpolate_trajectory(donor, c(0.5, 3, 12, 24, 36, 96))

  # Before the donor's first observation, and after its last.
  expect_true(all(is.na(out[1:3])))
  expect_true(is.na(out[[6]]))
  # Its own observations, and interpolation between them, are untouched.
  expect_equal(out[[4]], 8.5)
  expect_equal(out[[5]], 7.5)
  expect_equal(
    synpmx:::.interpolate_trajectory(donor, 30)[[1L]], 8.0
  )
})

test_that("a donor with one observation speaks only for that time", {
  # A single point carries no shape, so it says nothing about any other time.
  # It used to be broadcast across the whole target vector.
  donor <- list(time = 24, value = 8.5)
  out <- synpmx:::.interpolate_trajectory(donor, c(0.5, 24, 96))
  expect_equal(out, c(NA, 8.5, NA))
})
