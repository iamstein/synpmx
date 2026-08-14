load_nlmixr2_dataset <- function(name) {
  environment <- new.env(parent = emptyenv())
  utils::data(list = name, package = "nlmixr2data", envir = environment)
  get(name, envir = environment, inherits = FALSE)
}

integration_budget <- function(censoring = 0) {
  pmx_budget_allocation(.1, .15, .15, .1, .5 - censoring, censoring)
}

test_that("theo_md runs end to end with repeated dose-relative PK", {
  skip_if_not_installed("nlmixr2data")
  source <- load_nlmixr2_dataset("theo_md")
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", cmt = "CMT", covariates = "WT")
  endpoints <- list(cp = pmx_endpoint(
    alignment = "dose_relative", transform = "log", shape = "occasion",
    cmt = 2
  ))
  model <- .fit_private(
    source, roles, endpoints, 5, 0,
    pmx_bounds(c(0, 170), list(cp = c(0, 30)), amt = c(0, 500),
               covariates = list(WT = c(40, 130))),
    pmx_public_design(
      pmx_schema(source), dose_evid = 101, dose_cmt = 1
    ),
    pmx_contribution_limits(40, 8, 8, 30, 11), integration_budget(),
    backend = "public", public_source = TRUE
  )
  synthetic <- .generate_private(model, seed = 42)
  expect_true(validate_pmx(synthetic, roles, endpoints)$valid)
  expect_equal(length(unique(synthetic$ID)), length(unique(source$ID)))
  expect_true(all(vapply(split(synthetic$EVID, synthetic$ID),
                         function(x) sum(x != 0) == 7L, logical(1))))
  expect_equal(model$population$event$n_doses, 7L)
  expect_equal(model$population$event$dose_interval, 24, tolerance = 0.1)
  expect_true(all(vapply(split(synthetic$EVID, synthetic$ID), function(x) {
    sum(x == 0) == as.integer(round(model$population$event$observation_count))
  }, logical(1))))
  observation_occasions <- lapply(split(synthetic, synthetic$ID), function(subject) {
    doses <- sort(subject$TIME[subject$EVID != 0])
    observations <- subject$TIME[subject$EVID == 0]
    pmax(1L, findInterval(observations, doses))
  })
  expect_true(all(vapply(observation_occasions, function(x) {
    counts <- as.integer(table(factor(x, levels = 1:7)))
    all(counts[3:6] == 0L) && counts[7] == 11L &&
      (identical(counts[c(1, 2)], c(10L, 1L)) ||
       identical(counts[c(1, 2)], c(11L, 0L)))
  }, logical(1))))
  expect_equal(
    round(model$population$timing$cp$occasion_observation_count),
    c(10, 1, 0, 0, 0, 0, 11, 0)
  )
  expect_equal(
    model$population$timing$cp$occasion_presence_probability,
    c(1, 10 / 12, 0, 0, 0, 0, 1, 0)
  )
  fitted_sampling <- sampling_summary(model)
  expect_equal(fitted_sampling$sampling_probability[1:7],
               c(1, 10 / 12, 0, 0, 0, 0, 1))
  expect_equal(round(fitted_sampling$observations_if_sampled[1:7]),
               c(10, 1, 0, 0, 0, 0, 11))
  cp <- synthetic[synthetic$EVID == 0, ]
  directional_peaks <- unlist(lapply(split(synthetic, synthetic$ID), function(subject) {
    doses <- sort(subject$TIME[subject$EVID != 0])
    observations <- subject[subject$EVID == 0, , drop = FALSE]
    occasion <- pmax(1L, findInterval(observations$TIME, doses))
    vapply(split(observations, occasion), function(profile) {
      values <- profile$DV[order(profile$TIME)]
      direction <- sign(diff(values))
      direction <- direction[direction != 0]
      sum(diff(direction) < 0)
    }, integer(1))
  }))
  expect_true(all(directional_peaks <= 1L))
  first <- cp[cp$ID == unique(cp$ID)[1] & cp$TIME < 12, ]
  expect_gt(max(first$DV), first$DV[1L])
  source_vectors <- split(source$TIME, source$ID)
  synthetic_vectors <- split(synthetic$TIME, synthetic$ID)
  expect_false(any(vapply(synthetic_vectors, function(x) {
    any(vapply(source_vectors, identical, logical(1), y = x))
  }, logical(1))))
})

test_that("warfarin preserves lower-case endpoint-specific schema", {
  skip_if_not_installed("nlmixr2data")
  source <- load_nlmixr2_dataset("warfarin")
  roles <- pmx_roles(id = "id", time = "time", dv = "dv", amt = "amt",
                     evid = "evid", dvid = "dvid",
                     covariates = c("wt", "age", "sex"))
  endpoints <- list(
    cp = pmx_endpoint("cp", "dose_relative", "log", "occasion"),
    pca = pmx_endpoint("pca", "study_time", "identity", "global")
  )
  model <- .fit_private(
    source, roles, endpoints, 5, 0,
    pmx_bounds(c(0, 144), list(cp = c(0, 25), pca = c(0, 120)),
               amt = c(0, 200),
               covariates = list(wt = c(40, 150), age = c(18, 100))),
    pmx_public_design(
      pmx_schema(source), dose_evid = 1
    ),
    pmx_contribution_limits(30, 2, 2, c(cp = 20, pca = 12), 12),
    integration_budget(), backend = "public", public_source = TRUE
  )
  synthetic <- .generate_private(model, seed = 42)
  expect_true(validate_pmx(synthetic, roles, endpoints)$valid)
  expect_equal(length(unique(synthetic$id)), length(unique(source$id)))
  expect_identical(names(synthetic), names(source))
  expect_identical(levels(synthetic$dvid), c("cp", "pca"))
  expect_identical(levels(synthetic$sex), levels(source$sex))
  expect_setequal(unique(as.character(synthetic$dvid[synthetic$evid == 0])),
                  c("cp", "pca"))
  source_cp <- source[source$evid == 0 & source$dvid == "cp", ]
  synthetic_cp <- synthetic[synthetic$evid == 0 & synthetic$dvid == "cp", ]
  expect_lte(abs(
    nrow(synthetic_cp) / length(unique(synthetic$id)) -
      nrow(source_cp) / length(unique(source$id))
  ), 1)
  synthetic_cp_max <- vapply(split(synthetic_cp$time, synthetic_cp$id), max, numeric(1))
  expect_gte(stats::median(synthetic_cp_max), 72)
  expect_true(any(synthetic_cp$time > 24))
})

test_that("wbcSim creates coherent generalized infusion and recovery", {
  skip_if_not_installed("nlmixr2data")
  source <- load_nlmixr2_dataset("wbcSim")
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", cmt = "CMT", rate = "RATE",
                     covariates = c("V2I", "V1I", "CLI"))
  endpoints <- list(wbc = pmx_endpoint(
    alignment = "study_time", transform = "log", shape = "global", cmt = 3
  ))
  model <- .fit_private(
    source, roles, endpoints, 5, 0,
    pmx_bounds(c(0, 720), list(wbc = c(0, 30)), amt = c(-200, 200),
               rate = c(-200, 200),
               covariates = list(V2I = c(100, 1500), V1I = c(100, 1200),
                                 CLI = c(100, 800))),
    pmx_public_design(
      pmx_schema(source), dose_evid = 10101, dose_cmt = 1,
      endpoint_cmt = list(wbc = 3)
    ),
    pmx_contribution_limits(20, 2, 2, 12, 9), integration_budget(),
    backend = "public", public_source = TRUE
  )
  synthetic <- .generate_private(model, seed = 42)
  event <- synthetic$EVID != 0
  expect_true(validate_pmx(synthetic, roles, endpoints)$valid)
  expect_equal(length(unique(synthetic$ID)), length(unique(source$ID)))
  expect_equal(synthetic$AMT[event], synthetic$RATE[event], tolerance = 1e-8)
  expect_true(any(synthetic$AMT[event] > 0) && any(synthetic$AMT[event] < 0))
  expect_false(any(synthetic$TIME == 4580))
  first <- synthetic[synthetic$ID == unique(synthetic$ID)[1] & synthetic$EVID == 0, ]
  expect_lt(min(first$DV), first$DV[1L])
  expect_gt(first$DV[nrow(first)], min(first$DV))
})

# SIM-040, on a public dataset. `theo_md` is genuinely dosed by weight and is
# exactly the case inference cannot recover: the recorded mg/kg runs from 3.1 to
# 5.9, which is eight ratio levels for eleven distinct amounts, so detection
# declines. Pinned here because it is the public evidence for the claim the
# demo and the algorithm vignette both make.
test_that("theo_md: inference declines weight-based dosing, declaring it does not", {
  skip_if_not_installed("nlmixr2data")
  source <- load_nlmixr2_dataset("theo_md")
  roles <- function(declare) {
    args <- list(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                 evid = "EVID", cmt = "CMT", covariates = "WT")
    if (declare) args$dose_covariate <- "WT"
    do.call(pmx_roles, args)
  }
  per_kg <- function(x) {
    dosed <- x[x$EVID != 0, ]
    dosed$AMT / dosed$WT
  }
  # Each patient's whole course, as the thing M5 has to keep coherent.
  profiles <- function(x) {
    dosed <- x[x$EVID != 0, ]
    vapply(split(round(per_kg(x), 4), dosed$ID),
           function(r) paste(r, collapse = "/"), character(1))
  }

  inferred <- suppressWarnings(suppressMessages(
    synpmx_avatar(source, roles(FALSE), seed = 303)
  ))
  declared <- suppressWarnings(suppressMessages(
    synpmx_avatar(source, roles(TRUE), seed = 303)
  ))

  # Inference declines, and says so rather than leaving a blank.
  expect_true(is.na(attr(inferred, "pmx_settings")$dose_basis))
  expect_match(attr(inferred, "pmx_settings")$dose_basis_note, "WT")
  # So every avatar carries its anchor's milligrams over its own blended
  # weight -- a dose per kilogram the study never prescribed.
  expect_lt(mean(profiles(inferred) %in% profiles(source)), 0.5)

  expect_equal(attr(declared, "pmx_settings")$dose_basis, "WT")
  expect_true(attr(declared, "pmx_settings")$dose_basis_declared)
  # Declared, every avatar is on a course some real patient was given.
  expect_equal(mean(profiles(declared) %in% profiles(source)), 1)
  # And the cohort spans the source's own mg/kg range rather than a shrunken
  # one, which is the visible symptom in the demo vignette.
  expect_equal(range(per_kg(declared)), range(per_kg(source)),
               tolerance = 1e-8)
})

# SIM-051, on the dataset that is now its worked case. `nimoData` records ten
# weekly infusions at their actual times, so no two subjects share an opening
# past the first dose or two. Before SIM-053 this read B1b = 12 FAIL with the
# dosing intact; it now passes every check with the dosing gone instead. Both
# halves are pinned, because the point of the row is that the quiet one is the
# one carrying the finding.
test_that("nimoData: the guarantee is reached by shortening the dosing", {
  skip_if_not_installed("nlmixr2data")
  source <- load_nlmixr2_dataset("nimoData")
  roles <- pmx_roles(
    id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
    rate = "RATE", mdv = "MDV", tad = "TAD", occasion = "OCC",
    covariates = c("BSA", "AGE", "HGT"), keep = "DOS"
  )
  synthetic <- suppressWarnings(suppressMessages(
    synpmx_avatar(source, roles, seed = 606)
  ))
  settings <- attr(synthetic, "pmx_settings")
  doses_per_patient <- function(x) mean(table(x$ID[x$EVID != 0]))

  # Nothing fails, and B1b in particular no longer does.
  card <- synpmx_scorecard(source, synthetic, roles)
  expect_false(any(card$verdict == "FAIL"))
  expect_equal(settings$identifying_dose_schedules, 0L)

  # The price is the whole of the dosing, and A5b is the only row that says so.
  expect_equal(doses_per_patient(source), 10)
  expect_lt(doses_per_patient(synthetic), 3)
  expect_equal(settings$dose_truncated_fraction, 1)

  # Constructing the protocol grid the study was actually run on recovers it,
  # with the privacy rows unchanged. This is the fix the vignette walks through.
  nominal <- source
  nominal$NTIME <- (nominal$OCC - 1) * 168 +
    ifelse(nominal$EVID == 0, round(nominal$TAD / 24) * 24, 0)
  roles_nominal <- pmx_roles(
    id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
    rate = "RATE", mdv = "MDV", tad = "TAD", occasion = "OCC",
    nominal_time = "NTIME",
    covariates = c("BSA", "AGE", "HGT"), keep = "DOS"
  )
  fixed <- suppressWarnings(suppressMessages(
    synpmx_avatar(nominal, roles_nominal, seed = 606)
  ))
  expect_equal(doses_per_patient(fixed), 10)
  expect_equal(attr(fixed, "pmx_settings")$identifying_dose_schedules, 0L)
  expect_equal(attr(fixed, "pmx_settings")$dose_truncated_fraction, 0)
})

# Coverage the nlmixr2data five do not give us. None of theo_md, warfarin,
# wbcSim, nimoData or mavoglurant has a declared `nominal_time`, a `cens` role,
# more than two endpoints, or treatment arms to stratify on -- which is to say
# every feature a real study report exercises was untested on public data.
# `xgxr` is already a Suggests dependency and its datasets have all four.
test_that("xgxr case1_pkpd: nominal grid, two endpoints, six arms, 180 patients", {
  skip_if_not_installed("xgxr")
  source <- as.data.frame(
    get(utils::data(list = "case1_pkpd", package = "xgxr"))
  )
  roles <- pmx_roles(
    id = "ID", time = "TIME", dv = "LIDV", amt = "AMT", evid = "EVID",
    cmt = "CMT", dvid = "NAME", nominal_time = "NOMTIME",
    strata = c("TRTACT", "DOSE"), covariates = "WEIGHTB",
    keep = "STUDY"
  )
  expect_true(validate_pmx(source, roles)$valid)
  synthetic <- suppressWarnings(suppressMessages(
    synpmx_avatar(source, roles, seed = 1)
  ))
  settings <- attr(synthetic, "pmx_settings")

  expect_true(validate_pmx(synthetic, roles)$valid)
  # A declared nominal time is the case the whole coarsening design is aimed
  # at, and it should leave nothing exposed on the schedule axis.
  expect_equal(settings$time_grid, "nominal")
  expect_equal(settings$unique_schedule_n, 0L)
  expect_equal(settings$identifying_visit_sets, 0L)
  expect_equal(settings$identifying_dose_schedules, 0L)
  # Six arms, and an avatar never leaves the arm it was anchored in.
  expect_equal(length(unique(synthetic$TRTACT)),
               length(unique(source$TRTACT)))
})

test_that("xgxr mad: five observation endpoints, ordinal, count and binary", {
  skip_if_not_installed("xgxr")
  source <- as.data.frame(get(utils::data(list = "mad", package = "xgxr")))
  roles <- pmx_roles(
    id = "ID", time = "TIME", dv = "LIDV", amt = "AMT", evid = "EVID",
    cmt = "CMT", dvid = "NAME", mdv = "MDV", nominal_time = "NOMTIME",
    strata = c("TRTACT", "DOSE"),
    covariates = c("WEIGHTB", "SEX")
  )
  synthetic <- suppressWarnings(suppressMessages(
    synpmx_avatar(source, roles, seed = 2)
  ))
  expect_true(validate_pmx(synthetic, roles)$valid)
  # Nothing is lost: every endpoint in is an endpoint out. Endpoint loss is the
  # `SIM-036` failure mode and it is invisible without a multi-endpoint source.
  expect_setequal(unique(as.character(synthetic$NAME)),
                  unique(as.character(source$NAME)))
  expect_equal(attr(synthetic, "pmx_settings")$identifying_visit_sets, 0L)
})

test_that("pheno_sd: dosing that can only be masked by throwing it away", {
  skip_if_not_installed("nlmixr2data")
  source <- load_nlmixr2_dataset("pheno_sd")
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", covariates = c("WT", "APGR"))
  synthetic <- suppressWarnings(suppressMessages(
    synpmx_avatar(source, roles, seed = 3)
  ))
  settings <- attr(synthetic, "pmx_settings")
  expect_true(validate_pmx(synthetic, roles)$valid)
  # Both guarantees hold. The dose side used to leak here -- neonatal
  # phenobarbital is dosed to the individual infant, so 55 of 59 have a dose
  # schedule nobody shares -- and re-anchoring now finds the few who can be
  # built on safely instead of rejection-sampling twelve times and giving up.
  expect_equal(settings$identifying_visit_sets, 0L)
  expect_equal(settings$identifying_dose_schedules, 0L)

  # What that costs, pinned so it cannot quietly get worse or be forgotten.
  # No infant's complete schedule is reusable, but many share an OPENING: q12h
  # from time zero for several doses before the individual course diverges. So
  # about half the dosing survives -- 10 doses per patient in, about 5.5 out,
  # over roughly 38 of 56 regimens. It used to be 1.1 doses over 3 regimens,
  # because the truncation grid was the union of every dose time in the cohort
  # and this study has no single leading grid (`SIM-053`).
  doses_per_patient <- function(x) {
    mean(table(x$ID[x$EVID != 0]))
  }
  expect_gt(doses_per_patient(source), 9)
  expect_gt(doses_per_patient(synthetic), 4)
  expect_lt(doses_per_patient(synthetic), doses_per_patient(source))
  expect_gt(settings$dose_regimens_represented, 25L)
})

test_that("pheno_sd: no avatar wears one infant's opening", {
  # The truncation above is only legitimate because every schedule it emits is
  # one several infants could have produced. Checked on the source, where the
  # comparison is exact: generated times carry resampled deviations, so the
  # finished table cannot be compared to the grid the plan reasoned on.
  skip_if_not_installed("nlmixr2data")
  source <- load_nlmixr2_dataset("pheno_sd")
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", covariates = c("WT", "APGR"))
  coarsened <- synpmx:::.coarsen_source_time(source, roles)$source
  subjects <- as.character(synpmx:::.unique_in_order(coarsened$ID))
  dosed <- coarsened[coarsened$EVID != 0, ]
  times <- lapply(subjects, function(s) {
    sort(unique(as.numeric(dosed$TIME[as.character(dosed$ID) == s])))
  })
  plan <- synpmx:::.dose_truncation_plan(coarsened, roles, 2L)
  kept <- ifelse(is.na(plan$target), plan$depth, plan$target)

  key <- function(t) paste(sprintf("%.12g", t), collapse = ",")
  complete <- table(vapply(times, key, character(1)))
  opened <- table(unlist(lapply(times, function(t) {
    vapply(seq_along(t), function(d) key(t[seq_len(d)]), character(1))
  })))
  count <- function(tab, k) {
    j <- match(k, names(tab))
    if (is.na(j)) 0L else as.integer(tab[[j]])
  }

  safe <- vapply(seq_along(times), function(i) {
    k <- key(times[[i]][seq_len(kept[[i]])])
    n <- count(complete, k)
    n >= 2L || (n == 0L && count(opened, k) >= 2L)
  }, logical(1))
  expect_true(all(safe))
  # And the truncation is doing real work rather than passing everything back.
  expect_true(any(kept < lengths(times)))
})
