# What a fitted model releases, and how far one patient moves it. Everything
# here runs without a fitter except the two end-to-end studies at the bottom,
# which skip where nlmixr2 cannot build a model.

.warfarin_source <- function() {
  data <- as.data.frame(nlmixr2data::warfarin)
  data$ntime <- data$time
  data
}

.warfarin_roles <- function() {
  pmx_roles(id = "id", time = "time", nominal_time = "ntime", dv = "dv",
            amt = "amt", evid = "evid", dvid = "dvid",
            covariates = c("wt", "age", "sex"))
}

.stored_fits <- function() {
  files <- list.files(system.file("extdata", package = "synpmx"),
                      pattern = "model-fit[.]rds$", full.names = TRUE)
  stats::setNames(lapply(files, readRDS), basename(files))
}

# Every path in a nested list whose name says it holds random effects.
.eta_paths <- function(x, path = "fit") {
  if (!is.list(x) || is.data.frame(x)) return(character())
  found <- character()
  for (name in names(x)) {
    here <- paste0(path, "$", name)
    if (grepl("^eta", name)) found <- c(found, here)
    found <- c(found, .eta_paths(x[[name]], here))
  }
  found
}

test_that("no stored fit carries the individual random effects (REV-056)", {
  fits <- .stored_fits()
  skip_if(!length(fits), "no stored fits")
  for (name in names(fits)) {
    expect_identical(.eta_paths(unclass(fits[[name]])), character(),
                     info = name)
  }
})

test_that("the release holds only what generation reads (REV-057)", {
  fits <- .stored_fits()
  skip_if(!length(fits), "no stored fits")
  for (name in names(fits)) {
    release <- model_release(fits[[name]])
    expect_s3_class(release, "pmx_model_release")
    expect_s3_class(release, "pmx_fitted_model")
    expect_identical(.release_strays(release), character(), info = name)
    expect_null(release$candidates)
    expect_null(release$correlations)
    expect_null(release$timing)
    expect_null(release$design)
    expect_null(release$endpoints$signals)
    for (model in release$pk_models) {
      expect_null(model$movement$changes)
      expect_null(model$candidates)
    }
  }
  # The release is what is generated from: the same seed gives the same data.
  fit <- fits[["warfarin-model-fit.rds"]]
  skip_if(is.null(fit), "stored warfarin fit unavailable")
  a <- suppressWarnings(synpmx_model_generate(fit, seed = 3))
  b <- suppressWarnings(synpmx_model_generate(model_release(fit), seed = 3))
  attributes(a)$pmx_fitted_model <- attributes(b)$pmx_fitted_model <- NULL
  expect_identical(a, b)
  expect_error(model_candidates(model_release(fit)), "no candidate table")
})

test_that("generated data carries the release and no diagnostic (REV-057)", {
  stored <- system.file("extdata", "warfarin-model-fit.rds", package = "synpmx")
  skip_if(!nzchar(stored), "stored fit unavailable")
  fit <- readRDS(stored)
  synthetic <- suppressWarnings(synpmx_model_generate(fit, seed = 1))
  attached <- attr(synthetic, "pmx_fitted_model")
  expect_s3_class(attached, "pmx_model_release")
  expect_identical(.eta_paths(unclass(attached)), character())
  expect_null(attached$correlations)
  expect_null(attached$candidates)
  expect_null(attached$privacy$influence$change)
})

test_that("a factor subject ID leaves no label in the schema (REV-058)", {
  data <- .warfarin_source()
  labels <- sprintf("STUDY-%03d", data$id)
  data$id <- factor(labels)
  roles <- .warfarin_roles()
  schema <- .source_schema(data, roles, c("cp", "pca"),
                           rep("all", length(unique(data$id))))
  expect_null(schema$id_levels)
  expect_length(levels(schema$prototypes$id), 0L)
  ids <- .new_subject_ids(schema, 5L)
  expect_false(any(as.character(ids) %in% labels))
  # The PCA summary is built on the same schema, and no label reaches it.
  summary <- suppressWarnings(suppressMessages(
    synpmx_pca_summarize(data, roles, seed = 1)))
  text <- unlist(lapply(rapply(unclass(summary), function(x) {
    if (is.factor(x)) c(levels(x), as.character(x)) else as.character(x)
  }, classes = c("character", "factor"), how = "unlist"), identity))
  expect_false(any(labels %in% text))
})

test_that("a numeric ID keeps only an offset rounded up to a power of ten", {
  expect_equal(.id_offset(c(3, 32, 7)), 100)
  expect_equal(.id_offset(c(1001, 1042)), 10000)
  expect_equal(.id_offset(c(99, 12)), 100)
  expect_equal(.id_offset(letters), 0)
  schema <- list(id_class = "numeric", id_offset = .id_offset(1:32))
  expect_equal(.new_subject_ids(schema, 3L), c(101, 102, 103))
})

test_that("the emission floor is not any one patient's value (REV-059)", {
  expect_equal(vapply(c(0.3, 4.5, 7, 1, 25), .round_down_125, numeric(1)),
               c(0.2, 2, 5, 1, 20))
  data <- .warfarin_source()
  roles <- .warfarin_roles()
  floors <- .model_quantification_floor(data, roles, c("cp", "pca"), 3L)
  observed <- data$evid == 0 & data$dvid == "cp" & data$dv > 0
  smallest <- min(data$dv[observed])
  expect_false(isTRUE(all.equal(floors$cp, smallest / 2)))
  expect_equal(floors$cp, .round_down_125(floors$cp))
  # One patient reporting a value a thousand times lower than anyone else does
  # not move the floor: it is read at the third patient, not the first.
  low <- data
  first <- which(low$id == low$id[observed][1L] & low$dvid == "cp" &
                   low$evid == 0 & low$dv > 0)[1L]
  low$dv[first] <- smallest / 1000
  expect_equal(.model_quantification_floor(low, roles, "cp", 3L)$cp, floors$cp)
})

test_that("covariate summaries hold no median and resist one extreme patient", {
  data <- .warfarin_source()
  roles <- .warfarin_roles()
  clean <- .covariate_model(data, roles)
  expect_null(clean$wt$median)
  expect_null(clean$age$median)
  heavy <- data
  heavy$wt[heavy$id == heavy$id[1L]] <- 250
  shifted <- .covariate_model(heavy, roles)
  # Trimmed, the 250 kg patient barely moves the released distribution;
  # untrimmed, the SD on the log scale grew by more than half.
  expect_lt(abs(shifted$wt$meanlog - clean$wt$meanlog), 0.02)
  expect_lt(abs(shifted$wt$sdlog - clean$wt$sdlog), 0.02)
  first <- !duplicated(heavy$id)
  expect_gt(stats::sd(log(heavy$wt[first])) / stats::sd(log(data$wt[first])),
            1.5)
  # A level one patient holds is not released (REV-054).
  rare <- data
  rare$sex <- as.character(rare$sex)
  rare$sex[rare$id == rare$id[1L]] <- "unknown"
  levels_out <- .covariate_model(rare, roles)$sex$levels
  expect_false("unknown" %in% levels_out)
})

test_that("trimmed moments recover a normal spread and ignore the extremes", {
  set.seed(11)
  draws <- replicate(2000, .trimmed_moments(stats::rnorm(40))[["sd"]])
  expect_equal(mean(draws), 1, tolerance = 0.05)
  x <- c(stats::rnorm(39), 50)
  expect_lt(abs(.trimmed_moments(x)[["mean"]]), 0.5)
  expect_equal(.trimmed_moments(c(1, 2))[["mean"]], 1.5)
})

test_that("the allometric reference and planned doses are not one patient's value (REV-060)", {
  data <- .warfarin_source()
  roles <- pmx_roles(id = "id", time = "time", nominal_time = "ntime",
                     dv = "dv", amt = "amt", evid = "evid", dvid = "dvid",
                     covariates = "wt")
  weight <- .model_weight_covariate(data, roles)
  expect_equal(weight$reference, signif(weight$reference, 1L))
  # A mode only where several patients share the value.
  expect_equal(.shared_mode(c(100, 100, 100, 50), 3L), 100)
  individual <- c(101.3, 98.2, 140.1, 77.7)
  expect_equal(.shared_mode(individual, 3L), signif(mean(individual), 2L))
  expect_false(.shared_mode(individual, 3L) %in% individual)
})

test_that("a discrete level held by fewer than three patients is folded (REV-062)", {
  folded <- .fold_rare_levels(c(0, 0, 0, 1, 0, 1, 1, 2),
                              c("a", "b", "c", "d", "e", "f", "g", "h"), 3L)
  expect_equal(folded$levels, c(0, 1))
  expect_equal(folded$probability, c(5, 3) / 8)
  # Patients, not observations: one patient recorded twice is still one.
  twice <- .fold_rare_levels(c(0, 0, 0, 2, 2), c("a", "b", "c", "d", "d"), 3L)
  expect_equal(twice$levels, 0)
  expect_null(.fold_rare_levels(c(0, 1, 2), c("a", "b", "c"), 3L))
  expect_equal(.fold_rare_levels(c(0, 1, 2), c("a", "b", "c"), 1L)$levels,
               c(0, 1, 2))
})

test_that("single-patient influence fails a patient given a thousand times the dose (SIM-056)", {
  set.seed(7)
  etas <- data.frame(eta.cl = stats::rnorm(30, 0, 0.3),
                     eta.v = stats::rnorm(30, 0, 0.2))
  rownames(etas) <- paste0("S", seq_len(30))
  omega_of <- function(e) {
    out <- diag(colMeans(e^2))
    dimnames(out) <- list(c("cl", "v"), c("cl", "v"))
    out
  }
  clean <- .pk_influence("cp", etas, omega_of(etas), c(cl = 4, v = 40))
  expect_true(all(.influence_verdict(clean$change) == "pass"))
  # A thousand-fold exposure the record does not show lands in the random
  # effect on volume, and the variance that absorbs it is the patient's.
  etas$eta.v[7] <- log(1 / 1000)
  misdosed <- .pk_influence("cp", etas, omega_of(etas), c(cl = 4, v = 40))
  expect_equal(.influence_verdict(max(misdosed$change)), "FAIL")
  expect_equal(attr(misdosed, "who")[[which.max(misdosed$change)]], "S7")
  # The stored reading is per estimate; nothing in it is about a patient.
  expect_false(any(grepl("S7", unlist(misdosed))))
})

test_that("the privacy checks catch a fit built before these fixes", {
  stored <- system.file("extdata", "warfarin-model-fit.rds", package = "synpmx")
  skip_if(!nzchar(stored), "stored fit unavailable")
  fit <- readRDS(stored)
  checks <- model_privacy_checks(fit)
  expect_s3_class(checks, "pmx_privacy_checks")
  expect_identical(checks$check, paste0("P", 1:5))
  expect_true(all(checks$verdict[1:4] == "pass"))
  # The defects these checks exist for, put back by hand.
  old <- unclass(fit)
  old$quantification_floor <- list(cp = 0.3)
  old$covariates$wt$median <- 71.7
  old$schema$id_offset <- 32
  old$schema$prototypes$id <- factor(c("A", "B"))
  class(old) <- "pmx_fitted_model"
  flagged <- model_privacy_checks(old)
  expect_equal(flagged$verdict[flagged$check == "P2"], "FAIL")
  expect_equal(flagged$verdict[flagged$check == "P3"], "FAIL")
  stray <- model_release(fit)
  stray$pk_models[[1L]]$parameters$etas <- data.frame(eta.cl = 1:3)
  expect_match(paste(.release_strays(stray), collapse = " "), "etas")
})

test_that("scorecard E3 reads the release's influence verdict", {
  stored <- system.file("extdata", "warfarin-model-fit.rds", package = "synpmx")
  skip_if(!nzchar(stored), "stored fit unavailable")
  fit <- readRDS(stored)
  source <- .warfarin_source()
  roles <- .warfarin_roles()
  bad <- fit
  bad$privacy <- list(
    influence = .influence_row("PK", "cp: v", "between-subject SD", 0.7, 70,
                               "points", 32L),
    thresholds = .influence_thresholds)
  bad$privacy$influence$verdict <- "FAIL"
  synthetic <- suppressWarnings(synpmx_model_generate(bad, seed = 1))
  card <- as.data.frame(suppressMessages(
    synpmx_scorecard(source, synthetic, roles)))
  expect_equal(card$verdict[card$check == "E3"], "FAIL")
  expect_match(card$result[card$check == "E3"], "70")

  good <- fit
  good$privacy <- list(
    influence = .influence_row("PK", "cp: v", "between-subject SD", 0.2, 2,
                               "points", 32L),
    thresholds = .influence_thresholds)
  good$privacy$influence$verdict <- "pass"
  synthetic <- suppressWarnings(synpmx_model_generate(good, seed = 1))
  card <- as.data.frame(suppressMessages(
    synpmx_scorecard(source, synthetic, roles)))
  expect_equal(card$verdict[card$check == "E3"], "pass")
  # A passing release carries the threshold, not the number.
  expect_match(card$result[card$check == "E3"], "less than")
  expect_false(grepl("points", card$result[card$check == "E3"]))
})

test_that("E3 is not applicable to data no population model produced", {
  source <- private_fixture(24L)
  roles <- private_roles()
  synthetic <- suppressWarnings(suppressMessages(
    synpmx_avatar(source, roles, seed = 1)))
  plain <- as.data.frame(
    suppressMessages(synpmx_scorecard(source, synthetic, roles)))
  expect_equal(plain$verdict[plain$check == "E3"], "not applicable")
})

# End to end, with a population fitter --------------------------------------

.fitter_works <- local({
  answer <- NULL
  function() {
    if (!is.null(answer)) return(answer)
    answer <<- requireNamespace("nlmixr2est", quietly = TRUE) &&
      requireNamespace("rxode2", quietly = TRUE) &&
      !inherits(try(suppressMessages(suppressWarnings(
        rxode2::rxode2({ d / dt(central) <- -0.1 * central })
      )), silent = TRUE), "try-error")
    answer
  }
})

.privacy_study <- function(n = 24, misdose = NULL, seed = 5) {
  set.seed(seed)
  times <- c(0.5, 1, 2, 4, 8, 12, 24)
  pieces <- lapply(seq_len(n), function(subject) {
    p <- c(cl = 4 * exp(stats::rnorm(1, 0, 0.3)),
           v = 40 * exp(stats::rnorm(1, 0, 0.2)),
           ka = 1.1 * exp(stats::rnorm(1, 0, 0.3)))
    # The patient was given a thousand times the dose the record shows.
    given <- if (identical(subject, misdose)) 100000 else 100
    concentration <- .pk_profile(list(pk = "1cmt_oral"), times, given, 0, p)
    rows <- rbind(
      data.frame(TIME = 0, NTIME = 0, DV = 0, AMT = 100, EVID = 1L, CMT = 1L,
                 DVID = "cp", MDV = 1L),
      data.frame(TIME = times, NTIME = times,
                 DV = concentration * exp(stats::rnorm(length(times), 0, 0.1)),
                 AMT = 0, EVID = 0L, CMT = 2L, DVID = "cp", MDV = 0L))
    rows$ID <- sprintf("P%02d", subject)
    rows$WT <- 70 * exp(stats::rnorm(1, 0, 0.15))
    rows
  })
  out <- do.call(rbind, pieces)
  out$DVID <- factor(out$DVID)
  rownames(out) <- NULL
  out
}

.privacy_roles <- function() {
  pmx_roles(id = "ID", time = "TIME", nominal_time = "NTIME", dv = "DV",
            amt = "AMT", evid = "EVID", cmt = "CMT", dvid = "DVID",
            mdv = "MDV", covariates = "WT")
}

test_that("a fresh fit stores no per-patient quantity and names no patient", {
  skip_if_not(.fitter_works(), "no population fitter that can build a model")
  data <- .privacy_study()
  fit <- synpmx_model_estimate(data, .privacy_roles(), pk = "1cmt_oral",
                               quiet = TRUE, seed = 1)
  expect_identical(.eta_paths(unclass(fit)), character())
  # No patient identifier anywhere in the fit, as a value or a name.
  everything <- c(unlist(fit, use.names = FALSE),
                  names(unlist(fit)))
  expect_false(any(unique(data$ID) %in% everything))
  checks <- model_privacy_checks(fit)
  expect_true(all(checks$verdict[1:4] == "pass"))
  expect_false(identical(checks$verdict[checks$check == "P5"], "FAIL"))
  expect_true(is.data.frame(attr(checks, "influence")))
  expect_null(attr(fit$privacy$influence, "who"))
})

test_that("a patient misdosed a thousandfold fails P5 and is named at estimation", {
  skip_if_not(.fitter_works(), "no population fitter that can build a model")
  data <- .privacy_study(misdose = 9L)
  warnings <- character()
  fit <- withCallingHandlers(
    synpmx_model_estimate(data, .privacy_roles(), pk = "1cmt_oral",
                          quiet = TRUE, seed = 1),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("moves a released estimate", warnings)))
  expect_true(any(grepl("P09", warnings, fixed = TRUE)))
  checks <- model_privacy_checks(fit)
  expect_equal(checks$verdict[checks$check == "P5"], "FAIL")
  # The patient is named on the console and nowhere in the fit.
  expect_false(any(grepl("P09", unlist(fit), fixed = TRUE)))
})

test_that("the PCA summary's log offset is not one patient's value (REV-064)", {
  data <- .warfarin_source()
  summary <- suppressWarnings(suppressMessages(
    synpmx_pca_summarize(data, .warfarin_roles(), seed = 1)))
  offsets <- vapply(summary$basis$transforms, function(t) t$offset,
                    numeric(1))
  offsets <- offsets[offsets > 0]
  expect_true(length(offsets) > 0L)
  for (value in offsets) expect_equal(value, .round_down_125(value))
  observed <- data$evid == 0 & data$dvid == "cp" & data$dv > 0
  expect_false(isTRUE(all.equal(offsets[["cp"]], min(data$dv[observed]) / 2)))
})
