# Estimation. Everything here needs `nlmixr2` and a toolchain that can link a
# compiled model, so every test skips where either is missing. The parts of the
# generator that do not need a fitter are covered in `test-model-generate.R` and
# `test-model-design.R`, which is the whole reason the fitter arrives last.

# Compiling once and caching the answer: `nlmixr2` being installed is not the
# same as it being able to build a model, and the difference is a linker
# failure that reports itself as a missing C compiler.
.population_fitting_works <- local({
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

skip_without_fitter <- function() {
  testthat::skip_if_not(.population_fitting_works(),
                        "no population fitter that can build a model")
}

.estimate_roles <- function(...) {
  pmx_roles(id = "ID", time = "TIME", nominal_time = "NTIME", dv = "DV",
            amt = "AMT", evid = "EVID", cmt = "CMT", dvid = "DVID",
            mdv = "MDV", ...)
}

# A one-compartment oral study simulated from known parameters, so that what the
# search should recover is known rather than merely plausible.
.oral_study <- function(n = 30, cl = 4, v = 40, ka = 1.1, seed = 42) {
  set.seed(seed)
  times <- c(0.25, 0.5, 1, 2, 4, 8, 12, 24)
  pieces <- lapply(seq_len(n), function(subject) {
    p <- c(cl = cl * exp(stats::rnorm(1, 0, 0.3)),
           v = v * exp(stats::rnorm(1, 0, 0.2)),
           ka = ka * exp(stats::rnorm(1, 0, 0.3)))
    concentration <- .pk_profile(list(pk = "1cmt_oral"), times, 100, 0, p)
    rows <- rbind(
      data.frame(TIME = 0, NTIME = 0, DV = 0, AMT = 100, EVID = 1L, CMT = 1L,
                 DVID = "cp", MDV = 1L),
      data.frame(TIME = times, NTIME = times,
                 DV = concentration * (1 + stats::rnorm(length(times), 0, 0.1)),
                 AMT = 0, EVID = 0L, CMT = 2L, DVID = "cp", MDV = 0L)
    )
    rows$ID <- subject
    rows$WT <- 70 * exp(stats::rnorm(1, 0, 0.15))
    rows
  })
  out <- do.call(rbind, pieces)
  out$DVID <- factor(out$DVID)
  out$DV <- pmax(out$DV, 0.001)
  rownames(out) <- NULL
  out
}

# The default path performs exactly one population fit, which is the design
# point rather than an accident: a search costs a fit per candidate and testing
# allometry against AIC costs another. Fits are still slow enough that the ones
# the tests share are built once and reused.
.shared <- local({
  cache <- list()
  function(name, expression) {
    if (is.null(cache[[name]])) cache[[name]] <<- expression
    cache[[name]]
  }
})

.default_fit <- function() {
  .shared("default", synpmx_model_estimate(.oral_study(), .estimate_roles(),
                                           quiet = TRUE,
                                           covariate_effects = "none"))
}

test_that("the default fits one model and no more", {
  skip_without_fitter()
  # One row, because one model was fitted. Not a search that happened to have
  # one candidate: `pk` is what asks for a search.
  expect_identical(nrow(model_candidates(.default_fit())), 1L)
  expect_true(model_candidates(.default_fit())$converged)
})

test_that("the fit recovers the model and parameters it was simulated from", {
  skip_without_fitter()
  fit <- .default_fit()
  expect_s3_class(fit, "pmx_fitted_model")
  expect_identical(fit$structural, "1cmt_oral")
  expect_equal(unname(fit$parameters$fixed[["cl"]]), 4, tolerance = 0.3)
  expect_equal(unname(fit$parameters$fixed[["v"]]), 40, tolerance = 0.3)
  expect_identical(fit$endpoints$pk, "cp")
})

test_that("every model fitted is in the table, converged or not", {
  skip_without_fitter()
  fit <- .default_fit()
  table <- model_candidates(fit)
  expect_true(all(c("model", "converged", "aic", "note") %in% names(table)))
  expect_true(nrow(table) >= 1L)
  # A search that came down to one survivor must not look like a search that
  # had one candidate, so failures keep their row and their reason.
  expect_true(all(is.na(table$aic) | is.finite(table$aic)))
  expect_true(all(nzchar(table$note[!table$converged])))
  expect_identical(fit$structural, table$model[which.min(table$aic)])
})

test_that("`pk` naming several models is how a search is asked for", {
  skip_without_fitter()
  fit <- synpmx_model_estimate(.oral_study(), .estimate_roles(),
                               pk = c("1cmt_oral", "1cmt_iv"), quiet = TRUE,
                               covariate_effects = "none")
  expect_identical(nrow(model_candidates(fit)), 2L)
  expect_match(fit$design$reason, "searched over")
  expect_identical(fit$structural,
                   model_candidates(fit)$model[which.min(model_candidates(fit)$aic)])
})

test_that("`pk` forces one model", {
  skip_without_fitter()
  data <- .oral_study()
  fit <- synpmx_model_estimate(data, .estimate_roles(), pk = "1cmt_iv",
                               quiet = TRUE, covariate_effects = "none")
  expect_identical(fit$structural, "1cmt_iv")
  expect_identical(nrow(model_candidates(fit)), 1L)
  expect_match(fit$design$reason, "declared")
})

test_that("a two-compartment model is available by asking for it", {
  skip_without_fitter()
  data <- .oral_study()
  # Not fitted by default -- it costs about five times a one-compartment fit --
  # but nothing stops a caller who wants one.
  expect_false("2cmt_oral" %in% .default_fit()$candidates$model)
  fit <- synpmx_model_estimate(data, .estimate_roles(), pk = "2cmt_oral",
                               quiet = TRUE, covariate_effects = "none")
  expect_identical(fit$structural, "2cmt_oral")
  expect_true(all(c("q", "v2") %in% names(fit$parameters$fixed)))
  synthetic <- synpmx_model_generate(fit, n_subjects = 10, seed = 2)
  expect_true(validate_pmx(synthetic, .estimate_roles())$valid)
})

test_that("`pk` naming a model outside the closed-form set is refused", {
  data <- .oral_study()
  expect_error(
    synpmx_model_estimate(data, .estimate_roles(), pk = "3cmt_oral",
                          quiet = TRUE),
    "outside the closed-form set"
  )
})

test_that("allometric scaling is asserted, not tested, and costs no extra fit", {
  skip_without_fitter()
  auto <- .shared("auto", synpmx_model_estimate(
    .oral_study(), .estimate_roles(covariates = "WT"), quiet = TRUE,
    covariate_effects = "auto"))
  # Applied because a weight-like covariate is declared, with the standard
  # exponents. Testing it against a model without it would double the cost of
  # the only fit the default path performs.
  expect_identical(auto$covariate_effects$cl$covariate, "WT")
  expect_equal(auto$covariate_effects$cl$exponent, 0.75)
  expect_equal(auto$covariate_effects$v$exponent, 1)
  expect_identical(nrow(model_candidates(auto)), 1L)
  expect_length(.default_fit()$covariate_effects, 0L)
})

test_that("no weight-like covariate means no scaling, and still one fit", {
  skip_without_fitter()
  # `AGE` is declared and positive but is not a weight, and nothing here tries
  # to recognise a body weight from its values.
  data <- .oral_study()
  data$AGE <- 40
  fit <- synpmx_model_estimate(data, .estimate_roles(covariates = "AGE"),
                               quiet = TRUE)
  expect_length(fit$covariate_effects, 0L)
  expect_identical(nrow(model_candidates(fit)), 1L)
})

test_that("covariate_effects only takes the two documented values", {
  data <- .oral_study()
  expect_error(
    synpmx_model_estimate(data, .estimate_roles(), covariate_effects = "all",
                          quiet = TRUE),
    "\"auto\" or \"none\""
  )
})

test_that("the fitted model carries no per-subject quantity", {
  skip_without_fitter()
  fit <- .shared("auto", synpmx_model_estimate(
    .oral_study(), .estimate_roles(covariates = "WT"), quiet = TRUE,
    covariate_effects = "auto"))
  expect_null(fit$parameters$etas)
  expect_false(any(grepl("eta", names(fit$parameters))))
  # The random effects are read to report a correlation and then discarded.
  expect_true(is.null(fit$correlations) ||
                all(c("covariate", "parameter", "correlation") %in%
                      names(fit$correlations)))
})

test_that("estimation reads recorded times, not the nominal grid", {
  skip_without_fitter()
  data <- .oral_study()
  data$TIME <- ifelse(data$EVID == 0L, data$TIME * 2, data$TIME)
  fit <- synpmx_model_estimate(data, .estimate_roles(), pk = "1cmt_oral",
                               quiet = TRUE, covariate_effects = "none")
  # Doubling the recorded clock halves the apparent elimination rate. Had the
  # fit read `NTIME` the estimate would not move.
  expect_lt(fit$parameters$fixed[["cl"]], .default_fit()$parameters$fixed[["cl"]])
})

# How long it took is reported, because it is the number a caller weighs a
# rerun against.
test_that("the fit carries the time it took, per candidate and in total", {
  skip_without_fitter()
  fit <- .default_fit()
  table <- model_candidates(fit)
  expect_true("seconds" %in% names(table))
  expect_true(all(table$seconds > 0))
  expect_gte(fit$timing$total, fit$timing$fit)
  expect_gte(fit$timing$fit, sum(table$seconds) - 1e-6)
  expect_output(print(model_report(fit)), "time to fit")
  # Printing the object prints the report, so the wait is on both.
  expect_output(print(fit), "time to fit")
})

test_that("a whole study round-trips through estimate and generate", {
  skip_without_fitter()
  data <- .oral_study()
  roles <- .estimate_roles(covariates = "WT")
  synthetic <- synpmx_model(data, roles, n_subjects = 25, seed = 8,
                            pk = "1cmt_oral", quiet = TRUE)
  expect_true(validate_pmx(synthetic, roles)$valid)
  expect_length(unique(synthetic$ID), 25L)
  observed <- synthetic$DV[synthetic$EVID == 0L]
  source_observed <- data$DV[data$EVID == 0L]
  expect_equal(stats::median(observed), stats::median(source_observed),
               tolerance = 0.5)
  expect_s3_class(attr(synthetic, "pmx_fitted_model"), "pmx_fitted_model")
})

test_that("the accessors return what the object holds and nothing else", {
  skip_without_fitter()
  fit <- .default_fit()
  expect_setequal(names(model_parameters(fit)),
                  c("fixed", "omega", "residual"))
  expect_s3_class(model_report(fit), "pmx_model_report")
  out <- paste(utils::capture.output(print(model_report(fit))), collapse = " ")
  expect_match(out, "Estimated by nlmixr2")
  expect_match(out, "Summarized from the source, not estimated")
})

test_that("the accessors refuse anything that is not a fitted model", {
  expect_error(model_report(list()))
  expect_error(model_candidates(list()))
  expect_error(model_parameters(list()))
})

test_that("a baseline weight is recognised however the study spelled it", {
  # A weight column is rarely called `WT`. Refusing `WEIGHTB` means silently
  # fitting no covariate at all on a study that declared one, which is how
  # `xgxr::case1_pkpd` went through with no scaling.
  recognised <- function(name) {
    roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", evid = "EVID",
                       covariates = name)
    data <- data.frame(ID = 1, TIME = 0, DV = 1, EVID = 0)
    data[[name]] <- 70
    !is.null(.model_weight_covariate(data, roles))
  }
  for (name in c("WT", "wt", "WEIGHT", "WEIGHTB", "BW", "BWT", "WTBL",
                 "WEIGHTBL", "WT0", "BWEIGHT")) {
    expect_true(recognised(name), info = name)
  }
  # And the other direction matters more: scaling clearance by a height or a
  # cell count would be worse than scaling it by nothing.
  for (name in c("AGE", "HEIGHT", "SEX", "WBC", "EGFR", "CRCL")) {
    expect_false(recognised(name), info = name)
  }
})

test_that("a non-positive weight is not scaled on", {
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", evid = "EVID",
                     covariates = "WT")
  data <- data.frame(ID = 1:2, TIME = 0, DV = 1, EVID = 0, WT = c(70, 0))
  expect_null(.model_weight_covariate(data, roles))
})

test_that("the non-compartmental reading stays inside one dose interval", {
  # Pooling cycles is what made this wrong: sorting a multiple-dose subject's
  # samples by time after dose interleaves cycle 1 with cycle 85, and the area
  # and terminal slope of that sequence describe nothing. On `case1_pkpd` it
  # read a half-life of 371 hours.
  times <- c(0.25, 1, 2, 4, 8, 12)
  cycle <- function(start, scale) {
    data.frame(actual_tad = times, tad = times, dv = scale * exp(-0.1 * times),
               interval = start, endpoint = "cp",
               subject = "1", time = start * 24 + times,
               first_dose_amt = 100, ntime = start * 24 + times,
               first_dose_time = 0)
  }
  pooled <- rbind(cycle(1, 10), cycle(20, 30))
  one_interval <- .model_nca_subject(pooled[pooled$interval == 1, ])
  interleaved <- .model_nca_subject(pooled)
  # The pooled reading is not a profile: its peak sits in the wrong cycle.
  expect_equal(one_interval$cmax, 10 * exp(-0.1 * 0.25), tolerance = 1e-8)
  expect_gt(interleaved$cmax, one_interval$cmax)
  # And the estimator reads one interval, so it recovers the real slope.
  estimates <- .model_initial_estimates(pooled, "1cmt_iv", "cp")
  expect_equal(unname(estimates[["cl"]] / estimates[["v"]]), 0.1,
               tolerance = 0.05)
})

# SIM-065. A design with fewer than three samples in any one dose interval has
# no non-compartmental reading, and what stood in for one was a constant.
test_that("starting values come from the data where NCA cannot be read", {
  cl <- 0.006
  v <- 1.3
  dose_times <- c(0, seq(12, 120, by = 12))
  amounts <- c(25, rep(3.5, length(dose_times) - 1L))
  n <- 40
  sample_at <- seq(6, 130, length.out = n)
  data <- do.call(rbind, lapply(seq_len(n), function(i) {
    rbind(
      data.frame(TIME = dose_times, NTIME = dose_times, DV = NA_real_,
                 AMT = amounts, EVID = 1L, CMT = 1L, DVID = "cp", MDV = 1L,
                 ID = i),
      data.frame(TIME = sample_at[i], NTIME = sample_at[i],
                 DV = .pk_profile(list(pk = "1cmt_iv"), sample_at[i], amounts,
                                  dose_times, c(cl = cl, v = v)),
                 AMT = 0, EVID = 0L, CMT = 2L, DVID = "cp", MDV = 0L, ID = i)
    )
  }))
  data$DVID <- factor(data$DVID)
  roles <- .estimate_roles()
  observations <- .model_observations(data, roles)
  # Nobody has three samples in one interval, so every NCA quantity is missing.
  expect_true(all(is.na(unlist(.model_nca_subject(
    observations[observations$subject == "1", ])))))

  estimates <- .model_initial_estimates(observations, "1cmt_iv", "cp")
  # Within a factor of five of the truth, against the constants that stood here
  # before -- clearance 1 and volume 10, which are 170x and 8x out.
  expect_lt(abs(log(estimates[["cl"]] / cl)), log(5))
  expect_lt(abs(log(estimates[["v"]] / v)), log(5))
})

# SIM-061. The floor is half the smallest value the study reported, and only
# for an endpoint that lives on a positive scale.
test_that("the quantification floor is half the smallest reported value", {
  data <- .oral_study()
  roles <- .estimate_roles()
  observed <- data$EVID == 0L & !is.na(data$DV)
  endpoint <- unique(as.character(data$DVID[observed]))[1L]
  floor_value <- .model_quantification_floor(data, roles, endpoint)
  expect_equal(floor_value[[endpoint]],
               min(data$DV[observed & data$DVID == endpoint]) / 2)
})

# SIM-064. One non-positive reading is what an assay returns near its limit,
# not a statement that the endpoint reaches zero.
test_that("a lone non-positive value leaves the floor in place", {
  data <- .oral_study()
  roles <- .estimate_roles()
  observed <- which(data$EVID == 0L & !is.na(data$DV))
  endpoint <- as.character(data$DVID[observed[1L]])
  smallest <- min(data$DV[observed])
  data$DV[observed[1L]] <- -0.01
  floor_value <- .model_quantification_floor(data, roles, endpoint)
  expect_equal(floor_value[[endpoint]], smallest / 2)
})

# SIM-074. A study dosed both ways is one study, and the route is a property of
# each dose record. Before `adm`/`routes` existed, one nonzero rate anywhere
# made the whole cohort an infusion study and the subcutaneous patients were
# fitted with a model that has no absorption in it.
test_that("a declared administration column routes each dose", {
  data <- .oral_study(n = 20)
  data$ADM <- 2L
  half <- as.integer(data$ID) <= 10L
  data$ADM[half] <- 1L
  roles <- .estimate_roles(adm = "ADM",
                           routes = c("1" = "iv", "2" = "extravascular"))

  expect_setequal(.study_routes(data, roles), c("iv", "extravascular"))
  design <- .model_detect_design(data, roles, .model_observations(data, roles),
                                 "cp")
  expect_identical(design$route, "mixed")
  expect_identical(design$candidates, "1cmt_mixed")
  expect_match(design$reason, "declared through `adm` and `routes`")

  # The solver is told which compartment each dose enters: the depot for the
  # extravascular ones, central for the rest.
  table <- .model_estimation_data(data, roles, "cp")
  doses <- table[table$EVID != 0L, , drop = FALSE]
  routes <- .dose_routes(data, roles)[.dose_rows(data, roles)]
  expect_identical(doses$CMT[routes == "extravascular"],
                   rep(1L, sum(routes == "extravascular")))
  expect_identical(doses$CMT[routes == "iv"], rep(2L, sum(routes == "iv")))
  expect_true(all(table$CMT[table$EVID == 0L] == 2L))

  # One route declared is not a mixed study, and keeps its own model.
  single <- data
  single$ADM <- 2L
  expect_identical(.model_detect_design(single, roles,
                                        .model_observations(single, roles),
                                        "cp")$candidates, "1cmt_oral")
})

# Bioavailability scales the extravascular doses and not the intravenous ones,
# which is what it means and why a study dosed one way cannot identify it.
test_that("the mixed model routes each dose to its own form", {
  p <- c(cl = 4, v = 40, ka = 0.4, f = 0.7)
  at <- c(0.5, 2, 8, 24)
  iv <- .pk_profile(list(pk = "1cmt_mixed"), at, 100, 0, p, routes = "iv")
  ev <- .pk_profile(list(pk = "1cmt_mixed"), at, 100, 0, p,
                    routes = "extravascular")
  # The extravascular form is the oral one, bioavailability and all.
  expect_equal(ev, .pk_profile(list(pk = "1cmt_oral"), at, 100, 0, p))
  # The intravenous form is the bolus at the full dose, not 70% of it.
  expect_equal(iv, .pk_profile(list(pk = "1cmt_iv"), at, 100, 0,
                               replace(p, "f", 1)))
  # One patient given both receives both.
  expect_equal(.pk_profile(list(pk = "1cmt_mixed"), at, c(100, 100), c(0, 0), p,
                           routes = c("iv", "extravascular")), iv + ev)
})

test_that("the mixed model fits `f` on the depot and gives it no eta", {
  spec <- deparse(.model_nlmixr_function("1cmt_mixed",
                                         c(cl = 4, v = 40, ka = 0.4, f = 0.7),
                                         "prop", 0.2, NULL))
  expect_true(any(grepl("f(depot) <- f", spec, fixed = TRUE)))
  expect_true(any(grepl("f <- exp(tf)", spec, fixed = TRUE)))
  expect_false(any(grepl("eta.f", spec, fixed = TRUE)))
  expect_true(any(grepl("eta.ka", spec, fixed = TRUE)))
})

# SIM-069. There is no minimum number of observations. An endpoint measured
# once is a level, and a level is still something to generate -- the alternative
# is a synthetic study silently missing an endpoint the source has.
test_that("an endpoint measured once is generated as a level", {
  one <- data.frame(subject = "1", endpoint = "biomarker", dv = 42,
                    aligned = 0, stringsAsFactors = FALSE)
  shape <- .model_fit_pd(one, "biomarker")
  expect_identical(shape$pd, "constant")
  expect_equal(shape$typical[["baseline"]], 42)
  # One observation shows no scatter. `sd()` of one number is NA, and an NA
  # residual reaches generation as an NA observation.
  expect_equal(shape$residual$sd, 0)
  expect_equal(shape$baseline_cv, 0)
  expect_false(any(shape$candidates$converged &
                     shape$candidates$shape != "constant"))
})

# A shape needs something left over to be wrong about. Two points fit a line
# exactly, and `AIC` of a saturated fit is -Inf, so a line would win every
# two-point comparison and generate a curve through both with no residual.
test_that("a saturated shape is not a candidate", {
  two <- data.frame(subject = c("1", "2"), endpoint = "b", dv = c(5, 9),
                    aligned = c(0, 24), stringsAsFactors = FALSE)
  shape <- .model_fit_pd(two, "b")
  expect_identical(shape$pd, "constant")
  expect_false("linear" %in% shape$candidates$shape[shape$candidates$converged])

  # A third point at a third time gives the line a degree of freedom, and it
  # becomes a candidate rather than a certainty.
  three <- rbind(two, data.frame(subject = "3", endpoint = "b", dv = 7,
                                 aligned = 48, stringsAsFactors = FALSE))
  expect_true("linear" %in%
                .model_fit_pd(three, "b")$candidates$shape[
                  .model_fit_pd(three, "b")$candidates$converged])
})

# What is genuinely unfittable is an endpoint with no usable row at all: every
# value missing, or no time to fit against. The gates above make that hard to
# reach through `synpmx_model_estimate()` -- `nominal_time` is required on every
# observation before this runs, and an undosed subject keeps its own clock -- so
# the contract is pinned here and the call site reports it if anything ever
# produces one.
test_that("an endpoint with no usable observation is fitted to nothing", {
  none <- data.frame(subject = c("1", "2"), endpoint = "b", dv = c(5, 9),
                     aligned = c(NA_real_, NA_real_), stringsAsFactors = FALSE)
  expect_null(.model_fit_pd(none, "b"))
  missing_values <- data.frame(subject = "1", endpoint = "b", dv = NA_real_,
                               aligned = 0, stringsAsFactors = FALSE)
  expect_null(.model_fit_pd(missing_values, "b"))
})

# SIM-071. The fit table is ordered by the study's own subject order, not by how
# the identifiers sort as text. `ID` reaches the solver as character, so sorting
# on it puts subject 10 before subject 2 whenever a study numbers its patients,
# and `focei` answers a different subject order with a different optimum.
test_that("the fit table follows the study's subject order, not the text one", {
  data <- .oral_study(n = 12)
  roles <- .estimate_roles()
  table <- .model_estimation_data(data, roles, "cp")
  expect_identical(.unique_in_order(table$ID),
                   as.character(.unique_in_order(data$ID)))
  # Which is not what sorting the identifiers as text would give.
  expect_false(identical(.unique_in_order(table$ID),
                         sort(unique(as.character(data$ID)))))

  # The same study with its identifiers written as text is the same table.
  text <- data
  text$ID <- paste0("S", formatC(as.integer(data$ID), width = 3, flag = "0"))
  text_table <- .model_estimation_data(text, roles, "cp")
  expect_identical(text_table[, setdiff(names(text_table), "ID")],
                   table[, setdiff(names(table), "ID")])
})

# The wait, broken down. Every candidate is a separate compiled population fit,
# so a caller deciding whether to name `pk` needs to see which one cost the time.
test_that("the fit reports how long each fit took", {
  skip_without_fitter()
  fit <- .default_fit()
  expect_true(all(c("fit", "total", "candidates", "pd") %in%
                    names(fit$timing)))
  expect_identical(fit$timing$candidates$model, fit$candidates$model)
  expect_true(all(fit$timing$candidates$seconds >= 0))
  out <- paste(utils::capture.output(print(model_report(fit))),
               collapse = " ")
  expect_match(out, "nlmixr2: ")
})

# The arm of every subject, named by subject. `vapply()` names its result from a
# character input and leaves a numeric one unnamed, so a study whose `ID` is a
# number gave an unnamed vector and anything asking "which subjects are in this
# arm" got nothing -- silently, because an empty arm looks like a small one.
test_that("the subject-to-arm vector is named whatever type the ID is", {
  data <- .oral_study()
  data$ARM <- ifelse(as.integer(data$ID) %% 2L == 0L, "high", "low")
  roles <- .estimate_roles(strata = "ARM")

  numeric_ids <- .model_subject_arms(data, roles)
  expect_identical(names(numeric_ids),
                   as.character(.unique_in_order(data$ID)))

  character_data <- data
  character_data$ID <- paste0("S", data$ID)
  character_ids <- .model_subject_arms(character_data, roles)
  expect_identical(names(character_ids),
                   .unique_in_order(character_data$ID))
  expect_identical(unname(numeric_ids), unname(character_ids))
})

# SIM-060 option (a), as an option rather than the default. The pooled shape
# predicts one number for every arm, so a synthetic patient's dose never reaches
# their response; per arm it does, and asserts no dose-response form to get
# there.
test_that("a PD shape can be fitted per arm, and falls back where it cannot", {
  data <- .oral_study()
  # Two arms whose PD endpoint differs by construction: one flat, one rising.
  data$ARM <- ifelse(as.integer(data$ID) %% 2L == 0L, "high", "low")
  pd_rows <- data[data$EVID == 0L, , drop = FALSE]
  pd_rows$DVID <- "effect"
  pd_rows$DV <- ifelse(pd_rows$ARM == "high", 20 + 0.05 * pd_rows$TIME, 20)
  data <- rbind(data, pd_rows)
  data <- data[order(data$ID, data$TIME, data$EVID == 0L), , drop = FALSE]
  roles <- .estimate_roles(strata = "ARM")

  observations <- .model_observations(data, roles)
  group <- .model_subject_arms(data, roles)
  pooled <- .model_fit_pd(observations, "effect", NULL)
  arms <- .model_fit_pd_arms(observations, "effect", NULL, group, pooled)

  expect_named(arms, .unique_in_order(unname(group)))
  at <- max(observations$aligned)
  value <- function(shape) {
    .pd_profile(list(pd = shape$pd), at, numeric(), numeric(), shape$typical)
  }
  # The pooled shape answers one number to both arms; the per-arm shapes do not.
  expect_gt(value(arms[["high"]]), value(arms[["low"]]) + 1)
  expect_true(value(pooled) > value(arms[["low"]]) &&
                value(pooled) < value(arms[["high"]]))

  # An arm with too little of the endpoint to fit keeps the pooled shape rather
  # than losing the endpoint, so no arm is left with nothing to generate.
  thin <- observations[observations$endpoint != "effect" |
                         observations$subject %in%
                         names(group)[group == "low"][1L], , drop = FALSE]
  fallback <- .model_fit_pd_arms(thin, "effect", NULL, group, pooled)
  expect_identical(fallback[["high"]], pooled)
})

test_that("pd_by_arm is TRUE or FALSE", {
  expect_error(synpmx_model_estimate(.oral_study(), .estimate_roles(),
                                     pd_by_arm = "yes"),
               "`pd_by_arm` must be TRUE or FALSE")
})

# SIM-068. A censored row goes to the fitter as censored rather than as a value
# nobody measured. The convention is `nlmixr2`'s and `pmx_roles()`'s at once:
# `CENS` 1 with `DV` holding the limit, and `LIMIT` closing the interval.
test_that("censoring reaches the fit table instead of being imputed away", {
  data <- .oral_study()
  roles <- .estimate_roles(cens = "CENS")
  endpoint <- as.character(data$DVID[which(data$EVID == 0L)[1L]])
  observed <- which(data$EVID == 0L & !is.na(data$DV) & data$DVID == endpoint)
  limit <- stats::quantile(data$DV[observed], 0.2)
  data$CENS <- 0L
  data$CENS[observed[data$DV[observed] < limit]] <- 1L
  data$DV[data$CENS == 1L] <- limit

  table <- .model_estimation_data(data, roles, endpoint)
  expect_true("CENS" %in% names(table))
  expect_identical(sum(table$CENS != 0), sum(data$CENS != 0L))
  # Never on a dose record, which is not an observation of anything.
  expect_true(all(table$CENS[table$EVID != 0L] == 0))
  # And the value handed over is the study's own limit, not a draw below it.
  expect_true(all(table$DV[table$CENS != 0] == limit))

  # No `cens` role, no column: a study that declares none is handed the table
  # it was handed before.
  plain <- .model_estimation_data(data, .estimate_roles(), endpoint)
  expect_false("CENS" %in% names(plain))
})

test_that("a declared interval limit rides along with the flag", {
  data <- .oral_study()
  roles <- .estimate_roles(cens = "CENS", limit = "LIMIT")
  endpoint <- as.character(data$DVID[which(data$EVID == 0L)[1L]])
  observed <- which(data$EVID == 0L & !is.na(data$DV) & data$DVID == endpoint)
  data$CENS <- 0L
  data$LIMIT <- NA_real_
  data$CENS[observed[1:3]] <- 1L
  data$LIMIT[observed[1:3]] <- 0.01

  table <- .model_estimation_data(data, roles, endpoint)
  expect_identical(sum(!is.na(table$LIMIT)), 3L)
  expect_true(all(table$LIMIT[table$CENS != 0] == 0.01))
  # `LIMIT` is meaningless where nothing is censored, and stays empty there.
  expect_true(all(is.na(table$LIMIT[table$CENS == 0])))
})

# The floor sits beside the declared assay limit rather than under it: where
# the study says where its assay stopped, `.censor_latent()` puts that boundary
# back and a second floor beneath it would be counted and warned about while
# changing nothing.
test_that("an endpoint the study censors is given no floor", {
  data <- .oral_study()
  roles <- .estimate_roles(cens = "CENS")
  observed <- which(data$EVID == 0L & !is.na(data$DV))
  endpoint <- as.character(data$DVID[observed[1L]])
  data$CENS <- 0L
  data$CENS[observed[data$DV[observed] < stats::quantile(data$DV[observed],
                                                         0.05)]] <- 1L
  expect_length(.model_quantification_floor(data, roles, endpoint), 0L)
  # An endpoint the same study leaves uncensored still gets one.
  data$CENS <- 0L
  expect_length(.model_quantification_floor(data, roles, endpoint), 1L)
})

# SIM-064. The same reading decides the residual model: `nimoData` reports one
# negative concentration in 321, and treating that as evidence that the
# endpoint reaches zero fitted an additive residual of 1.46 to values whose
# median is 3.
test_that("one non-positive reading does not make the scale include zero", {
  nimo_shaped <- c(-0.2319, seq(0.2647, 10, length.out = 320))
  expect_equal(.model_assay_floor(nimo_shaped), 0.2647 / 2)
  # An endpoint whose scale really includes zero says so in many rows.
  expect_null(.model_assay_floor(c(rep(0, 20), seq(0.3, 10, length.out = 100))))
})

test_that("an endpoint that routinely reports zero is given no floor", {
  data <- .oral_study()
  roles <- .estimate_roles()
  observed <- which(data$EVID == 0L & !is.na(data$DV))
  endpoint <- as.character(data$DVID[observed[1L]])
  data$DV[observed[seq_len(ceiling(0.2 * length(observed)))]] <- 0
  expect_length(.model_quantification_floor(data, roles, endpoint), 0L)
})

# SIM-062. The PD residual is what is left around each subject's own curve, so
# a study where subjects share a shape at different levels does not come back
# as scatter.
test_that("the PD residual excludes between-subject level", {
  set.seed(4)
  times <- c(0, 12, 24, 48, 72)
  observations <- do.call(rbind, lapply(seq_len(24), function(i) {
    level <- 100 * exp(stats::rnorm(1, 0, 0.30))   # a wide spread of levels
    data.frame(subject = as.character(i), endpoint = "resp", aligned = times,
               dv = level * exp(-0.02 * times) + stats::rnorm(length(times), 0, 1))
  }))
  fit <- .model_fit_pd(observations, "resp")
  # The level spread belongs to the baseline term, and the residual is the
  # measurement noise the data was built with rather than the whole spread.
  expect_gt(fit$baseline_cv, 0.2)
  expect_lt(fit$residual$sd, 0.25 * stats::sd(observations$dv))
})

test_that("a PD shape that fails to converge stays in the candidate table", {
  set.seed(5)
  times <- c(0, 12, 24, 48, 72)
  observations <- do.call(rbind, lapply(seq_len(12), function(i) {
    data.frame(subject = as.character(i), endpoint = "resp", aligned = times,
               dv = 50 + stats::rnorm(length(times), 0, 0.5))
  }))
  fit <- .model_fit_pd(observations, "resp")
  expect_true(all(c("shape", "converged", "aic", "note") %in%
                    names(fit$candidates)))
  expect_setequal(fit$candidates$shape,
                  c("constant", "linear", "exponential"))
})

# SIM-077 / REV-050. Bioavailability has to reach the extravascular doses and
# no others. Under `linCmt()` it reached the intravenous ones instead, so `f`
# came back as its own reciprocal while the fit's own predictions tracked the
# data -- a wrong answer that looked like a right one.
#
# The true value is deliberately far from 0.7, which is where `f` starts. The
# claim this replaces used a true value equal to the starting value and would
# have passed against a parameter nothing ever touched.

test_that("the mixed model is written with states, not linCmt", {
  spec <- .model_nlmixr_function("1cmt_mixed", c(cl = 4, v = 40, ka = 0.4,
                                                 f = 0.7), "prop", 0.2)
  text <- paste(deparse(spec), collapse = "\n")
  expect_match(text, "d/dt(depot)", fixed = TRUE)
  expect_match(text, "d/dt(central)", fixed = TRUE)
  expect_match(text, "f(depot) <- f", fixed = TRUE)
  expect_false(grepl("linCmt", text, fixed = TRUE))

  # A single-route study has no `f` to place and keeps the closed form.
  oral <- .model_nlmixr_function("1cmt_oral", c(cl = 4, v = 40, ka = 0.4),
                                 "prop", 0.2)
  expect_match(paste(deparse(oral), collapse = "\n"), "linCmt", fixed = TRUE)
})

test_that("the mixed model's compartments carry the dose the data sends them", {
  # Evaluated rather than fitted: at known parameters the two routes have to
  # come out where `.pk_profile()` puts them, which is the thing `linCmt()` got
  # wrong without ever looking wrong.
  skip_if_not_installed("nlmixr2est")
  p <- c(cl = 4, v = 40, ka = 0.4, f = 0.35)
  at <- c(0.5, 2, 8, 24)
  expected_ev <- .pk_profile(list(pk = "1cmt_mixed"), at, 100, 0, p,
                             routes = "extravascular")
  expected_iv <- .pk_profile(list(pk = "1cmt_mixed"), at, 100, 0, p,
                             routes = "iv")
  # Bioavailability separates them: without it the two would differ only in
  # shape, and the inverted binding would pass this check.
  expect_lt(expected_ev[[1L]], expected_iv[[1L]])
  expect_equal(unname(expected_ev[[1L]] /
                        .pk_profile(list(pk = "1cmt_mixed"), at, 100, 0,
                                    c(p[c("cl", "v", "ka")], f = 1),
                                    routes = "extravascular")[[1L]]),
               0.35, tolerance = 1e-8)
})

# `max_fit_subjects`: the population model is fitted to a subset drawn in
# proportion to the arms, and everything else reads the whole study.

test_that("the cap cannot sit below the floor", {
  expect_error(synpmx_model_estimate(.oral_study(), .estimate_roles(),
                                     max_fit_subjects = 10L, min_subjects = 20L),
               "below `min_subjects`")
})

test_that("a capped fit reports the count and still generates", {
  skip_if_not_installed("nlmixr2est")
  data <- .oral_study(n = 30)
  fit <- suppressWarnings(suppressMessages(synpmx_model_estimate(
    data, .estimate_roles(), max_fit_subjects = 12L, min_subjects = 10L,
    seed = 2, quiet = TRUE)))
  expect_equal(fit$fit_subjects$fitted, 12L)
  expect_equal(fit$fit_subjects$of, 30L)
  expect_equal(fit$settings$max_fit_subjects, 12L)
  expect_match(paste(capture.output(print(model_report(fit))), collapse = "\n"),
               "12 of 30 patients")
  # The apparatus read the whole study.
  expect_equal(fit$n_source, 30L)
  synthetic <- suppressWarnings(synpmx_model_generate(fit, n_subjects = 10L, seed = 1))
  expect_true(validate_pmx(synthetic, .estimate_roles())$valid)
})
