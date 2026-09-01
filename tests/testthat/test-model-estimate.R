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
