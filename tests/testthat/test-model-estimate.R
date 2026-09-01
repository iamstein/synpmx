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

# Fitting is slow -- every candidate compiles -- so the fits the tests share are
# built once and reused. A test that needs its own fit says so by building one.
.shared <- local({
  cache <- list()
  function(name, expression) {
    if (is.null(cache[[name]])) cache[[name]] <<- expression
    cache[[name]]
  }
})

.searched_fit <- function() {
  .shared("searched", synpmx_model_estimate(.oral_study(), .estimate_roles(),
                                            quiet = TRUE,
                                            covariate_effects = "none"))
}

test_that("the search recovers the model and parameters it was simulated from", {
  skip_without_fitter()
  fit <- .searched_fit()
  expect_s3_class(fit, "pmx_fitted_model")
  expect_identical(fit$structural, "1cmt_oral")
  expect_equal(unname(fit$parameters$fixed[["cl"]]), 4, tolerance = 0.3)
  expect_equal(unname(fit$parameters$fixed[["v"]]), 40, tolerance = 0.3)
  expect_identical(fit$endpoints$pk, "cp")
})

test_that("every candidate the design admitted is in the table, converged or not", {
  skip_without_fitter()
  fit <- .searched_fit()
  table <- model_candidates(fit)
  expect_true(all(c("model", "converged", "aic", "note") %in% names(table)))
  expect_true(nrow(table) >= 1L)
  # A search that came down to one survivor must not look like a search that
  # had one candidate, so failures keep their row and their reason.
  expect_true(all(is.na(table$aic) | is.finite(table$aic)))
  expect_true(all(nzchar(table$note[!table$converged])))
  expect_identical(fit$structural, table$model[which.min(table$aic)])
})

test_that("`pk` forces a model and skips the search", {
  skip_without_fitter()
  data <- .oral_study()
  fit <- synpmx_model_estimate(data, .estimate_roles(), pk = "1cmt_oral",
                               quiet = TRUE, covariate_effects = "none")
  expect_identical(fit$structural, "1cmt_oral")
  expect_identical(nrow(model_candidates(fit)), 1L)
  expect_match(fit$design$reason, "declared")
})

test_that("`pk` naming a model outside the closed-form set is refused", {
  data <- .oral_study()
  expect_error(
    synpmx_model_estimate(data, .estimate_roles(), pk = "3cmt_oral",
                          quiet = TRUE),
    "should be one of"
  )
})

test_that("allometric scaling is kept only where it improves AIC", {
  skip_without_fitter()
  auto <- .shared("auto", synpmx_model_estimate(
    .oral_study(), .estimate_roles(covariates = "WT"), pk = "1cmt_oral",
    quiet = TRUE, covariate_effects = "auto"))
  none <- .searched_fit()
  expect_length(none$covariate_effects, 0L)
  # Weight is simulated independent of clearance here, so allometry should not
  # earn its two parameters. Either answer is legitimate; what must hold is that
  # the effect is recorded when it is kept and absent when it is not.
  if (length(auto$covariate_effects)) {
    expect_identical(auto$covariate_effects$cl$covariate, "WT")
    expect_equal(auto$covariate_effects$cl$exponent, 0.75)
  } else {
    expect_identical(auto$structural, none$structural)
  }
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
    .oral_study(), .estimate_roles(covariates = "WT"), pk = "1cmt_oral",
    quiet = TRUE, covariate_effects = "auto"))
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
  expect_lt(fit$parameters$fixed[["cl"]], .searched_fit()$parameters$fixed[["cl"]])
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
  fit <- .searched_fit()
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
