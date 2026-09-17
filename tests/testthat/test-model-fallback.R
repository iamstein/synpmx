# SIM-087: acceptance is independent of a finite AIC, and fallback must not
# silently replace a model explicitly requested by the caller.
.fallback_search <- function(failed = character(), candidates = c("1cmt_oral", "2cmt_oral"),
                              fallback = TRUE, errors = character()) {
  tried <- character()
  local_mocked_bindings(
    .model_initial_estimates = function(...) c(cl = 4, v = 40, ka = 1),
    .model_nlmixr_function = function(candidate, ...) candidate,
    .model_run_fit = function(spec, ...) {
      tried <<- c(tried, spec)
      if (spec %in% errors) stop("solver failed")
      list(model = spec)
    },
    .model_fit_aic = function(fit) if (grepl("^1cmt", fit$model)) 1 else 10,
    .model_assess_fit = function(fit, ...) list(
      converged = TRUE, accepted = !fit$model %in% failed,
      note = if (fit$model %in% failed) "population simulation failed" else ""),
    .package = "synpmx"
  )
  answer <- .model_fit_candidates(data.frame(DV = 1:3), candidates, NULL,
    "cp", "prop", "focei", TRUE, fallback = fallback)
  list(answer = answer, tried = tried)
}

test_that("an accepted two-compartment fit skips its sibling", {
  x <- .fallback_search()
  expect_identical(x$tried, "2cmt_oral")
  expect_identical(x$answer$selected, "2cmt_oral")
  expect_match(x$answer$reason, "passed acceptance")
})

test_that("generation failure and solver errors both trigger fallback", {
  x <- .fallback_search(failed = "2cmt_oral")
  expect_identical(x$tried, c("2cmt_oral", "1cmt_oral"))
  expect_identical(x$answer$selected, "1cmt_oral")
  expect_match(x$answer$reason, "population simulation failed")
  expect_true(x$answer$table$converged[[1]])
  expect_false(x$answer$table$accepted[[1]])
  y <- .fallback_search(errors = "2cmt_oral")
  expect_identical(y$answer$selected, "1cmt_oral")
  expect_false(y$answer$table$converged[[1]])
  expect_match(y$answer$reason, "solver failed")
  expect_error(.fallback_search(failed = c("1cmt_oral", "2cmt_oral")),
                "No candidate model passed")
})

test_that("explicit models are never replaced and explicit sets compare AIC", {
  expect_error(.fallback_search(failed = "2cmt_oral", candidates = "2cmt_oral",
                                 fallback = FALSE), "No candidate model passed")
  x <- .fallback_search(fallback = FALSE)
  expect_identical(x$tried, c("1cmt_oral", "2cmt_oral"))
  expect_identical(x$answer$selected, "1cmt_oral")
})

test_that("ambiguous routes fall back separately before comparison", {
  x <- .fallback_search(failed = "2cmt_iv",
    candidates = c("1cmt_iv", "1cmt_oral", "2cmt_iv", "2cmt_oral"))
  expect_identical(x$tried, c("2cmt_iv", "2cmt_oral", "1cmt_iv"))
  expect_identical(x$answer$selected, "1cmt_iv")
  expect_match(x$answer$reason, "lowest AIC among accepted routes")
})

test_that("missing AIC cannot silently decide an explicit comparison", {
  table <- data.frame(model = c("1cmt_oral", "2cmt_oral"),
                       accepted = TRUE, aic = c(1, NA_real_), note = "")
  expect_error(.model_select_candidate(table), "AIC is unavailable")
  expect_identical(.model_select_candidate(table[2, ])$selected, "2cmt_oral")
})

.check_fixture <- function() {
  p <- list(fixed = c(cl = 4, v = 40, ka = 1),
    omega = diag(c(0.04, 0.04, 0.04)),
    residual = list(kind = "proportional", cv = 0.1))
  dimnames(p$omega) <- list(names(p$fixed), names(p$fixed))
  times <- c(0.25, 0.5, 1, 2, 4, 8, 12, 24)
  one <- data.frame(ID = 1, TIME = c(0, times), EVID = c(1L, rep(0L, 8)),
    AMT = c(100, rep(0, 8)), DV = c(NA, .pk_profile(list(pk = "1cmt_oral"),
                                                 times, 100, 0, p$fixed)))
  data <- do.call(rbind, lapply(1:20, function(i) { one$ID <- i; one }))
  list(parameters = p, data = data)
}

test_that("optimizer failure is detected even when a fit has finite AIC", {
  expect_null(.model_check_convergence(list(convergence = 0), "focei"))
  for (code in list(1, 8, NA_real_, NULL)) {
    expect_match(.model_check_convergence(list(convergence = code, AIC = 1,
                   message = "false convergence"), "focei"), "successful convergence")
  }
  expect_match(.model_check_convergence(list(convergence = 0,
    message = "false convergence (8)"), "focei"), "successful convergence")
})

test_that("SAEM acceptance reads the trajectory rather than optimizer codes", {
  stable <- data.frame(iter = 1:100, tcl = rep(1, 100), tv = rep(3, 100))
  expect_null(.model_check_convergence(list(parHistData = stable), "saem"))
  stable$tcl[91:100] <- 10
  expect_match(.model_check_convergence(list(parHistData = stable), "saem"), "drifting")
  expect_match(.model_check_convergence(list(), "saem"), "trajectory")
})

test_that("invalid fixed effects, residuals and covariance cannot reach generation", {
  p <- .check_fixture()$parameters
  expect_null(.model_check_parameters(p))
  bad <- p; bad$fixed[1] <- Inf
  expect_match(.model_check_parameters(bad), "fixed effects")
  bad <- p; bad$residual$cv <- -1
  expect_match(.model_check_parameters(bad), "residual")
  bad <- p; bad$omega[1, 1] <- -1
  expect_match(.model_check_parameters(bad), "covariance")
  bad <- p; bad$omega[1, 2] <- 2
  expect_match(.model_check_parameters(bad), "symmetric")
})

test_that("generation checks use new random effects and preserve the RNG state", {
  x <- .check_fixture()
  set.seed(7); before <- .Random.seed
  check <- .model_check_generation(x$parameters, "1cmt_oral", x$data)
  expect_null(check$failure)
  expect_identical(.Random.seed, before)
  expect_identical(check, .model_check_generation(x$parameters, "1cmt_oral", x$data))
  # A plausible typical curve cannot hide an explosive population covariance.
  x$parameters$omega <- x$parameters$omega * 1000000
  expect_match(.model_check_generation(x$parameters, "1cmt_oral", x$data)$failure,
                "invalid parameters|non-finite|100")
  expect_identical(.Random.seed, before)
})

test_that("moderate discrepancies warn while gross level failures reject", {
  x <- .check_fixture()
  x$data$DV <- x$data$DV / 3
  check <- .model_check_generation(x$parameters, "1cmt_oral", x$data)
  expect_null(check$failure)
  expect_match(check$note, "median ratio")
  x$data$DV <- x$data$DV / 1000
  expect_match(.model_check_generation(x$parameters, "1cmt_oral", x$data)$failure,
                "100")
  x$data$CENS <- ifelse(x$data$EVID == 0, 1, 0)
  check <- .model_check_generation(x$parameters, "1cmt_oral", x$data)
  expect_match(check$failure, "100")
  expect_match(check$note, "omitted for censored")
  x$data$CENS <- -x$data$CENS
  check <- .model_check_generation(x$parameters, "1cmt_oral", x$data)
  expect_null(check$failure)
  expect_match(check$note, "source maximum unknown")
})

test_that("mixed doses and infusions use the same profiles as generation", {
  x <- .check_fixture()
  x$parameters$fixed <- c(x$parameters$fixed, f = 0.5)
  x$data$RATE <- ifelse(x$data$EVID == 1, 50, 0)
  x$data$CMT <- 2L
  check <- .model_check_generation(x$parameters, "1cmt_mixed", x$data)
  expect_null(check$failure)
})

test_that("missing observed designs cannot silently pass generation checks", {
  x <- .check_fixture()
  expect_match(.model_check_generation(x$parameters, "1cmt_oral",
    x$data[FALSE, ])$failure, "no observed design")
})

test_that("the generation screen applies the declared allometry", {
  x <- .check_fixture()
  x$parameters$omega[] <- 0
  x$parameters$residual$cv <- 0
  x$data$WT <- 140
  p <- x$parameters$fixed
  p[c("cl", "v")] <- p[c("cl", "v")] * 2^c(0.75, 1)
  observed <- x$data$EVID == 0
  x$data$DV[observed] <- .pk_profile(list(pk = "1cmt_oral"),
    x$data$TIME[observed], 100, 0, p)
  check <- .model_check_generation(x$parameters, "1cmt_oral", x$data,
                                    list(covariate = "WT", reference = 70))
  expect_null(check$failure)
  expect_equal(check$median_ratio, 1)
  expect_equal(check$spread_ratio, 1)
})

test_that("overflowing diagnostic ratios still count as gross discrepancies", {
  x <- .check_fixture()
  x$parameters$residual$cv <- 0
  x$data$DV[x$data$EVID == 0] <- 1e-200
  x$data$DV[x$data$TIME == 24] <- 1e308
  local_mocked_bindings(.pk_profile = function(model, time, ...) rep(1e200, length(time)),
                        .package = "synpmx")
  expect_match(.model_check_generation(x$parameters, "1cmt_oral", x$data)$failure,
                "generated median")
})


test_that("left censoring cannot hide explosive residual population draws", {
  # A SAEM trial on the public case1 study returned an enormous residual with
  # a flat trajectory. Both convergence and finite-parameter checks passed;
  # the source's left censoring must not disable the upper-tail rejection.
  x <- .check_fixture()
  x$data$CENS <- ifelse(x$data$TIME == 24, 1, 0)
  x$parameters$residual$cv <- 1e50
  expect_match(.model_check_generation(x$parameters, "1cmt_oral", x$data)$failure,
                "100 times")
})


test_that("partly left-censored data still reject collapsed concentrations", {
  x <- .check_fixture()
  x$data$CENS <- ifelse(x$data$TIME == 24, 1, 0)
  x$parameters$fixed[c("cl", "v")] <- x$parameters$fixed[c("cl", "v")] * 1e40
  x$parameters$residual$cv <- 1e20
  check <- .model_check_generation(x$parameters, "1cmt_oral", x$data)
  expect_match(check$failure, "source median bounds")
  expect_gt(check$source_median_bounds[["lower"]], 0)
})

test_that("left, right and interval limits give conservative median bounds", {
  x <- data.frame(DV = c(2, 4, 6, 8, 3), CENS = c(1, -1, 1, -1, 0),
                   LIMIT = c(NA, NA, 1, 10, 400))
  b <- .model_observation_bounds(x)
  expect_equal(b$lower, c(0, 4, 1, 8, 3))
  expect_equal(b$upper, c(2, Inf, 6, 10, 3))
  expect_identical(b$observed, c(FALSE, FALSE, FALSE, FALSE, TRUE))
})

# SIM-090: optimizer diagnostics and generation acceptance are separate.
test_that("false convergence alone warns and keeps a usable fit", {
  x <- .check_fixture()
  local_mocked_bindings(.model_read_fit = function(...) x$parameters,
                        .package = "synpmx")
  fit <- list(convergence = 1, message = "false convergence (8)")
  assess <- function(f) .model_assess_fit(f, "1cmt_oral", x$parameters$fixed,
                                         "prop", "focei", x$data)
  result <- assess(fit)
  expect_false(result$converged)
  expect_true(result$accepted)
  expect_match(result$note, "Warning: the optimizer stalled before convergence was confirmed. Check that the synthetic data reasonably reproduces the source data's patterns and variability.", fixed = TRUE)
  fit$message <- "iteration limit reached"
  expect_false(assess(fit)$accepted)
  fit$message <- "false convergence (8); evaluation limit reached"
  expect_false(assess(fit)$accepted)
  fit$message <- "false convergence (8)"
  x$parameters$fixed[1] <- Inf
  expect_false(assess(fit)$accepted)
})
