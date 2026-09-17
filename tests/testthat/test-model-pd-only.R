# SIM-088: declaring a response must bypass PK inference and compartment fitting.
.pd_only_fixture <- function(doses = TRUE) {
  observations <- expand.grid(ID = 1:24, TIME = c(0, 1, 2, 4, 8))
  observations$NTIME <- observations$TIME
  observations$DV <- 10 + observations$ID / 24 + 0.5 * observations$TIME
  observations$EVID <- 0L
  observations$AMT <- 0
  observations$CMT <- 3L
  observations$ARM <- ifelse(observations$ID <= 12, "A", "B")
  if (doses) {
    dose <- observations[observations$TIME == 0, ]
    dose$DV <- NA_real_; dose$EVID <- 1L; dose$AMT <- 100; dose$CMT <- 1L
    observations <- rbind(dose, observations)
  }
  list(data = observations, roles = pmx_roles(id = "ID", time = "TIME",
    nominal_time = "NTIME", dv = "DV", evid = "EVID", amt = "AMT",
    cmt = "CMT", strata = "ARM"))
}

test_that("a declared PD-only study never invokes the population fitter", {
  local_mocked_bindings(.model_run_fit = function(...) stop("PK fitter called"),
                        .package = "synpmx")
  x <- .pd_only_fixture()
  fit <- synpmx_model_estimate(x$data, x$roles, endpoint_roles = c(pd = "DV"),
                              seed = 12, quiet = TRUE)
  expect_identical(fit$endpoints$pk, character())
  expect_identical(fit$endpoints$pd, "DV")
  expect_null(fit$structural)
  expect_null(fit$parameters)
  expect_length(fit$pk_models, 0)
  expect_equal(nrow(model_candidates(fit)), 0)
  expect_identical(model_parameters(fit)$pd, fit$pd)
  report <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(report, "PD-only")
  expect_false(grepl("Estimated by nlmixr2", report, fixed = TRUE))
  synthetic <- synpmx_model_generate(fit, seed = 12)
  expect_true(all(is.finite(synthetic$DV[synthetic$EVID == 0])))
  expect_true(all(synthetic$CMT[synthetic$EVID == 0] == 3))
  expect_true(any(synthetic$EVID != 0))
  expect_equal(synthetic, synpmx_model_generate(fit, seed = 12))
  direct <- synpmx_model(x$data, x$roles, endpoint_roles = c(pd = "DV"),
                         seed = 12, quiet = TRUE)
  expect_length(attr(direct, "pmx_fitted_model")$pk_models, 0)
  expect_true(all(is.finite(direct$DV[direct$EVID == 0])))
})

test_that("PD-only fitting and generation also work without dose records", {
  x <- .pd_only_fixture(FALSE)
  fit <- synpmx_model_estimate(x$data, x$roles,
    endpoint_roles = list(pk = character()), pd_by_arm = TRUE, quiet = TRUE)
  expect_setequal(names(fit$pd$DV$arms), c("A", "B"))
  synthetic <- synpmx_model_generate(fit, seed = 7)
  expect_gt(nrow(synthetic), 0)
  expect_true(all(synthetic$EVID == 0))
  expect_true(all(is.finite(synthetic$DV)))
})

test_that("response declarations validate names and conflicting PK requests", {
  x <- .pd_only_fixture()
  expect_error(synpmx_model_estimate(x$data, x$roles,
    endpoint_roles = c(pd = "missing"), quiet = TRUE), "PD endpoint.*not in")
  expect_error(synpmx_model_estimate(x$data, x$roles,
    endpoint_roles = c(pk = "DV", pd = "DV"), quiet = TRUE), "both PK and PD")
  expect_error(synpmx_model_estimate(x$data, x$roles,
    endpoint_roles = c(pd = "DV"), pk = "1cmt_oral", quiet = TRUE), "PD-only")
  expect_error(synpmx_model_estimate(x$data, x$roles,
    endpoint_roles = c(pd = "DV"), start_param = c(cl = 1), quiet = TRUE), "PD-only")
})

test_that("discrete PD endpoints use their visit frequencies without a PK fit", {
  x <- .pd_only_fixture(FALSE)
  x$data$DV <- as.integer(x$data$ID %% 2)
  fit <- synpmx_model_estimate(x$data, x$roles, endpoint_roles = c(pd = "DV"),
                              quiet = TRUE)
  expect_length(fit$pk_models, 0)
  expect_identical(fit$endpoints$discrete, "DV")
  synthetic <- synpmx_model_generate(fit, seed = 4)
  expect_true(all(synthetic$DV %in% c(0, 1)))
})

test_that("the public white-cell response is inferred as PD", {
  skip_if_not_installed("nlmixr2data")
  x <- as.data.frame(nlmixr2data::wbcSim); x$NTIME <- x$TIME
  roles <- pmx_roles(id = "ID", time = "TIME", nominal_time = "NTIME", dv = "DV",
                    amt = "AMT", evid = "EVID", cmt = "CMT", rate = "RATE")
  local_mocked_bindings(.model_run_fit = function(...) stop("PK fitter called"),
                        .package = "synpmx")
  fit <- synpmx_model_estimate(x, roles, seed = 1, quiet = TRUE)
  expect_length(fit$pk_models, 0)
  expect_identical(names(fit$pd), "DV")
  synthetic <- synpmx_model_generate(fit, seed = 1)
  expect_true(all(synthetic$CMT[synthetic$EVID == 0] == 3))
  expect_true(all(is.finite(synthetic$DV[synthetic$EVID == 0])))
})

# SIM-089: first-dose baselines are observations even at the dose timestamp.
test_that("positive first-dose baselines infer PD without invoking PK", {
  x <- .pd_only_fixture()
  local_mocked_bindings(.model_run_fit = function(...) stop("PK fitter called"),
                        .package = "synpmx")
  fit <- synpmx_model_estimate(x$data, x$roles, quiet = TRUE)
  expect_length(fit$pk_models, 0)
  expect_identical(fit$endpoints$decided_by, "inferred")
  generated <- synpmx_model_generate(fit, seed = 1)
  expect_true(all(is.finite(generated$DV[generated$EVID == 0])))
  # Missing doses likewise cannot justify a compartment fit.
  x <- .pd_only_fixture(FALSE)
  expect_length(synpmx_model_estimate(x$data, x$roles, quiet = TRUE)$pk_models, 0)
})

test_that("baseline inference uses actual dose times and permits a PK override", {
  x <- .pd_only_fixture()
  x$data$NTIME[x$data$EVID == 1] <- -1
  obs <- .model_observations(x$data, x$roles)
  expect_false(.model_endpoint_signals(x$data, x$roles, obs)$post_dose)
  expect_identical(.model_classify_endpoints(x$data, x$roles, obs,
                    c(pk = "DV"))$pk, "DV")
  # A zero baseline is not evidence of PD, nor is a later pre-dose trough.
  x$data$DV[x$data$TIME == 0 & x$data$EVID == 0] <- 0
  obs <- .model_observations(x$data, x$roles)
  expect_true(.model_endpoint_signals(x$data, x$roles, obs)$post_dose)
})

test_that("baseline inference preserves public PK endpoints beside responses", {
  skip_if_not_installed("nlmixr2data")
  warfarin <- as.data.frame(nlmixr2data::warfarin)
  warfarin$ntime <- warfarin$time
  roles <- pmx_roles(id = "id", time = "time", nominal_time = "ntime", dv = "dv",
                    amt = "amt", evid = "evid", dvid = "dvid")
  classified <- .model_classify_endpoints(warfarin, roles,
                                           .model_observations(warfarin, roles))
  expect_identical(classified$pk, "cp")
  expect_identical(classified$pd, "pca")
  roles <- pmx_roles(id = "ID", time = "TIME", nominal_time = "NTIME", dv = "DV",
                    amt = "AMT", evid = "EVID", dvid = "NAME", cmt = "CMT",
                    addl = "ADDL", ii = "II", cens = "CENS")
  source <- pmx_expand_doses(onc_sim, roles)
  classified <- .model_classify_endpoints(source, roles,
                                           .model_observations(source, roles))
  expect_identical(classified$pk, "Everolimus trough")
  expect_identical(classified$pd, "SLD")
})
