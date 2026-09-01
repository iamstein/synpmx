# The PMX model generator's fitted-model object and its input gates. One
# deterministic test per gate, on fixtures, base R only -- no suggested package,
# and nothing here fits anything, because the gates and the object are built
# before any estimation exists.

.model_roles <- function(...) {
  pmx_roles(id = "ID", time = "TIME", nominal_time = "NTIME", dv = "DV",
            amt = "AMT", evid = "EVID", cmt = "CMT", dvid = "DVID",
            mdv = "MDV", ...)
}

test_that("the cohort floor refuses below min_subjects and names the alternatives", {
  expect_error(.model_require_subjects(19L, 20L), "at least 20 subjects")
  expect_error(.model_require_subjects(19L, 20L), "synpmx_pca\\(\\)")
  expect_true(.model_require_subjects(20L, 20L))
  # Higher than PCA's floor of 10, deliberately: a cohort PCA will summarize is
  # not necessarily one this generator will fit.
  expect_error(.model_require_subjects(12L, 20L))
})

test_that("min_subjects itself must be one positive integer", {
  expect_error(.model_require_subjects(30L, 0L), "min_subjects")
  expect_error(.model_require_subjects(30L, c(5L, 10L)), "min_subjects")
})

test_that("an undeclared nominal grid is refused, and says why", {
  data <- pmx_simulated_fixture(24)
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", cmt = "CMT", dvid = "DVID", mdv = "MDV")
  expect_error(.model_require_nominal_time(data, roles),
               "requires `nominal_time`")
  expect_error(.model_require_nominal_time(data, roles),
               "statement about the protocol")
  expect_true(.model_require_nominal_time(data, .model_roles()))
})

test_that("a nominal grid with holes in it is refused", {
  data <- pmx_simulated_fixture(24)
  observed <- which(data$EVID == 0 & !is.na(data$DV))
  data$NTIME[observed[1:3]] <- NA_real_
  expect_error(.model_require_nominal_time(data, .model_roles()),
               "missing on 3 of")
})

test_that("time coverage counts distinct nominal times after a dose", {
  data <- pmx_simulated_fixture(24)
  roles <- .model_roles()
  # The fixture doses at 0 and 12 and samples the same offsets in each
  # interval, so the count is the protocol's slots and not the row count.
  expect_identical(.model_time_coverage(data, roles), 7L)
  expect_true(.model_require_time_coverage(data, roles, 6L))
  expect_error(.model_require_time_coverage(data, roles, 8L),
               "8 distinct nominal times")
})

test_that("a trough-only study is refused and pointed at the other generators", {
  data <- pmx_simulated_fixture(24)
  # Keep only the pre-dose samples: two nominal times after dose survive, which
  # is what a study sampled at troughs alone looks like.
  keep <- data$EVID != 0 | data$NTIME %in% c(0, 12)
  data <- data[keep, , drop = FALSE]
  roles <- .model_roles()
  expect_lt(.model_time_coverage(data, roles), 6L)
  expect_error(.model_require_time_coverage(data, roles, 6L),
               "starting values")
  expect_error(.model_require_time_coverage(data, roles, 6L),
               "synpmx_avatar\\(\\)")
})

test_that("the arm floor is the shared one, and names the caller", {
  group <- c(rep("a", 10), rep("b", 2))
  expect_error(.require_arms(group, 3L, "synpmx_model_estimate()"),
               "`synpmx_model_estimate\\(\\)` needs at least 3 patients")
  expect_error(.require_arms(group, 3L, "synpmx_pca()"),
               "`synpmx_pca\\(\\)` needs at least 3 patients")
  expect_true(.require_arms(group, 2L, "synpmx_model_estimate()"))
})

test_that("synpmx_pca_summarize() still refuses a short arm after the lift", {
  data <- pmx_simulated_fixture(24)
  data$ARM <- ifelse(data$ID <= 22, "a", "b")
  roles <- .model_roles(strata = "ARM")
  expect_error(synpmx_pca_summarize(data, roles),
               "`synpmx_pca\\(\\)` needs at least 3 patients")
})

# The privacy gate. `REV-042`: a model estimated from the confidential data
# entering the differentially private path as a public structural input would
# consume no budget for information it took from that data.

.fitted_stub <- function() {
  structure(list(structural = "1cmt_iv"), class = "pmx_fitted_model")
}

test_that("a fitted model is not a structural model", {
  expect_false(inherits(.fitted_stub(), "pmx_structural_model"))
})

test_that("the private path refuses a fitted model and names the reason", {
  fit <- .fitted_stub()
  expect_error(.reject_fitted_model(fit, "model", "synpmx_calibrated()"),
               "no privacy budget")
  expect_error(pmx_prior(fit, source = "literature"), "no privacy budget")
  expect_error(
    synpmx_prior(fit, design = NULL),
    "no privacy budget"
  )
  # And the reason is what comes back, rather than a generic type error.
  expect_error(pmx_prior(fit, source = "literature"), "pmx_fitted_model")
})

test_that("the gate passes anything that is not a fitted model through", {
  expect_true(.reject_fitted_model(c(0.5, 2), "range", "pmx_prior()"))
  expect_s3_class(pmx_prior(c(0.5, 2), source = "literature"), "pmx_prior")
})

# The object.

.fitted_fixture <- function(structural = "1cmt_iv", ...) {
  args <- list(
    structural = structural,
    candidates = data.frame(
      model = c("1cmt_iv", "2cmt_iv"), converged = c(TRUE, FALSE),
      aic = c(120.4, NA_real_), note = c("", "boundary"),
      stringsAsFactors = FALSE
    ),
    parameters = list(
      fixed = c(cl = 3.1, v = 21.0),
      omega = matrix(c(0.09, 0, 0, 0.04), 2, 2,
                     dimnames = list(c("cl", "v"), c("cl", "v"))),
      residual = list(kind = "proportional", cv = 0.18)
    ),
    endpoints = list(pk = "cp", pd = c(pd = "linear")),
    arms = list(arms = "all", sizes = c(all = 24L)),
    dosing = list(), visits = list(), schema = list(),
    roles = .model_roles(),
    settings = list(min_subjects = 20L, min_arm_patients = 3L),
    n_source = 24L
  )
  do.call(.pmx_fitted_model, utils::modifyList(args, list(...)))
}

test_that("a well-formed fitted model is constructed and prints", {
  fit <- .fitted_fixture()
  expect_s3_class(fit, "pmx_fitted_model")
  out <- paste(utils::capture.output(print(fit)), collapse = " ")
  expect_match(out, "1cmt_iv")
  expect_match(out, "24 patients")
  # The out-of-scope statement prints with the object, because its contents
  # look exactly like a real population analysis and will be read as one.
  expect_match(out, "not estimates to report")
})

test_that("the constructor refuses a structural model it cannot simulate", {
  expect_error(.fitted_fixture(structural = "1cmt_michaelis"), "must be one of")
})

test_that("the constructor refuses a selection no candidate supports", {
  expect_error(.fitted_fixture(structural = "2cmt_iv"),
               "not among the converged candidates")
})

test_that("the constructor refuses fixed effects the model needs and lacks", {
  expect_error(
    .fitted_fixture(structural = "1cmt_oral",
                    candidates = data.frame(model = "1cmt_oral",
                                            converged = TRUE, aic = 1,
                                            note = "",
                                            stringsAsFactors = FALSE)),
    "missing: ka"
  )
})

test_that("every random effect needs a fixed effect of the same name", {
  omega <- matrix(0.09, 1, 1, dimnames = list("ka", "ka"))
  expect_error(
    .fitted_fixture(parameters = list(fixed = c(cl = 3.1, v = 21),
                                      omega = omega,
                                      residual = list(kind = "proportional",
                                                      cv = 0.18))),
    "needs a fixed effect"
  )
})

test_that("the object carries no per-subject quantity", {
  # Empirical Bayes estimates are per-subject quantities, so an object carrying
  # them would be a description of each real patient in the study. Generation
  # draws random effects from the covariance matrix instead, which is why the
  # object has nowhere to put them.
  fit <- .fitted_fixture()
  expect_false(any(grepl("eta|eb|individual|subject",
                         names(fit$parameters), ignore.case = TRUE)))
  expect_length(fit$parameters$fixed, 2L)
})
