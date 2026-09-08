# The PMX model generator's fitted-model object and its input gates. One
# deterministic test per gate, on fixtures, base R only -- no suggested package,
# and nothing here fits anything, because the gates and the object are built
# before any estimation exists.

.model_roles <- function(...) {
  pmx_roles(id = "ID", time = "TIME", nominal_time = "NTIME", dv = "DV",
            amt = "AMT", evid = "EVID", cmt = "CMT", dvid = "DVID",
            mdv = "MDV", ...)
}

test_that("a cohort under min_subjects warns, fits, and names what it costs", {
  # The message is wrapped for reading, so a phrase can span a line break.
  # Match the squished text rather than the layout.
  note <- squish(tryCatch(.model_note_subjects(19L, 20L),
                          warning = conditionMessage))
  expect_match(note, "fitting a population model to 19 subjects")
  expect_match(note, "the scorecard asks")
  expect_match(note, "synpmx_pca\\(\\)")
  expect_warning(.model_note_subjects(19L, 20L))
  expect_silent(.model_note_subjects(20L, 20L))
  # Higher than PCA's floor of 10, deliberately: a cohort PCA will summarize is
  # not necessarily one this generator will fit well.
  expect_warning(.model_note_subjects(12L, 20L))
})

test_that("min_subjects itself must be one positive integer", {
  expect_error(.model_note_subjects(30L, 0L), "min_subjects")
  expect_error(.model_note_subjects(30L, c(5L, 10L)), "min_subjects")
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
  # Observations alone: dropping them is a real option and is offered.
  expect_error(.model_require_nominal_time(data, .model_roles()),
               "drop those rows")
})

# The advice has to depend on what the rows are. "Drop those rows" applied to
# the dose records takes the dosing out of the study, and the run then fails a
# gate later with a message about the sampling schedule -- which is how a
# missing nominal time on a dose record turns into a confusing error about
# observations.
test_that("dose records missing a nominal time are named, and not dropped", {
  data <- pmx_simulated_fixture(24)
  data$NTIME[data$EVID != 0] <- NA_real_
  expect_error(.model_require_nominal_time(data, .model_roles()),
               "of them dose records")
  expect_error(.model_require_nominal_time(data, .model_roles()),
               "takes the dosing out of the study")
  expect_false(grepl("drop those rows", tryCatch(
    .model_require_nominal_time(data, .model_roles()),
    error = function(e) conditionMessage(e))))
})

test_that("time coverage counts distinct nominal times after a dose", {
  data <- pmx_simulated_fixture(24)
  roles <- .model_roles()
  # The fixture doses at 0 and 12 and samples the same offsets in each
  # interval, so the count is the protocol's slots and not the row count.
  expect_identical(.model_time_coverage(data, roles), 7L)
  expect_silent(.model_require_time_coverage(data, roles, 6L))
  expect_warning(.model_require_time_coverage(data, roles, 8L),
                 "below `min_time_bins` = 8")
})

test_that("a trough-only study is warned about rather than refused", {
  data <- pmx_simulated_fixture(24)
  # Keep only the pre-dose samples: two nominal times after dose survive, which
  # is what a study sampled at troughs alone looks like.
  keep <- data$EVID != 0 | data$NTIME %in% c(0, 12)
  data <- data[keep, , drop = FALSE]
  roles <- .model_roles()
  expect_lt(.model_time_coverage(data, roles), 6L)
  expect_warning(.model_require_time_coverage(data, roles, 6L),
                 "starting values")
  expect_warning(.model_require_time_coverage(data, roles, 6L),
                 "synpmx_avatar\\(\\)")
})

# Nothing to fit is still an error, and the message says which of the two role
# columns the count of zero was really about.
test_that("no post-dose observation is refused, naming what the roles read", {
  data <- pmx_simulated_fixture(24)
  roles <- .model_roles()

  no_observations <- data[data$EVID != 0, , drop = FALSE]
  expect_identical(.model_time_coverage(no_observations, roles), 0L)
  expect_error(.model_require_time_coverage(no_observations, roles, 6L),
               "none of the .* rows is an observation")
  # Each branch opens with its own finding. "No observation recorded after a
  # dose" in front of a study with no dose records sends the reader to look at
  # sampling times that are fine.
  expect_error(.model_require_time_coverage(no_observations, roles, 6L),
               "found nothing to fit")

  no_doses <- data[data$EVID == 0, , drop = FALSE]
  expect_identical(.model_time_coverage(no_doses, roles), 0L)
  expect_error(.model_require_time_coverage(no_doses, roles, 6L),
               "found no dose records")
  # The message shows what the column held, so a dataset that marks its doses
  # some other way -- or lost them to a filter -- can be read off the error.
  expect_error(.model_require_time_coverage(no_doses, roles, 6L),
               "`EVID` holds 0 \\([0-9]+ rows\\)")
  expect_error(.model_require_time_coverage(no_doses, roles, 6L),
               "0 of [0-9]+ rows have a positive `AMT`")
  amt_only <- no_doses
  amt_only$AMT[seq(1, nrow(amt_only), by = 10)] <- 100
  expect_error(.model_require_time_coverage(amt_only, roles, 6L),
               "and [1-9][0-9]* of [0-9]+ rows have a positive `AMT`")

  # Doses and observations that never meet in one subject.
  split_ids <- data
  split_ids[[roles$id]] <- ifelse(data$EVID == 0, data$ID, data$ID + 1000L)
  expect_identical(.model_time_coverage(split_ids, roles), 0L)
  expect_error(.model_require_time_coverage(split_ids, roles, 6L),
               "no subject holds an observation after a dose")
  expect_error(.model_require_time_coverage(split_ids, roles, 6L),
               "found no observation recorded after a dose")
})

test_that("a short arm is warned about and dropped, naming the caller", {
  group <- c(rep("a", 10), rep("b", 2))
  expect_warning(keep <- .drop_short_arms(group, 3L,
                                          "synpmx_model_estimate()"),
                 "`synpmx_model_estimate\\(\\)` dropped 2 patient\\(s\\)")
  expect_identical(keep, group == "a")
  expect_warning(.drop_short_arms(group, 3L, "synpmx_pca()"),
                 "`synpmx_pca\\(\\)` dropped 2 patient\\(s\\) in 1 arm")
  expect_silent(keep <- .drop_short_arms(group, 2L,
                                         "synpmx_model_estimate()"))
  expect_true(all(keep))
})

test_that("synpmx_pca_summarize() drops a short arm rather than refusing", {
  data <- pmx_simulated_fixture(24)
  data$ARM <- ifelse(data$ID <= 22, "a", "b")
  roles <- .model_roles(strata = "ARM")
  expect_warning(trial_summary <- synpmx_pca_summarize(data, roles),
                 "`synpmx_pca\\(\\)` dropped 2 patient\\(s\\)")
  expect_identical(names(trial_summary$arms$sizes), "a")
  expect_identical(trial_summary$n_source, 22L)
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

# SIM-075. A fit that reports an objective and an AIC has not necessarily
# estimated anything. Where the optimizer takes no effective step, `nlmixr2`
# returns the starting values and every number downstream is a starting value
# wearing the costume of an estimate.

test_that("a fit that never moved is detected", {
  fixed <- c(cl = 6.172, v = 55.61, ka = 1.366, f = 0.7)
  omega <- diag(c(0.1, 0.1, 0.1))
  rownames(omega) <- colnames(omega) <- c("cl", "v", "ka")
  start <- c(cl = 6.172, v = 55.61, ka = 1.366, f = 0.7)

  movement <- .model_fit_movement(fixed, omega, start)
  expect_false(movement$moved)
  expect_true(all(movement$changes$relative_change == 0))
  # Every parameter is accounted for, the between-subject terms included: the
  # study that found this had three etas all sitting on sqrt(0.1).
  expect_setequal(movement$changes$parameter,
                  c("cl", "v", "ka", "f", "omega.cl", "omega.v", "omega.ka"))
})

test_that("a fit that moved anywhere is not flagged", {
  start <- c(cl = 4, v = 40, ka = 0.4, f = 0.7)
  omega <- diag(c(0.1, 0.1, 0.1))
  rownames(omega) <- colnames(omega) <- c("cl", "v", "ka")

  # One fixed effect moving is enough; the gate asks whether the optimizer
  # moved at all, not whether it converged tightly.
  moved_fixed <- .model_fit_movement(c(cl = 5, v = 40, ka = 0.4, f = 0.7),
                                     omega, start)
  expect_true(moved_fixed$moved)

  # So is one between-subject term moving while the fixed effects sit still.
  omega_moved <- omega
  diag(omega_moved) <- c(0.25, 0.1, 0.1)
  expect_true(.model_fit_movement(start, omega_moved, start)$moved)
})

test_that("the unmoved warning says so in words nobody skims past", {
  fixed <- c(cl = 4, v = 40)
  omega <- diag(c(0.1, 0.1))
  rownames(omega) <- colnames(omega) <- c("cl", "v")
  movement <- .model_fit_movement(fixed, omega, c(cl = 4, v = 40))

  note <- squish(tryCatch(.model_warn_unmoved("1cmt_iv", movement),
                          warning = conditionMessage))
  expect_match(note, "THE FIT DID NOT MOVE")
  expect_match(note, "starting values rather than estimates")
  expect_match(note, "Do not generate from this fit")
  # It names the parameters, because "did not move" is a claim the reader has
  # to be able to check.
  expect_match(note, "cl")
})

# `max_fit_subjects`: the population model is fitted to a subset drawn in
# proportion to the arms, and everything else reads the whole study.

test_that("the fit subset is proportional to the arms and seeded", {
  subjects <- as.character(1:100)
  arms <- rep(c("A", "B", "C"), c(60, 30, 10))
  drawn <- .model_fit_subset(subjects, subjects, arms, 20L, seed = 4)
  expect_length(drawn, 20L)
  expect_equal(as.vector(table(arms[match(drawn, subjects)])), c(12L, 6L, 2L))
  expect_identical(drawn, .model_fit_subset(subjects, subjects, arms, 20L, seed = 4))
  expect_false(identical(drawn, .model_fit_subset(subjects, subjects, arms, 20L, seed = 5)))
  # Order is the source's, so the empirical Bayes estimates line up.
  expect_identical(drawn, subjects[subjects %in% drawn])
  # An arm with one fitted patient keeps that patient.
  tiny <- c(rep("A", 98), "B", "C")
  small <- .model_fit_subset(subjects, subjects, tiny, 10L, seed = 1)
  expect_true(all(c("99", "100") %in% small))
  # Below the cap nothing is drawn.
  expect_null(.model_fit_subset(subjects[1:15], subjects, arms, 20L, seed = 1))
})
