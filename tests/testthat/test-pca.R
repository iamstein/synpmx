# The PCA generator. `pca-algorithm.Rmd` is the prose these pin.

# `pmx_generated_roles()` names an assigned-dose column the fixture does not
# carry, so the roles are spelled out against the columns that exist.
pca_fixture <- function(n = 60) {
  list(
    data = pmx_simulated_fixture(n),
    roles = pmx_roles(
      id = "ID", time = "TIME", nominal_time = "NTIME", tad = "TAD",
      occasion = "OCC", dv = "DV", amt = "AMT", rate = "RATE", evid = "EVID",
      cmt = "CMT", dvid = "DVID", mdv = "MDV", cens = "CENS",
      covariates = c("WT", "AGE", "SEX")
    )
  )
}

test_that("synpmx_pca generates a legal table at the requested size", {
  fixture <- pca_fixture(60)
  synthetic <- synpmx_pca(fixture$data, fixture$roles, seed = 1)
  expect_true(validate_pmx(synthetic, fixture$roles)$valid)
  expect_equal(
    length(unique(synthetic[[fixture$roles$id]])),
    length(unique(fixture$data[[fixture$roles$id]]))
  )
  expect_setequal(names(synthetic),
                  intersect(names(fixture$data), names(synthetic)))
  expect_true(all(names(synthetic) %in% names(fixture$data)))
})

test_that("the same seed generates the same table", {
  fixture <- pca_fixture(40)
  a <- synpmx_pca(fixture$data, fixture$roles, seed = 7)
  b <- synpmx_pca(fixture$data, fixture$roles, seed = 7)
  expect_equal(a, b)
})

test_that("no generated subject copies a real subject's values or times", {
  fixture <- pca_fixture(60)
  synthetic <- synpmx_pca(fixture$data, fixture$roles, seed = 3)
  card <- synpmx_scorecard(fixture$data, synthetic, fixture$roles)
  copies <- as.data.frame(card)
  for (check in c("B4a", "B4b")) {
    row <- copies[copies$check == check, , drop = FALSE]
    expect_identical(row$verdict, "pass", info = check)
  }
})

# The guard the owner asked for: a late grid cell held by one or two patients
# describes those patients and nobody else, so it is dropped rather than filled
# from the median for everyone.
test_that("a grid cell held by too few patients is dropped", {
  fixture <- pca_fixture(60)
  data <- fixture$data
  subjects <- unique(data$ID)
  late <- max(data$NTIME)
  # Only two subjects reach the last visit.
  drop <- data$NTIME == late & data$EVID == 0 & !(data$ID %in% subjects[1:2])
  sparse <- data[!drop, , drop = FALSE]

  strict <- attr(synpmx_pca(sparse, fixture$roles, seed = 5,
                            min_column_patients = 3L), "pmx_pca_fit")
  loose <- attr(synpmx_pca(sparse, fixture$roles, seed = 5,
                           min_column_patients = 2L), "pmx_pca_fit")
  held <- function(fit) {
    vapply(fit$members[fit$kinds == "endpoint_cell"],
           function(cell) as.numeric(cell$patients), numeric(1))
  }
  expect_true(all(held(strict) >= 3L))
  expect_true(any(held(loose) == 2L))
  expect_gt(length(held(loose)), length(held(strict)))
})

# The basis is not capped the way AVATAR's distance metric is: every nominal
# time a large enough share of the cohort reached becomes a feature, where
# `.common_grid()` would have quantiled them down to fifteen.
test_that("the grid is the whole nominal grid, not the AVATAR profile width", {
  fixture <- pca_fixture(80)
  synthetic <- synpmx_pca(fixture$data, fixture$roles, seed = 11,
                          min_column_patients = 3L)
  fit <- attr(synthetic, "pmx_pca_fit")
  cells <- fit$members[fit$kinds == "endpoint_cell"]
  for (endpoint in unique(vapply(cells, function(c) c$endpoint, character(1)))) {
    modelled <- sum(vapply(cells, function(c) identical(c$endpoint, endpoint),
                           logical(1)))
    available <- length(fit$grids[[endpoint]])
    expect_equal(modelled, available, info = endpoint)
  }
})

test_that("components are capped at a fifth of the cohort", {
  fixture <- pca_fixture(30)
  fit <- attr(synpmx_pca(fixture$data, fixture$roles, seed = 2),
              "pmx_pca_fit")
  expect_lte(fit$k, 6L)
})

test_that("the report inventories every released quantity", {
  fixture <- pca_fixture(60)
  report <- pmx_pca_report(synpmx_pca(fixture$data, fixture$roles, seed = 1))
  expect_s3_class(report, "pmx_pca_report")
  expect_true(all(c("quantity", "what", "numbers", "min_patients") %in%
                    names(report)))
  expect_true(all(report$min_patients[!is.na(report$min_patients)] >= 1))
})

test_that("components report one loading per retained cell and component", {
  fixture <- pca_fixture(60)
  synthetic <- synpmx_pca(fixture$data, fixture$roles, seed = 1)
  fit <- attr(synthetic, "pmx_pca_fit")
  components <- pmx_pca_components(synthetic)
  expect_equal(nrow(components),
               sum(fit$kinds == "endpoint_cell") * fit$k)
  expect_equal(nrow(attr(components, "variance_explained")), fit$k)
})

test_that("a source too small to fit a basis is refused rather than fitted", {
  fixture <- pca_fixture(8)
  expect_error(synpmx_pca(fixture$data, fixture$roles, seed = 1),
               "at least 10 subjects")
})

# `dose_term = "log"` spends one coefficient instead of one mean per arm. Both
# must run; which one is faithful is a property of the study, not of the code.
test_that("both dose terms produce a legal table", {
  fixture <- pca_fixture(60)
  for (term in c("factor", "log")) {
    synthetic <- synpmx_pca(fixture$data, fixture$roles, seed = 4,
                            dose_term = term)
    expect_true(validate_pmx(synthetic, fixture$roles)$valid, info = term)
  }
})

# The owner's rule, 2026-08-26: an arm of one or two patients is not modelled.
test_that("an arm below the minimum is refused rather than generated from", {
  fixture <- pca_fixture(60)
  data <- fixture$data
  data$ARM <- ifelse(data$ID %in% unique(data$ID)[1:2], "rare", "main")
  roles <- pmx_roles(
    id = "ID", time = "TIME", nominal_time = "NTIME", tad = "TAD",
    occasion = "OCC", dv = "DV", amt = "AMT", rate = "RATE", evid = "EVID",
    cmt = "CMT", dvid = "DVID", mdv = "MDV", cens = "CENS",
    covariates = c("WT", "AGE", "SEX"), strata = "ARM"
  )
  expect_error(synpmx_pca(data, roles, seed = 1),
               "at least 3 patients in every arm")
  expect_true(validate_pmx(
    synpmx_pca(data, roles, seed = 1, min_arm_patients = 2L), roles
  )$valid)
})

# The nominal grid is the model's axis, so it is declared rather than inferred.
test_that("a missing nominal_time is refused rather than derived", {
  fixture <- pca_fixture(60)
  roles <- pmx_roles(
    id = "ID", time = "TIME", tad = "TAD", occasion = "OCC", dv = "DV",
    amt = "AMT", rate = "RATE", evid = "EVID", cmt = "CMT", dvid = "DVID",
    mdv = "MDV", cens = "CENS", covariates = c("WT", "AGE", "SEX")
  )
  expect_error(synpmx_pca(fixture$data, roles, seed = 1),
               "requires `nominal_time`")

  gappy <- fixture$data
  gappy$NTIME[gappy$EVID == 0][1:5] <- NA_real_
  expect_error(synpmx_pca(gappy, fixture$roles, seed = 1),
               "`nominal_time` is missing on")
})

# Dose rows and observation rows are placed on the same declared grid, so the
# event structure that makes a dataset fittable survives.
test_that("dose and observation timing survive generation", {
  fixture <- pca_fixture(60)
  synthetic <- synpmx_pca(fixture$data, fixture$roles, seed = 9)
  roles <- fixture$roles

  spans <- function(data) {
    doses <- data[data[[roles$evid]] != 0, , drop = FALSE]
    observed <- data[.observation_rows(data, roles, require_present = TRUE), ,
                     drop = FALSE]
    by_subject <- function(frame, column, fun) {
      vapply(split(frame[[column]], frame[[roles$id]]), fun, numeric(1))
    }
    list(
      doses = mean(table(doses[[roles$id]])),
      last_dose = stats::median(by_subject(doses, roles$time, max)),
      max_obs = stats::median(by_subject(observed, roles$time, max))
    )
  }
  source_spans <- spans(fixture$data)
  synthetic_spans <- spans(synthetic)
  expect_equal(synthetic_spans$doses, source_spans$doses, tolerance = 0.05)
  expect_equal(synthetic_spans$last_dose, source_spans$last_dose)
  expect_equal(synthetic_spans$max_obs, source_spans$max_obs)

  # Every observation sits on a nominal time the source used.
  observed <- synthetic[.observation_rows(synthetic, roles,
                                          require_present = TRUE), ]
  expect_true(all(observed[[roles$time]] %in% fixture$data[[roles$nominal_time]]))
})

# The owner's requirement, 2026-08-26: generation is a function of the summaries
# and nothing else. This is the check that says so rather than the comment.
test_that("generation uses the model alone and never the source rows", {
  fixture <- pca_fixture(60)
  synthetic <- synpmx_pca(fixture$data, fixture$roles, seed = 21)
  model <- attr(synthetic, "pmx_pca_model")
  expect_s3_class(model, "pmx_pca_model")

  # The generator takes the model and a count. There is no argument through
  # which a patient row could reach it.
  expect_identical(names(formals(synpmx:::.pca_generate)),
                   c("model", "n_subjects"))

  # Regenerating from the model alone reproduces the table exactly.
  again <- .with_local_seed(21, synpmx:::.pca_generate(model, 60L))
  expect_equal(as.data.frame(synthetic), as.data.frame(again),
               ignore_attr = TRUE)
})

test_that("the model holds no per-patient rows", {
  fixture <- pca_fixture(60)
  model <- attr(synpmx_pca(fixture$data, fixture$roles, seed = 1),
                "pmx_pca_model")
  n_source <- length(unique(fixture$data$ID))
  lengths_in <- function(x) {
    if (is.list(x) && !is.data.frame(x)) {
      return(unlist(lapply(x, lengths_in), use.names = FALSE))
    }
    if (is.data.frame(x)) return(nrow(x))
    length(x)
  }
  # Every stored vector is a summary: nothing is one value per source patient.
  expect_false(any(lengths_in(model$dosing) == n_source))
  expect_equal(nrow(model$schema$prototypes[[1L]]), NULL)
  expect_true(all(vapply(model$schema$prototypes, length, integer(1)) == 0L))
})

test_that("the dosing model is the arm's shared schedule, and says so", {
  fixture <- pca_fixture(60)
  synthetic <- synpmx_pca(fixture$data, fixture$roles, seed = 1)
  dosing <- pmx_pca_dosing(synthetic)
  expect_true(all(c("arm", "dose", "time", "amt", "share", "patients",
                    "distinct") %in% names(dosing)))
  expect_true(all(dosing$share > 0 & dosing$share <= 1))

  # The generated dose times are exactly the model's, arm by arm.
  roles <- fixture$roles
  doses <- synthetic[synthetic[[roles$evid]] != 0, , drop = FALSE]
  expect_setequal(unique(doses[[roles$time]]), unique(dosing$time))
})

test_that("the visit model covers every modelled cell in every arm", {
  fixture <- pca_fixture(60)
  synthetic <- synpmx_pca(fixture$data, fixture$roles, seed = 1)
  model <- attr(synthetic, "pmx_pca_model")
  visits <- pmx_pca_visits(synthetic)
  cells <- sum(model$basis$kinds == "endpoint_cell")
  expect_equal(nrow(visits), cells * length(model$arms$arms))
  expect_true(all(visits$probability >= 0 & visits$probability <= 1))
})

# A placebo arm records its administrations with `AMT = 0`. `.dose_rows()` drops
# those whenever any positive amount exists, which left the arm with no dosing
# events at all; the schedule is read from every dosing event instead.
test_that("a zero-amount arm keeps its dosing events", {
  fixture <- pca_fixture(60)
  data <- fixture$data
  placebo <- unique(data$ID)[1:20]
  data$AMT[data$ID %in% placebo & data$EVID != 0] <- 0
  data$ARM <- ifelse(data$ID %in% placebo, "placebo", "active")
  roles <- pmx_roles(
    id = "ID", time = "TIME", nominal_time = "NTIME", tad = "TAD",
    occasion = "OCC", dv = "DV", amt = "AMT", rate = "RATE", evid = "EVID",
    cmt = "CMT", dvid = "DVID", mdv = "MDV", cens = "CENS",
    covariates = c("WT", "AGE", "SEX"), strata = "ARM"
  )
  synthetic <- synpmx_pca(data, roles, seed = 2)

  expect_gt(sum(synthetic$EVID != 0 & synthetic$ARM == "placebo"), 0)
  expect_equal(sum(synthetic$EVID != 0), sum(data$EVID != 0))
  dosing <- pmx_pca_dosing(synthetic)
  expect_setequal(unique(dosing$arm), c("placebo", "active"))
  expect_true(all(dosing$amt[dosing$arm == "placebo"] == 0))
})

# The two stages are the public interface: summarize, look, then generate.
test_that("summarize and generate compose to the same table as synpmx_pca", {
  fixture <- pca_fixture(60)
  one_call <- synpmx_pca(fixture$data, fixture$roles, seed = 31)
  model <- synpmx_pca_summarize(fixture$data, fixture$roles)
  two_calls <- synpmx_pca_generate(model, seed = 31)
  expect_equal(as.data.frame(one_call), as.data.frame(two_calls),
               ignore_attr = TRUE)
})

test_that("generate takes only a model and a count", {
  expect_identical(names(formals(synpmx_pca_generate)),
                   c("model", "n_subjects", "seed"))
  fixture <- pca_fixture(60)
  model <- synpmx_pca_summarize(fixture$data, fixture$roles)
  expect_s3_class(model, "pmx_pca_model")
  expect_equal(model$n_source, 60L)

  # It defaults to the size it was fitted on, and honours a different one.
  expect_equal(
    length(unique(synpmx_pca_generate(model, seed = 1)[[fixture$roles$id]])),
    60L
  )
  expect_equal(
    length(unique(
      synpmx_pca_generate(model, n_subjects = 25L, seed = 1)[[fixture$roles$id]]
    )),
    25L
  )
})

test_that("the model prints what it holds", {
  fixture <- pca_fixture(60)
  model <- synpmx_pca_summarize(fixture$data, fixture$roles)
  printed <- paste(utils::capture.output(print(model)), collapse = "\n")
  expect_match(printed, "60 patients")
  expect_match(printed, "components")
  expect_match(printed, "dosing")
})

test_that("settings passed through synpmx_pca reach the summary", {
  fixture <- pca_fixture(60)
  synthetic <- synpmx_pca(fixture$data, fixture$roles, seed = 1,
                          n_components = 2L, dose_term = "log")
  model <- attr(synthetic, "pmx_pca_model")
  expect_equal(model$basis$k, 2L)
  expect_equal(model$settings$dose_term, "log")
})

# Values below the assay limit are reported on the boundary with CENS = 1, as
# the source reports them. Without this an arm entirely below quantification
# comes back as a spread of small numbers rather than the flat line recorded.
test_that("censoring is reapplied to generated values", {
  fixture <- pca_fixture(60)
  data <- fixture$data
  roles <- fixture$roles
  observed <- data$EVID == 0 & !is.na(data$DV) & data$DVID == "cp"
  lloq <- stats::quantile(data$DV[observed], 0.4, names = FALSE)
  below <- observed & data$DV < lloq
  data$DV[below] <- lloq
  data$CENS[below] <- 1L

  synthetic <- synpmx_pca(data, roles, seed = 41)
  generated <- synthetic$EVID == 0 & !is.na(synthetic$DV) &
    synthetic$DVID == "cp"

  # Something was censored, and everything censored sits on the boundary.
  expect_gt(sum(synthetic$CENS[generated] == 1L), 0L)
  on_limit <- synthetic$DV[generated & synthetic$CENS == 1L]
  expect_true(all(abs(on_limit - lloq) < 1e-8))
  # Nothing is emitted below the limit at all.
  expect_gte(min(synthetic$DV[generated]), lloq - 1e-8)
  # The censored fraction tracks the source rather than collapsing or vanishing.
  expect_equal(mean(synthetic$CENS[generated] == 1L),
               mean(data$CENS[observed] == 1L), tolerance = 0.15)
})

test_that("the assay limit is in the model and in the report", {
  fixture <- pca_fixture(60)
  data <- fixture$data
  observed <- data$EVID == 0 & !is.na(data$DV) & data$DVID == "cp"
  lloq <- stats::quantile(data$DV[observed], 0.4, names = FALSE)
  data$DV[observed & data$DV < lloq] <- lloq
  data$CENS[observed & data$DV <= lloq] <- 1L

  model <- synpmx_pca_summarize(data, fixture$roles, seed = 1)
  expect_equal(model$schema$censoring$cp$left, lloq)
  expect_true("assay limits" %in% pmx_pca_report(model)$quantity)
})

# Summarizing is stochastic only through the censoring imputation, so it takes a
# seed of its own and honours it.
test_that("summarize is reproducible under its seed", {
  fixture <- pca_fixture(60)
  a <- synpmx_pca_summarize(fixture$data, fixture$roles, seed = 5)
  b <- synpmx_pca_summarize(fixture$data, fixture$roles, seed = 5)
  expect_equal(a$basis$rotation, b$basis$rotation)
  expect_equal(a$scores, b$scores)
})
