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

  basis_of <- function(minimum) {
    attr(synpmx_pca(sparse, fixture$roles, seed = 5,
                    min_column_patients = minimum), "pmx_trial_summary")$basis
  }
  strict <- basis_of(3L)
  loose <- basis_of(2L)
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
  fit <- attr(synthetic, "pmx_trial_summary")$basis
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
              "pmx_trial_summary")$basis
  expect_lte(fit$k, 6L)
})

test_that("the report inventories every released quantity", {
  fixture <- pca_fixture(60)
  report <- pca_report(synpmx_pca(fixture$data, fixture$roles, seed = 1))
  expect_s3_class(report, "pca_report")
  expect_true(all(c("quantity", "what", "numbers", "min_patients") %in%
                    names(report)))
  expect_true(all(report$min_patients[!is.na(report$min_patients)] >= 1))
})

test_that("components report one loading per feature and component", {
  fixture <- pca_fixture(60)
  synthetic <- synpmx_pca(fixture$data, fixture$roles, seed = 1)
  fit <- attr(synthetic, "pmx_trial_summary")$basis
  components <- pca_components(synthetic)
  # Every feature, covariates included: a component loading on both a covariate
  # and an endpoint is what carries their relationship into the output, and a
  # cell-only table cannot show it.
  expect_equal(nrow(components), length(fit$columns) * fit$k)
  expect_true(any(components$kind != "endpoint_cell"))
  expect_true(all(is.na(components$time[components$kind != "endpoint_cell"])))
  expect_equal(nrow(attr(components, "variance_explained")), fit$k)

  # Each component's squared loadings sum to one, so the mass per block reads
  # as a share.
  mass <- tapply(components$loading^2, components$component, sum)
  expect_true(all(abs(mass - 1) < 1e-8))
})

test_that("features and scores cover the grid and the score model", {
  fixture <- pca_fixture(60)
  trial_summary <- synpmx_pca_summarize(fixture$data, fixture$roles, seed = 1)
  features <- pca_features(trial_summary)
  scores <- pca_scores(trial_summary)

  expect_equal(nrow(features), length(trial_summary$basis$columns))
  expect_true(all(c("feature", "kind", "endpoint", "time", "covariate",
                    "level", "patients", "center", "scale", "transform") %in%
                    names(features)))
  expect_true(all(features$scale > 0))

  expect_equal(nrow(scores),
               length(trial_summary$arms$arms) * trial_summary$basis$k)
  expect_true(all(scores$sd > 0))

  # Both dose terms report the same shape, so two runs can be compared.
  logged <- pca_scores(synpmx_pca_summarize(fixture$data, fixture$roles,
                                            seed = 1, dose_term = "log"))
  expect_identical(names(logged), names(scores))
  expect_equal(nrow(logged), nrow(scores))
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

# The nominal grid is the trial_summary's axis, so it is declared rather than inferred.
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
test_that("generation uses the trial_summary alone and never the source rows", {
  fixture <- pca_fixture(60)
  synthetic <- synpmx_pca(fixture$data, fixture$roles, seed = 21)
  trial_summary <- attr(synthetic, "pmx_trial_summary")
  expect_s3_class(trial_summary, "pmx_trial_summary")

  # The generator takes the trial_summary and a count. There is no argument through
  # which a patient row could reach it.
  expect_identical(names(formals(synpmx:::.pca_generate)),
                   c("trial_summary", "n_subjects"))

  # Regenerating from the trial_summary alone reproduces the table exactly.
  again <- .with_local_seed(21, synpmx:::.pca_generate(trial_summary, 60L))
  expect_equal(as.data.frame(synthetic), as.data.frame(again),
               ignore_attr = TRUE)
})

test_that("the trial_summary holds no per-patient rows", {
  fixture <- pca_fixture(60)
  trial_summary <- attr(synpmx_pca(fixture$data, fixture$roles, seed = 1),
                "pmx_trial_summary")
  n_source <- length(unique(fixture$data$ID))
  lengths_in <- function(x) {
    if (is.list(x) && !is.data.frame(x)) {
      return(unlist(lapply(x, lengths_in), use.names = FALSE))
    }
    if (is.data.frame(x)) return(nrow(x))
    length(x)
  }
  # Every stored vector is a summary: nothing is one value per source patient.
  expect_false(any(lengths_in(trial_summary$dosing) == n_source))
  expect_equal(nrow(trial_summary$schema$prototypes[[1L]]), NULL)
  expect_true(all(vapply(trial_summary$schema$prototypes, length, integer(1)) == 0L))
})

test_that("the dosing trial_summary is the arm's planned schedule, and says so", {
  fixture <- pca_fixture(60)
  synthetic <- synpmx_pca(fixture$data, fixture$roles, seed = 1)
  dosing <- pca_dosing(synthetic)
  rates <- pca_dose_rates(synthetic)
  expect_true(all(c("arm", "cycle", "time", "planned_amt") %in% names(dosing)))
  expect_true(all(c("arm", "planned_cycles", "levels", "discontinuation",
                    "interruption", "reduction", "patients", "source_doses",
                    "distinct") %in% names(rates)))
  for (column in c("discontinuation", "interruption", "reduction")) {
    expect_true(all(rates[[column]] >= 0 & rates[[column]] <= 1),
                info = column)
  }

  # The generated dose times come from the planned grid and nowhere else.
  roles <- fixture$roles
  doses <- synthetic[synthetic[[roles$evid]] != 0, , drop = FALSE]
  expect_true(all(doses[[roles$time]] %in% dosing$time))
})

test_that("the visit trial_summary covers every modelled cell in every arm", {
  fixture <- pca_fixture(60)
  synthetic <- synpmx_pca(fixture$data, fixture$roles, seed = 1)
  trial_summary <- attr(synthetic, "pmx_trial_summary")
  visits <- pca_visits(synthetic)
  cells <- sum(trial_summary$basis$kinds == "endpoint_cell")
  expect_equal(nrow(visits), cells * length(trial_summary$arms$arms))
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
  dosing <- pca_dosing(synthetic)
  expect_setequal(unique(dosing$arm), c("placebo", "active"))
  expect_true(all(dosing$amt[dosing$arm == "placebo"] == 0))
})

# The two stages are the public interface: summarize, look, then generate.
test_that("summarize and generate compose to the same table as synpmx_pca", {
  fixture <- pca_fixture(60)
  one_call <- synpmx_pca(fixture$data, fixture$roles, seed = 31)
  trial_summary <- synpmx_pca_summarize(fixture$data, fixture$roles)
  two_calls <- synpmx_pca_generate(trial_summary, seed = 31)
  expect_equal(as.data.frame(one_call), as.data.frame(two_calls),
               ignore_attr = TRUE)
})

test_that("generate takes only a trial_summary and a count", {
  expect_identical(names(formals(synpmx_pca_generate)),
                   c("trial_summary", "n_subjects", "seed"))
  fixture <- pca_fixture(60)
  trial_summary <- synpmx_pca_summarize(fixture$data, fixture$roles)
  expect_s3_class(trial_summary, "pmx_trial_summary")
  expect_equal(trial_summary$n_source, 60L)

  # It defaults to the size it was fitted on, and honours a different one.
  expect_equal(
    length(unique(synpmx_pca_generate(trial_summary, seed = 1)[[fixture$roles$id]])),
    60L
  )
  expect_equal(
    length(unique(
      synpmx_pca_generate(trial_summary, n_subjects = 25L, seed = 1)[[fixture$roles$id]]
    )),
    25L
  )
})

test_that("the trial_summary prints what it holds", {
  fixture <- pca_fixture(60)
  trial_summary <- synpmx_pca_summarize(fixture$data, fixture$roles)
  printed <- paste(utils::capture.output(print(trial_summary)), collapse = "\n")
  expect_match(printed, "60 patients")
  expect_match(printed, "components")
  expect_match(printed, "dosing")
})

test_that("settings passed through synpmx_pca reach the summary", {
  fixture <- pca_fixture(60)
  synthetic <- synpmx_pca(fixture$data, fixture$roles, seed = 1,
                          n_components = 2L, dose_term = "log")
  trial_summary <- attr(synthetic, "pmx_trial_summary")
  expect_equal(trial_summary$basis$k, 2L)
  expect_equal(trial_summary$settings$dose_term, "log")
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

test_that("the assay limit is in the trial_summary and in the report", {
  fixture <- pca_fixture(60)
  data <- fixture$data
  observed <- data$EVID == 0 & !is.na(data$DV) & data$DVID == "cp"
  lloq <- stats::quantile(data$DV[observed], 0.4, names = FALSE)
  data$DV[observed & data$DV < lloq] <- lloq
  data$CENS[observed & data$DV <= lloq] <- 1L

  trial_summary <- synpmx_pca_summarize(data, fixture$roles, seed = 1)
  expect_equal(trial_summary$schema$censoring$cp$left, lloq)
  expect_true("assay limits" %in% pca_report(trial_summary)$quantity)
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

# The log transform's offset comes from the reported values, not the imputed
# ones. `.impute_censored()` draws uniformly on the raw scale, so a value
# censored at 0.05 can land at 5e-05 and drag the offset three orders of
# magnitude down. The log range then runs far below anything observed, the
# Gaussian score residual is fitted on that inflated spread, and exponentiating
# it puts detectable concentrations into an arm that is entirely below the limit.
test_that("censored imputation cannot widen the log transform", {
  fixture <- pca_fixture(60)
  data <- fixture$data
  roles <- fixture$roles
  observed <- data$EVID == 0 & !is.na(data$DV) & data$DVID == "cp"

  # A late window nobody is above the limit in, and a high limit so most of the
  # endpoint is censored.
  lloq <- stats::quantile(data$DV[observed], 0.6, names = FALSE)
  late <- observed & data$NTIME >= stats::median(data$NTIME[observed])
  data$DV[observed & data$DV < lloq] <- lloq
  data$CENS[observed & data$DV <= lloq] <- 1L
  data$DV[late] <- lloq
  data$CENS[late] <- 1L

  trial_summary <- synpmx_pca_summarize(data, roles, seed = 3)
  # The offset is half the smallest positive REPORTED value, so it is on the
  # scale of the limit rather than on the scale of a uniform draw below it.
  expect_gte(trial_summary$basis$transforms$cp$offset, lloq / 4)

  synthetic <- synpmx_pca_generate(trial_summary, seed = 3)
  generated_late <- synthetic$EVID == 0 & !is.na(synthetic$DV) &
    synthetic$DVID == "cp" &
    synthetic$NTIME >= stats::median(data$NTIME[observed])
  # A window the whole cohort was below the limit in comes back below it.
  expect_gt(mean(synthetic$CENS[generated_late] == 1L), 0.8)
})

# An oncology-shaped source: a planned q21d schedule that almost nobody
# completes, with reductions to 75% and 50% of the starting dose, skipped
# cycles, and discontinuation. The variability in the schedule is what such a
# dataset is usually for, so it has to survive generation.
oncology_fixture <- function(n = 80, reduction = 0.10, interruption = 0.12,
                             discontinuation = 0.09, seed = 11) {
  set.seed(seed)
  cycles <- (0:11) * 21 * 24
  subject <- function(id, arm) {
    start <- if (arm == "A") 100 else 200
    level <- 1L
    rows <- list()
    for (i in seq_along(cycles)) {
      if (level < 3L && stats::runif(1) < reduction) level <- level + 1L
      if (i == 1L || stats::runif(1) >= interruption) {
        rows[[length(rows) + 1L]] <- data.frame(
          TIME = cycles[i], NTIME = cycles[i],
          AMT = start * c(1, 0.75, 0.5)[level], EVID = 1,
          DV = NA_real_, DVID = NA_character_
        )
      }
      if (i > 1L && stats::runif(1) < discontinuation) break
    }
    doses <- do.call(rbind, rows)
    seen <- cycles[cycles <= max(doses$TIME)]
    observations <- data.frame(
      TIME = seen + 24, NTIME = seen + 24, AMT = 0, EVID = 0,
      DV = exp(stats::rnorm(length(seen), log(start / 50), 0.4)), DVID = "cp"
    )
    out <- rbind(doses, observations)
    out$ID <- id
    out$ARM <- arm
    out[order(out$TIME, out$EVID == 0), , drop = FALSE]
  }
  do.call(rbind, lapply(seq_len(n), function(i) {
    subject(i, if (i <= n / 2) "A" else "B")
  }))
}

oncology_roles <- function() {
  pmx_roles(id = "ID", time = "TIME", nominal_time = "NTIME", dv = "DV",
            amt = "AMT", evid = "EVID", dvid = "DVID", strata = "ARM")
}

test_that("dose reductions and interruptions are recovered as rates", {
  data <- oncology_fixture()
  trial_summary <- synpmx_pca_summarize(data, oncology_roles(), seed = 1)
  rates <- pca_dose_rates(trial_summary)

  # The ladder is the one the fixture used, and nothing else: reductions are
  # within-patient decreases, so ordinary between-patient amount differences
  # cannot invent a level.
  for (arm in rates$arm) {
    levels <- as.numeric(strsplit(rates$levels[rates$arm == arm], ", ")[[1L]])
    expect_equal(levels, c(1, 0.75, 0.5), tolerance = 0.02, info = arm)
  }
  expect_true(all(abs(rates$reduction - 0.10) < 0.06))
  expect_true(all(abs(rates$interruption - 0.12) < 0.06))
  expect_true(all(abs(rates$discontinuation - 0.09) < 0.06))
})

test_that("schedule variability survives generation", {
  data <- oncology_fixture()
  roles <- oncology_roles()
  synthetic <- synpmx_pca(data, roles, seed = 1)
  expect_true(validate_pmx(synthetic, roles)$valid)

  key <- function(d) {
    doses <- d[d$EVID == 1, , drop = FALSE]
    length(unique(vapply(split(doses, doses$ID), function(p) {
      paste(sprintf("%.6g", p$TIME), sprintf("%.6g", p$AMT), collapse = "|")
    }, character(1))))
  }
  # The whole point: not one schedule per arm.
  expect_gt(key(synthetic), 20L)

  per_patient <- function(d) {
    doses <- d[d$EVID == 1, , drop = FALSE]
    counts <- vapply(split(doses, doses$ID), nrow, integer(1))
    reduced <- vapply(split(doses, doses$ID),
                      function(p) length(unique(p$AMT)) > 1L, logical(1))
    c(doses = mean(counts), spread = stats::sd(counts),
      reduced = mean(reduced))
  }
  source_shape <- per_patient(data)
  synthetic_shape <- per_patient(synthetic)
  expect_equal(synthetic_shape[["doses"]], source_shape[["doses"]],
               tolerance = 0.2)
  expect_equal(synthetic_shape[["spread"]], source_shape[["spread"]],
               tolerance = 0.3)
  expect_gt(synthetic_shape[["reduced"]], 0.25)

  # Nothing is generated past the last cycle the arm shares.
  planned <- pca_dosing(synthetic)
  expect_true(all(synthetic$TIME[synthetic$EVID == 1] %in% planned$time))
})

# The degenerate case, and it is the common one: a study with a fixed schedule
# must come through unchanged rather than acquiring dropout it never had.
test_that("a fixed schedule stays fixed", {
  fixture <- pca_fixture(60)
  synthetic <- synpmx_pca(fixture$data, fixture$roles, seed = 1)
  rates <- pca_dose_rates(synthetic)

  expect_true(all(rates$discontinuation == 0))
  expect_true(all(rates$interruption == 0))
  expect_true(all(rates$reduction == 0))
  expect_true(all(rates$levels == "1"))

  source_doses <- mean(table(fixture$data$ID[fixture$data$EVID != 0]))
  synthetic_doses <- mean(table(synthetic$ID[synthetic$EVID != 0]))
  expect_equal(synthetic_doses, source_doses)
})

# `pca_component_effect()` is the components on the scale the study reported in.
# It runs the generator's own inversion, so what it says a component does is
# what a drawn subject actually gets rather than a second implementation.
test_that("a component's effect is the model's own inversion of its loadings", {
  fixture <- pca_fixture(60)
  trial_summary <- synpmx_pca_summarize(fixture$data, fixture$roles, seed = 1)
  basis <- trial_summary$basis
  effect <- pca_component_effect(trial_summary, sds = 1)

  expect_setequal(unique(effect$score_sd), c(-1, 0, 1))
  expect_equal(nrow(effect), 3L * basis$k * length(basis$columns))
  expect_setequal(unique(effect$component), paste0("PC", seq_len(basis$k)))

  # The centre is score zero on every component, so it cannot depend on which
  # component is being displaced.
  centre <- subset(effect, score_sd == 0)
  by_component <- split(centre$value, centre$component)
  for (values in by_component[-1L]) expect_equal(values, by_component[[1L]])

  # And the centre is the inverted column mean, not something recomputed here.
  reference <- by_component[[1L]]
  expected <- basis$centers
  for (endpoint in names(basis$transforms)) {
    cells <- !is.na(effect$endpoint[seq_along(expected)]) &
      effect$endpoint[seq_along(expected)] == endpoint
    expected[cells] <- synpmx:::.inverse_dv(expected[cells],
                                            basis$transforms[[endpoint]])
  }
  expect_equal(reference, unname(expected))
})

test_that("a component's effect moves the features its loadings sit on", {
  fixture <- pca_fixture(60)
  trial_summary <- synpmx_pca_summarize(fixture$data, fixture$roles, seed = 1)
  basis <- trial_summary$basis
  effect <- pca_component_effect(trial_summary, sds = 2)

  # Both inverse transforms are increasing, so on the reported scale a feature
  # must move in its loading's direction and cannot move at all where the
  # loading is zero. That holds whichever transform a cell carries, which is
  # why it is the property asserted rather than a rank on the modelling scale.
  first <- subset(effect, component == "PC1")
  moved <- subset(first, score_sd == 2)$value - subset(first, score_sd == 0)$value
  loading <- basis$rotation[, 1L]
  expect_equal(sign(moved), sign(loading), ignore_attr = TRUE)

  # On the modelling scale the step is linear in `sds`, which is what makes the
  # argument readable as a distance along the component. Checked on the
  # identity-transform cells, where the reported scale is the modelling scale.
  one <- pca_component_effect(trial_summary, sds = 1)
  centre <- subset(one, component == "PC1" & score_sd == 0)$value
  step_one <- subset(one, component == "PC1" & score_sd == 1)$value
  identity_cells <- is.na(first$endpoint[first$score_sd == 0])
  expect_true(any(identity_cells))
  expect_equal(
    (subset(first, score_sd == 2)$value - centre)[identity_cells],
    2 * (step_one - centre)[identity_cells]
  )

  expect_error(pca_component_effect(trial_summary, sds = 0), "positive")
  expect_error(pca_component_effect(trial_summary, sds = c(1, 2)), "single")
})
