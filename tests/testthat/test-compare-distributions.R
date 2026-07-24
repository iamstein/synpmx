# compare_pmx_distributions(): per-covariate and per-endpoint summary tables,
# source vs synthetic. private_fixture() carries two endpoints (cp, pd), two
# numeric covariates (WT, AGE), and one categorical (SEX).

test_that("the summary covers every endpoint and covariate for both datasets", {
  source <- private_fixture()
  roles <- private_roles()
  synthetic <- suppressWarnings(synpmx_avatar(source, roles, seed = 11))

  summary <- compare_pmx_distributions(source, synthetic, roles)
  expect_s3_class(summary, "pmx_distribution_summary")

  # Both endpoints, both datasets.
  expect_setequal(unique(summary$endpoints$variable), c("cp", "pd"))
  expect_setequal(unique(summary$endpoints$dataset), c("source", "synthetic"))
  expect_equal(nrow(summary$endpoints), 4L)
  expect_true(all(c("n", "n_subjects", "mean", "sd", "min", "q25",
                    "median", "q75", "max") %in% names(summary$endpoints)))

  # Numeric covariates split from the categorical one.
  expect_setequal(unique(summary$covariates_numeric$variable), c("WT", "AGE"))
  expect_equal(unique(summary$covariates_categorical$variable), "SEX")
  expect_setequal(unique(summary$covariates_categorical$level),
                  c("female", "male"))
})

test_that("summary statistics match a hand computation on the source", {
  source <- private_fixture()
  roles <- private_roles()

  summary <- compare_pmx_distributions(source, roles = roles)

  # WT is one baseline value per subject: 55 + 4 * subject for subjects 1..8.
  wt_row <- summary$covariates_numeric[
    summary$covariates_numeric$variable == "WT" &
      summary$covariates_numeric$dataset == "source", ]
  wt <- 55 + 4 * seq_len(8L)
  expect_equal(wt_row$n, 8L)
  expect_equal(wt_row$mean, mean(wt))
  expect_equal(wt_row$median, stats::median(wt))
  expect_equal(wt_row$max, max(wt))

  # cp endpoint: DV on observation rows with a present, non-event value.
  selected <- source$EVID == 0 & source$MDV == 0 & !is.na(source$DV) &
    as.character(source$DVID) == "cp"
  cp_dv <- source$DV[selected]
  cp_row <- summary$endpoints[summary$endpoints$variable == "cp", ]
  expect_equal(cp_row$n, length(cp_dv))
  expect_equal(cp_row$mean, mean(cp_dv))
  expect_equal(cp_row$n_subjects, length(unique(source$ID[selected])))
})

test_that("categorical proportions sum to one within each variable/dataset", {
  source <- private_fixture()
  roles <- private_roles()
  synthetic <- suppressWarnings(synpmx_avatar(source, roles, seed = 12))

  cat <- compare_pmx_distributions(source, synthetic, roles)$covariates_categorical
  by_group <- split(cat$proportion, interaction(cat$variable, cat$dataset,
                                                 drop = TRUE))
  for (props in by_group) expect_equal(sum(props), 1)
  # Counts are per subject, so they total the subject count, not the row count.
  totals <- tapply(cat$n, interaction(cat$variable, cat$dataset, drop = TRUE),
                   sum)
  expect_true(all(totals == length(unique(source$ID))))
})

test_that("synthetic = NULL summarizes the source alone", {
  source <- private_fixture()
  roles <- private_roles()

  summary <- compare_pmx_distributions(source, roles = roles)
  expect_equal(unique(summary$endpoints$dataset), "source")
  expect_equal(unique(summary$covariates_numeric$dataset), "source")
})

test_that("a dataset with no covariates yields NULL covariate tables", {
  source <- private_fixture()
  roles <- private_roles()
  roles$covariates <- NULL

  summary <- compare_pmx_distributions(source, roles = roles)
  expect_null(summary$covariates_numeric)
  expect_null(summary$covariates_categorical)
  expect_false(is.null(summary$endpoints))
})

test_that("every table is marked restricted, and print returns invisibly", {
  source <- private_fixture()
  roles <- private_roles()
  synthetic <- suppressWarnings(synpmx_avatar(source, roles, seed = 13))

  summary <- compare_pmx_distributions(source, synthetic, roles)
  expect_equal(attr(summary$endpoints, "release_status"),
               "restricted_not_releasable")
  expect_equal(attr(summary$covariates_numeric, "release_status"),
               "restricted_not_releasable")

  expect_output(out <- withVisible(print(summary)), "distribution summary")
  expect_false(out$visible)
  expect_identical(out$value, summary)
})
