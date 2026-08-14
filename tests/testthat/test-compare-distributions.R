# compare_pmx_distributions(): per-covariate and per-endpoint distributions,
# source vs synthetic, drawn by default and tabulated on request.
# private_fixture() carries two endpoints (cp, pd), two numeric covariates
# (WT, AGE), and one categorical (SEX). Its numeric covariates hold one distinct
# value per subject, so `private_fixture(n)` is also the knob that walks a
# covariate across the eight-value bars/density boundary.

test_that("the summary covers every endpoint and covariate for both datasets", {
  source <- private_fixture()
  roles <- private_roles()
  synthetic <- suppressWarnings(synpmx_avatar(source, roles, seed = 11))

  summary <- compare_pmx_distributions(source, synthetic, roles,
                                      output = "tables")
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

  summary <- compare_pmx_distributions(source, roles = roles,
                                      output = "tables")

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

  cat <- compare_pmx_distributions(source, synthetic, roles,
                                 output = "tables")$covariates_categorical
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

  summary <- compare_pmx_distributions(source, roles = roles,
                                      output = "tables")
  expect_equal(unique(summary$endpoints$dataset), "source")
  expect_equal(unique(summary$covariates_numeric$dataset), "source")
})

test_that("a dataset with no covariates yields NULL covariate tables", {
  source <- private_fixture()
  roles <- private_roles()
  roles$covariates <- NULL

  summary <- compare_pmx_distributions(source, roles = roles,
                                      output = "tables")
  expect_null(summary$covariates_numeric)
  expect_null(summary$covariates_categorical)
  expect_false(is.null(summary$endpoints))
})

test_that("every table is marked restricted, and print returns invisibly", {
  source <- private_fixture()
  roles <- private_roles()
  synthetic <- suppressWarnings(synpmx_avatar(source, roles, seed = 13))

  summary <- compare_pmx_distributions(source, synthetic, roles,
                                      output = "tables")
  expect_equal(attr(summary$endpoints, "release_status"),
               "restricted_not_releasable")
  expect_equal(attr(summary$covariates_numeric, "release_status"),
               "restricted_not_releasable")

  expect_output(out <- withVisible(print(summary)), "distribution summary")
  expect_false(out$visible)
  expect_identical(out$value, summary)
})

# The figure ------------------------------------------------------------------

test_that("the default draws a figure covering every endpoint and covariate", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  source <- private_fixture()
  roles <- private_roles()
  synthetic <- suppressWarnings(synpmx_avatar(source, roles, seed = 21))

  figure <- compare_pmx_distributions(source, synthetic, roles)
  expect_s3_class(figure, "ggplot")
  expect_equal(attr(figure, "release_status"), "restricted_not_releasable")

  panels <- .distribution_panel_data(source, synthetic, roles)
  expect_setequal(unique(panels$variable), c("cp", "pd", "WT", "AGE", "SEX"))
  # Endpoints first, then continuous covariates, then categorical ones, so the
  # figure reads in the order the tables are printed in.
  expect_equal(unique(panels$section),
               c("endpoint", "continuous", "categorical"))
  expect_setequal(unique(panels$dataset), c("source", "synthetic"))
})

test_that("eight distinct values is the boundary between bars and a density", {
  roles <- private_roles()
  kind_of <- function(n, variable) {
    panels <- .distribution_panel_data(private_fixture(n), NULL, roles)
    unique(panels$kind[panels$variable == variable])
  }
  # WT is 55 + 4 * subject, so `n` subjects give exactly `n` distinct values.
  expect_equal(kind_of(8L, "WT"), "bars")
  expect_equal(kind_of(9L, "WT"), "density")
  # A categorical covariate is bars whatever its level count.
  expect_equal(kind_of(9L, "SEX"), "bars")
})

test_that("bar levels are shared, ordered numerically, and sum to one", {
  roles <- private_roles()
  panels <- .distribution_panel_data(private_fixture(4L), NULL, roles)
  wt <- panels[panels$variable == "WT", ]
  expect_equal(wt$kind[[1L]], "bars")
  # 59, 63, 67, 71 -- numeric order, not the 59/63/67/71 alphabetical accident
  # that would put 100 before 20 on a study with wider values.
  expect_equal(wt$level, format(55 + 4 * seq_len(4L), trim = TRUE))
  expect_equal(sum(wt$proportion), 1)
})

test_that("a wide positive spread takes a log10 axis and a narrow one does not", {
  roles <- private_roles()
  source <- private_fixture(12L)
  # Three orders of magnitude across twelve subjects: enough distinct values to
  # be a density, and a wide enough span to be drawn on a log axis.
  source$WT <- rep(10^seq(0, 3, length.out = 12L), each = 16L)
  panels <- .distribution_panel_data(source, NULL, roles)

  expect_true(unique(panels$log10[panels$variable == "WT"]))
  # AGE is 26..37 on the same fixture: positive, but nothing like as wide.
  expect_false(unique(panels$log10[panels$variable == "AGE"]))
  # The values stay in their own units -- the panel gets scale_x_log10(), so
  # the axis reads 1/10/100 rather than 0/1/2.
  expect_equal(max(panels$value[panels$variable == "WT"]), 1000)
})

test_that("synthetic = NULL draws the source alone", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  source <- private_fixture()
  roles <- private_roles()

  panels <- .distribution_panel_data(source, NULL, roles)
  expect_equal(unique(panels$dataset), "source")
  expect_s3_class(compare_pmx_distributions(source, roles = roles), "ggplot")
})

test_that("the height grows a row of panels at a time and has a floor", {
  roles <- private_roles()
  source <- private_fixture()

  # Five panels: two endpoints and three covariates, three across, two rows.
  expect_equal(compare_pmx_distributions_height(source, roles),
               0.6 + 2.1 * 2)
  # One covariate is one row, which the floor takes over.
  roles_one <- roles
  roles_one$covariates <- "WT"
  expect_equal(compare_pmx_distributions_height(source, roles_one), 3)
  expect_gt(compare_pmx_distributions_height(source, roles),
            compare_pmx_distributions_height(source, roles_one))
})

test_that("a missing plotting package names the way through", {
  expect_error(
    with_mocked_bindings(
      compare_pmx_distributions(private_fixture(), roles = private_roles()),
      requireNamespace = function(...) FALSE,
      .package = "base"
    ),
    'output = "tables"',
    fixed = TRUE
  )
})
