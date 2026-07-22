test_that("small groups and k larger than available donors use warnings and fallbacks", {
  roles <- fixture_roles()
  one <- pmx_fixture(n = 1L)
  expect_warning(
    mock <- mock_pmx(one, roles, n_subjects = 2, seed = 2, k = 20),
    "anchor"
  )
  expect_equal(length(unique(mock$ID)), 2L)
  expect_true(validate_pmx(mock, roles)$valid)

  two <- pmx_fixture(n = 2L)
  expect_warning(
    mock_pmx(two, roles, n_subjects = 2, seed = 2, k = 20),
    "fewer than two"
  )
})

test_that("zero-variance and duplicated profiles are handled safely", {
  source <- pmx_fixture(duplicated_profiles = TRUE)
  roles <- fixture_roles()
  expect_warning(
    mock <- mock_pmx(
      source, roles, n_subjects = 3, seed = 72,
      subject_noise_sd = 0, residual_noise_sd = 0
    ),
    "zero neighbor distances"
  )
  expect_true(all(is.finite(mock$DV)))
  expect_true(validate_pmx(mock, roles)$valid)
})

test_that("clearly different repeat-dose gaps form separate compatible groups", {
  source <- pmx_fixture(n = 8L)
  late_subject <- source$ID > 4L
  second_occasion <- late_subject & source$TIME >= 12
  source$TIME[second_occasion] <- source$TIME[second_occasion] + 12
  roles <- fixture_roles()

  mock <- mock_pmx(source, roles, n_subjects = 2, seed = 5, k = 2)
  expect_equal(attr(mock, "pmx_settings")$compatible_event_groups, 2L)
  expect_true(validate_pmx(mock, roles)$valid)
})

test_that("generator options fail clearly", {
  source <- pmx_fixture()
  roles <- fixture_roles()
  expect_error(mock_pmx(source, roles, n_subjects = 0), "positive integer")
  expect_error(mock_pmx(source, roles, k = 0), "positive integer")
  expect_error(mock_pmx(source, roles, pca_variance = 2), "in \\(0, 1\\]")
  expect_error(mock_pmx(source, roles, residual_phi = 1), "strictly between")
  expect_error(mock_pmx(source, roles, seed = -1), "integer from 0")
  expect_error(mock_pmx(source, roles, event_method = "pca"), "template")
  expect_error(mock_pmx(source, roles, dv_method = "gan"), "avatar_blend")
})

test_that("comparison reports structure and exploratory plots", {
  source <- pmx_fixture()
  roles <- fixture_roles()
  mock <- mock_pmx(source, roles, n_subjects = 3, seed = 90)
  comparison <- compare_pmx(source, mock, roles)

  expect_s3_class(comparison, "pmx_comparison")
  expect_equal(comparison$summary$dataset, c("source", "mock"))
  expect_true(all(comparison$column_classes$matches))
  expect_true(comparison$validation$mock$valid)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_named(comparison$plots, c("overlay", "faceted"))
    expect_s3_class(comparison$plots$overlay, "ggplot")
  }
})
