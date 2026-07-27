# REV-025: `k` sets how many real patients are blended, but `max_donor_weight`
# is what bounds any one patient's contribution to an avatar. The cap must hold
# for every donor, not only the largest -- capping just the maximum leaves the
# runner-up above the cap once the excess is redistributed.

test_that("no donor exceeds the cap, not just the largest one", {
  # Capping only which.max() would return (0.30, 0.35, 0.21, 0.084, 0.056):
  # the stated maximum, violated by the donor the redistribution created.
  weights <- synpmx:::.cap_weights(c(0.50, 0.25, 0.15, 0.06, 0.04), 0.30)

  expect_true(all(weights <= 0.30 + 1e-9))
  expect_equal(sum(weights), 1, tolerance = 1e-12)
  expect_equal(weights, c(0.30, 0.30, 0.24, 0.096, 0.064), tolerance = 1e-9)
})

test_that("a cap no weight vector can satisfy relaxes to uniform", {
  # Three donors cannot all sit at or below 0.30; 1/3 each is the flattest
  # blend available, and is returned instead of an error or a violated cap.
  expect_equal(synpmx:::.cap_weights(c(0.99, 0.005, 0.005), 0.30),
               rep(1 / 3, 3), tolerance = 1e-9)
  expect_equal(synpmx:::.cap_weights(c(0.60, 0.40), 0.30),
               c(0.5, 0.5), tolerance = 1e-9)
})

test_that("weights already under the cap are left alone", {
  weights <- c(0.28, 0.24, 0.20, 0.16, 0.12)

  expect_equal(synpmx:::.cap_weights(weights, 0.30), weights, tolerance = 1e-12)
})

test_that("the cap holds through the sampled weights, over many draws", {
  set.seed(11)
  for (trial in seq_len(200)) {
    weights <- synpmx:::.randomized_weights(runif(5, 0.1, 5), max_weight = 0.30)
    expect_true(max(weights) <= 0.30 + 1e-9)
    expect_equal(sum(weights), 1, tolerance = 1e-12)
  }
})

test_that("a single donor is unaffected: the cap cannot invent a second", {
  expect_equal(synpmx:::.randomized_weights(0.5), 1)
  expect_equal(synpmx:::.cap_weights(1, 0.30), 1)
})

test_that("synpmx_avatar records the cap and the effective donor count", {
  skip_if_not_installed("nlmixr2data")
  data("theo_md", package = "nlmixr2data", envir = environment())
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", cmt = "CMT", covariates = "WT")

  settings <- attr(
    suppressWarnings(synpmx_avatar(theo_md, roles, seed = 101)), "pmx_settings"
  )

  expect_equal(settings$max_donor_weight, 0.30)
  # 1 / sum(w^2) cannot exceed k, and a 0.30 cap must keep it above the 1/0.30
  # a maximally concentrated capped vector would give.
  expect_lte(settings$mean_effective_donors, settings$k)
  expect_gte(settings$min_effective_donors, 1 / 0.30)
})

test_that("max_donor_weight is validated", {
  skip_if_not_installed("nlmixr2data")
  data("theo_md", package = "nlmixr2data", envir = environment())
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", cmt = "CMT", covariates = "WT")

  for (bad in list(0, -0.1, 1.5, NA_real_, c(0.3, 0.4))) {
    expect_error(
      synpmx_avatar(theo_md, roles, max_donor_weight = bad),
      "max_donor_weight"
    )
  }
})
