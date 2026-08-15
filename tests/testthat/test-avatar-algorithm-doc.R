# `vignettes/avatar-algorithm.Rmd` is the specification of `synpmx_avatar()`,
# and it prints numbers: the default donor count, the ceiling on one donor's
# share, how many levels make an endpoint ordinal. Nothing re-reads that prose
# when a default moves, so this file is what keeps the two in step -- it asserts
# from the code every number the vignette states.
#
# If this fails, the vignette is describing a package that no longer behaves
# that way. Fix the package or update the vignette -- do not relax the test.
#
# Four numbers the vignette states are NOT pinned here, because they are inline
# literals with no name to assert on. Pinning them means promoting them to named
# constants or asserting on behavior, and neither is done yet:
#
#   - the screen's reference quantile, the 90th, inside
#     `.structural_outlier_anchors()`; its `mult = 2` is pinned below, so
#     "twice the cohort's 90th percentile" is half covered.
#   - the covariate log-scale rule, `max / median > 3`, in
#     `.synthesize_covariates()`;
#   - the covariate noise fallback `max(0.05|mu|, 0.01)`, same function;
#   - the endpoint noise fallback `max(0.1|median|, 0.01)`, in
#     `.synthesize_trajectories()`.
#
# That is the known gap. Do not read this file as covering them.

# Every one of these appears as a number in the vignette's prose, which is the
# only reason it is here: a default the vignette does not name does not need
# pinning, and a default it names cannot be allowed to drift silently.
test_that("the defaults the vignette states are the defaults in the signature", {
  defaults <- formals(synpmx_avatar)

  # Step 8: k donors, and Step 6's route floor of `k + 1`.
  expect_equal(defaults$k, 5)
  # Step 9: the ceiling on any one donor's share of the avatar.
  expect_equal(defaults$max_donor_weight, 0.50)
  # Step 4: the retained-variance target for the PCA.
  expect_equal(defaults$pca_variance, 0.90)
  # Steps 10 and 11: the subject-level and residual noise scales, and the AR(1)
  # correlation between consecutive residuals.
  expect_equal(defaults$subject_noise_sd, 0.15)
  expect_equal(defaults$residual_noise_sd, 0.05)
  expect_equal(defaults$residual_phi, 0.6)
  # Step 7: the floor on how many patients must share a visit set, a shape, or
  # a dose schedule before an avatar may be given it. The vignette works M4's
  # example through "with a floor of two".
  expect_equal(defaults$min_pattern_share, 2L)

  # The three flags the vignette names by their default state.
  expect_true(defaults$screen)
  expect_true(defaults$coarsen_time)
  expect_true(defaults$preserve_strata_balance)

  # Step 8's table marks "drop" as the default, which here means first in the
  # `match.arg()` vector rather than a scalar default.
  expect_equal(eval(defaults$on_donor_shortfall),
               c("drop", "noise", "error"))
})

# These four are named constants or formals of internals, so the vignette's
# numbers can be read straight off them.
test_that("the internal constants the vignette states hold their values", {
  # Step 3: "at most twelve distinct levels" is `ordinal`, more is `integer`.
  expect_equal(synpmx:::.endpoint_ordinal_max_levels, 12L)
  # Step 6: "Strata holding fewer than three source patients are deliberately
  # not balanced".
  expect_equal(synpmx:::.strata_balance_floor, 3L)
  # Step 4: "a common grid of at most fifteen times".
  expect_equal(formals(synpmx:::.common_grid)$max_points, 15L)
  # Step 12: inference "requires the dose-to-covariate ratio to collapse onto a
  # handful of levels within 2%".
  expect_equal(formals(synpmx:::.detect_dose_basis)$tolerance, 0.02)
  # Step 6, M2: "exceeds twice the cohort's 90th percentile". The multiplier is
  # a formal; the quantile is not, and is listed in the gap above.
  expect_equal(formals(synpmx:::.structural_outlier_anchors)$mult, 2)
})
