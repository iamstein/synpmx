# README.md is hand-maintained: there is no README.Rmd knitting it, so nothing
# re-runs its example when behavior changes. This test is what keeps the
# example honest. It runs exactly the code shown under "A first synthetic
# dataset" and asserts exactly the output printed beneath it.
#
# If this fails, the README is telling readers something the package no longer
# does. Fix the package or update the README — do not relax the test.

test_that("the README example produces the output the README shows", {
  skip_if_not_installed("nlmixr2data")

  data("theo_md", package = "nlmixr2data", envir = environment())

  roles <- pmx_roles(
    id = "ID", time = "TIME", dv = "DV", amt = "AMT",
    evid = "EVID", cmt = "CMT", covariates = "WT"
  )

  synthetic <- suppressWarnings(synpmx_avatar(theo_md, roles, seed = 101))

  #> [1] TRUE
  expect_true(validate_pmx(synthetic, roles)$valid)

  # head(synthetic, 4), as printed in README.md.
  expect_equal(names(synthetic),
               c("ID", "TIME", "DV", "AMT", "EVID", "CMT", "WT"))

  head4 <- head(synthetic, 4)
  expect_equal(as.character(head4$ID), rep("13", 4))
  expect_equal(head4$TIME, c(0.00, 0.00, 0.30, 0.63), tolerance = 1e-8)
  expect_equal(head4$DV,
               c(0.00000000, 0.02403989, 10.00272146, 11.89216329),
               tolerance = 1e-6)
  expect_equal(head4$AMT, c(267.84, 0.00, 0.00, 0.00), tolerance = 1e-6)
  expect_equal(head4$EVID, c(101L, 0L, 0L, 0L))
  expect_equal(head4$CMT, c(1L, 2L, 2L, 2L))
  expect_equal(head4$WT, rep(85.25496, 4), tolerance = 1e-5)
})
