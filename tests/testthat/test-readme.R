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
  expect_equal(head4$TIME, c(0.0000000, 0.0000000, 0.2816667, 0.5416667),
               tolerance = 1e-6)
  expect_equal(head4$DV,
               c(0.000000, 0.000000, 2.861275, 3.600730),
               tolerance = 1e-6)
  expect_equal(head4$AMT, c(267.84, 0.00, 0.00, 0.00), tolerance = 1e-6)
  expect_equal(head4$EVID, c(101L, 0L, 0L, 0L))
  expect_equal(head4$CMT, c(1L, 2L, 2L, 2L))
  expect_equal(head4$WT, rep(69.80615, 4), tolerance = 1e-5)
})

# The second and third blocks, "Declaring every column" and "The masking
# options, and their defaults". These exist to show every role and every masking
# argument at once, so the test's job is to prove that the full declaration is
# still accepted and still produces the printed table -- an argument renamed or
# a role dropped shows up here rather than in a reader's session.
test_that("the README full-declaration example runs as shown", {
  study <- pmx_simulated_fixture(24)
  study$YTYPE <- ifelse(study$DVID == "cp", 1L, 2L)
  study$NAME <- as.character(study$DVID)
  study$DVID <- NULL
  study$TRTN <- ifelse(study$ID %% 2L == 1L, 1L, 2L)
  study$TRT <- ifelse(study$TRTN == 1L, "100 mg QD", "200 mg QD")
  study$AMT[study$EVID != 0] <- 100 * study$TRTN[study$EVID != 0]
  study$STUDYID <- "EXAMPLE-001"
  bloq <- study$EVID == 0 & study$NAME == "cp" & study$DV < 1.2
  study$DV[bloq] <- 1.2
  study$CENS[bloq] <- 1L
  study$LIMIT <- ifelse(bloq, 0, NA_real_)

  roles <- pmx_roles(
    id                 = "ID",
    time               = "TIME",
    dv                 = "DV",
    evid               = "EVID",
    amt                = "AMT",
    rate               = "RATE",
    cmt                = "CMT",
    dvid               = c("YTYPE", "NAME"),
    mdv                = "MDV",
    nominal_time       = "NTIME",
    tad                = "TAD",
    occasion           = "OCC",
    cens               = "CENS",
    limit              = "LIMIT",
    covariates         = c("WT", "AGE", "SEX"),
    subject_properties = c("TRT", "TRTN"),
    keep               = "STUDYID"
  )

  #> [1] TRUE
  expect_true(validate_pmx(study, roles)$valid)

  # The README shows no warning and no dropped-column message under this call,
  # which is only true because every column is declared and every arm clears the
  # donor floor. If that stops holding, the README is showing a clean run that
  # readers will not get.
  expect_silent(
    synthetic <- synpmx_avatar(
      study, roles,
      n_subjects         = NULL,
      seed               = 2026,
      k                  = 5,
      max_donor_weight   = 0.50,
      on_donor_shortfall = "drop",
      screen             = TRUE,
      coarsen_time       = TRUE,
      min_pattern_share  = 2L,
      subject_noise_sd   = 0.15,
      residual_noise_sd  = 0.05,
      residual_phi       = 0.6,
      time_jitter        = 0,
      pca_variance       = 0.90
    )
  )

  #> [1] TRUE
  expect_true(validate_pmx(synthetic, roles)$valid)

  #> [1] 384  21
  expect_equal(dim(synthetic), c(384L, 21L))

  # head(synthetic[, ...], 6), as printed in README.md.
  shown <- head(
    synthetic[, c("ID", "TIME", "NTIME", "OCC", "NAME", "DV", "CENS", "TRT")], 6
  )
  expect_equal(as.character(shown$ID), rep("25", 6))
  expect_equal(shown$TIME, c(0, 0, 0.25, 1, 2, 4), tolerance = 1e-6)
  expect_equal(shown$NTIME, c(0, 0, 0.25, 1, 2, 4), tolerance = 1e-6)
  expect_equal(shown$OCC, rep(1L, 6))
  expect_equal(shown$NAME, c("cp", "pd", "cp", "cp", "cp", "pd"))
  expect_equal(shown$DV,
               c(0.000000, 90.823217, 1.200000, 8.249096, 3.766242, 78.636255),
               tolerance = 1e-6)
  expect_equal(shown$CENS, c(0L, 0L, 1L, 0L, 0L, 0L))
  expect_equal(shown$TRT, rep("100 mg QD", 6))
})
