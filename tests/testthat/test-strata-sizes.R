# C3 reports "how many arms kept their size" as one count. The question a user
# actually has is which arm moved and by how much, and that needs both sides.

sizes_source <- function(n = 60L, arms = c("A", "B", "C")) {
  source <- pmx_simulated_fixture(n)
  ids <- unique(as.character(source$ID))
  source$ARM <- rep(arms, length.out = length(ids))[
    match(as.character(source$ID), ids)
  ]
  source
}

sizes_roles <- function(...) {
  pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
            cmt = "CMT", dvid = "DVID", covariates = "WT", ...)
}

sizes_synthetic <- function(source, roles, ...) {
  suppressWarnings(suppressMessages(synpmx_avatar(source, roles, ...)))
}

test_that("both sides are reported for every stratum level", {
  source <- sizes_source()
  roles <- sizes_roles(strata = "ARM")
  synthetic <- sizes_synthetic(source, roles, seed = 3)

  sizes <- compare_pmx_strata_sizes(source, synthetic, roles)

  expect_identical(sort(sizes$level), c("A", "B", "C"))
  expect_identical(sizes$source_patients, c(20L, 20L, 20L))
  expect_identical(sum(sizes$synthetic_patients), 60L)
  # Every arm is well above the floor, so every arm is held exactly.
  expect_true(all(sizes$balanced))
  expect_identical(sizes$synthetic_patients, sizes$source_patients)
})

test_that("`expected` scales to the synthetic cohort rather than the source", {
  # Without this, asking for a third of the cohort reads as every arm drifting.
  source <- sizes_source()
  roles <- sizes_roles(strata = "ARM")
  synthetic <- sizes_synthetic(source, roles, n_subjects = 20, seed = 3)

  sizes <- compare_pmx_strata_sizes(source, synthetic, roles)

  expect_equal(sizes$expected, rep(20 / 3, 3), tolerance = 1e-8)
  # Largest-remainder rounding, so each arm lands on one side of it or the other.
  expect_false(any(synpmx:::.strata_sizes_drifted(sizes)))
})

test_that("a stratum under the floor is marked as free to vary", {
  source <- sizes_source(arms = c("A", "B", "C"))
  ids <- unique(as.character(source$ID))
  source$ARM[as.character(source$ID) %in% ids[1:2]] <- "TINY"
  roles <- sizes_roles(strata = "ARM")
  synthetic <- sizes_synthetic(source, roles, seed = 3)

  sizes <- compare_pmx_strata_sizes(source, synthetic, roles)

  expect_false(sizes$balanced[sizes$level == "TINY"])
  expect_true(all(sizes$balanced[sizes$level != "TINY"]))
})

test_that("joint cells are reported, and the marginals do not claim a floor", {
  # The floor applies to the joint cell, so a 20-patient arm can still be split
  # into cells that all fall below it -- which is how an arm moves a long way
  # without any single-column reading explaining it.
  source <- sizes_source()
  ids <- unique(as.character(source$ID))
  source$COHORT <- rep(paste0("c", 1:7), length.out = length(ids))[
    match(as.character(source$ID), ids)
  ]
  roles <- sizes_roles(strata = c("ARM", "COHORT"))
  synthetic <- sizes_synthetic(source, roles, seed = 3)

  sizes <- compare_pmx_strata_sizes(source, synthetic, roles)

  expect_true(all(c("ARM", "COHORT", "ARM x COHORT") %in% sizes$column))
  expect_true(all(is.na(sizes$balanced[sizes$column %in% c("ARM", "COHORT")])))
  expect_false(anyNA(sizes$balanced[sizes$column == "ARM x COHORT"]))
})

test_that("no strata means no rows rather than an error", {
  source <- sizes_source()
  roles <- sizes_roles()
  synthetic <- sizes_synthetic(source, roles, seed = 3)

  sizes <- compare_pmx_strata_sizes(source, synthetic, roles)

  expect_identical(nrow(sizes), 0L)
  expect_output(print(sizes), "nothing to size")
})

test_that("C3 sends the reader to the helper", {
  source <- sizes_source()
  roles <- sizes_roles(strata = "ARM")
  synthetic <- sizes_synthetic(source, roles, seed = 3)

  card <- synpmx_scorecard(source, synthetic, roles)

  expect_identical(card$explore[card$check == "C3"],
                   "compare_pmx_strata_sizes(source, synthetic, roles)")
})
