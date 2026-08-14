# C3 reports "how many arms kept their endpoints" as one count. The question a
# user then has is which arm lost which endpoint, and that needs both sides, one
# row per arm and endpoint.
#
# Counts are patients, not rows. A truncated sampling schedule moves the row
# count of an arm that kept every endpoint and every patient, so a row count
# here would report noise and miss the loss.

endpoints_source <- function(n = 40L, arms = c("A", "B")) {
  source <- pmx_simulated_fixture(n)
  ids <- unique(as.character(source$ID))
  source$ARM <- rep(arms, length.out = length(ids))[
    match(as.character(source$ID), ids)
  ]
  source
}

endpoints_roles <- function(...) {
  pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
            cmt = "CMT", dvid = "DVID", covariates = "WT", ...)
}

endpoints_synthetic <- function(source, roles, ...) {
  suppressWarnings(suppressMessages(synpmx_avatar(source, roles, ...)))
}

test_that("every arm and endpoint pair is reported from both sides", {
  source <- endpoints_source()
  roles <- endpoints_roles(strata = "ARM")
  synthetic <- endpoints_synthetic(source, roles, seed = 3)

  held <- compare_pmx_strata_endpoints(source, synthetic, roles)

  expect_s3_class(held, "pmx_strata_endpoints")
  expect_identical(sort(unique(held$level)), c("A", "B"))
  expect_identical(nrow(held),
                   length(unique(held$level)) * length(unique(held$endpoint)))
  expect_true(all(held$source_patients > 0L))
  expect_true(all(held$synthetic_patients > 0L))
  expect_identical(attr(held, "release_status"), "restricted_not_releasable")
})

test_that("patients are counted, not rows", {
  source <- endpoints_source(n = 20L)
  roles <- endpoints_roles(strata = "ARM")
  synthetic <- endpoints_synthetic(source, roles, seed = 5)

  held <- compare_pmx_strata_endpoints(source, synthetic, roles)

  # One arm holds half the cohort, and no cell may exceed that however many
  # samples each patient contributed.
  expect_true(all(held$source_patients <= 10L))
  observed <- source[source$EVID == 0, ]
  expect_gt(nrow(observed), sum(held$source_patients))
})

test_that("an endpoint an arm lost is the row that stands out", {
  source <- endpoints_source()
  roles <- endpoints_roles(strata = "ARM")
  synthetic <- endpoints_synthetic(source, roles, seed = 4)

  endpoint <- sort(unique(as.character(synthetic$DVID)))[[1L]]
  lost <- synthetic[!(synthetic$ARM == "A" & synthetic$DVID == endpoint &
                        synthetic$EVID == 0), ]

  held <- compare_pmx_strata_endpoints(source, lost, roles)
  gone <- held[held$source_patients > 0L & held$synthetic_patients == 0L, ]

  expect_identical(nrow(gone), 1L)
  expect_identical(gone$level, "A")
  expect_identical(gone$endpoint, endpoint)
  expect_output(print(held), "LOST BY THIS ARM")
})

test_that("no declared strata gives an empty table rather than an error", {
  source <- endpoints_source()
  roles <- endpoints_roles()
  synthetic <- endpoints_synthetic(source, roles, seed = 3)

  held <- compare_pmx_strata_endpoints(source, synthetic, roles)

  expect_identical(nrow(held), 0L)
  expect_output(print(held), "nothing to compare")
})
