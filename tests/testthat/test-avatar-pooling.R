# REV-025: AVATAR must blend a floor of `k` real patients into every synthetic
# subject, borrowing the nearest donors across dose/schedule groups when a
# subject's own group is too small, and alerting loudly only when the source is
# smaller than the floor.

mk_profiles <- function(sig, coord) {
  list(
    subjects = seq_along(sig),
    coordinates = matrix(coord, ncol = 1L,
                         dimnames = list(as.character(seq_along(sig)), "pc1")),
    signatures = sig
  )
}

test_that("a singleton group borrows the nearest k donors from other groups", {
  # Subject 6 (signature B, coord 10) is alone in its group.
  p <- mk_profiles(c("A", "A", "A", "A", "A", "B"), c(0, 1, 2, 3, 4, 10))
  w <- synpmx:::.warning_collector()

  donors <- synpmx:::.select_donors(6L, p, k = 5L, w)

  expect_length(donors$indices, 5L)
  expect_false(6L %in% donors$indices)          # never its own donor
  # nearest to coord 10 among {1..5} (coords 0..4) is 5, 4, 3, 2, 1
  expect_equal(donors$indices, c(5L, 4L, 3L, 2L, 1L))
  expect_true(any(grepl("borrowed", w$messages)))
})

test_that("a well-populated group uses same-schedule donors, nearest first", {
  p <- mk_profiles(c("A", "A", "A", "A", "A", "B"), c(0, 1, 2, 3, 4, 10))
  w <- synpmx:::.warning_collector()

  donors <- synpmx:::.select_donors(1L, p, k = 3L, w)

  # subject 1 (coord 0); same group {2,3,4,5} coords 1,2,3,4 -> nearest 2,3,4
  expect_equal(donors$indices, c(2L, 3L, 4L))
  expect_length(w$messages, 0L)                 # nothing borrowed
})

test_that("an undersized group is topped up from the nearest other subjects", {
  # group A = {1,2}; anchor 1, k = 4 -> take 2, borrow 2 nearest from B.
  p <- mk_profiles(c("A", "A", "B", "B", "B", "B"), c(0, 1, 5, 6, 7, 20))
  w <- synpmx:::.warning_collector()

  donors <- synpmx:::.select_donors(1L, p, k = 4L, w)

  # same-group {2}; then nearest others to coord 0 among coords 5,6,7,20 -> 3,4,5
  expect_equal(donors$indices, c(2L, 3L, 4L, 5L))
  expect_true(any(grepl("borrowed", w$messages)))
})

test_that("a source smaller than the floor triggers a loud alert", {
  # 3 subjects, default k = 5 -> at most 2 donors, below the floor.
  msgs <- testthat::capture_messages(
    suppressWarnings(synpmx_avatar(private_fixture(3L), private_roles(), seed = 1))
  )
  expect_true(any(grepl("SYNPMX ALERT", msgs)))
})

test_that("a unique-dose subject no longer collapses to a sole donor", {
  one <- function(id, dose) data.frame(
    ID = id, TIME = c(0, 0.5, 1, 2, 4),
    DV = c(0, dose * 0.01, dose * 0.02, dose * 0.015, dose * 0.008),
    AMT = c(dose, 0, 0, 0, 0), EVID = c(1L, 0L, 0L, 0L, 0L),
    CMT = c(1L, 2L, 2L, 2L, 2L), WT = 70
  )
  # Five subjects at dose 100 and one unique subject at dose 999 (own signature).
  src <- do.call(rbind, c(lapply(1:5, one, dose = 100), list(one(6L, 999))))
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", cmt = "CMT", covariates = "WT")

  msgs <- testthat::capture_messages(
    syn <- suppressWarnings(synpmx_avatar(src, roles, n_subjects = 12L, seed = 3))
  )
  # 6 subjects means 5 donors are available, which meets the floor: no alert.
  expect_false(any(grepl("SYNPMX ALERT", msgs)))
  expect_true(validate_pmx(syn, roles)$valid)
})
