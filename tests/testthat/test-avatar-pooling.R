# REV-025: AVATAR must blend a floor of `k` real patients into every synthetic
# subject, borrowing the nearest donors across dose/schedule groups when a
# subject's own group is too small, and alerting loudly only when the source is
# smaller than the floor.

mk_profiles <- function(sig, coord, routes = rep("oral", length(sig))) {
  list(
    subjects = seq_along(sig),
    coordinates = matrix(coord, ncol = 1L,
                         dimnames = list(as.character(seq_along(sig)), "pc1")),
    signatures = sig,
    routes = routes
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

test_that("route is never crossed, even when it starves the donor pool", {
  # Subject 1 is the only IV subject; subjects 2-6 are oral and much nearer in
  # profile space. Borrowing them would reach k, and it must not happen.
  p <- mk_profiles(
    sig = c("A", "B", "B", "B", "B", "B"),
    coord = c(0, 1, 2, 3, 4, 5),
    routes = c("iv", rep("oral", 5))
  )
  w <- synpmx:::.warning_collector()

  donors <- synpmx:::.select_donors(1L, p, k = 5L, w)

  # No legal donor exists, so the anchor stands alone rather than borrowing
  # across routes. The caller is what drops such anchors before generation.
  expect_equal(donors$indices, 1L)
  expect_equal(donors$weights, 1)
})

test_that("borrowing stays inside the anchor's route", {
  # Two IV subjects (1, 2) sit far apart in profile space; oral subjects 3-6 are
  # nearer to anchor 1 than subject 2 is, and still must not be borrowed.
  p <- mk_profiles(
    sig = c("A", "B", "C", "C", "C", "C"),
    coord = c(0, 20, 1, 2, 3, 4),
    routes = c("iv", "iv", rep("oral", 4))
  )
  w <- synpmx:::.warning_collector()

  donors <- synpmx:::.select_donors(1L, p, k = 5L, w)

  expect_equal(donors$indices, 2L)   # the one route-compatible subject
  expect_false(any(c(3L, 4L, 5L, 6L) %in% donors$indices))
})

test_that("a route arm below the floor is dropped from the anchor pool", {
  one <- function(id, dose, rate) data.frame(
    ID = id, TIME = c(0, 0.5, 1, 2, 4),
    DV = c(0, dose * 0.01, dose * 0.02, dose * 0.015, dose * 0.008),
    AMT = c(dose, 0, 0, 0, 0), RATE = c(rate, 0, 0, 0, 0),
    EVID = c(1L, 0L, 0L, 0L, 0L), CMT = c(1L, 2L, 2L, 2L, 2L), WT = 70
  )
  # Six oral subjects meet the floor; two infusion subjects cannot.
  src <- do.call(rbind, c(
    lapply(1:6, one, dose = 100, rate = 0),
    lapply(7:8, one, dose = 100, rate = 50)
  ))
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", cmt = "CMT", rate = "RATE",
                     covariates = "WT")

  msgs <- testthat::capture_messages(
    syn <- suppressWarnings(
      synpmx_avatar(src, roles, n_subjects = 20L, seed = 5)
    )
  )

  expect_true(any(grepl("SYNPMX ALERT", msgs)))
  expect_true(any(grepl("never blended across routes", msgs)))
  # No avatar carries an infusion, because no infusion subject could anchor one.
  expect_true(all(syn$RATE == 0))
})

test_that("on_donor_shortfall keeps, drops, or refuses the starved arm", {
  one <- function(id, rate) data.frame(
    ID = id, TIME = c(0, 0.5, 1, 2, 4),
    DV = c(0, 1, 2, 1.5, 0.8) + id / 100,
    AMT = c(100, 0, 0, 0, 0), RATE = c(rate, 0, 0, 0, 0),
    EVID = c(1L, 0L, 0L, 0L, 0L), CMT = c(1L, 2L, 2L, 2L, 2L), WT = 70
  )
  src <- do.call(rbind, c(lapply(1:6, one, rate = 0),
                          lapply(7:8, one, rate = 50)))
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", cmt = "CMT", rate = "RATE",
                     covariates = "WT")

  # "noise" keeps the starved arm, so infusion avatars appear.
  kept_msgs <- testthat::capture_messages(
    kept <- suppressWarnings(synpmx_avatar(
      src, roles, n_subjects = 20L, seed = 5, on_donor_shortfall = "noise"
    ))
  )
  expect_true(any(kept$RATE != 0))
  expect_true(any(grepl("individually identifying", kept_msgs)))

  # "error" refuses, and names both alternatives so the message is actionable.
  expect_error(
    synpmx_avatar(src, roles, n_subjects = 20L, seed = 5,
                  on_donor_shortfall = "error"),
    "on_donor_shortfall"
  )
  expect_error(
    synpmx_avatar(src, roles, n_subjects = 20L, seed = 5,
                  on_donor_shortfall = "error"),
    "\"noise\""
  )

  # The default drops, and says how to keep them instead.
  dropped_msgs <- testthat::capture_messages(
    dropped <- suppressWarnings(
      synpmx_avatar(src, roles, n_subjects = 20L, seed = 5)
    )
  )
  expect_true(all(dropped$RATE == 0))
  expect_true(any(grepl("on_donor_shortfall = \\\"noise\\\"", dropped_msgs)))

  expect_equal(attr(dropped, "pmx_settings")$on_donor_shortfall, "drop")
  expect_error(
    synpmx_avatar(src, roles, on_donor_shortfall = "sometimes"),
    "should be one of"
  )
})

test_that("every arm below the floor generates anyway rather than emptying", {
  # No arm can reach the floor, so "drop" has nothing left to keep; generation
  # proceeds as if "noise" rather than failing.
  one <- function(id, rate) data.frame(
    ID = id, TIME = c(0, 1, 2), DV = c(0, 2, 1) + id / 100,
    AMT = c(100, 0, 0), RATE = c(rate, 0, 0), EVID = c(1L, 0L, 0L),
    CMT = c(1L, 2L, 2L), WT = 70
  )
  src <- do.call(rbind, c(lapply(1:2, one, rate = 0),
                          lapply(3:4, one, rate = 50)))
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", cmt = "CMT", rate = "RATE",
                     covariates = "WT")

  msgs <- testthat::capture_messages(
    syn <- suppressWarnings(synpmx_avatar(src, roles, n_subjects = 6L, seed = 2))
  )

  expect_true(any(grepl("nothing to generate", msgs)))
  expect_true(validate_pmx(syn, roles)$valid)
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
