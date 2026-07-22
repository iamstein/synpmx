test_that("complete repeated-dose templates and tied row order are coherent", {
  source <- pmx_fixture()
  roles <- fixture_roles()
  mock <- mock_pmx(source, roles, n_subjects = 6, seed = 19)
  controls <- c("TIME", "AMT", "RATE", "EVID", "CMT", "MDV")
  source_controls <- lapply(split(source, source$ID), function(subject) {
    unname(subject[, controls, drop = FALSE])
  })

  for (id in unique(mock$ID)) {
    subject <- mock[mock$ID == id, ]
    expect_true(all(diff(subject$TIME) >= 0))
    expect_identical(subject$EVID[subject$TIME == 0], c(1L, 0L))
    expect_identical(subject$EVID[subject$TIME == 12], c(1L, 0L))
    expect_true(all(subject$DV[subject$EVID != 0] == 0))
    expect_true(all(subject$AMT[subject$EVID == 0] == 0))
    expect_true(all(subject$RATE == 0))
    expect_true(all(subject$CMT[subject$EVID != 0] == 1L))
    expect_true(all(subject$MDV[subject$EVID != 0] == 1L))
    candidate <- unname(subject[, controls, drop = FALSE])
    matches_source <- vapply(source_controls, function(source_control) {
      isTRUE(all.equal(candidate, source_control, check.attributes = FALSE))
    }, logical(1))
    expect_true(any(matches_source))
  }
})

test_that("observation DVs are finite, constrained, and not copied trajectories", {
  source <- pmx_fixture()
  roles <- fixture_roles()
  mock <- mock_pmx(source, roles, n_subjects = 6, seed = 29)
  mock_observed <- mock$EVID == 0 & mock$MDV == 0
  expect_true(all(is.finite(mock$DV[mock_observed])))
  expect_true(all(mock$DV[mock_observed] >= 0))

  source_trajectories <- split(source$DV[source$EVID == 0],
                               source$ID[source$EVID == 0])
  mock_trajectories <- split(mock$DV[mock$EVID == 0],
                             mock$ID[mock$EVID == 0])
  copied <- vapply(mock_trajectories, function(trajectory) {
    any(vapply(source_trajectories, identical, logical(1), y = trajectory))
  }, logical(1))
  expect_false(any(copied))
})

test_that("multiple endpoints retain factor semantics and stay separated", {
  source <- multidvid_fixture()
  roles <- multidvid_roles()
  mock <- mock_pmx(source, roles, n_subjects = 5, seed = 51)
  observed <- mock$evid == 0

  expect_identical(levels(mock$dvid), c("cp", "pd"))
  expect_setequal(unique(as.character(mock$dvid[observed])), c("cp", "pd"))
  expect_true(all(is.finite(mock$dv[observed])))
  expect_true(mean(mock$dv[observed & mock$dvid == "pd"]) >
                mean(mock$dv[observed & mock$dvid == "cp"]))
  expect_true(validate_pmx(mock, roles)$valid)
})

test_that("missing observation patterns remain attached to templates", {
  source <- pmx_fixture(missing_dv = TRUE)
  roles <- fixture_roles()
  mock <- mock_pmx(source, roles, n_subjects = 12, seed = 4)
  missing <- is.na(mock$DV)
  expect_true(all(mock$MDV[missing] == 1L))
  expect_true(all(mock$DV[mock$MDV == 0] >= 0))
  expect_true(validate_pmx(mock, roles)$valid)
})

test_that("coherent time jitter retains ties and ordering", {
  source <- pmx_fixture()
  roles <- fixture_roles()
  mock <- mock_pmx(source, roles, n_subjects = 3, seed = 15,
                   time_jitter = 0.02)
  for (id in unique(mock$ID)) {
    subject <- mock[mock$ID == id, ]
    expect_true(all(diff(subject$TIME) >= 0))
    tied_groups <- table(subject$TIME)
    expect_equal(sum(tied_groups == 2L), 2L)
  }
})
