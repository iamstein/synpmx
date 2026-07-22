load_nlmixr2_dataset <- function(name) {
  environment <- new.env(parent = emptyenv())
  utils::data(list = name, package = "nlmixr2data", envir = environment)
  get(name, envir = environment, inherits = FALSE)
}

test_that("theo_md runs end to end", {
  skip_if_not_installed("nlmixr2data")
  source <- load_nlmixr2_dataset("theo_md")
  roles <- pmx_roles(
    "ID", "TIME", "DV", "AMT", "EVID", "CMT", covariates = "WT"
  )
  mock <- mock_pmx(source, roles, n_subjects = 4, seed = 42)
  comparison <- compare_pmx(source, mock, roles)

  expect_true(validate_pmx(mock, roles)$valid)
  expect_equal(length(unique(mock$ID)), 4L)
  expect_true(all(mock$DV[mock$EVID == 0] >= 0))
  expect_true(all(vapply(split(mock$EVID, mock$ID),
                         function(x) sum(x != 0) == 7L, logical(1))))
  expect_true(comparison$validation$mock$valid)
})

test_that("warfarin runs endpoints separately end to end", {
  skip_if_not_installed("nlmixr2data")
  source <- load_nlmixr2_dataset("warfarin")
  roles <- pmx_roles(
    "id", "time", "dv", "amt", "evid", dvid = "dvid",
    covariates = c("wt", "age", "sex")
  )
  mock <- mock_pmx(source, roles, n_subjects = 4, seed = 42)
  comparison <- compare_pmx(source, mock, roles)

  expect_true(validate_pmx(mock, roles)$valid)
  expect_identical(levels(mock$dvid), levels(source$dvid))
  expect_identical(levels(mock$sex), levels(source$sex))
  expect_setequal(unique(as.character(mock$dvid[mock$evid == 0])),
                  c("cp", "pca"))
  expect_true(all(mock$dv[mock$evid == 0] >= 0))
  expect_true(comparison$validation$mock$valid)
})

test_that("wbcSim retains infusion records and delayed trajectories end to end", {
  skip_if_not_installed("nlmixr2data")
  source <- load_nlmixr2_dataset("wbcSim")
  roles <- pmx_roles(
    "ID", "TIME", "DV", "AMT", "EVID", "CMT", rate = "RATE",
    covariates = c("V2I", "V1I", "CLI")
  )
  mock <- mock_pmx(source, roles, n_subjects = 4, seed = 42)
  comparison <- compare_pmx(source, mock, roles)
  event <- mock$EVID != 0
  observation <- mock$EVID == 0

  expect_true(validate_pmx(mock, roles)$valid)
  expect_true(all(mock$AMT[event] == mock$RATE[event]))
  expect_true(any(mock$AMT[event] > 0))
  expect_true(any(mock$AMT[event] < 0))
  expect_true(all(mock$CMT[event] == 1L))
  expect_true(all(mock$CMT[observation] == 3L))
  expect_true(all(is.finite(mock$DV[observation])))
  expect_true(all(mock$DV[observation] >= 0))
  expect_true(comparison$validation$mock$valid)
})
