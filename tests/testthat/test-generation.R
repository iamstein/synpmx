test_that("generation is reproducible post-processing and leaves accounting unchanged", {
  model <- fit_public_fixture()
  accounting <- model$privacy$accounting
  set.seed(987)
  state <- .Random.seed
  first <- generate_pmx(model, 5, seed = 123)
  expect_identical(.Random.seed, state)
  second <- generate_pmx(model, 5, seed = 123)
  different <- generate_pmx(model, 5, seed = 124)
  expect_identical(first, second)
  expect_false(identical(first, different))
  expect_identical(model$privacy$accounting, accounting)
})

test_that("generation defaults to the released fitted cohort size", {
  source <- private_fixture(12L)
  model <- fit_public_fixture(source)
  synthetic <- generate_pmx(model, seed = 125)
  expect_equal(model$population$private_subject_count, 12)
  expect_equal(length(unique(synthetic$ID)), length(unique(source$ID)))
})

test_that("timing-cell selection retains high-probability late cells", {
  expect_identical(.select_timing_cells(c(0, 0, 1, 1), 2L), 3:4)
})

test_that("schema, classes, factors, new IDs, and covariates are coherent", {
  source <- private_fixture()
  synthetic <- generate_pmx(fit_public_fixture(source), 6, seed = 22)
  expect_identical(names(synthetic), names(source))
  expect_identical(vapply(synthetic, class, character(1)),
                   vapply(source, class, character(1)))
  expect_identical(levels(synthetic$DVID), levels(source$DVID))
  expect_identical(levels(synthetic$SEX), levels(source$SEX))
  expect_type(synthetic$ID, "integer")
  expect_false(any(unique(synthetic$ID) %in% source$ID))
  for (column in private_roles()$covariates) {
    expect_true(all(vapply(split(synthetic[[column]], synthetic$ID),
                           function(x) length(unique(x)) == 1L, logical(1))))
  }
  expect_true(validate_pmx(synthetic, private_roles(), private_endpoints())$valid)
})

test_that("subject properties remain coherent with conditioned regimens", {
  source <- private_fixture(12L)
  source$ARM <- as.integer(ifelse(source$ID %% 2L, 1L, 2L))
  event <- source$EVID != 0 & source$AMT > 0
  source$AMT[event] <- ifelse(source$ARM[event] == 1L, 50, 150)
  source$DOSE <- ave(
    source$AMT, interaction(source$ID, source$OCC, drop = TRUE),
    FUN = function(value) max(value)
  )
  roles <- pmx_roles(
    id = "ID", time = "TIME", nominal_time = "NTIME", tad = "TAD",
    occasion = "OCC", dv = "DV", amt = "AMT", evid = "EVID",
    cmt = "CMT", dvid = "DVID", mdv = "MDV", rate = "RATE",
    cens = "CENS", limit = "LIMIT", assigned_dose = "DOSE",
    covariates = c("WT", "AGE", "SEX"), subject_properties = "ARM"
  )
  design <- pmx_public_design(
    pmx_schema(source), dose_evid = 1, dose_cmt = 1,
    endpoint_grids = lapply(private_endpoints(), `[[`, "grid"),
    endpoint_cmt = list(cp = 2, pd = 3),
    category_levels = list(ARM = c(1, 2)), time_jitter_sd = .01
  )
  model <- suppressWarnings(fit_private_pmx(
    source, roles, private_endpoints(), 5, 0,
    pmx_bounds(
      c(0, 24), list(cp = c(0, 20), pd = c(0, 120)),
      amt = c(0, 200), rate = c(-200, 200),
      covariates = list(WT = c(40, 120), AGE = c(18, 90)),
      limit = list(cp = c(0, 20), pd = c(0, 120))
    ),
    design, private_limits(), private_budget(),
    backend = "public", public_source = TRUE
  ))
  summary <- subject_property_summary(model)
  expect_equal(summary$ARM, 1:2)
  expect_equal(summary$probability, c(.5, .5))
  expect_equal(summary$dose_amount, c(50, 150), tolerance = 1e-8)
  event_entry <- model$privacy$accounting$entries[
    model$privacy$accounting$entries$query == "event_and_regimen", ]
  expect_equal(event_entry$sensitivity, 10)
  expect_equal(event_entry$dimensions, 20)

  synthetic <- generate_pmx(model, 60, seed = 2201)
  expect_true(validate_pmx(synthetic, roles, private_endpoints())$valid)
  property_by_subject <- vapply(
    split(synthetic$ARM, synthetic$ID), function(value) unique(value)[1L], integer(1)
  )
  amount_by_subject <- vapply(split(synthetic, synthetic$ID), function(subject) {
    unique(subject$AMT[subject$EVID != 0 & subject$AMT > 0])[1L]
  }, numeric(1))
  expect_equal(
    amount_by_subject,
    ifelse(property_by_subject == 1L, 50, 150), tolerance = 1e-8
  )
  expect_true(all(synthetic$DOSE == ave(
    synthetic$DOSE, interaction(synthetic$ID, synthetic$OCC, drop = TRUE),
    FUN = function(value) value[1L]
  )))
  positive_event <- synthetic$EVID != 0 & synthetic$AMT > 0
  expect_equal(synthetic$DOSE[positive_event], synthetic$AMT[positive_event])
})

test_that("dose-relative PK repeats while study-time PD remains global", {
  synthetic <- generate_pmx(fit_public_fixture(), 8, seed = 51)
  cp <- synthetic[synthetic$EVID == 0 & synthetic$DVID == "cp", ]
  by_occasion <- split(cp, interaction(cp$ID, cp$OCC, drop = TRUE))
  expect_true(all(vapply(by_occasion, function(x) {
    peak <- which.max(x$DV)
    peak > 1L && peak < nrow(x) && x$DV[peak] > x$DV[1L] &&
      x$DV[peak] > x$DV[nrow(x)]
  }, logical(1))))

  pd <- synthetic[synthetic$EVID == 0 & synthetic$DVID == "pd", ]
  early <- mean(pd$DV[abs(pd$NTIME - 0) < 1e-8])
  mid <- mean(pd$DV[abs(pd$NTIME - 12) < 1e-8])
  expect_lt(mid, early)
  expect_equal(length(unique(pd$NTIME)), length(private_endpoints()$pd$grid))
})

test_that("TAD, occasions, tied ordering, doses, and times are coherent", {
  model <- fit_public_fixture()
  synthetic <- generate_pmx(model, 4, seed = 91)
  expected_observations <- as.integer(round(
    model$population$event$observation_count
  ))
  expect_true(all(vapply(split(synthetic$EVID, synthetic$ID), function(x) {
    sum(x == 0) == expected_observations
  }, logical(1))))
  for (id in unique(synthetic$ID)) {
    subject <- synthetic[synthetic$ID == id, ]
    expect_true(all(diff(subject$TIME) >= 0))
    expect_true(all(subject$TAD >= 0))
    dose_time <- subject$TIME[subject$EVID != 0 & subject$AMT > 0]
    observations <- subject$EVID == 0
    occasion <- subject$OCC[observations]
    expected <- subject$TIME[observations] - dose_time[occasion]
    expect_equal(subject$TAD[observations], pmax(0, expected), tolerance = 1e-8)
    at_zero <- which(subject$TIME == 0)
    if (length(at_zero) > 1L) expect_true(subject$EVID[at_zero[1L]] != 0)
  }
  observations <- synthetic$EVID == 0
  tied_blocks <- split(
    synthetic$TIME[observations],
    interaction(synthetic$ID[observations], synthetic$NTIME[observations], drop = TRUE)
  )
  expect_true(all(vapply(tied_blocks, function(x) length(unique(x)) == 1L,
                         logical(1))))
  expect_true(any(abs(synthetic$TIME[observations] - synthetic$NTIME[observations]) >
                    1e-10))
})

test_that("occasion and hybrid alignments execute separately", {
  endpoints <- private_endpoints()
  endpoints$cp$alignment <- "occasion"
  endpoints$pd$alignment <- "hybrid"
  source <- private_fixture()
  model <- suppressWarnings(fit_private_pmx(
    source, private_roles(), endpoints, 5, 0, private_bounds(),
    private_design(source), private_limits(), private_budget(),
    backend = "public", public_source = TRUE
  ))
  synthetic <- generate_pmx(model, 3, 11)
  expect_true(validate_pmx(synthetic, private_roles(), endpoints)$valid)
  expect_setequal(unique(as.character(synthetic$DVID[synthetic$EVID == 0])),
                  c("cp", "pd"))
  cp <- synthetic[synthetic$EVID == 0 & synthetic$DVID == "cp", ]
  expect_true(all(vapply(split(cp, interaction(cp$ID, cp$OCC)), function(x) {
    peak <- which.max(x$DV)
    peak > 1L && peak < nrow(x)
  }, logical(1))))
  expect_true(length(model$population$trajectories$pd$local_grid) > 1L)
  expect_equal(
    sort(unique(synthetic$NTIME[synthetic$EVID == 0 & synthetic$DVID == "pd"])),
    endpoints$pd$grid
  )
})
