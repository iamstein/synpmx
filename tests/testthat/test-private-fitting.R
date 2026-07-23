test_that("privacy configuration fails closed", {
  source <- private_fixture()
  args <- list(
    data = source, roles = private_roles(), endpoints = private_endpoints(),
    epsilon = 5, delta = 0, bounds = private_bounds(),
    public_design = private_design(source),
    contribution_limits = private_limits(),
    budget_allocation = private_budget()
  )
  if (dp_backend_status()$available) {
    production <- suppressWarnings(do.call(fit_private_pmx, args))
    expect_true(production$privacy$formal_dp)
    expect_identical(production$privacy$backend$name, "OpenDP")
  } else {
    expect_error(do.call(fit_private_pmx, args), "backend.*unavailable|OpenDP")
  }
  args$backend <- "public"
  expect_error(do.call(fit_private_pmx, args), "public_source = TRUE")
  args$public_source <- TRUE
  args$epsilon <- 0
  expect_error(do.call(fit_private_pmx, args), "epsilon")
  args$epsilon <- 5
  args$delta <- 1e-6
  expect_error(do.call(fit_private_pmx, args), "delta_justification")
  expect_false("seed" %in% names(formals(fit_private_pmx)))
})

test_that("contributions are bounded before all private aggregate groups", {
  source <- private_fixture(2)
  extra <- source[source$ID == 1L, , drop = FALSE]
  extra$TIME <- extra$TIME + 24
  extra$NTIME <- extra$TIME
  extra$OCC <- extra$OCC + 2L
  source <- rbind(source, extra, extra)
  limits <- private_limits(max_rows = 20, max_cells = 6)
  endpoints <- synpmx:::.normalize_endpoints(
    private_endpoints(), private_roles(), private_bounds(),
    private_design(source), limits
  )
  bounded <- synpmx:::.bound_subject_contributions(
    source, private_roles(), endpoints, private_bounds(), limits
  )
  expect_true(all(vapply(bounded, function(x) nrow(x$data) <= 20L, logical(1))))
  expect_true(all(vapply(bounded, function(x) length(x$doses) <= 4L, logical(1))))
  expect_true(all(vapply(bounded, function(x) {
    all(vapply(x$observations, length, integer(1)) <= 16L)
  }, logical(1))))

  explicit <- private_fixture(1)
  explicit$OCC[3L] <- 99L
  explicit_bounded <- synpmx:::.bound_subject_contributions(
    explicit, private_roles(), endpoints, private_bounds(), limits
  )
  expect_false(3L %in% explicit_bounded[[1L]]$observations$cp)
})

test_that("dose-occasion sampling density is learned in the timing release", {
  model <- fit_public_fixture()
  expect_equal(
    model$population$timing$cp$occasion_presence_probability,
    c(1, 1, 0, 0)
  )
  expect_equal(
    model$population$timing$cp$occasion_observation_count,
    c(4, 4, 0, 0)
  )
  expect_length(
    model$population$timing$pd$occasion_observation_count,
    0L
  )
  expect_equal(
    length(model$population$timing$cp$grid_probability),
    length(private_endpoints()$cp$grid)
  )
  summary <- sampling_summary(model)
  expect_equal(summary$sampling_probability, c(1, 1))
  expect_equal(summary$observations_if_sampled, c(4, 4))
})

test_that("unoccupied trajectory cells do not create artificial troughs", {
  filled <- synpmx:::.fill_unoccupied_curve(
    value_unit = c(0.9, NA, 0.3), presence = c(1, 0, 1),
    grid = c(0, 1, 3)
  )
  expect_equal(filled, c(0.9, 0.7, 0.3))
})

test_that("regimen and sampling are inferred when no schedules are supplied", {
  source <- private_fixture()
  endpoints <- private_endpoints()
  endpoints <- lapply(endpoints, function(endpoint) {
    endpoint$grid <- NULL
    endpoint
  })
  design <- pmx_public_design(
    schema = pmx_schema(source), dose_evid = 1, dose_cmt = 1,
    endpoint_cmt = list(cp = 2, pd = 3), time_jitter_sd = 0.01
  )
  model <- suppressWarnings(fit_private_pmx(
    source, private_roles(), endpoints, 5, 0, private_bounds(), design,
    private_limits(), private_budget(), backend = "public", public_source = TRUE
  ))

  expect_null(model$public$design$dose_times)
  expect_null(model$public$design$dose_interval)
  expect_null(model$public$design$n_doses)
  expect_length(model$public$design$endpoint_grids, 0L)
  expect_length(model$public$design$endpoint_occasion_grids, 0L)
  expect_equal(model$population$event$n_doses, 2L)
  expect_equal(model$population$event$dose_interval, 12)
  expect_equal(model$population$event$dose_amount, 104.5)
  expect_true(model$public$endpoints$cp$grid_automatic)
  expect_true(model$public$endpoints$pd$grid_automatic)
  expect_true(all(diff(model$public$endpoints$cp$grid) > 0))

  synthetic <- generate_pmx(model, seed = 381)
  expect_equal(length(unique(synthetic$ID)), length(unique(source$ID)))
  expect_true(all(vapply(split(synthetic$EVID, synthetic$ID), function(x) {
    sum(x != 0) == 2L
  }, logical(1))))
  expect_true(validate_pmx(synthetic, private_roles(), endpoints)$valid)
})

test_that("rare dense sampling is not treated as sparse sampling for everyone", {
  source <- private_fixture()
  remove <- source$ID > 2L & source$EVID == 0 &
    source$DVID == "cp" & source$OCC == 2L
  model <- fit_public_fixture(source[!remove, , drop = FALSE])
  cp <- sampling_summary(model)
  second <- cp[cp$endpoint == "cp" & cp$occasion == 2L, ]
  expect_equal(second$sampling_probability, 0.25)
  expect_equal(second$observations_if_sampled, 4)
  expect_equal(second$expected_observations, 1)

  synthetic <- generate_pmx(model, 200, seed = 141)
  sampled <- vapply(split(synthetic, synthetic$ID), function(subject) {
    any(subject$EVID == 0 & subject$DVID == "cp" & subject$OCC == 2L)
  }, logical(1))
  expect_gt(mean(sampled), 0.15)
  expect_lt(mean(sampled), 0.35)
})

test_that("numeric domains are clipped before every aggregate", {
  extreme <- private_fixture()
  observed_cp <- extreme$EVID == 0 & extreme$DVID == "cp"
  observed_pd <- extreme$EVID == 0 & extreme$DVID == "pd"
  dose <- extreme$EVID != 0
  extreme$DV[observed_cp] <- 1e6
  extreme$DV[observed_pd] <- 1e6
  extreme$AMT[dose] <- 1e6
  extreme$RATE[dose] <- 1e6
  extreme$WT <- 1e6
  extreme$AGE <- 1e6

  clipped <- extreme
  clipped$DV[observed_cp] <- 20
  clipped$DV[observed_pd] <- 120
  clipped$AMT[dose] <- 200
  clipped$RATE[dose] <- 200
  clipped$WT <- 120
  clipped$AGE <- 90

  extreme_model <- fit_public_fixture(extreme)
  clipped_model <- fit_public_fixture(clipped)
  expect_identical(extreme_model$population, clipped_model$population)
})

test_that("accounting, ledger, and serialization contain no raw payload", {
  set.seed(1301)
  rng_state <- .Random.seed
  model <- fit_public_fixture()
  expect_identical(.Random.seed, rng_state)
  validation <- validate_private_model(model)
  expect_true(validation$valid)
  report <- privacy_report(model)
  expect_false(report$formal_dp)
  expect_match(report$delta_justification, "asserted public")
  expect_match(report$qualification, "No privacy guarantee")
  expect_lte(report$accounting$realized_epsilon, report$epsilon)
  expect_lte(report$accounting$realized_delta, report$delta)
  expect_equal(report$adjacency, "add-or-remove one complete subject")
  expect_true(nzchar(report$release_id))
  banned <- c("raw_rows", "raw_data", "source_ids", "subject_profiles",
              "event_templates", "raw_residuals", "unnoised_aggregates")
  expect_length(intersect(
    tolower(synpmx:::.recursive_names(model)), banned
  ), 0L)
  expect_false(any(grepl(
    "seed|unnoised|privacy_noise",
    tolower(synpmx:::.recursive_names(model))
  )))
  expect_null(model$source_ids)
  expect_null(model$raw_data)

  second <- fit_public_fixture()
  expect_false(identical(model$ledger$release_id, second$ledger$release_id))

  forged <- model
  forged$privacy$formal_dp <- TRUE
  forged$privacy$backend$validated <- TRUE
  expect_false(validate_private_model(forged)$valid)

  overspent <- model
  overspent$privacy$accounting$realized_epsilon <- 6
  expect_false(validate_private_model(overspent)$valid)
})

test_that("serialized models omit source identifier values, including levels", {
  source <- private_fixture()
  source$ID <- factor(
    paste0("source-secret-id-", source$ID),
    levels = paste0("source-secret-id-", seq_len(8L))
  )
  model <- fit_public_fixture(source)
  id_schema <- model$public$design$schema$columns[[1L]]
  expect_identical(id_schema$levels, character())
  serialized <- rawToChar(serialize(model, NULL, ascii = TRUE))
  expect_false(grepl("source-secret-id-", serialized, fixed = TRUE))
  synthetic <- generate_pmx(model, 4, 77)
  expect_true(is.factor(synthetic$ID))
  expect_true(all(grepl("^syn_", as.character(synthetic$ID))))
})

test_that("missing optional event amount or rate values remain bounded", {
  source <- private_fixture()
  dose <- source$EVID != 0
  source$AMT[which(dose)[1L]] <- NA_real_
  source$RATE[which(dose)[2L]] <- NA_real_
  model <- fit_public_fixture(source)
  expect_true(validate_private_model(model)$valid)
  expect_true(validate_pmx(
    generate_pmx(model, 4, 88), private_roles(), private_endpoints()
  )$valid)
})

test_that("neighboring inputs differ by one complete subject", {
  full <- private_fixture(8)
  neighbor <- full[full$ID != 8L, , drop = FALSE]
  first <- fit_public_fixture(full)
  second <- fit_public_fixture(neighbor)
  expect_false(identical(first$population, second$population))
  expect_equal(first$privacy$adjacency, "add-or-remove one complete subject")
})

test_that("direct identifiers and unmodeled datetimes fail", {
  source <- private_fixture()
  source$PATIENT_NAME <- "example"
  design <- private_design(source)
  expect_error(
    fit_private_pmx(
      source, private_roles(), private_endpoints(), 5, 0,
      private_bounds(), design, private_limits(), private_budget(),
      backend = "public", public_source = TRUE
    ),
    "direct identifier"
  )

  source <- private_fixture()
  source$contact_email <- "example@example.invalid"
  design <- private_design(source)
  expect_error(
    fit_private_pmx(
      source, private_roles(), private_endpoints(), 5, 0,
      private_bounds(), design, private_limits(), private_budget(),
      backend = "public", public_source = TRUE
    ),
    "direct identifier"
  )

  source <- private_fixture()
  source$COLLECTION_DATE <- as.Date("2025-01-01")
  roles <- private_roles()
  design <- private_design(source)
  expect_error(
    fit_private_pmx(
      source, roles, private_endpoints(), 5, 0, private_bounds(), design,
      private_limits(), private_budget(), backend = "public",
      public_source = TRUE
    ), "Date/POSIX"
  )
})

test_that("OpenDP adapter tests fail closed or pass canonically", {
  status <- dp_backend_status()
  if (!status$available) {
    expect_error(run_dp_backend_tests(), "unavailable")
  } else {
    result <- run_dp_backend_tests()
    expect_true(result$passed)
  }
})
