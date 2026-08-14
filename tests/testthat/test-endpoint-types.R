# SIM-045. A PMX table keeps every endpoint in one numeric DV column, so
# restoring the *column's* class restores nothing about an endpoint: on
# `xgxr::mad` a 0/1 endpoint came back as 600 distinct values spanning -0.13 to
# 1.08, an ordinal 1/2/3 endpoint reached 4.69, and an integer count came back
# in fractions. Blending is a weighted mean and the noise terms are continuous,
# so this is what the mechanism does unless the values are put back on their
# scale at emit.
#
# These tests pin four separate things, because the earlier defects in this
# package all slipped through where one part was verified and the whole question
# was asked of nobody: what the detector decides, what generation emits, that
# declaring overrides both, and that a continuous endpoint is untouched.

discrete_fixture <- function(n = 24L) {
  pieces <- lapply(seq_len(n), function(subject) {
    time <- c(0, 0, 0, 0, 2, 2, 2, 6, 6, 6, 12, 12, 12)
    endpoint <- c("dose", "conc", "score", "flag",
                  "conc", "score", "flag",
                  "conc", "score", "flag",
                  "conc", "score", "flag")
    evid <- c(1L, rep(0L, 12L))
    # `conc` is continuous, `score` is a 1-3 scale, `flag` is 0/1, and every
    # one of them lives in the same DV column.
    shape <- c(0, 5.5, 2, 1, 8.25, 3, 1, 4.5, 2, 0, 1.75, 1, 0)
    dv <- shape
    dv[endpoint == "conc"] <- shape[endpoint == "conc"] *
      (1 + (subject - n / 2) * 0.02)
    dv[endpoint == "score"] <- pmin(3, pmax(1, shape[endpoint == "score"] +
                                              (subject %% 3L) - 1L))
    dv[endpoint == "flag"] <- as.numeric((shape[endpoint == "flag"] +
                                            subject) %% 2L == 0L)
    dv[evid != 0] <- 0
    data.frame(
      ID = as.integer(subject), TIME = time, DV = dv,
      AMT = ifelse(evid != 0, 100, 0), EVID = evid,
      DVID = endpoint, WT = 60 + subject,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

discrete_roles <- function(endpoint_types = NULL) {
  pmx_roles(
    id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
    dvid = "DVID", covariates = "WT", endpoint_types = endpoint_types
  )
}

observed_values <- function(data, roles, endpoint) {
  rows <- data$EVID == 0 & data$DVID == endpoint & !is.na(data$DV)
  data$DV[rows]
}

test_that("the detector separates binary, ordinal and continuous endpoints", {
  types <- as.data.frame(pmx_endpoint_types(discrete_fixture(),
                                            discrete_roles()))
  expect_identical(types$type[types$endpoint == "flag"], "binary")
  expect_identical(types$type[types$endpoint == "score"], "ordinal")
  expect_identical(types$type[types$endpoint == "conc"], "continuous")
  expect_identical(types$levels[types$endpoint == "score"], "1, 2, 3")
  expect_identical(types$decided_by[types$endpoint == "flag"], "inferred")
})

test_that("many whole-number levels are `integer`, not a scale to snap onto", {
  data <- discrete_fixture()
  # A percentage recorded without decimals, which is warfarin's `pca`: whole
  # numbers, but far too many levels to be a scale.
  data$DV[data$DVID == "score"] <- rep(
    seq(9, 100, length.out = 4L) + rep(c(0, 1, 2, 3, 4, 5), each = 4L)[1:4],
    length.out = sum(data$DVID == "score")
  )
  data$DV[data$DVID == "score"] <- round(data$DV[data$DVID == "score"]) +
    rep(seq_len(20L), length.out = sum(data$DVID == "score"))
  types <- as.data.frame(pmx_endpoint_types(data, discrete_roles()))
  expect_identical(types$type[types$endpoint == "score"], "integer")
  expect_identical(types$levels[types$endpoint == "score"], "--")
})

test_that("too few observations is left continuous rather than guessed at", {
  data <- discrete_fixture(n = 2L)
  types <- as.data.frame(pmx_endpoint_types(data, discrete_roles()))
  expect_identical(types$type[types$endpoint == "flag"], "continuous")
  expect_match(types$reason[types$endpoint == "flag"], "too few to call")
})

test_that("generated discrete endpoints stay on their source scale", {
  data <- discrete_fixture()
  roles <- discrete_roles()
  synthetic <- suppressWarnings(synpmx_avatar(data, roles, seed = 11))

  flag <- observed_values(synthetic, roles, "flag")
  score <- observed_values(synthetic, roles, "score")
  conc <- observed_values(synthetic, roles, "conc")

  expect_true(all(flag %in% c(0, 1)))
  expect_true(all(score %in% c(1, 2, 3)))
  # The continuous endpoint must NOT be quantized as a side effect, and the
  # blend has to still be doing its job: identical values everywhere would pass
  # the two checks above for the wrong reason.
  expect_gt(length(unique(conc)), 10L)
  expect_false(all(abs(conc - round(conc)) < 1e-8))

  settings <- attr(synthetic, "pmx_settings")
  expect_identical(settings$endpoint_type_violations, 0L)
  expect_identical(settings$endpoint_value_types$flag$type, "binary")
  expect_identical(settings$endpoint_value_types$score$levels, c(1, 2, 3))
})

test_that("without the repair the same run leaves the source scale", {
  # Pins the defect rather than only the fix: with every endpoint declared
  # continuous, generation reproduces `SIM-045` exactly, so a regression that
  # quietly stops snapping cannot pass this file.
  data <- discrete_fixture()
  roles <- discrete_roles(endpoint_types = c(flag = "continuous",
                                             score = "continuous"))
  synthetic <- suppressWarnings(synpmx_avatar(data, roles, seed = 11))
  flag <- observed_values(synthetic, roles, "flag")
  expect_false(all(flag %in% c(0, 1)))
})

test_that("a declared type overrides the inference in both directions", {
  data <- discrete_fixture()
  # `conc` is continuous by measurement; declaring it integer rounds it.
  roles <- discrete_roles(endpoint_types = c(conc = "count"))
  types <- as.data.frame(pmx_endpoint_types(data, roles))
  expect_identical(types$type[types$endpoint == "conc"], "integer")
  expect_identical(types$decided_by[types$endpoint == "conc"], "declared")

  synthetic <- suppressWarnings(synpmx_avatar(data, roles, seed = 11))
  conc <- observed_values(synthetic, roles, "conc")
  expect_true(all(abs(conc - round(conc)) < 1e-8))
})

test_that("an impossible declaration is refused with the alternative named", {
  data <- discrete_fixture()
  expect_error(
    pmx_endpoint_types(data, discrete_roles(endpoint_types = c(conc = "binary"))),
    "declared binary"
  )
  expect_error(
    pmx_endpoint_types(data, discrete_roles(endpoint_types = c(conc = "ordinal"))),
    "distinct observed values"
  )
  expect_error(
    pmx_endpoint_types(data, discrete_roles(endpoint_types = c(nope = "binary"))),
    "not present in the data"
  )
  expect_error(pmx_roles(id = "ID", time = "TIME", dv = "DV", evid = "EVID",
                         endpoint_types = c(flag = "categorical")),
               "must be one of")
  expect_error(pmx_roles(id = "ID", time = "TIME", dv = "DV", evid = "EVID",
                         endpoint_types = "binary"),
               "must be named by endpoint")
})

test_that("`endpoint_types` is not treated as a column name", {
  data <- discrete_fixture()
  roles <- discrete_roles(endpoint_types = c(flag = "binary"))
  # `.assert_roles()` and the retained-column allowlist both unlist the roles
  # object; an entry keyed by endpoint rather than by column must not reach
  # either, or every run errors with "role columns not found".
  expect_silent(synpmx:::.assert_roles(data, roles))
  expect_false("binary" %in% synpmx:::.retained_role_columns(roles))
})

test_that("A6 is on every card, discrete endpoints or not", {
  data <- discrete_fixture()
  roles <- discrete_roles()
  synthetic <- suppressWarnings(synpmx_avatar(data, roles, seed = 11))
  card <- as.data.frame(synpmx_scorecard(data, synthetic, roles))
  expect_identical(card$result[card$check == "A6"], "2 of 2")
  expect_identical(card$verdict[card$check == "A6"], "pass")

  continuous_only <- data[data$DVID %in% c("dose", "conc"), ]
  continuous_roles <- discrete_roles()
  continuous_synthetic <- suppressWarnings(
    synpmx_avatar(continuous_only, continuous_roles, seed = 11)
  )
  continuous_card <- as.data.frame(
    synpmx_scorecard(continuous_only, continuous_synthetic, continuous_roles)
  )
  # Present and passing, not absent. A card whose rows depend on the study
  # cannot be compared across studies, and an absent row reads as one that
  # passed when it means the question was never asked.
  expect_identical(continuous_card$result[continuous_card$check == "A6"],
                   "no discrete endpoint")
  expect_identical(continuous_card$verdict[continuous_card$check == "A6"],
                   "pass")
})

test_that("the scorecard fails when a discrete endpoint left its scale", {
  data <- discrete_fixture()
  roles <- discrete_roles()
  synthetic <- suppressWarnings(synpmx_avatar(data, roles, seed = 11))
  damaged <- synthetic
  target <- which(damaged$EVID == 0 & damaged$DVID == "flag")[[1L]]
  damaged$DV[target] <- 0.4
  attr(damaged, "pmx_settings") <- attr(synthetic, "pmx_settings")
  card <- as.data.frame(synpmx_scorecard(data, damaged, roles))
  expect_identical(card$result[card$check == "A6"], "1 of 2")
  expect_identical(card$verdict[card$check == "A6"], "FAIL")
})
