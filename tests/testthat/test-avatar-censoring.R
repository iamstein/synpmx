# AVATAR censoring reconstruction.
#
# Before this existed, a declared `cens` role was accepted for schema validation
# but CENS was copied from the anchor template while DV was independently
# blended, so the two disagreed and the characteristic spike of identical values
# at the assay limit was smeared into distinct numbers. What is tested here is
# that DV, CENS, and LIMIT are now reconstructed together from one latent value.

cens_roles <- function(limit = "LIMIT") {
  pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
            dvid = "DVID", cmt = "CMT", mdv = "MDV", cens = "CENS",
            limit = limit, covariates = "WT")
}

# A left-censored source in the Monolix convention: DV sits at the limit on a
# censored row, LIMIT is absent.
left_censored_source <- function(n = 30, quantile = 0.25) {
  source <- pmx_simulated_fixture(n)
  source <- source[source$DVID == "cp", ]
  observed <- source$EVID == 0 & !is.na(source$DV)
  lloq <- as.numeric(stats::quantile(source$DV[observed], quantile,
                                     na.rm = TRUE))
  source$CENS <- ifelse(observed & source$DV < lloq, 1L, 0L)
  source$LIMIT <- NA_real_
  source$DV[source$CENS == 1L] <- lloq
  list(data = source, lloq = lloq)
}

test_that("censored synthetic rows report DV exactly at the limit", {
  fixture <- left_censored_source()
  synthetic <- suppressWarnings(
    synpmx_avatar(fixture$data, cens_roles(), n_subjects = 30, seed = 7)
  )
  observed <- synthetic[synthetic$EVID == 0 & !is.na(synthetic$DV), ]
  censored <- observed$CENS == 1

  expect_true(any(censored))
  # The source reports one value on every censored row; so must the output.
  expect_length(unique(round(observed$DV[censored], 9)), 1L)
  expect_equal(unique(observed$DV[censored]), fixture$lloq)
})

test_that("CENS and DV agree on every generated row", {
  fixture <- left_censored_source()
  synthetic <- suppressWarnings(
    synpmx_avatar(fixture$data, cens_roles(), n_subjects = 30, seed = 7)
  )
  observed <- synthetic[synthetic$EVID == 0 & !is.na(synthetic$DV), ]
  tolerance <- 1e-9
  # A flagged row above the limit, or an unflagged row below it, is exactly the
  # incoherence this feature exists to remove.
  expect_equal(sum(observed$CENS == 1 & observed$DV > fixture$lloq + tolerance), 0L)
  expect_equal(sum(observed$CENS == 0 & observed$DV < fixture$lloq - tolerance), 0L)
})

test_that("the censored fraction stays in the neighbourhood of the source", {
  # Blending shrinks variance toward the middle and the imputation draws
  # uniformly below the limit, so the rate is reproduced approximately rather
  # than exactly. A wide band still catches a mechanism that has stopped
  # censoring at all, or that censors everything.
  fixture <- left_censored_source(quantile = 0.25)
  synthetic <- suppressWarnings(
    synpmx_avatar(fixture$data, cens_roles(), n_subjects = 30, seed = 7)
  )
  observed <- synthetic[synthetic$EVID == 0 & !is.na(synthetic$DV), ]
  rate <- mean(observed$CENS == 1)
  expect_gt(rate, 0.10)
  expect_lt(rate, 0.45)
})

test_that("an uncensored source is untouched by the censoring path", {
  source <- pmx_simulated_fixture(20)
  source <- source[source$DVID == "cp", ]
  source$CENS <- 0L
  source$LIMIT <- NA_real_
  synthetic <- suppressWarnings(
    synpmx_avatar(source, cens_roles(), n_subjects = 20, seed = 3)
  )
  observed <- synthetic[synthetic$EVID == 0 & !is.na(synthetic$DV), ]
  # No boundary can be read from a source with nothing censored, so nothing is
  # imposed on the output.
  expect_true(all(observed$CENS == 0))
})

test_that("right-censored rows are reconstructed at the upper boundary", {
  source <- pmx_simulated_fixture(30)
  source <- source[source$DVID == "cp", ]
  observed <- source$EVID == 0 & !is.na(source$DV)
  uloq <- as.numeric(stats::quantile(source$DV[observed], 0.85, na.rm = TRUE))
  source$CENS <- ifelse(observed & source$DV > uloq, -1L, 0L)
  source$LIMIT <- NA_real_
  source$DV[source$CENS == -1L] <- uloq

  synthetic <- suppressWarnings(
    synpmx_avatar(source, cens_roles(), n_subjects = 30, seed = 11)
  )
  rows <- synthetic[synthetic$EVID == 0 & !is.na(synthetic$DV), ]
  censored <- rows$CENS == -1
  expect_true(any(censored))
  expect_length(unique(round(rows$DV[censored], 9)), 1L)
  expect_equal(unique(rows$DV[censored]), uloq)
  # Nothing above the boundary escapes the flag.
  expect_equal(sum(rows$CENS == 0 & rows$DV > uloq + 1e-9), 0L)
})

test_that("a declared cens role with no readable boundary warns", {
  # CENS flagged but DV missing on those rows: nothing can be recovered, and
  # silently carrying the flag through is what used to produce wrong data.
  source <- pmx_simulated_fixture(20)
  source <- source[source$DVID == "cp", ]
  observed <- source$EVID == 0 & !is.na(source$DV)
  lloq <- as.numeric(stats::quantile(source$DV[observed], 0.2, na.rm = TRUE))
  source$CENS <- ifelse(observed & source$DV < lloq, 1L, 0L)
  source$LIMIT <- NA_real_
  source$DV[source$CENS == 1L] <- NA_real_
  source$MDV[source$CENS == 1L] <- 1L

  expect_warning(
    synpmx_avatar(source, cens_roles(), n_subjects = 20, seed = 5),
    "no censoring boundary could be read"
  )
})

test_that("imputation keeps censored donors from dragging the blend up", {
  # Every censored source row reports the limit. Blending those substituted
  # values directly would put a floor under the synthetic data that the real
  # study does not have, so the imputed values must sit below the limit.
  fixture <- left_censored_source()
  imputed <- .impute_censored(fixture$data, cens_roles())
  was_censored <- fixture$data$CENS == 1L
  expect_true(all(imputed$DV[was_censored] < fixture$lloq))
  expect_true(all(imputed$DV[was_censored] >= 0))
  # Uncensored rows are left exactly as they were.
  expect_identical(imputed$DV[!was_censored], fixture$data$DV[!was_censored])
})

test_that("interval-censored rows carry both ends of the interval", {
  source <- pmx_simulated_fixture(30)
  source <- source[source$DVID == "cp", ]
  observed <- source$EVID == 0 & !is.na(source$DV)
  lower <- as.numeric(stats::quantile(source$DV[observed], 0.10, na.rm = TRUE))
  upper <- as.numeric(stats::quantile(source$DV[observed], 0.30, na.rm = TRUE))
  inside <- observed & source$DV > lower & source$DV < upper
  source$CENS <- ifelse(inside, 1L, 0L)
  source$LIMIT <- ifelse(inside, lower, NA_real_)
  source$DV[inside] <- upper

  synthetic <- suppressWarnings(
    synpmx_avatar(source, cens_roles(), n_subjects = 30, seed = 13)
  )
  rows <- synthetic[synthetic$EVID == 0 & !is.na(synthetic$DV), ]
  censored <- rows$CENS == 1
  expect_true(any(censored))
  # DV reports the upper end, LIMIT the lower: an interval, not a point.
  expect_equal(unique(rows$DV[censored]), upper)
  expect_equal(unique(rows$LIMIT[censored]), lower)
  expect_true(all(is.na(rows$LIMIT[!censored])))
})

# REV-021: validate_pmx() must reject data whose CENS flag disagrees with DV.

test_that("validate_pmx flags a censored value above an uncensored one", {
  source <- pmx_simulated_fixture(20)
  source <- source[source$DVID == "cp", ]
  observed <- source$EVID == 0 & !is.na(source$DV)
  # Flag a quarter of rows at random without touching DV: exactly the old
  # copied-from-anchor behaviour, where CENS and DV were produced independently.
  source$CENS <- 0L
  source$LIMIT <- NA_real_
  flagged <- sample(which(observed), max(1L, floor(sum(observed) * 0.25)))
  source$CENS[flagged] <- 1L

  report <- validate_pmx(source, cens_roles())
  expect_false(report$valid)
  coherence <- report$checks$status[report$checks$check == "cens_dv_coherence"]
  expect_identical(coherence, "error")
})

test_that("validate_pmx accepts coherently censored data", {
  fixture <- left_censored_source()
  # The source itself is coherent: every censored row sits at the limit, every
  # uncensored row above it.
  report <- validate_pmx(fixture$data, cens_roles())
  coherence <- report$checks$status[report$checks$check == "cens_dv_coherence"]
  expect_identical(coherence, "pass")
})

test_that("the coherence check does not fire without a cens role", {
  source <- pmx_simulated_fixture(15)
  source <- source[source$DVID == "cp", ]
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", dvid = "DVID", cmt = "CMT", mdv = "MDV")
  report <- validate_pmx(source, roles)
  expect_false("cens_dv_coherence" %in% report$checks$check)
  expect_true(report$valid)
})

test_that("right-censoring coherence is checked in the mirror direction", {
  source <- pmx_simulated_fixture(20)
  source <- source[source$DVID == "cp", ]
  observed <- source$EVID == 0 & !is.na(source$DV)
  source$CENS <- 0L
  source$LIMIT <- NA_real_
  # Flag a low-valued row as right-censored: a right-censored value should be an
  # upper measurement, so one below ordinary values is incoherent.
  lowest <- which(observed)[which.min(source$DV[observed])]
  source$CENS[lowest] <- -1L

  report <- validate_pmx(source, cens_roles())
  coherence <- report$checks$status[report$checks$check == "cens_dv_coherence"]
  expect_identical(coherence, "error")
})
