# REV-051. A regimen written as one record plus `ADDL`/`II` has to reach every
# reader as the doses it stands for, and leave as the record it arrived as.

compact_roles <- function(...) {
  pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
            addl = "ADDL", ii = "II", ...)
}
compact_study <- function() data.frame(
  ID = c(1, 1, 1, 2, 2),
  TIME = c(0, 12, 36, 0, 6),
  DV = c(NA, 3.1, 2.4, NA, 4.0),
  AMT = c(100, 0, 0, 50, 0),
  EVID = c(1L, 0L, 0L, 1L, 0L),
  ADDL = c(2L, 0L, 0L, 0L, 0L),
  II = c(24, 0, 0, 0, 0)
)

test_that("expansion writes every implied dose out, in order", {
  out <- pmx_expand_doses(compact_study(), compact_roles())
  doses <- out[out$EVID == 1L, ]
  expect_equal(nrow(out), 7L)
  expect_equal(doses$TIME[doses$ID == 1], c(0, 24, 48))
  expect_equal(doses$AMT[doses$ID == 1], rep(100, 3))
  expect_true(all(out$ADDL == 0L) && all(out$II == 0))
  expect_true(is.integer(out$ADDL))
  # Time-ordered within subject, dose before observation at a shared time.
  expect_equal(out$TIME[out$ID == 1], c(0, 12, 24, 36, 48))
  expect_equal(out$TIME[out$ID == 2], c(0, 6))
})

test_that("expansion moves a declared nominal time with the dose", {
  study <- compact_study()
  study$NTIME <- study$TIME
  out <- pmx_expand_doses(study, compact_roles(nominal_time = "NTIME"))
  doses <- out[out$EVID == 1L & out$ID == 1, ]
  expect_equal(doses$NTIME, c(0, 24, 48))
})

test_that("expansion refuses an interval it cannot place doses on", {
  study <- compact_study()
  study$II[1L] <- 0
  expect_error(pmx_expand_doses(study, compact_roles()), "positive interval")
  study <- compact_study()
  study$ADDL[1L] <- 1.5
  expect_error(pmx_expand_doses(study, compact_roles()), "whole number")
})

test_that("neither function touches a study that declares no compression", {
  plain <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID")
  study <- compact_study()
  expect_identical(pmx_expand_doses(study, plain), study)
  expect_identical(pmx_compress_doses(study, plain), study)
  expect_error(pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                         evid = "EVID", addl = "ADDL") |>
                 pmx_expand_doses(data = study),
               "declared together")
})

test_that("compression is the exact inverse of expansion", {
  roles <- compact_roles()
  study <- compact_study()
  back <- pmx_compress_doses(pmx_expand_doses(study, roles), roles)
  expect_equal(nrow(back), nrow(study))
  expect_equal(back$ADDL, study$ADDL)
  expect_equal(back$II, study$II)
  expect_equal(back$TIME, study$TIME)
})

test_that("compression folds onc_sim back to itself", {
  roles <- pmx_roles(id = "ID", time = "TIME", nominal_time = "NTIME",
                     dv = "DV", amt = "AMT", evid = "EVID", cmt = "CMT",
                     dvid = "NAME", addl = "ADDL", ii = "II")
  expanded <- pmx_expand_doses(onc_sim, roles)
  # Every dose is now a row: 237 records stood for 71,180 doses.
  expect_equal(sum(expanded$EVID == 1L),
               sum(onc_sim$ADDL[onc_sim$EVID == 1L] + 1L))
  back <- pmx_compress_doses(expanded, roles)
  expect_equal(nrow(back), nrow(onc_sim))
  expect_equal(back$ADDL, onc_sim$ADDL)
  expect_equal(back$II, onc_sim$II)
  expect_equal(back$AMT, onc_sim$AMT)
})

test_that("compression leaves what expansion could not have produced", {
  roles <- compact_roles()
  # A run ends where the interval breaks. Two doses are always one interval
  # apart, so the first pair folds; the drifted third dose stays its own record
  # because `ADDL` could not place it.
  drift <- data.frame(ID = 1, TIME = c(0, 24, 48.1), DV = NA, AMT = 100,
                      EVID = 1L, ADDL = 0L, II = 0)
  out <- pmx_compress_doses(drift, roles)
  expect_equal(nrow(out), 2L)
  expect_equal(out$ADDL, c(1L, 0L))
  expect_equal(out$TIME, c(0, 48.1))
  # A dose change is a different record and starts its own run.
  change <- data.frame(ID = 1, TIME = c(0, 24, 48, 72), DV = NA,
                       AMT = c(100, 100, 50, 50), EVID = 1L, ADDL = 0L, II = 0)
  out <- pmx_compress_doses(change, roles)
  expect_equal(nrow(out), 2L)
  expect_equal(out$AMT, c(100, 50))
  expect_equal(out$ADDL, c(1L, 1L))
  # Both clocks have to agree with the interval: a regular recorded time over
  # a nominal grid that drifts breaks the run where the nominal gap does.
  nominal <- data.frame(ID = 1, TIME = c(0, 24, 48), NTIME = c(0, 24, 47),
                        DV = NA, AMT = 100, EVID = 1L, ADDL = 0L, II = 0)
  out <- pmx_compress_doses(nominal, compact_roles(nominal_time = "NTIME"))
  expect_equal(nrow(out), 2L)
  expect_equal(out$ADDL, c(1L, 0L))
})

test_that("compression creates the columns a generator never wrote", {
  roles <- compact_roles()
  bare <- data.frame(ID = 1, TIME = c(0, 24, 48, 12), DV = c(NA, NA, NA, 2),
                     AMT = c(100, 100, 100, 0), EVID = c(1L, 1L, 1L, 0L))
  out <- pmx_compress_doses(bare, roles)
  expect_true(all(c("ADDL", "II") %in% names(out)))
  expect_equal(nrow(out), 2L)
  expect_equal(out$ADDL[out$EVID == 1L], 2L)
  expect_equal(out$II[out$EVID == 1L], 24)
  # Beside the amount, where a reader looks for them.
  expect_equal(match("ADDL", names(out)), match("AMT", names(out)) + 1L)
})

test_that("the dose skeleton sees every dose once expanded", {
  roles <- compact_roles()
  study <- compact_study()
  expect_equal(sum(.dose_rows(study, roles)), 2L)
  expect_equal(sum(.dose_rows(pmx_expand_doses(study, roles), roles)), 4L)
})

test_that("validate_pmx checks a declared TAD against the expanded doses", {
  study <- compact_study()
  study$TAD <- c(0, 12, 12, 0, 6)         # 36 h is 12 h after the dose at 24
  report <- validate_pmx(study, compact_roles(tad = "TAD"))
  row <- report$checks[report$checks$check == "tad_agreement", ]
  expect_equal(row$status, "pass")
  study$TAD[3L] <- 36                     # the un-expanded reading
  report <- validate_pmx(study, compact_roles(tad = "TAD"))
  row <- report$checks[report$checks$check == "tad_agreement", ]
  expect_equal(row$status, "warning")
  expect_false(grepl("not checked", row$message))
})

test_that("the avatar returns a compressed study from a compressed one", {
  one <- function(id) rbind(
    data.frame(ID = id, TIME = 0, DV = NA_real_, AMT = 10, EVID = 1L,
               CMT = 1L, ADDL = 4L, II = 1, WT = 70 + id),
    data.frame(ID = id, TIME = c(0.5, 1, 2, 3, 4, 5),
               DV = c(2, 3, 3.6, 3.9, 4, 4.1) * (1 + id / 50),
               AMT = 0, EVID = 0L, CMT = 2L, ADDL = 0L, II = 0, WT = 70 + id))
  source <- do.call(rbind, lapply(1:12, one))
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", cmt = "CMT", addl = "ADDL", ii = "II",
                     covariates = "WT")
  synthetic <- suppressWarnings(suppressMessages(
    synpmx_avatar(source, roles, n_subjects = 12L, seed = 3)))
  doses <- synthetic[synthetic$EVID == 1L, ]
  # One record per subject standing for five doses, as the source wrote it.
  expect_equal(nrow(doses), 12L)
  expect_true(all(doses$ADDL == 4L))
  expect_true(all(doses$II == 1))
  expect_true(validate_pmx(synthetic, roles)$valid)
})

test_that("the scorecard counts doses, not dose records", {
  one <- function(id, addl, ii) rbind(
    data.frame(ID = id, TIME = 0, DV = NA_real_, AMT = 10, EVID = 1L, CMT = 1L,
               ADDL = addl, II = ii, WT = 70 + id),
    data.frame(ID = id, TIME = c(1, 2, 3, 4, 5),
               DV = c(2, 3, 3.6, 3.9, 4) * (1 + id / 60),
               AMT = 0, EVID = 0L, CMT = 2L, ADDL = 0L, II = 0, WT = 70 + id))
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", cmt = "CMT", addl = "ADDL", ii = "II",
                     covariates = "WT")
  # Same doses on both sides, written in different numbers of records: one
  # block of five against five blocks of one.
  source <- do.call(rbind, lapply(1:10, one, addl = 4L, ii = 1))
  spread <- do.call(rbind, lapply(1:10, function(id) {
    rows <- one(id, addl = 0L, ii = 0)
    doses <- rows[rep(1L, 5L), ]
    doses$TIME <- 0:4
    rbind(doses, rows[rows$EVID == 0L, ])
  }))
  a5b <- function(a, b) {
    card <- suppressMessages(synpmx_scorecard(a, b, roles))
    as.data.frame(card)$result[card$check == "A5b"]
  }
  # Five doses per patient either way, so the row reads no change at all.
  expect_match(a5b(source, spread), "^5 -> 5$")
  # And it is not fooled the other way round.
  expect_match(a5b(spread, source), "^5 -> 5$")
})
