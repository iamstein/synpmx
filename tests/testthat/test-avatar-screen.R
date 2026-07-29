# synpmx_avatar(screen = TRUE) is the default: a source subject whose follow-up
# or dose count is more than twice the cohort's 90th percentile is not used as
# an anchor, so no avatar inherits its extreme skeleton. A source with no such
# outlier is unaffected.

scr_sub <- function(id, obs) {
  data.frame(
    ID = as.character(id), TIME = c(0, obs),
    DV = c(0, 5 * exp(-0.1 * obs)),
    AMT = c(100, rep(0, length(obs))),
    EVID = c(1L, rep(0L, length(obs))),
    CMT = c(1L, rep(2L, length(obs))),
    stringsAsFactors = FALSE
  )
}
scr_roles <- function() {
  pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
            evid = "EVID", cmt = "CMT")
}
obs_max <- function(d) max(d$TIME[d$EVID == 0])

test_that("a gross follow-up outlier is never anchored by default", {
  src <- rbind(
    do.call(rbind, lapply(1:8, scr_sub, obs = c(0.5, 1, 2, 4))),
    scr_sub(9, obs = c(0.5, 1, 2, 100))          # extreme long follow-up
  )
  roles <- scr_roles()

  # `min_pattern_share = 1` keeps attendance sampling out of this test. Sampling
  # would also suppress the 100-hour structure -- the outlier's pattern is held
  # by nobody else, so it never enters the pool -- and then `screen = FALSE`
  # would look identical to `screen = TRUE` for reasons that have nothing to do
  # with the screen. The two mechanisms overlap here by design; this test is
  # about one of them.
  screened <- suppressWarnings(
    synpmx_avatar(src, roles, n_subjects = 40, seed = 1, min_pattern_share = 1)
  )
  raw <- suppressWarnings(
    synpmx_avatar(src, roles, n_subjects = 40, seed = 1, screen = FALSE,
                  min_pattern_share = 1)
  )

  expect_lte(obs_max(screened), 4)        # the 100-hour structure never appears
  expect_gt(obs_max(raw), 4)              # unscreened, it is inherited
  expect_false(identical(screened, raw))
  expect_equal(length(unique(screened$ID)), 40L)   # count preserved
})

test_that("a source with no structural outlier is unaffected by screening", {
  src <- do.call(rbind, lapply(1:8, scr_sub, obs = c(0.5, 1, 2, 4)))
  roles <- scr_roles()

  expect_identical(
    suppressWarnings(synpmx_avatar(src, roles, n_subjects = 20, seed = 3)),
    suppressWarnings(synpmx_avatar(src, roles, n_subjects = 20, seed = 3,
                                   screen = FALSE))
  )
})

test_that("screening leaves dose-magnitude variation alone", {
  # Weight-based-style dosing: every subject a slightly different dose. None of
  # this is a follow-up or dose-count outlier, so screening changes nothing.
  src <- do.call(rbind, Map(function(i, d) {
    s <- scr_sub(i, obs = c(0.5, 1, 2, 4)); s$AMT[1] <- d; s
  }, 1:8, seq(300, 370, length.out = 8)))
  roles <- scr_roles()

  expect_identical(
    suppressWarnings(synpmx_avatar(src, roles, n_subjects = 20, seed = 5)),
    suppressWarnings(synpmx_avatar(src, roles, n_subjects = 20, seed = 5,
                                   screen = FALSE))
  )
})
