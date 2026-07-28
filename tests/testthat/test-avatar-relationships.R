# AVATAR preserves covariate, dose, and treatment-arm relationships as a *side
# effect*: `.synthesize_covariates()` and `.synthesize_trajectories()` are handed
# the same donors and the same weights, and donors are ranked in a profile space
# built from the covariates and the DV trajectory together. Nothing declares
# that behaviour, so nothing would complain if a future change broke it --- for
# example by drawing covariates separately, or by reweighting donors per
# endpoint, which has been floated to stop a densely-sampled endpoint dominating
# donor selection.
#
# These tests pin the behaviour loosely. They are deliberately not tight: the
# retained fraction is a tendency, not a specification, and measured values move
# with cohort size and signal-to-noise (see the worked article
# `vignettes/articles/example-avatar-PKPD-covariate-treatment-effect.Rmd`).
# What they catch is the relationship being *destroyed*, not it drifting.

pkpd_source <- function(n_per_arm = 30, seed = 2026) {
  set.seed(seed)
  one <- function(id, arm, dose) {
    wt <- rnorm(1, 70, 12)
    cl <- 5 * (wt / 70)^0.75 * exp(rnorm(1, 0, 0.20))
    v <- 40 * (wt / 70) * exp(rnorm(1, 0, 0.15))
    ka <- 1.2
    tt <- c(0.5, 1, 2, 4, 8, 12, 24)
    cp <- dose / v * ka / (ka - cl / v) * (exp(-cl / v * tt) - exp(-ka * tt))
    rbind(
      data.frame(ID = id, TIME = 0, DV = NA_real_, AMT = dose, EVID = 1L,
                 CMT = 1L, WT = wt, ARM = arm),
      data.frame(ID = id, TIME = tt, DV = cp * exp(rnorm(length(tt), 0, 0.08)),
                 AMT = 0, EVID = 0L, CMT = 2L, WT = wt, ARM = arm)
    )
  }
  out <- do.call(rbind, c(
    lapply(seq_len(n_per_arm), one, arm = "low", dose = 300),
    lapply(n_per_arm + seq_len(n_per_arm), one, arm = "high", dose = 600)
  ))
  out <- out[order(out$ID, out$TIME, out$EVID == 0L), ]
  rownames(out) <- NULL
  out
}

pkpd_roles <- function() {
  pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
            cmt = "CMT", covariates = "WT", keep = "ARM")
}

subject_auc <- function(d) {
  obs <- d$EVID == 0 & !is.na(d$DV)
  pk <- d[obs, ]
  ids <- as.character(pk$ID)
  auc <- tapply(seq_len(nrow(pk)), ids, function(i) {
    t <- pk$TIME[i]; y <- pk$DV[i]; o <- order(t)
    sum(diff(t[o]) * (head(y[o], -1) + tail(y[o], -1)) / 2)
  })
  list(auc = auc,
       wt = tapply(pk$WT, ids, function(v) v[1]),
       arm = tapply(as.character(pk$ARM), ids, function(v) v[1]))
}

test_that("a covariate-exposure relationship survives blending", {
  # The weight effect must be estimated *within* dose, not pooled across arms:
  # AUC differs twofold by arm, which swamps the covariate signal and makes a
  # pooled slope almost pure dose. Conditioning on ARM isolates the effect the
  # allometric term actually put in.
  wt_slope <- function(d) {
    s <- subject_auc(d)
    coef(lm(log(s$auc) ~ log(s$wt) + s$arm))[[2]]
  }
  src <- pkpd_source()
  source_slope <- wt_slope(src)

  # A single run is genuinely noisy -- the worked article measures a run-to-run
  # SD around half the slope itself, enough for one draw to invert the sign of
  # a weak relationship. Averaging a few runs tests the tendency, which is what
  # this behaviour actually is.
  synthetic_slope <- mean(vapply(c(11L, 12L, 13L, 14L, 15L), function(sd_) {
    wt_slope(suppressWarnings(synpmx_avatar(src, pkpd_roles(), seed = sd_)))
  }, numeric(1)))

  expect_lt(source_slope, 0)
  expect_lt(synthetic_slope, 0)
  # Wide bounds on purpose. The article measures roughly 143% of the source
  # slope on a larger cohort -- amplified, because WT helps choose the donors --
  # so this is a destruction check, not an equality check.
  expect_gt(synthetic_slope / source_slope, 0.4)
  expect_lt(synthetic_slope / source_slope, 2.5)
})

test_that("the dose-exposure separation between arms survives", {
  src <- pkpd_source()
  s <- subject_auc(src)
  source_ratio <- mean(s$auc[s$arm == "high"]) / mean(s$auc[s$arm == "low"])

  syn <- suppressWarnings(synpmx_avatar(src, pkpd_roles(), seed = 12))
  y <- subject_auc(syn)
  synthetic_ratio <- mean(y$auc[y$arm == "high"]) / mean(y$auc[y$arm == "low"])

  expect_gt(source_ratio, 1.5)
  # Dose is part of the event signature, so donors never cross arms and this
  # ratio is preserved closely rather than merely surviving.
  expect_gt(synthetic_ratio / source_ratio, 0.75)
  expect_lt(synthetic_ratio / source_ratio, 1.25)
})

test_that("a kept treatment label stays coherent with the dose it implies", {
  src <- pkpd_source()
  syn <- suppressWarnings(synpmx_avatar(src, pkpd_roles(), seed = 13))
  doses <- syn[syn$EVID == 1L, c("ARM", "AMT")]

  # ARM is copied from the anchor while DV comes from the donors, so this is
  # only true because the event signature keeps donors inside the anchor's dose
  # group. If that ever changes, a "high" label will arrive on "low" data.
  expect_setequal(unique(doses$AMT[doses$ARM == "low"]), 300)
  expect_setequal(unique(doses$AMT[doses$ARM == "high"]), 600)
})

test_that("covariates and endpoints are filled from one donor draw", {
  # The structural guarantee under all of the above: one call to
  # `.select_donors()` per synthetic subject, whose result is passed to both
  # the covariate and the trajectory synthesiser.
  body_text <- paste(deparse(body(synpmx_avatar)), collapse = " ")
  expect_match(body_text, "donors <- \\.select_donors")
  expect_match(body_text, "\\.synthesize_covariates\\(")
  expect_match(body_text, "\\.synthesize_trajectories\\(")
  # Both must receive the same objects; a per-endpoint or per-column redraw
  # would break every relationship test above.
  expect_equal(lengths(regmatches(body_text, gregexpr("\\.select_donors", body_text))), 1L)
})
