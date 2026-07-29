# synpmx_avatar(coarsen_time = TRUE) is the default. It exists to close SIM-014
# ("Exact source timing vectors could be copied") on the AVATAR engine, where the
# gate had only ever been enforced against the structural/DP path even though
# AVATAR is the mode that actually copies an event skeleton verbatim.
#
# The mechanism is a two-step that only works in one order: snap every subject
# onto a shared visit grid, which is many-to-one and *destroys* the per-subject
# deviation, and only then resample deviations pooled across the cohort. Adding
# noise to the original time instead -- what `time_jitter` does -- leaves the
# original recoverable, because `.offset_unique_times()` clamps every time inside
# its own Voronoi cell. These tests pin the order, not just the outcome.

crs_sub <- function(id, times, offset) {
  data.frame(
    ID = as.character(id),
    TIME = c(0, times + offset),
    NTIME = c(0, times),
    DV = c(0, 5 * exp(-0.15 * times)),
    AMT = c(100, rep(0, length(times))),
    EVID = c(1L, rep(0L, length(times))),
    CMT = c(1L, rep(2L, length(times))),
    WT = 70 + id,
    stringsAsFactors = FALSE
  )
}
crs_roles <- function(nominal = NULL) {
  pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
            cmt = "CMT", covariates = "WT", nominal_time = nominal)
}
# Twelve subjects on one protocol schedule, each recorded at their own actual
# times: the regime where every event signature is unique.
crs_source <- function(n = 12L, sd = 0.05, seed = 4) {
  set.seed(seed)
  nominal <- c(0.5, 1, 2, 4, 8)
  do.call(rbind, lapply(seq_len(n), function(id) {
    crs_sub(id, nominal, stats::rnorm(length(nominal), sd = sd))
  }))
}

timing_vector_copied <- function(source, synthetic, roles) {
  src <- split(as.numeric(source[[roles$time]]), source[[roles$id]])
  gen <- split(as.numeric(synthetic[[roles$time]]), synthetic[[roles$id]])
  any(vapply(gen, function(g) {
    any(vapply(src, function(s) isTRUE(all.equal(g, s)), logical(1)))
  }, logical(1)))
}

test_that("actual recorded times leave every subject alone in its class", {
  report <- skeleton_uniqueness(crs_source(), crs_roles())
  # The premise of the whole mechanism: under actual times the schedule token is
  # unique per subject, so the verbatim skeleton copy is identifying.
  expect_equal(attr(report, "n_alone"), 12L)
  expect_true(all(report$obs_time_class == 1L))
  # The event signature does *not* include observation times, so it is blind to
  # this: one dose at a shared amount puts all twelve in one signature class.
  expect_equal(attr(report, "n_alone_signature"), 0L)
  # ... while the observation count is shared by everyone. That asymmetry is why
  # coarsening and the outlier screen do different jobs.
  expect_true(all(report$n_obs_class == 12L))
})

test_that("coarsening closes SIM-014 on AVATAR and leaving it off does not", {
  source <- crs_source()
  roles <- crs_roles()

  coarsened <- suppressWarnings(
    synpmx_avatar(source, roles, n_subjects = 20, seed = 1)
  )
  verbatim <- suppressWarnings(
    synpmx_avatar(source, roles, n_subjects = 20, seed = 1,
                  coarsen_time = FALSE)
  )

  expect_false(timing_vector_copied(source, coarsened, roles))
  # The defect this closes. If this ever passes, the skeleton stopped being
  # copied verbatim and the mechanism can be reconsidered.
  expect_true(timing_vector_copied(source, verbatim, roles))
  expect_true(validate_pmx(coarsened, roles)$valid)
})

test_that("time_jitter cannot substitute, at any magnitude", {
  source <- crs_source()
  roles <- crs_roles()
  # The Voronoi clamp holds each time within half a gap of the source value, so
  # jitter saturates instead of hiding the schedule. A generated vector stops
  # being byte-identical, but stays nearest to the subject it was copied from.
  worst_departure <- function(jitter_sd) {
    jittered <- suppressWarnings(
      synpmx_avatar(source, roles, n_subjects = 12, seed = 1,
                    coarsen_time = FALSE, time_jitter = jitter_sd)
    )
    src <- split(as.numeric(source$TIME), source$ID)
    gen <- split(as.numeric(jittered$TIME), jittered$ID)
    distances <- vapply(gen, function(g) {
      matches <- Filter(function(s) length(s) == length(g), src)
      if (!length(matches)) return(NA_real_)
      min(vapply(matches, function(s) max(abs(g - s)), numeric(1)))
    }, numeric(1))
    max(distances[is.finite(distances)])
  }

  # Every time is clamped into its own Voronoi cell, so a jittered vector can
  # never depart its source by more than half the widest gap in the schedule --
  # here 4 hours between the last two visits, so 2. The standard deviation is
  # therefore irrelevant past that point: a hundred and a hundred thousand give
  # the same bound. Whatever `time_jitter` is set to, the source schedule stays
  # recoverable, which is why it cannot be the privacy mechanism.
  bound <- max(diff(sort(unique(source$TIME)))) / 2 + 0.25
  expect_lt(worst_departure(100), bound)
  expect_lt(worst_departure(1e5), bound)
})

test_that("subjects sharing a schedule become exchangeable before generation", {
  source <- crs_source()
  roles <- crs_roles()
  coarsened <- synpmx:::.coarsen_source_time(source, roles)

  # The property the order of operations buys: after snapping, a subject's time
  # vector no longer depends on which subject it was. Anything downstream that
  # copies a skeleton is copying the class, not the patient.
  vectors <- tapply(as.numeric(coarsened$source$TIME),
                    as.character(coarsened$source$ID),
                    function(v) paste(format(v, digits = 12), collapse = ","))
  expect_equal(length(unique(vectors)), 1L)
  # And the removed deviations are kept, so realism can be restored.
  expect_true(stats::sd(coarsened$deviations) > 0)
})

test_that("a declared nominal_time is used as the grid, and TIME stays distinct", {
  source <- crs_source()
  roles <- crs_roles(nominal = "NTIME")
  synthetic <- suppressWarnings(
    synpmx_avatar(source, roles, n_subjects = 12, seed = 2)
  )
  settings <- attr(synthetic, "pmx_settings")
  expect_identical(settings$time_grid, "nominal")
  expect_true(settings$time_deviation_sd > 0)
  # Workflow coverage: code that reconciles actual against nominal time still has
  # a nonzero deviation to reconcile.
  expect_false(isTRUE(all.equal(synthetic$TIME, synthetic$NTIME)))
  expect_true(validate_pmx(synthetic, roles)$valid)
})

test_that("a source already on its nominal grid is generated unchanged", {
  # No deviation to remove means none to restore, and the RNG stream must not
  # move either -- the same guarantee `screen` gives a cohort with no outlier.
  source <- crs_source(sd = 0)
  roles <- crs_roles()
  on <- suppressWarnings(synpmx_avatar(source, roles, n_subjects = 10, seed = 5))
  off <- suppressWarnings(
    synpmx_avatar(source, roles, n_subjects = 10, seed = 5,
                  coarsen_time = FALSE)
  )
  # Data identical; only the recorded settings differ, by design.
  strip <- function(x) {
    x <- as.data.frame(x)
    attributes(x) <- attributes(x)[c("names", "row.names", "class")]
    x
  }
  expect_equal(strip(on), strip(off))
})

test_that("a structural outlier does not distort the derived grid", {
  # A fixed cluster count spends a cluster on the outlier's lone far-out visit
  # and under-resolves the schedule everyone else shares, manufacturing
  # deviation from an already-gridded cohort. The gap threshold does not.
  source <- rbind(crs_source(n = 8L, sd = 0),
                  crs_sub(99, c(0.5, 1, 2, 4, 500), 0))
  coarsened <- synpmx:::.coarsen_source_time(source, crs_roles())
  expect_equal(coarsened$deviations, numeric())
})

test_that("TAD is recomputed rather than inherited from the anchor", {
  source <- crs_source()
  source$TAD <- source$TIME
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", cmt = "CMT", covariates = "WT",
                     nominal_time = "NTIME", tad = "TAD")
  synthetic <- suppressWarnings(
    synpmx_avatar(source, roles, n_subjects = 10, seed = 3)
  )
  # TAD describes elapsed time since the last dose, so carrying the anchor's
  # value over would describe a schedule the avatar no longer has.
  # The dose row itself is offset by the resampled deviation, so TAD is measured
  # from the avatar's own dose time -- which is the point of recomputing it.
  by_subject <- split(synthetic, as.character(synthetic$ID))
  for (one in by_subject) {
    dose_time <- one$TIME[one$EVID != 0][1L]
    expect_equal(one$TAD, pmax(one$TIME - dose_time, 0), tolerance = 1e-8)
  }
  expect_true(validate_pmx(synthetic, roles)$valid)
})

test_that("coarsen_time is validated and recorded", {
  source <- crs_source()
  roles <- crs_roles()
  expect_error(
    synpmx_avatar(source, roles, coarsen_time = "yes"),
    "`coarsen_time` must be TRUE or FALSE"
  )
  settings <- attr(
    suppressWarnings(synpmx_avatar(source, roles, n_subjects = 5, seed = 1)),
    "pmx_settings"
  )
  expect_true(settings$coarsen_time)
  expect_identical(settings$time_grid, "derived")
})
