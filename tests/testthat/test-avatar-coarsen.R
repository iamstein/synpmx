# synpmx_avatar(coarsen_time = TRUE) is the default. It exists to close SIM-014
# ("Exact source timing vectors could be copied") on the AVATAR engine, where the
# gate had only ever been enforced against the structural/DP path even though
# AVATAR is the mode that actually copies an event skeleton verbatim.
#
# The mechanism is a two-step that only works in one order: snap every subject
# onto a shared visit grid, which is many-to-one and *destroys* the per-subject
# deviation, and only then resample deviations pooled across the cohort. Adding
# noise to the original time instead would leave the original recoverable,
# because `.offset_unique_times()` clamps every time inside its own Voronoi
# cell. These tests pin the order, not just the outcome.

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

test_that("actual recorded times give every patient a unique schedule", {
  report <- skeleton_uniqueness(crs_source(), crs_roles())
  # The premise of the whole mechanism: under actual times the schedule token is
  # unique per subject, so the verbatim skeleton copy is identifying.
  expect_equal(attr(report, "n_unique_schedule"), 12L)
  expect_true(all(report$n_share_schedule == 1L))
  # The event signature does *not* include observation times, so it is blind to
  # this: one dose at a shared amount puts all twelve in one signature class.
  expect_equal(attr(report, "n_unique_dose_signature"), 0L)
  # ... while the observation count is shared by everyone. That asymmetry is why
  # coarsening and the outlier screen do different jobs.
  expect_true(all(report$n_share_obs_count == 12L))
})

test_that("the printed near-miss sentence survives an even patient count", {
  # An even number of patients alone on their schedule makes `median()` average
  # the middle pair, so the two medians are doubles and the `%d` that formatted
  # them was an error rather than a number: "invalid format '%d'; use format
  # %f, %e, %g or %a for numeric objects". Twelve patients, all unique.
  report <- skeleton_uniqueness(crs_source(), crs_roles())
  expect_equal(nrow(report) %% 2L, 0L)
  expect_output(print(report), "not necessarily far apart")
})

test_that("coarsening closes SIM-014 on AVATAR and leaving it off does not", {
  source <- crs_source()
  roles <- crs_roles()

  # `min_pattern_share = 1` keeps attendance sampling out of both arms. It
  # rewrites which visits an avatar has, so with it on neither arm reproduces a
  # source timing vector and the contrast this test exists to draw disappears
  # for reasons that have nothing to do with coarsening.
  coarsened <- suppressWarnings(
    synpmx_avatar(source, roles, n_subjects = 20, seed = 1,
                  min_pattern_share = 1)
  )
  verbatim <- suppressWarnings(
    synpmx_avatar(source, roles, n_subjects = 20, seed = 1,
                  coarsen_time = FALSE, min_pattern_share = 1)
  )

  expect_false(timing_vector_copied(source, coarsened, roles))
  # The defect this closes. If this ever passes, the skeleton stopped being
  # copied verbatim and the mechanism can be reconsidered.
  expect_true(timing_vector_copied(source, verbatim, roles))
  expect_true(validate_pmx(coarsened, roles)$valid)
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

# REV-031 -- `.derive_time_grid()` refuses any merge that would put one subject
# in a cell twice, but `.snap_to_grid()` assigns by nearest centre and knew
# nothing about that, so a time near a boundary could land in a neighbouring
# cell and undo it. The doubled visit then became its own attendance cell and
# was drawn by many avatars.
test_that("coarsening never collapses two of one subject's own visits", {
  # Two visits close enough that a nearest-centre snap pulls them together,
  # against a cohort whose spacing makes the surrounding merges legal.
  source <- do.call(rbind, lapply(1:10, function(subject) {
    drift <- (subject - 5) * 0.05
    time <- c(0, 1 + drift, 1.45 + drift, 4 + drift, 8 + drift)
    evid <- c(1L, 0L, 0L, 0L, 0L)
    data.frame(
      ID = subject, TIME = time,
      DV = c(0, 8, 7.5, 4, 1) * (1 + 0.05 * sin(subject)),
      AMT = ifelse(evid != 0, 100, 0), EVID = evid,
      CMT = c(1L, 2L, 2L, 2L, 2L), WT = 70 + subject,
      stringsAsFactors = FALSE
    )
  }))
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", cmt = "CMT", covariates = "WT")

  coarsened <- synpmx:::.coarsen_source_time(source, roles)$source
  observed <- coarsened$EVID == 0
  key <- paste(coarsened$ID[observed], coarsened$TIME[observed])
  expect_false(any(duplicated(key)))

  # Every subject keeps all five of its rows, and its observation times stay
  # distinct: the repair keeps recorded times rather than dropping a visit.
  per_subject <- table(coarsened$ID)
  expect_true(all(per_subject == 5L))
})

test_that("genuine duplicate records are left alone by the repair", {
  # Two observations recorded at the *same* time are not a collapse the grid
  # caused, so coarsening must not invent a difference between them.
  source <- do.call(rbind, lapply(1:10, function(subject) {
    evid <- c(1L, 0L, 0L, 0L)
    data.frame(
      ID = subject, TIME = c(0, 1, 1, 4),
      DV = c(0, 8, 8.1, 4) * (1 + 0.05 * sin(subject)),
      AMT = ifelse(evid != 0, 100, 0), EVID = evid,
      CMT = c(1L, 2L, 2L, 2L), WT = 70 + subject, stringsAsFactors = FALSE
    )
  }))
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID", cmt = "CMT", covariates = "WT")

  coarsened <- synpmx:::.coarsen_source_time(source, roles)$source
  observed <- coarsened$EVID == 0
  key <- paste(coarsened$ID[observed], coarsened$TIME[observed])
  expect_equal(sum(duplicated(key)), 10L)
})

# A rename once left `sum(exposure$alone)` reading a column that no longer
# existed, so the reported exposure silently became 0 while the alerts still
# fired off the right number. Pin the settings against the screen they come
# from, so the two cannot drift apart again.
test_that("reported exposure agrees with skeleton_uniqueness on the coarsened source", {
  source <- crs_source()
  roles <- crs_roles()
  synthetic <- suppressWarnings(suppressMessages(
    synpmx_avatar(source, roles, seed = 4)
  ))
  settings <- attr(synthetic, "pmx_settings")

  # `coarsen_time = TRUE` is the public way to score the same grid the
  # generator builds, so this reads the settings against the argument a user
  # would reach for rather than against an internal.
  screen <- skeleton_uniqueness(source, roles, coarsen_time = TRUE)

  # A rename is only caught if the columns are named, not merely summed: a
  # missing column sums to 0 and would match a settings value of 0.
  expect_true(all(c("unique_schedule", "n_share_rarest_time") %in% names(screen)))
  expect_equal(settings$unique_schedule_n, sum(screen$unique_schedule))
  expect_equal(settings$unique_obs_time_n,
               sum(screen$n_share_rarest_time == 1L, na.rm = TRUE))
  # The split is exhaustive: a unique schedule is either a unique sample time
  # or a unique set of visits, never neither and never both.
  expect_equal(settings$unique_schedule_n,
               settings$unique_obs_time_n + settings$unique_visit_set_n)
  expect_gte(settings$unique_visit_set_n, 0L)
})

# SIM-038. `format(x, digits = 12)` lays out a whole vector at once, so the same
# time keyed differently depending on what else that patient was sampled at.
# Both fixtures below turn on exactly that: one patient carries a fractional
# time and the other does not, while the visits they share are identical.
sim038_source <- function() {
  # Two patients attend the same four visits. The first also has one PK sample
  # at a fractional time, which is what used to change how its *other* times
  # were formatted.
  one <- function(id, extra) {
    time <- sort(c(0, 24, 48, 72, extra))
    data.frame(
      ID = id, TIME = time, DV = c(NA, seq_along(time[-1L])),
      AMT = ifelse(time == 0, 100, 0),
      EVID = ifelse(time == 0, 1L, 0L),
      CMT = ifelse(time == 0, 1L, 2L), WT = 70 + id,
      stringsAsFactors = FALSE
    )
  }
  rbind(one(1L, 1.2142857142857), one(2L, numeric()), one(3L, numeric()))
}

sim038_roles <- function() {
  pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
            cmt = "CMT", covariates = "WT")
}

test_that("a time shared by other patients is never scored as a one-off", {
  screen <- skeleton_uniqueness(sim038_source(), sim038_roles())
  by_id <- split(screen, screen$subject_id)

  # Patients 2 and 3 hold only times that all three patients share, so neither
  # has a rarest time of 1. Before the fix, patient 1's fractional sample made
  # its copy of hour 24 key differently, splitting hour 24 into two keys of one
  # and two holders and reporting a one-off time for everybody.
  expect_equal(by_id[["2"]]$n_share_rarest_time, 3L)
  expect_equal(by_id[["3"]]$n_share_rarest_time, 3L)
  # Patient 1 genuinely is alone at its extra sample.
  expect_equal(by_id[["1"]]$n_share_rarest_time, 1L)
  expect_true(by_id[["1"]]$unique_schedule)
  expect_equal(by_id[["1"]]$why_unique, "one-off observation time")
})

test_that("patients attending the same visits share one attendance key", {
  source <- sim038_source()
  roles <- sim038_roles()
  rows <- split(seq_len(nrow(source)), source$ID)
  key <- function(id) {
    synpmx:::.attendance_key(source, roles,
                             seq_len(nrow(source)) %in% rows[[id]])
  }
  # The point of the fix: identical visit sets must produce identical keys, or
  # every pattern has one holder, `min_pattern_share` discards it, and the
  # avatar keeps its anchor's own absences.
  expect_identical(key("2"), key("3"))
  expect_false(identical(key("1"), key("2")))
})
