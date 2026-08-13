.validate_generator_options <- function(n_subjects, source_n, event_method,
                                        dv_method, k, pca_variance,
                                        subject_noise_sd, residual_noise_sd,
                                        residual_phi, time_jitter,
                                        max_donor_weight, coarsen_time,
                                        min_pattern_share) {
  if (is.null(n_subjects)) n_subjects <- source_n
  if (length(n_subjects) != 1L || is.na(n_subjects) ||
      n_subjects < 1 || n_subjects != as.integer(n_subjects)) {
    stop("`n_subjects` must be NULL or one positive integer.", call. = FALSE)
  }
  if (!identical(event_method, "template")) {
    stop("The prototype supports only `event_method = \"template\"`.",
         call. = FALSE)
  }
  if (!identical(dv_method, "avatar_blend")) {
    stop("The prototype supports only `dv_method = \"avatar_blend\"`.",
         call. = FALSE)
  }
  if (length(k) != 1L || is.na(k) || k < 1 || k != as.integer(k)) {
    stop("`k` must be one positive integer.", call. = FALSE)
  }
  if (length(pca_variance) != 1L || is.na(pca_variance) ||
      pca_variance <= 0 || pca_variance > 1) {
    stop("`pca_variance` must be in (0, 1].", call. = FALSE)
  }
  for (argument in c("subject_noise_sd", "residual_noise_sd", "time_jitter")) {
    value <- get(argument)
    if (length(value) != 1L || is.na(value) || !is.finite(value) || value < 0) {
      stop("`", argument, "` must be one finite nonnegative number.",
           call. = FALSE)
    }
  }
  if (length(residual_phi) != 1L || is.na(residual_phi) ||
      !is.finite(residual_phi) || abs(residual_phi) >= 1) {
    stop("`residual_phi` must be finite and strictly between -1 and 1.",
         call. = FALSE)
  }
  if (length(max_donor_weight) != 1L || is.na(max_donor_weight) ||
      !is.finite(max_donor_weight) || max_donor_weight <= 0 ||
      max_donor_weight > 1) {
    stop("`max_donor_weight` must be one number in (0, 1].", call. = FALSE)
  }
  if (length(coarsen_time) != 1L || is.na(coarsen_time) ||
      !is.logical(coarsen_time)) {
    stop("`coarsen_time` must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(min_pattern_share) != 1L || is.na(min_pattern_share) ||
      !is.finite(min_pattern_share) || min_pattern_share < 1 ||
      min_pattern_share != as.integer(min_pattern_share)) {
    stop("`min_pattern_share` must be one integer of 1 or more.", call. = FALSE)
  }
  as.integer(n_subjects)
}

# Move each *distinct* time in a skeleton by an offset, coherently: every row
# recorded at one time moves together, so a dose and its pre-dose sample stay
# tied and the row order can never invert. Offsets are clamped into each time's
# own Voronoi cell, which keeps visits from crossing or colliding.
#
# The clamp is why this cannot be a privacy control. Whatever the offsets are,
# every time stays within half a gap of where it started, so an offset applied
# to a *source* time leaves the source time recoverable. `coarsen_time` works
# because it moves the centre to a grid point shared by the whole equivalence
# class first; only then is the offset harmless. Order matters, and
# `.coarsen_source_time()` is what establishes it.
.offset_unique_times <- function(skeleton, roles, offsets) {
  time <- skeleton[[roles$time]]
  unique_time <- sort(unique(time))
  if (!length(unique_time)) return(skeleton)
  moved <- unique_time + offsets(length(unique_time))
  if (length(unique_time) > 1L) {
    midpoint <- (unique_time[-1L] + unique_time[-length(unique_time)]) / 2
    first_lower <- if (min(unique_time) >= 0) 0 else -Inf
    lower <- c(first_lower, midpoint)
    # The last cell is bounded by half the final gap, not left open. An open
    # upper bound let a large offset stretch follow-up arbitrarily, which is
    # exactly the structurally extreme output `screen = TRUE` exists to prevent
    # -- and it let one visit escape the "no time moves more than half a gap"
    # property every other visit obeys.
    last <- length(unique_time)
    upper <- c(midpoint, unique_time[last] + (unique_time[last] -
                                              unique_time[last - 1L]) / 2)
    moved <- pmin(pmax(moved, lower), upper)
  } else {
    moved <- max(0, moved)
  }
  skeleton[[roles$time]] <- moved[match(time, unique_time)]
  skeleton
}

.jitter_skeleton_time <- function(skeleton, roles, time_jitter) {
  if (time_jitter == 0) return(skeleton)
  .offset_unique_times(skeleton, roles, function(n) {
    stats::rnorm(n, sd = time_jitter)
  })
}

# Nearest-grid-point snap. `findInterval()` keeps this O(n log k) and monotone,
# so a subject's times can never reorder.
.snap_to_grid <- function(x, grid) {
  grid <- sort(unique(grid[is.finite(grid)]))
  if (!length(grid)) return(x)
  if (length(grid) == 1L) return(rep(grid, length(x)))
  index <- findInterval(x, grid, all.inside = TRUE)
  lower <- grid[index]
  upper <- grid[index + 1L]
  ifelse(abs(x - lower) <= abs(upper - x), lower, upper)
}

# The visit grid, when no `nominal_time` role names one. Best-effort: it
# recovers a shared schedule where the data contains one, and declines to invent
# one where it does not. `nominal_time` is the reliable path, and
# `scripts/measure_skeleton_uniqueness.R` is how you tell which case you are in.
#
# Pooled times are merged into visit cells agglomeratively, closest pair first,
# and a merge is taken only when the resulting cell holds no subject twice. That
# constraint is the definition of "too coarse": collapsing two of one subject's
# own visits into a single time destroys a real visit rather than hiding a
# schedule.
#
# The merge decision is *local*, which is the whole point. Two earlier rules
# used a single global threshold -- a fraction of the smallest within-subject
# gap, then the largest threshold satisfying the constraint everywhere -- and
# both are hostage to the tightest pair of samples anywhere in the study. One
# subject with two draws ten minutes apart at hour 1 drove the global threshold
# to ten minutes and left the hour-24 visits, spread across a whole day, sitting
# unmerged. Deciding each boundary on its own pins only the cells around that
# tight pair, so a dense early phase and a sparse late one can be gridded at the
# resolution each actually has.
.derive_time_grid <- function(times, subject_index) {
  keep <- is.finite(times)
  times <- times[keep]
  subject_index <- subject_index[keep]
  distinct <- sort(unique(times))
  n <- length(distinct)
  if (n < 2L) return(distinct)

  position <- match(times, distinct)
  members <- vector("list", n)
  observed_at <- split(subject_index, position)
  for (name in names(observed_at)) {
    members[[as.integer(name)]] <- unique(observed_at[[name]])
  }

  # Disjoint subject sets are not enough on their own. An outlier's lone
  # observation at hour 500 shares no subject with the hour-4 visit, so nothing
  # would stop the two merging and dragging every hour-4 observation to their
  # common mean. A cell is meant to be one visit, and two visits of a typical
  # subject sit a typical spacing apart, so no boundary wider than that can be
  # internal to a visit.
  #
  # The median, deliberately, not the minimum. A minimum is hostage to the
  # tightest pair of samples anywhere in the study -- the same fragility that
  # made a global threshold useless -- and this guard does not need to be tight.
  # Distinguishing adjacent visits is the subject-disjointness check's job; all
  # this has to do is refuse the absurd.
  spacing <- unlist(lapply(split(times, subject_index), function(v) {
    diff(sort(unique(v)))
  }), use.names = FALSE)
  spacing <- spacing[is.finite(spacing) & spacing > 0]
  widest <- if (length(spacing)) stats::median(spacing) else Inf

  # Union-find over cells. Merges are only ever between neighbours, so cells
  # stay contiguous runs of `distinct` and the roots need no ordering.
  parent <- seq_len(n)
  root <- function(i) {
    while (parent[i] != i) i <- parent[i]
    i
  }

  gaps <- diff(distinct)
  for (boundary in order(gaps)) {
    left <- root(boundary)
    right <- root(boundary + 1L)
    if (left == right) next
    # Blocked boundaries stay blocked: cells only grow, so two that share a
    # subject now will still share one later. No need to revisit.
    if (length(intersect(members[[left]], members[[right]]))) next
    if (gaps[boundary] > widest) next
    parent[right] <- left
    members[[left]] <- c(members[[left]], members[[right]])
    # `members[[right]] <- NULL` would *delete* the element and shift every
    # index after it; the cell is unreachable now, so empty it in place.
    members[right] <- list(NULL)
  }

  cell <- vapply(seq_len(n), root, integer(1))
  if (length(unique(cell)) == n) return(distinct)
  # Centre each visit on the mean of the times actually recorded there, so a
  # visit sampled by many patients is not pulled by a sparse neighbour.
  as.numeric(tapply(times, cell[position], mean))
}

# Collapse every subject onto a shared visit grid, and keep the deviations that
# were removed so they can be resampled back as a pool.
#
# This is what makes an avatar's schedule non-identifying. `.event_signature()`
# keys on dose gaps and amounts, so under actual recorded times almost every
# subject is alone in its own equivalence class (`skeleton_uniqueness()` shows
# this directly) and the verbatim skeleton copy at generation reproduces one real
# patient's visit schedule exactly. Snapping is many-to-one: it destroys the
# per-subject deviation rather than perturbing it, so subjects that shared a
# protocol schedule become genuinely exchangeable before any avatar is built.
#
# The removed deviations are pooled across the whole cohort and resampled
# independently afterwards. A single deviation value is weakly identifying in the
# same sense a resampled covariate value is (`design/METHOD_DISCUSSION.md` s4),
# while the vector an avatar receives is a fresh combination. Resampling them is
# what keeps `TIME` distinct from `NTIME` in the output, so workflow code that
# reconciles the two still has something to reconcile.
.coarsen_source_time <- function(source, roles) {
  time <- suppressWarnings(as.numeric(source[[roles$time]]))
  finite <- is.finite(time)
  none <- list(source = source, deviations = numeric(), grid = "none",
               grid_derived_rows = 0L, grid_rows = 0L)
  if (!any(finite)) return(none)

  target <- rep(NA_real_, length(time))
  kind <- "nominal"
  if (!is.null(roles$nominal_time)) {
    nominal <- suppressWarnings(as.numeric(source[[roles$nominal_time]]))
    usable <- finite & is.finite(nominal)
    target[usable] <- nominal[usable]
  }
  outstanding <- finite & !is.finite(target)
  if (any(outstanding)) {
    kind <- if (all(outstanding[finite])) "derived" else "mixed"
    grid <- .derive_time_grid(time[outstanding],
                              as.character(source[[roles$id]])[outstanding])
    target[outstanding] <- .snap_to_grid(time[outstanding], grid)
  }
  if (!any(is.finite(target))) return(none)

  # `.derive_time_grid()` refuses any merge that would put one subject into a
  # cell twice, because collapsing two of that subject's own visits deletes a
  # real visit rather than hiding a schedule. `.snap_to_grid()` then assigns by
  # *nearest centre*, which knows nothing about that constraint, so a time near
  # a boundary can land in a neighbouring cell and undo it (`REV-031`). The
  # damage is not the collision itself but what happens downstream: the doubled
  # visit becomes a distinct attendance cell, enters the pool, and is drawn by
  # many avatars -- on `warfarin`, 4 collisions in the source became 35 in the
  # output. Any subject whose own visits would collide keeps its recorded times
  # instead, which is the outcome the merge had already decided on.
  #
  # Only the derived grid is repaired. A declared `nominal_time` legitimately
  # places two samples in one protocol slot -- two draws on the same study day
  # are one visit by definition -- and that is the user's statement about their
  # own study, not an artefact to undo.
  if (any(outstanding)) {
    endpoint <- .endpoint(source, roles)
    key <- paste(as.character(source[[roles$id]]), endpoint,
                 ifelse(.dose_rows(source, roles), "<dose>", "<obs>"),
                 format(target, digits = 12, trim = TRUE))
    collides <- outstanding & (duplicated(key) | duplicated(key, fromLast = TRUE))
    # Only a genuine collapse counts: rows already recorded at the same time
    # were never distinct visits and are left alone.
    original <- paste(as.character(source[[roles$id]]), endpoint,
                      format(time, digits = 12, trim = TRUE))
    spurious <- collides &
      !(duplicated(original) | duplicated(original, fromLast = TRUE))
    target[spurious] <- NA_real_
  }
  if (!any(is.finite(target))) return(none)

  deviations <- (time - target)[finite & is.finite(target)]
  deviations <- deviations[is.finite(deviations)]
  # A source already on its nominal grid has nothing to remove and nothing to
  # restore. Emptying the pool rather than resampling zeros keeps the RNG stream
  # untouched too, so such a source generates byte-identical output whether
  # `coarsen_time` is on or off -- the same guarantee `screen` gives a cohort
  # with no structural outlier.
  if (!any(abs(deviations) > sqrt(.Machine$double.eps))) {
    deviations <- numeric()
  }
  snapped <- time
  snapped[finite & is.finite(target)] <- target[finite & is.finite(target)]
  source[[roles$time]] <- snapped
  list(source = source, deviations = deviations, grid = kind,
       # "mixed" on its own tells a caller who declared `nominal_time` that
       # something fell through, but not how much, and the answer decides
       # whether it matters. These two say how many rows were snapped to an
       # inferred grid rather than to the declared one.
       grid_derived_rows = sum(outstanding & is.finite(target)),
       grid_rows = sum(finite & is.finite(target)))
}

# The stratum a subject was assigned to: the combination of every
# `strata` column, or one shared stratum when none is declared.
# Protocol properties -- what dose level the arm received, which visits the
# schedule called for -- are constant within a stratum and differ between them,
# so it is the right grouping for both the dose basis below and the attendance
# pattern pool. It is deliberately *not* a blending barrier; only `.route_key()`
# is.
.subject_strata <- function(data, roles) {
  if (!length(roles$strata)) {
    return(rep("all", nrow(data)))
  }
  do.call(paste, c(lapply(roles$strata, function(column) {
    as.character(data[[column]])
  }), list(sep = "\r")))
}

.relative_spread <- function(values) {
  values <- values[is.finite(values)]
  if (length(values) < 2L) return(0)
  centre <- mean(values)
  if (!is.finite(centre) || abs(centre) < .Machine$double.eps) return(Inf)
  stats::sd(values) / abs(centre)
}

# Is the dose a fixed multiple of a baseline covariate -- mg/kg, mg/m^2 -- within
# each assigned stratum?
#
# This matters twice over. An avatar's `AMT` is copied verbatim from its anchor
# while its covariates are blended from donors, so under weight-based dosing the
# avatar's own mg/kg is wrong: source subjects all sit at exactly 5 mg/kg and
# avatars come out anywhere from 4.4 to 5.3, violating the protocol every
# generated subject claims to follow. And the copied amount is one real
# patient's real dose, which under proportional dosing discloses that patient's
# weight exactly -- the reason `skeleton_uniqueness()` reports nearly every
# subject as unique on dose (`REV-025` found the same thing from the other end,
# where it made every subject its own donor group).
#
# The fix is the decomposition that made time coarsening work. The multiplier is
# a protocol property shared by the whole stratum; the covariate is the
# individual part, and it is already blended. So keep the multiplier, recompute
# the amount from the avatar's own covariate, and the dose becomes both coherent
# and non-identifying.
#
# Detection is deliberately conservative, because recomputing amounts on a study
# that is not dose-proportional would be worse than leaving them alone. Dividing
# by the covariate has to *collapse* the variation: a handful of ratios standing
# in for many distinct amounts. See `.dose_ratio_levels()` for why the levels are
# clustered rather than averaged within a declared group.
#
# Detection is silent when it fails, and "silent" and "there is nothing here"
# look identical from the outside -- a reader of a run report cannot tell
# whether their weight-based study was recognised or whether the test refused
# it. So every exit records a one-line reason in `note`, and the run report
# prints it. `NULL` still means "leave the amounts alone"; the note says why.
.dose_basis_none <- function(note) list(covariate = NULL, note = note)

# Declared beats inferred. `.detect_dose_basis()` has to *guess* that a study is
# dose-proportional, and it guesses conservatively: the dose-to-covariate ratio
# has to collapse onto a handful of protocol levels. A study that rounds doses
# to vial sizes, or escalates within a patient by ratios that are not quite
# equal, fails that test and has its amounts left alone -- correctly, because
# rewriting amounts on a study that is not proportional would be worse. But a
# caller who *knows* their study is weight-based should not have to argue with
# an inference engine, so `dose_covariate` says it outright.
#
# The declared path also keeps each dose row's OWN ratio rather than snapping to
# clustered levels, so intra-patient escalation survives exactly: three doses at
# 1, 2 and 4 mg/kg stay at 1, 2 and 4 mg/kg even though no clustering step would
# have recovered those three levels from rounded milligrams.
.dose_basis_for <- function(source, roles) {
  declared <- roles$dose_covariate
  if (is.null(declared)) return(.detect_dose_basis(source, roles))
  if (is.null(roles$amt)) {
    stop("`dose_covariate` is declared but no `amt` role is, so there is no ",
         "amount to recompute.", call. = FALSE)
  }
  if (!declared %in% names(source)) {
    stop("`dose_covariate` names a column that is not in the data: `",
         declared, "`.", call. = FALSE)
  }
  amount <- suppressWarnings(as.numeric(source[[roles$amt]]))
  values <- suppressWarnings(as.numeric(source[[declared]]))
  usable <- .dose_rows(source, roles) & is.finite(amount) & amount > 0
  if (!any(usable)) {
    return(.dose_basis_none(sprintf(
      "`dose_covariate = \"%s\"` is declared, but no dose row carries a positive amount",
      declared)))
  }
  bad <- usable & (!is.finite(values) | values <= 0)
  if (any(bad)) {
    stop(sprintf(paste("`dose_covariate` (`%s`) is missing or non-positive on",
                       "%d dose row(s). It must be positive wherever a dose is",
                       "given, since the amount is divided by it."),
                 declared, sum(bad)), call. = FALSE)
  }
  ratio <- amount[usable] / values[usable]
  list(
    covariate = declared,
    # NULL levels is what marks the declared path for `.dose_levels_for()`:
    # hold the anchor's own per-row ratio instead of snapping to a level.
    levels = NULL,
    note = sprintf(paste("declared with `dose_covariate`, so no inference was",
                         "run. Each avatar's amount is rebuilt from its own",
                         "`%s` at its anchor's dose per unit (%s to %s across",
                         "the source)"),
                   declared, signif(min(ratio), 4), signif(max(ratio), 4))
  )
}

.detect_dose_basis <- function(source, roles, tolerance = 0.02,
                               max_levels = 10L) {
  if (is.null(roles$amt)) {
    return(.dose_basis_none("no `amt` role is declared, so there is no amount to explain"))
  }
  if (!length(roles$covariates)) {
    return(.dose_basis_none("no `covariates` are declared, so there is nothing to test the amounts against"))
  }
  amount <- suppressWarnings(as.numeric(source[[roles$amt]]))
  usable <- .dose_rows(source, roles) & is.finite(amount) & amount > 0
  if (sum(usable) < 3L) {
    return(.dose_basis_none(sprintf(
      "only %d dose row(s) carry a positive amount; at least 3 are needed",
      sum(usable))))
  }
  amount <- amount[usable]
  distinct_amounts <- length(unique(signif(amount, 8)))
  # Flat dosing: nothing varies, so nothing to explain and nothing identifying
  # about the amount.
  if (distinct_amounts < 2L) {
    return(.dose_basis_none(paste(
      "every dose row has the same amount, so nothing needed explaining and",
      "the amount identifies nobody")))
  }

  declined <- character()
  for (covariate in roles$covariates) {
    if (is.factor(source[[covariate]])) {
      declined <- c(declined, paste0(covariate, " (categorical)"))
      next
    }
    values <- suppressWarnings(as.numeric(source[[covariate]]))
    if (!is.numeric(values)) {
      declined <- c(declined, paste0(covariate, " (not numeric)"))
      next
    }
    values <- values[usable]
    if (any(!is.finite(values) | values <= 0)) {
      declined <- c(declined, paste0(covariate, " (missing or non-positive on dose rows)"))
      next
    }

    levels <- .dose_ratio_levels(amount / values, tolerance)
    if (is.null(levels)) {
      declined <- c(declined, paste0(covariate, " (ratios do not cluster)"))
      next
    }
    # Ratios must be far more concentrated than the amounts they came from.
    # A protocol has a handful of dose levels and many patients, so genuine
    # proportional dosing collapses dozens of distinct amounts onto a few
    # ratios; a study where dose is unrelated to the covariate produces about as
    # many ratios as amounts and is refused here.
    if (length(levels) > max_levels || length(levels) * 2L > sum(usable) ||
        distinct_amounts < 2L * length(levels)) {
      declined <- c(declined, sprintf(
        "%s (%d ratio levels for %d distinct amounts -- too many to be a protocol)",
        covariate, length(levels), distinct_amounts))
      next
    }
    return(list(covariate = covariate, levels = levels, note = sprintf(
      "the %d distinct dose amounts are a fixed multiple of `%s`, at %d protocol level(s)",
      distinct_amounts, covariate, length(levels))))
  }
  .dose_basis_none(sprintf(
    paste("the %d distinct dose amounts are not a fixed multiple of any",
          "declared covariate: %s"),
    distinct_amounts, paste(declined, collapse = "; ")))
}

# Group the observed dose-to-covariate ratios into protocol levels: sort, and
# cut wherever consecutive ratios differ by more than `tolerance` in relative
# terms. Returns the level values, or NULL when a cluster is too loose to be one
# level.
#
# Clustering rather than averaging within a declared stratum is what lets this
# see *intra-patient* escalation. A stratum is constant within subject, so it
# cannot tell a subject's 1 mg/kg dose from that same subject's later 2 mg/kg
# dose, and the pooled ratio is not constant, so the stratum-based test failed
# closed on exactly the design most likely to need it. A ratio does not care
# which subject or which occasion it came from.
.dose_ratio_levels <- function(ratio, tolerance) {
  ratio <- ratio[is.finite(ratio) & ratio > 0]
  if (!length(ratio)) return(NULL)
  sorted <- sort(ratio)
  if (length(sorted) == 1L) return(sorted)
  midpoint <- (sorted[-1L] + sorted[-length(sorted)]) / 2
  cluster <- cumsum(c(TRUE, abs(diff(sorted)) / midpoint > tolerance))
  levels <- vapply(split(sorted, cluster), mean, numeric(1))
  spread <- vapply(split(sorted, cluster), .relative_spread, numeric(1))
  if (any(!is.finite(spread) | spread > tolerance)) return(NULL)
  unname(levels)
}

# Which protocol level each dose row sits at, read from the anchor *before* its
# covariates are blended away. Kept separate from application because the level
# is a property of the anchor and the amount is a property of the avatar.
.dose_levels_for <- function(skeleton, roles, basis) {
  if (is.null(basis)) return(NULL)
  covariate <- suppressWarnings(as.numeric(
    .first_present(skeleton[[basis$covariate]])
  ))
  if (!is.finite(covariate) || covariate <= 0) return(NULL)
  amount <- suppressWarnings(as.numeric(skeleton[[roles$amt]]))
  dosed <- .dose_rows(skeleton, roles) & is.finite(amount) & amount > 0
  if (!any(dosed)) return(NULL)
  level <- rep(NA_real_, length(amount))
  ratio <- amount[dosed] / covariate
  level[dosed] <- if (is.null(basis$levels)) {
    # Declared basis: this row's own dose per unit, so a titration or an
    # escalation is carried through as it stands rather than rounded to the
    # nearest level some clustering step happened to find.
    ratio
  } else {
    basis$levels[max.col(-abs(outer(ratio, basis$levels, "-")),
                         ties.method = "first")]
  }
  level
}

# Recompute the dose from the avatar's own blended covariate, holding the level
# the anchor was dosed at. Any declared RATE scales with the amount so the
# infusion duration the protocol specified is preserved.
.apply_dose_basis <- function(skeleton, roles, basis, level) {
  if (is.null(basis) || is.null(level)) return(skeleton)
  covariate <- suppressWarnings(as.numeric(
    .first_present(skeleton[[basis$covariate]])
  ))
  if (!is.finite(covariate) || covariate <= 0) return(skeleton)
  amount <- suppressWarnings(as.numeric(skeleton[[roles$amt]]))
  if (length(level) != length(amount)) return(skeleton)
  change <- is.finite(level) & is.finite(amount) & amount > 0
  if (!any(change)) return(skeleton)

  target <- level[change] * covariate
  if (!is.null(roles$rate)) {
    rate <- suppressWarnings(as.numeric(skeleton[[roles$rate]]))
    scalable <- change & is.finite(rate) & rate > 0
    if (any(scalable)) {
      rate[scalable] <- rate[scalable] *
        (level[scalable] * covariate / amount[scalable])
      skeleton[[roles$rate]] <- rate
    }
  }
  amount[change] <- target
  skeleton[[roles$amt]] <- amount
  skeleton
}

# Attendance patterns ---------------------------------------------------------
#
# Once `coarsen_time` has put every subject on a shared visit grid, what is left
# of a subject's schedule is a bitmap: which of the shared visits they actually
# attended, per endpoint. On `warfarin` the grid holds 15 visits and 32 subjects
# produce 14 distinct bitmaps -- one covering 18 subjects, the rest held by one
# or two people each. Every column is shared, so no single time is identifying;
# the *combination of absences* is. A patient who missed weeks 2 and 3 is
# singled out by a fingerprint made of gaps.
#
# No grid can fix that, at any resolution. Merging two columns to blur a bitmap
# would delete a real visit from everybody who attended both, which is exactly
# what the disjointness constraint in `.derive_time_grid()` refuses. The grid
# decides where the columns are; it has no say in which cells hold a 1.
#
# So sample the bitmap instead of copying the anchor's, drawing only from
# patterns that at least `min_pattern_share` subjects hold. No avatar then wears
# a one-of-a-kind attendance pattern, and -- unlike dropping the exposed
# subjects -- nobody leaves the cohort: a subject with a rare pattern still
# contributes their measurements as a donor, only their distinctive absences
# stop being reproduced.
#
# Whole patterns are sampled, not individual bits. Drawing each visit
# independently would preserve the marginal attendance rate but destroy the
# correlation that makes an ending an ending -- once a subject stops they are
# gone, and independent draws produce implausible attend/miss/attend sequences.
# Copying a real pattern that at least `min_pattern_share` subjects share keeps
# that structure exactly, and sampling frequency-weighted keeps the per-visit
# attendance rates close to the source.
#
# Dose events are untouched. All the residual exposure is in observations, and
# sampling dose events risks emitting a regimen the protocol never permitted --
# the objection that ruled out randomized dose dropping.

.attendance_key <- function(data, roles, rows) {
  observed <- rows & .observation_rows(data, roles, require_present = TRUE)
  if (!any(observed)) return("")
  endpoint <- .endpoint(data, roles)[observed]
  time <- suppressWarnings(as.numeric(data[[roles$time]][observed]))
  # `.time_key()` and not `format(digits = 12)`: this key is built one subject
  # at a time and then compared *across* subjects, and `format()` lays out a
  # whole vector at once. Hour 12 keyed as "12" for a subject sampled only on
  # the hour and as "12.00000000000" for a subject who also has a 1.2142857
  # sample, so two patients who attended exactly the same visits produced
  # different keys. Every such pattern then had a single holder, was discarded
  # by `min_pattern_share`, and its avatar fell back to the anchor's own set of
  # visits -- the outcome the mechanism exists to prevent.
  cells <- paste0(endpoint, "@", .time_key(time))
  # A subject can hold two observation rows with the same endpoint at the same
  # time -- a duplicate record, or two compartments read at one visit where no
  # `dvid` tells them apart. Collapsing them to one cell meant rebuilding the
  # pattern returned fewer rows than the subject had, silently (`SIM-036`).
  # Numbering repeats past the first keeps them distinct; a lone occurrence
  # keeps its bare cell, so a source without repeats produces the keys it always
  # did.
  repeated <- stats::ave(seq_along(cells), cells, FUN = seq_along)
  extra <- repeated > 1L
  cells[extra] <- paste0(cells[extra], "#", repeated[extra])
  paste(sort(unique(cells)), collapse = ";")
}

# The time a cell names, with any repeat marker removed.
.cell_time <- function(cells) {
  suppressWarnings(as.numeric(sub("#.*$", "", sub("^.*@", "", cells))))
}

# Which occurrence at that endpoint and time a cell names. 1 unless marked.
.cell_repeat <- function(cells) {
  marked <- grepl("#", cells, fixed = TRUE)
  out <- rep(1L, length(cells))
  out[marked] <- as.integer(sub("^.*#", "", cells[marked]))
  out
}

# The pool an avatar may draw from, per stratum. Patterns are kept only when
# enough subjects share them; a stratum with nothing shared widely enough
# contributes no pool and its anchors keep their own pattern.
.attendance_pool <- function(source, roles, profiles, min_pattern_share) {
  if (min_pattern_share <= 1L) return(NULL)
  subject_rows <- profiles$subject_rows
  keys <- vapply(subject_rows, function(rows) {
    .attendance_key(source, roles, rows)
  }, character(1))
  strata <- vapply(subject_rows, function(rows) {
    .subject_strata(source, roles)[which(rows)[1L]]
  }, character(1))
  # An avatar must not be given a pattern naming endpoints its anchor does not
  # have, so the pool is separated by endpoint set as well as by stratum.
  endpoints <- vapply(subject_rows, function(rows) {
    observed <- rows & .observation_rows(source, roles, require_present = TRUE)
    paste(sort(unique(.endpoint(source, roles)[observed])), collapse = ",")
  }, character(1))
  group <- paste(strata, endpoints, sep = "\r")

  pool <- list()
  total_patterns <- 0L
  dropped_patterns <- 0L
  dropped_subjects <- 0L
  for (name in unique(group)) {
    member <- group == name
    counts <- table(keys[member])
    counts <- counts[names(counts) != ""]
    if (!length(counts)) next
    total_patterns <- total_patterns + length(counts)
    kept <- counts >= min_pattern_share

    # The grid this group's subjects were observed on, in time order, and each
    # distinct pattern expressed as a shape over it. `shape` is what lets a
    # pattern that only one patient holds still contribute: two patients who each
    # missed one visit are the same shape even when they missed different visits,
    # so together they clear a floor that neither clears alone.
    cells <- .attendance_cells(names(counts))
    droppable <- .droppable_cells(cells, names(counts))
    max_depth <- .shared_depth(cells, names(counts), as.integer(counts),
                               min_pattern_share)
    shapes <- vapply(names(counts), function(key) {
      .attendance_shape(key, cells, droppable)
    }, character(1))
    shape_counts <- tapply(as.integer(counts), shapes, sum)
    shape_kept <- shape_counts[shape_counts >= min_pattern_share]

    # A third, coarser abstraction, reached only when the second cleared
    # nothing: the KIND of missingness alone, ignoring how many visits it cost.
    #
    # Staggered discontinuation is the case that needs it, and it is the normal
    # shape of a real study: twenty-one patients who each stopped at a different
    # visit share no exact pattern AND no (kind, count) shape, because each
    # count has exactly one holder. The group then got no pool at all and every
    # avatar in it kept its anchor's own visit set -- one real patient's
    # absences, copied. "These patients dropped out" is plainly a property
    # twenty-one of them share, so it is enough to build from; the placement
    # search picks a depth that is not any single patient's.
    # A pattern survives if it can be reproduced itself, or if its shape can.
    carried <- kept | (shapes %in% names(shape_kept))
    if (!length(shape_kept)) {
      kinds <- sub("\\|.*$", "", shapes)
      kind_counts <- tapply(as.integer(counts), kinds, sum)
      kind_kept <- kind_counts[kind_counts >= min_pattern_share]
      if (length(kind_kept)) {
        # Request the typical depth for that kind; `.place_attendance()` moves
        # it when the exact depth belongs to somebody.
        typical <- vapply(names(kind_kept), function(k) {
          max(1L, as.integer(round(stats::median(
            as.integer(sub("^.*\\|", "", shapes[kinds == k]))
          ))))
        }, integer(1))
        shape_kept <- kind_kept
        names(shape_kept) <- paste0(names(kind_kept), "|", typical)
        carried <- kinds %in% names(kind_kept)
        # `kept` is necessarily all FALSE here -- an exact pattern clearing the
        # floor would have carried its own shape through above -- so renaming
        # the shapes invalidates no `exact` mapping.
        shapes <- rep(NA_character_, length(shapes))
      }
    }

    # What is lost is what none of the three abstractions can carry.
    lost <- !carried
    dropped_patterns <- dropped_patterns + sum(lost)
    dropped_subjects <- dropped_subjects + sum(as.integer(counts[lost]))
    if (!length(shape_kept)) next

    pool[[name]] <- list(
      cells = cells,
      droppable = droppable,
      max_depth = max_depth,
      shapes = names(shape_kept), shape_weight = as.numeric(shape_kept),
      # Exact patterns that clear the floor on their own, indexed by shape, so a
      # real pattern is preferred wherever reusing one is safe.
      exact = split(names(counts)[kept], shapes[kept]),
      exact_weight = split(as.numeric(counts[kept]), shapes[kept]),
      # Patterns too rare to reuse. A generated placement that lands on one of
      # these is rejected, so "nothing was copied" does not quietly weaken into
      # "something was reproduced by accident".
      rare = names(counts)[!kept],
      # The widely held set to fall back ON when a draw finds nothing legal and
      # the anchor's own set is one of the rare ones. Copying THAT would be the
      # disclosure the mechanism exists to prevent, while copying a set six
      # patients share discloses nothing, so the substitution is strictly better
      # than the alternative and needs no permission.
      common = if (any(kept)) {
        names(counts)[kept][[which.max(as.integer(counts[kept]))]]
      } else NA_character_
    )
  }
  # Which SUBJECTS hold a visit set too rare to be reused. Computed for every
  # group, including the ones that got no pool, because the question "would
  # copying this anchor's set disclose anything?" has to be answerable even
  # where there is nothing to replace it with.
  identifying <- logical(length(keys))
  for (name in unique(group)) {
    member <- group == name
    shared <- table(keys[member])
    identifying[member] <- nzchar(keys[member]) &
      as.integer(shared[keys[member]]) < min_pattern_share
  }
  # A second pool keyed on the ENDPOINT SET alone, ignoring the stratum. Groups
  # are (stratum x endpoint set), and splitting by stratum is what makes a small
  # arm hold too few patients to share anything. Donors already cross strata to
  # reach the floor -- only route of administration is absolute -- so a visit set
  # may cross one too. The endpoint set may not: an avatar must never be handed a
  # set naming an endpoint its anchor does not have.
  endpoint_of <- function(name) sub("^.*\r", "", name)
  endpoint_pool <- list()
  for (set in unique(endpoint_of(names(pool)))) {
    candidates <- pool[endpoint_of(names(pool)) == set]
    # The one with the most widely held set is the safest to borrow from.
    best <- which.max(vapply(candidates, function(p) {
      if (is.na(p$common)) 0 else max(p$shape_weight)
    }, numeric(1)))
    endpoint_pool[[set]] <- candidates[[best]]
  }

  list(pool = pool, endpoint_pool = endpoint_pool, endpoint_of = endpoint_of,
       group = group, keys = keys, identifying = identifying,
       total_patterns = total_patterns, dropped_patterns = dropped_patterns,
       dropped_subjects = dropped_subjects)
}

# Which cells a generated placement is allowed to drop.
#
# An endpoint where every patient in the group holds the identical set of times
# has no missingness to reproduce -- it is complete by construction. Placing a
# miss on one of its cells does not mask anything; it invents missingness the study
# does not have. Measured on a two-endpoint fixture, a biomarker drawn at all
# ten visits for all twenty-one patients came out of generation with seven to
# ten, because the shape was placed over the POOLED grid in time order and took
# whatever cells sat late, biomarker and PK alike.
#
# So placements are confined to the endpoints that actually vary. Exact patterns
# reused from the pool are untouched by this: they are real patterns, and a real
# pattern never drops a cell of a complete endpoint anyway.
.droppable_cells <- function(cells, keys) {
  endpoint <- sub("@.*$", "", cells)
  held <- strsplit(keys, ";", fixed = TRUE)
  varies <- vapply(unique(endpoint), function(name) {
    mine <- cells[endpoint == name]
    sets <- vapply(held, function(h) sum(mine %in% h), integer(1))
    length(unique(sets)) > 1L ||
      any(vapply(held, function(h) !all(mine %in% h), logical(1)) &
          vapply(held, function(h) any(mine %in% h), logical(1)))
  }, logical(1))
  endpoint %in% names(which(varies))
}

# The union of every cell any subject in the group was observed at, in time
# order. Time first, so "trailing" means late in the study -- follow-up ending --
# rather
# than late in an alphabet.
.attendance_cells <- function(keys) {
  cells <- unique(unlist(strsplit(keys, ";", fixed = TRUE), use.names = FALSE))
  cells <- cells[nzchar(cells)]
  cells[order(.cell_time(cells), cells)]
}

# How deep into the study a *generated* visit set is allowed to reach: the last
# cell that `min_pattern_share` subjects were still being observed at.
#
# A shape says how many visits were missed and how they were arranged, not where
# on the calendar the survivors sat, and for `trailing` and `block` that is
# harmless -- both leave the tail of the grid attended, so neither can push an
# avatar past where its own kind of patient stops. `scattered` can. It is free
# to hold any cells at all, and the grid is the union over the whole group, so
# one patient followed far longer than the rest lends their tail to everybody.
#
# That is `SIM-047`, seen on `wbcSim`: the grid runs to 4580 h because a single
# subject was followed that long, and an avatar holding a two-visit `scattered`
# shape was placed on cells 501 h and 3910 h -- a follow-up length no real
# patient of that shape has, and a fresh signal that *somebody* was followed to
# 3910 h. Nothing else caught it: B1a and B4a ask whether a real patient's set
# was copied, and this set is precisely one no real patient holds.
#
# The floor is `min_pattern_share` and the reading is the usual one -- a
# generated set may reach a depth several real patients reach, and no further.
.shared_depth <- function(cells, keys, counts, min_pattern_share) {
  held <- strsplit(keys, ";", fixed = TRUE)
  depth <- vapply(held, function(h) {
    reached <- which(cells %in% h)
    if (!length(reached)) 0L else max(reached)
  }, integer(1))
  # Subjects, not distinct patterns: `counts` is how many subjects hold each.
  reach <- vapply(seq_along(cells), function(d) sum(counts[depth >= d]),
                  integer(1))
  deep_enough <- which(reach >= min_pattern_share)
  if (!length(deep_enough)) length(cells) else max(deep_enough)
}

# Describe a pattern by *how* it departs from full attendance rather than by
# which cells it holds: the number of missed visits and their arrangement.
#
#   complete   attended everything
#   trailing   every miss is at the end -- a discontinuation, or follow-up
#              that has not reached those visits yet
#   block      the misses are contiguous but not terminal -- an interruption
#   scattered  anything else
#
# This is the whole point of the fallback. The exact-cell definition cannot see
# that "missed one visit, early" and "missed one visit, late" are the same kind
# of event, so each is a group of one and both are discarded. As a shape they are
# one group of two.
.attendance_shape <- function(key, cells, droppable = NULL) {
  if (!is.null(droppable)) cells <- cells[droppable]
  held <- strsplit(key, ";", fixed = TRUE)[[1L]]
  missed <- which(!(cells %in% held))
  if (!length(missed)) return("complete|0")
  n <- length(cells)
  kind <- if (all(missed > n - length(missed))) {
    "trailing"
  } else if (all(diff(missed) == 1L)) {
    "block"
  } else {
    "scattered"
  }
  paste0(kind, "|", length(missed))
}

# The visit set an avatar anchored here would be given, and whether it would be
# identifying. Split out of the generation loop so the loop can ASK before it
# builds, and pick a different anchor when the answer is bad.
# Whose DOSE schedule nobody else shares.
#
# The attendance mechanism deliberately never touches dose events: resampling
# them risks emitting a regimen the protocol never permitted, which would be a
# worse error than the one it fixes. That leaves the dose side unguarded, and it
# is exposed in exactly the same way the observation side is -- the anchor's
# event skeleton is copied verbatim, so a patient whose set of dose times nobody
# else shares hands it to every avatar built on them.
#
# The answer is therefore not to rewrite the regimen but to decline to build on
# it, the same rule the route floor and the visit-set machinery already apply.
# Measured after coarsening, since that is the grid generation runs on: `theo_md`
# and `warfarin` and `mavoglurant` have one shared dose schedule each and are
# unaffected, `wbcSim` has 3 of 45 alone, and `nimoData` has 12 of 12 -- a study
# whose dosing genuinely cannot be masked, which the caller needs told.
# Keys for the end-to-end checks, computed on the shared visit grid rather than
# on recorded times.
#
# This is the whole reason the checks are trustworthy. Generation resamples the
# per-visit deviations coarsening removed, so an avatar sitting in exactly the
# same grid cell as a real patient does NOT hold the same numeric time -- and a
# check comparing the numbers therefore reports "nothing identifying" on a study
# where every schedule is identifying, which is worse than no check. Snapping
# both sides back onto the grid asks the question that matters: is this avatar
# in the same set of cells as some real patient?
# Dosing that stopped early, as a shape.
#
# A patient who stopped dosing after four of ten doses holds a schedule that is
# a PREFIX of the full one. That is the dose-side analogue of a trailing
# attendance shape, and it is the one dose edit that is protocol-valid: a
# regimen truncated at a real dose time is a regimen the study actually
# permitted, whereas moving or inventing dose times is not. So a truncation
# depth held by one patient can be exchanged for one that is held by several,
# or by nobody at all, exactly as `.miss_counts()` walks a visit-set depth.
#
# It cannot rescue every cohort, and the arithmetic says which. Nineteen
# patients on three doses, one on two and one on one has depths {1:1, 2:1,
# 3:19}: depths 1 and 2 have a single holder each and there is no free depth
# between them, so "stopped after one dose" cannot be represented without
# pointing at the person who did. Ten stopping depths across forty patients
# leaves plenty free, which is the case this exists for.
.dose_truncation_plan <- function(source, roles, min_pattern_share) {
  subjects <- as.character(.unique_in_order(source[[roles$id]]))
  empty <- list(depth = rep(NA_integer_, length(subjects)),
                target = rep(NA_integer_, length(subjects)),
                cells = numeric())
  if (!length(subjects) || is.null(roles$amt)) return(empty)
  time <- suppressWarnings(as.numeric(source[[roles$time]]))
  id <- as.character(source[[roles$id]])
  dosed <- .dose_rows(source, roles) & is.finite(time)
  if (!any(dosed)) return(empty)

  times <- lapply(subjects, function(s) sort(unique(time[dosed & id == s])))
  cells <- sort(unique(unlist(times)))
  # A schedule is a truncation only when it is the first `d` dose times of the
  # cohort's union. Anything else -- a missed middle dose, a delay -- is not a
  # prefix and is left to the anchor-level rule.
  depth <- vapply(times, function(t) {
    d <- length(t)
    if (!d || d > length(cells)) return(NA_integer_)
    if (isTRUE(all.equal(t, cells[seq_len(d)]))) as.integer(d) else NA_integer_
  }, integer(1))

  held <- table(depth[!is.na(depth)])
  count_at <- function(d) {
    if (is.na(d) || !as.character(d) %in% names(held)) 0L else
      as.integer(held[[as.character(d)]])
  }
  # Safe means shared by enough patients, or held by none at all: a depth
  # nobody stopped at identifies nobody.
  safe <- function(d) {
    d >= 1L && d <= length(cells) &&
      (count_at(d) == 0L || count_at(d) >= as.integer(min_pattern_share))
  }
  target <- vapply(seq_along(subjects), function(i) {
    d <- depth[[i]]
    if (is.na(d) || safe(d)) return(NA_integer_)   # nothing to change
    # Downward only. Truncation drops dose rows; it cannot invent one, so a
    # deeper target silently leaves the schedule where it was -- which put an
    # avatar back on the very depth that one patient alone had used.
    for (candidate in rev(seq_len(max(d - 1L, 0L)))) {
      if (safe(candidate)) return(as.integer(candidate))
    }
    NA_integer_
  }, integer(1))
  list(depth = depth, target = target, cells = cells)
}

# Drop the dose rows past `depth`, so the avatar stops dosing where the plan
# says rather than where its anchor did. Observation rows are untouched: what
# happens after dosing stops is the attendance mechanism's business.
.apply_dose_truncation <- function(skeleton, roles, cells, depth) {
  if (is.na(depth) || !length(cells)) return(skeleton)
  time <- suppressWarnings(as.numeric(skeleton[[roles$time]]))
  dosed <- .dose_rows(skeleton, roles) & is.finite(time)
  if (!any(dosed)) return(skeleton)
  cutoff <- cells[[min(depth, length(cells))]]
  drop <- dosed & time > cutoff + sqrt(.Machine$double.eps)
  # Never drop every dose: an avatar with no dosing at all is not a truncated
  # regimen, it is a different kind of record.
  if (!any(drop) || !any(dosed & !drop)) return(skeleton)
  skeleton[!drop, , drop = FALSE]
}

.grid_snapped_keys <- function(data, roles, grid, rows) {
  time <- suppressWarnings(as.numeric(data[[roles$time]]))
  keep <- rows & is.finite(time)
  if (!any(keep)) return(character())
  snapped <- .snap_to_grid(time[keep], grid)
  vapply(split(snapped, as.character(data[[roles$id]])[keep]),
         function(x) paste(.time_key(sort(unique(x))), collapse = ","),
         character(1))
}

# Each patient's dose schedule as a text key, on whatever grid the data is on.
.dose_schedule_keys <- function(source, roles) {
  subjects <- .unique_in_order(source[[roles$id]])
  time <- suppressWarnings(as.numeric(source[[roles$time]]))
  id <- as.character(source[[roles$id]])
  dosed <- .dose_rows(source, roles) & is.finite(time)
  keys <- vapply(split(time[dosed], id[dosed]),
                 function(t) paste(.time_key(sort(unique(t))), collapse = ","),
                 character(1))[as.character(subjects)]
  keys[is.na(keys)] <- ""
  stats::setNames(keys, as.character(subjects))
}

.dose_schedule_sharing <- function(source, roles, min_pattern_share) {
  keys <- .dose_schedule_keys(source, roles)
  if (!length(keys)) return(logical())
  shared <- table(keys[nzchar(keys)])
  out <- rep(FALSE, length(keys))
  present <- nzchar(keys)
  out[present] <- as.integer(shared[keys[present]]) < as.integer(min_pattern_share)
  stats::setNames(out, names(keys))
}

.attendance_choice <- function(anchor, attendance) {
  none <- list(key = NA_character_, generated = FALSE, shifted = FALSE,
               substituted = FALSE,
               identifying = isTRUE(attendance$identifying[[anchor]]))
  # The anchor's own schedule group first; failing that, any group measuring
  # the same endpoints, which is what rescues a small arm. The endpoint set may
  # never be crossed -- an avatar must not be handed a set naming an endpoint
  # its anchor does not have.
  available <- attendance$pool[[attendance$group[[anchor]]]]
  if (is.null(available)) {
    available <- attendance$endpoint_pool[[
      attendance$endpoint_of(attendance$group[[anchor]])
    ]]
  }
  if (is.null(available)) return(none)

  drawn <- .draw_attendance(available)
  if (!is.na(drawn$key)) {
    return(list(key = drawn$key, generated = drawn$generated,
                shifted = isTRUE(drawn$shifted), substituted = FALSE,
                identifying = FALSE))
  }
  # Nothing legal was placed. Keeping the anchor's own set is fine unless that
  # set is one nobody else has, in which case the group's most widely held set
  # goes in instead: less faithful to this avatar, and it discloses nothing.
  if (none$identifying && !is.na(available$common)) {
    return(list(key = available$common, generated = FALSE, shifted = FALSE,
                substituted = TRUE, identifying = FALSE))
  }
  none
}

.attendance_key_from <- function(cells, missed) {
  # `cells[-integer(0)]` is `cells[0]`, not `cells`, so "nothing was missed" has
  # to be said outright or a complete pattern comes back as the empty key.
  if (!length(missed)) return(paste(sort(cells), collapse = ";"))
  paste(sort(cells[-missed]), collapse = ";")
}

# Every way of placing `n_miss` misses of this kind over `n` cells, in the order
# they should be tried. `trailing` yields exactly one, which is the fact the
# caller is built around.
.attendance_placements <- function(kind, n_miss, n, tries) {
  if (n_miss < 1L || n_miss >= n) return(list())
  switch(
    kind,
    trailing = list(seq.int(n - n_miss + 1L, n)),
    block = {
      # Contiguous but not terminal, so a block stays distinguishable from
      # an ending. Enumerated and shuffled rather than sampled with replacement:
      # the same start was previously drawn several times out of 24 tries while
      # other legal starts went untried.
      last_start <- n - n_miss
      if (last_start < 1L) return(list())
      lapply(sample(seq_len(last_start)), function(start) {
        seq.int(start, start + n_miss - 1L)
      })
    },
    scattered = {
      # choose(n, n_miss) is far too large to enumerate, so this one stays a
      # sample. A draw that came out contiguous or terminal is a different
      # shape, so it is discarded rather than returned.
      out <- list()
      for (attempt in seq_len(tries)) {
        missed <- sort(sample.int(n, n_miss))
        if (all(diff(missed) == 1L) || all(missed > n - n_miss)) next
        out[[length(out) + 1L]] <- missed
      }
      out
    },
    list()
  )
}

# The requested number of misses first, then outward: k, k+1, k-1, k+2, ...
.miss_counts <- function(wanted, n) {
  outward <- as.vector(rbind(wanted + seq_len(n), wanted - seq_len(n)))
  candidates <- c(wanted, outward)
  unique(candidates[candidates >= 1L & candidates < n])
}

# Place a shape's misses on the grid without copying anyone's arrangement, and
# refuse any placement that happens to reproduce a pattern too rare to have been
# reusable. Without that rejection the guarantee would weaken from "no synthetic
# patient carries a schedule unique to a real one" to merely "nothing was copied
# on purpose", and an attacker cannot tell those apart.
#
# The number of misses is allowed to move if, and only if, nothing at the
# requested number is legal. `trailing` is why: an ending after exactly k visits
# has ONE possible placement, and if any real patient dropped out at that depth
# then that placement is rare and there is nothing else at that depth to try.
# The old loop retried the identical placement 24 times, failed, and the caller
# then let the avatar keep its ANCHOR's visit set -- copying one real patient's
# absences exactly, which is worse on every axis than an ending a visit deeper.
# Dropout is the commonest kind of missingness, so this was the normal outcome
# rather than an edge case: 86% of avatars on a real 21-patient study.
.place_attendance <- function(cells, shape, rare, tries = 24L,
                              droppable = NULL,
                              max_depth = length(cells)) {
  parts <- strsplit(shape, "|", fixed = TRUE)[[1L]]
  kind <- parts[[1L]]
  wanted <- as.integer(parts[[2L]])
  # Positions are chosen among the cells that may be dropped, then mapped back
  # onto the full grid, so a complete endpoint's cells are never candidates.
  index <- if (is.null(droppable)) seq_along(cells) else which(droppable)
  # Droppable cells past the depth real patients share (`.shared_depth()`).
  # They are not candidates for placement -- they are missed outright, and the
  # arrangement is worked out over what is left. Non-droppable cells beyond the
  # cap stay attended: a complete endpoint is complete by construction and
  # dropping one of its cells invents missingness the study does not have.
  beyond <- index[index > max_depth]
  index <- index[index <= max_depth]
  n <- length(index)
  # The deepest set a placement may return, and so what "attended everything"
  # can mean here.
  complete <- .attendance_key_from(cells, beyond)
  if (kind == "complete" || wanted < 1L) {
    return(list(key = complete, shifted = length(beyond) > 0L))
  }
  # `wanted` counts misses over every droppable cell. The capped ones are now
  # misses too, so what is left to arrange is the remainder -- which keeps the
  # number of visits the avatar ends up with unchanged, and only moves them
  # inside the study's shared follow-up.
  wanted <- max(wanted - length(beyond), 0L)

  for (n_miss in .miss_counts(wanted, n)) {
    for (position in .attendance_placements(kind, n_miss, n, tries)) {
      key <- .attendance_key_from(cells, c(beyond, index[position]))
      if (!(key %in% rare)) {
        return(list(key = key,
                    shifted = !identical(n_miss, wanted) ||
                      length(beyond) > 0L))
      }
    }
  }
  # Last resort before the caller falls back to the anchor's own visit set:
  # attend everything. Legal unless exactly one real patient was the sole
  # complete attender, in which case it is that patient's pattern.
  if (!(complete %in% rare)) return(list(key = complete, shifted = TRUE))
  list(key = NA_character_, shifted = FALSE)
}

# Draw one attendance pattern for an avatar. A real pattern is used wherever
# reusing one is safe; a placement is generated only for shapes whose every
# realisation is too rare to reuse. The shape is drawn first either way, so the
# cohort's mix of complete, interrupted, and dropped-out subjects is preserved
# even when the individual arrangements are not.
.draw_attendance <- function(available) {
  index <- sample.int(length(available$shapes), 1L, prob = available$shape_weight)
  shape <- available$shapes[[index]]
  exact <- available$exact[[shape]]
  if (length(exact)) {
    weight <- available$exact_weight[[shape]]
    return(list(key = exact[[sample.int(length(exact), 1L, prob = weight)]],
                generated = FALSE, shifted = FALSE))
  }
  placed <- .place_attendance(available$cells, shape, available$rare,
                              droppable = available$droppable,
                              max_depth = available$max_depth %||%
                                length(available$cells))
  list(key = placed$key, generated = !is.na(placed$key),
       shifted = isTRUE(placed$shifted))
}

# Rebuild the skeleton's observation rows to match a sampled pattern, keeping
# every dose row exactly as the anchor had it. Rows for a visit the anchor did
# not attend are cloned from one of its own observation rows for that endpoint,
# so compartment, MDV, and every carried column stay coherent; the DV itself is
# filled by the donor blend afterwards, exactly as for a copied skeleton.
#
# The clone is the anchor's *nearest* observation of that endpoint, not its
# first. Only `time` is rewritten here, so every other row-varying column --
# nominal time, occasion, a per-row `keep` label -- arrives from whichever visit
# it was attached to. Cloning one template row per endpoint stamped the first
# visit's metadata onto all of them, which is `SIM-035`.
.apply_attendance <- function(skeleton, roles, key) {
  observed <- .observation_rows(skeleton, roles, require_present = TRUE)
  if (!any(observed) || !nzchar(key)) return(skeleton)
  wanted <- strsplit(key, ";", fixed = TRUE)[[1L]]
  wanted_endpoint <- sub("@.*$", "", wanted)
  wanted_time <- .cell_time(wanted)
  wanted_repeat <- .cell_repeat(wanted)
  endpoint <- .endpoint(skeleton, roles)
  time <- suppressWarnings(as.numeric(skeleton[[roles$time]]))

  chosen <- vapply(seq_along(wanted), function(index) {
    candidate <- which(observed & endpoint == wanted_endpoint[[index]])
    if (!length(candidate)) return(NA_integer_)
    distance <- abs(time[candidate] - wanted_time[[index]])
    distance[!is.finite(distance)] <- Inf
    # The n-th occurrence at a visit takes the n-th nearest row, so two cells
    # that name the same endpoint and time clone two different source rows
    # rather than the same one twice.
    rank <- min(wanted_repeat[[index]], length(candidate))
    candidate[[order(distance)[[rank]]]]
  }, integer(1))
  # NULL, not the untouched skeleton. Returning the skeleton meant the caller
  # believed a visit set had been applied when none had, so the avatar quietly
  # kept its ANCHOR's set -- and if that set was one nobody shared, the run
  # emitted it while reporting the pattern as sampled. Reachable whenever the
  # drawn key names an endpoint this anchor does not have, which is exactly what
  # borrowing a set from another schedule group can do.
  if (anyNA(chosen)) return(NULL)

  built <- skeleton[chosen, , drop = FALSE]
  built[[roles$time]] <- wanted_time
  # A source sitting on its nominal grid -- what `coarsen_time` produces from a
  # declared `nominal_time` -- must have the visit label move with the visit, or
  # `time` and the nominal column disagree in the output and any step that
  # reconciles the two sees a contradiction the source never had.
  if (!is.null(roles$nominal_time)) {
    nominal <- suppressWarnings(as.numeric(skeleton[[roles$nominal_time]]))
    gap <- abs(time[observed] - nominal[observed])
    if (length(gap) && !anyNA(gap) &&
        all(gap < sqrt(.Machine$double.eps))) {
      built[[roles$nominal_time]] <- wanted_time
    }
  }
  rbind(skeleton[!observed, , drop = FALSE], built)
}

# TAD is time since the most recent dose, so it goes stale the moment TIME
# moves. Recomputed from the generated skeleton rather than carried over from
# the anchor, where it would describe a schedule the avatar no longer has.
.recompute_tad <- function(skeleton, roles) {
  if (is.null(roles$tad)) return(skeleton)
  derived <- .derived_tad(skeleton, roles)
  if (all(is.na(derived))) return(skeleton)
  skeleton[[roles$tad]] <- ifelse(is.na(derived), 0, derived)
  skeleton
}

.subject_value <- function(data, rows, column) {
  .first_present(data[[column]][rows])
}

.synthesize_covariates <- function(skeleton, data, roles, donor_indices,
                                   weights, profiles, subject_noise_sd) {
  for (covariate in roles$covariates) {
    template <- data[[covariate]]
    values <- lapply(donor_indices, function(donor) {
      .subject_value(data, profiles$subject_rows[[donor]], covariate)
    })

    if (is.numeric(template) && !is.factor(template)) {
      numeric_values <- vapply(values, function(value) {
        if (!length(value) || is.na(value)) NA_real_ else as.numeric(value)
      }, numeric(1))
      okay <- is.finite(numeric_values)
      if (!any(okay)) {
        value <- .first_present(skeleton[[covariate]])
      } else {
        available <- numeric_values[okay]
        available_weights <- weights[okay] / sum(weights[okay])
        positive <- all(available > 0)
        skewed <- positive && length(available) > 1L &&
          max(available) / max(stats::median(available), .Machine$double.eps) > 3
        working <- if (skewed) log(available) else available
        center <- sum(working * available_weights)
        spread <- if (length(working) > 1L) stats::sd(working) else 0
        if (!is.finite(spread) || spread <= 0) {
          spread <- max(abs(center) * 0.05, 0.01)
        }
        value <- center + stats::rnorm(1L, sd = subject_noise_sd * spread)
        if (skewed) value <- exp(value)
        if (positive) value <- max(value, sqrt(.Machine$double.eps))
      }
    } else {
      character_values <- vapply(values, function(value) {
        if (!length(value) || is.na(value)) NA_character_ else as.character(value)
      }, character(1))
      okay <- !is.na(character_values)
      if (!any(okay)) {
        value <- .first_present(skeleton[[covariate]])
      } else {
        available_weights <- weights[okay] / sum(weights[okay])
        value <- sample(character_values[okay], 1L, prob = available_weights)
        if (is.logical(template)) value <- identical(value, "TRUE")
      }
    }
    skeleton[[covariate]][] <- value
  }
  skeleton
}

# Some roles only mean something to the differentially private engines: they
# carry accounting semantics AVATAR has no equivalent for. Rather than silently
# ignore them -- which is how a user ends up believing a treatment arm was
# handled when it was not -- reject them and point at the role that does the job.
# `cmt` says where a dose goes; `dvid` says which endpoint a measurement is.
# Nothing infers the second from the first, so a dataset whose observations sit
# in several compartments with no `dvid` has every measurement treated as one
# endpoint: one value transform across both scales, one censoring boundary, and
# one attendance cell per visit no matter how many endpoints were read at it.
# That last one used to delete rows outright (`SIM-036`). Refusing is better
# than guessing, because `CMT` is not reliably the endpoint -- two endpoints can
# be read from one compartment, and plenty of datasets put every observation in
# `CMT = 2` and separate endpoints only by `DVID`.
.require_endpoint_key <- function(data, roles) {
  if (!is.null(roles$dvid) || is.null(roles$cmt)) return(invisible(TRUE))
  observed <- .observation_rows(data, roles, require_present = TRUE)
  if (!any(observed)) return(invisible(TRUE))
  compartments <- unique(as.character(data[[roles$cmt]][observed]))
  compartments <- sort(compartments[!is.na(compartments)])
  if (length(compartments) < 2L) return(invisible(TRUE))
  stop("Observations occupy ", length(compartments), " compartments (",
       paste(compartments, collapse = ", "), ") but no `dvid` is declared, so ",
       "every measurement would be treated as one endpoint.\n",
       "  Declare the endpoint key. One column may serve as both, which is ",
       "what a NONMEM `CMT` usually is:\n",
       "    pmx_roles(..., cmt = \"", roles$cmt, "\", dvid = \"", roles$cmt,
       "\")", call. = FALSE)
}

# Say out loud what an absent `dvid` means. Where `cmt` is declared the check
# above can refuse, but with no endpoint key at all there is nothing to test
# against: two endpoints measured at disjoint times are indistinguishable from
# one endpoint measured at all of them, and the package cannot tell. The note is
# the only defense left, so it states the consequence rather than merely the
# fact. Repeated subject-and-time observations are reported when present because
# they are the one visible hint -- they are also perfectly ordinary as replicate
# measurements, which is why this informs rather than refuses.
.note_single_endpoint <- function(data, roles) {
  if (!is.null(roles$dvid)) return(invisible(TRUE))
  observed <- .observation_rows(data, roles, require_present = TRUE)
  if (!any(observed)) return(invisible(TRUE))
  lines <- paste0(
    "synpmx_avatar(): no `dvid` declared, so every observation is treated as ",
    "one endpoint.\n  Correct for a single-endpoint study; declare `dvid` if ",
    "this one has more."
  )
  repeated <- sum(duplicated(paste(
    as.character(data[[roles$id]][observed]),
    format(suppressWarnings(as.numeric(data[[roles$time]][observed])),
           digits = 12, trim = TRUE)
  )))
  if (repeated) {
    lines <- paste0(lines, "\n  ", repeated, " observation row(s) share a ",
                    "subject and time with another; that is ordinary for ",
                    "replicate\n  measurements and expected if two endpoints ",
                    "are being read at one visit.")
  }
  message(lines)
  invisible(TRUE)
}

.reject_dp_only_roles <- function(roles) {
  guidance <- c(
    assigned_dose = "carry the column with `keep`",
    exclude = "simply leave the column undeclared; it is then dropped"
  )
  present <- names(guidance)[
    vapply(names(guidance), function(r) !is.null(roles[[r]]), logical(1))
  ]
  if (length(present)) {
    lines <- vapply(present, function(r) {
      paste0("  `", r, "` is a differential-privacy role; synpmx_avatar() ",
             "does not use it. Instead, ", guidance[[r]], ".")
    }, character(1))
    stop("Roles that do not apply to AVATAR were declared:\n",
         paste(lines, collapse = "\n"), call. = FALSE)
  }
}

# Censoring -------------------------------------------------------------------
#
# AVATAR makes no formal privacy guarantee, so it may read the censoring
# boundary straight off the source data. A differentially private engine may
# not: there the limit is a declared public input, because a boundary inferred
# from confidential records is itself a release. That asymmetry is why this
# lives here rather than being shared with `fit_private`'s path.

# One endpoint's censoring boundaries, recovered from how the source reports
# censored rows under the Monolix convention: DV sits at the boundary, and
# LIMIT carries the other end of an interval when there is one.
.source_censoring <- function(data, roles, endpoint_name) {
  if (is.null(roles$cens)) return(NULL)
  observed <- .observation_rows(data, roles, require_present = TRUE) &
    .endpoint(data, roles) == endpoint_name
  if (!any(observed)) return(NULL)
  cens <- suppressWarnings(as.numeric(as.character(data[[roles$cens]][observed])))
  dv <- suppressWarnings(as.numeric(data[[roles$dv]][observed]))
  limit <- if (is.null(roles$limit)) rep(NA_real_, sum(observed)) else
    suppressWarnings(as.numeric(data[[roles$limit]][observed]))

  out <- list()
  left <- is.finite(cens) & cens == 1 & is.finite(dv) & !is.finite(limit)
  interval <- is.finite(cens) & cens == 1 & is.finite(dv) & is.finite(limit)
  right <- is.finite(cens) & cens == -1 & is.finite(dv)
  # A study with more than one assay limit collapses to the most conservative
  # boundary rather than inventing a per-row rule the source cannot support.
  if (any(left)) out$left <- max(dv[left])
  if (any(right)) out$right <- min(dv[right])
  if (any(interval)) {
    out$interval <- c(min(limit[interval]), max(dv[interval]))
  }
  if (!length(out)) NULL else out
}

# Replace reported values on censored rows with a latent draw inside the
# censoring region, so that blending sees plausible values rather than a stack
# of identical boundary substitutions. Without this every censored donor drags
# the blend toward the limit, and the synthetic data inherits a floor the real
# study does not have.
.impute_censored <- function(data, roles) {
  if (is.null(roles$cens)) return(data)
  observed <- .observation_rows(data, roles, require_present = TRUE)
  if (!any(observed)) return(data)
  endpoint <- .endpoint(data, roles)
  cens <- suppressWarnings(as.numeric(as.character(data[[roles$cens]])))
  dv <- suppressWarnings(as.numeric(data[[roles$dv]]))
  limit <- if (is.null(roles$limit)) rep(NA_real_, nrow(data)) else
    suppressWarnings(as.numeric(data[[roles$limit]]))

  for (name in unique(endpoint[observed])) {
    rows <- which(observed & endpoint == name & is.finite(cens) & cens != 0)
    if (!length(rows)) next
    for (i in rows) {
      if (cens[i] == 1 && is.finite(limit[i])) {
        # Interval-censored: the value lies between LIMIT and DV.
        dv[i] <- stats::runif(1L, min(limit[i], dv[i]), max(limit[i], dv[i]))
      } else if (cens[i] == 1) {
        # Left-censored: uniform below the limit. A uniform draw rather than a
        # fixed LLOQ/2 avoids replacing one artificial spike with another.
        dv[i] <- stats::runif(1L, 0, dv[i])
      } else if (cens[i] == -1) {
        # Right-censored: no upper bound is reported, so extend by the same
        # relative width a left-censored draw would span.
        dv[i] <- dv[i] * stats::runif(1L, 1, 2)
      }
    }
  }
  data[[roles$dv]] <- dv
  data
}

.donor_trajectory <- function(data, roles, rows, endpoint_name, transform) {
  subject_data <- data[rows, , drop = FALSE]
  observed <- .observation_rows(subject_data, roles, require_present = TRUE)
  endpoint <- .endpoint(subject_data, roles)
  selected <- observed & endpoint == endpoint_name
  time <- .aligned_time(subject_data, roles)[selected]
  value <- .transform_dv(subject_data[[roles$dv]][selected], transform)
  okay <- is.finite(time) & is.finite(value)
  list(time = time[okay], value = value[okay])
}

# A donor's contribution to the anchor's observation times, and NA everywhere
# the donor was not actually followed.
#
# `rule = 1` is the whole point. Interpolating *between* a donor's own
# observations is a claim about a patient who was measured either side of the
# gap; extending *beyond* them is a claim about a patient nobody measured, and
# the two must not be confused.
#
# This used to fall back to rescaling the anchor's time range onto the donor's,
# which reads as "the same shape on a stretched clock". That is wrong whenever
# the missing region is a distinct kinetic phase rather than a shifted window,
# and pharmacometric data is full of that case: on `warfarin`, 19 of 32 subjects
# have no `cp` sample before 24 h, so a donor's 24 h *elimination* concentration
# was served up as an anchor's 0.5 h *absorption* concentration. The generated
# median at 0.5 h was 8.8 against a source median of 0.0, the absorption limb
# was replaced by a flat plateau, and the warped segment rejoined the real data
# with a visible upward step at 24 h (`SIM-046`).
#
# Returning NA instead hands the decision to `.weighted_available()`, which
# blends over the donors that *do* cover that time and renormalizes their
# weights. A row no donor covers falls to the caller's dataset-median fallback,
# which says so out loud -- the correct answer when a study simply never
# measured anybody at that time, and one the caller can act on.
.interpolate_trajectory <- function(trajectory, target_time) {
  time <- trajectory$time
  value <- trajectory$value
  if (!length(time)) return(rep(NA_real_, length(target_time)))
  # One observation carries no shape, so it says nothing about any other time.
  # It is still the donor's value *at* that time, which interpolation would
  # give back, so keep it there and nowhere else.
  if (length(unique(time)) == 1L) {
    return(ifelse(target_time == time[[1L]], mean(value), NA_real_))
  }
  stats::approx(time, value, xout = target_time, ties = mean, rule = 1)$y
}

# What the cohort did at each time, for the rows no donor can speak to.
#
# The old fallback was one number -- the endpoint's median over every
# observation at every time -- which is the one quantity guaranteed to have the
# wrong shape. On `warfarin` it put 6.4 at 0.5 h, where the four subjects
# actually sampled then all read 0.0, because it averaged the absorption limb
# together with five days of elimination.
#
# A median *per time* has the right shape and is still a cohort marginal rather
# than any patient's value. Times backed by a single observation are dropped
# for that second reason: a median over one subject is that subject. The floor
# is `min_pattern_share`'s default of 2 and for the same reason -- "at least two
# real patients stand behind this" is the weakest claim worth making.
.cohort_median_trajectory <- function(data, roles, endpoint_name, transform,
                                      source_observed, source_endpoint,
                                      floor = 2L) {
  selected <- source_observed & source_endpoint == endpoint_name
  time <- .aligned_time(data, roles)[selected]
  value <- .transform_dv(data[[roles$dv]][selected], transform)
  okay <- is.finite(time) & is.finite(value)
  time <- time[okay]
  value <- value[okay]
  if (!length(value)) {
    stop("No usable DV values exist for endpoint `", endpoint_name, "`.",
         call. = FALSE)
  }
  held <- table(time)
  shared <- as.numeric(names(held))[held >= floor]
  # Every time point held by one subject only. Nothing time-resolved can be
  # said without naming somebody, so this degrades to the old flat median --
  # wrong in shape, but it is the only cohort-level answer left.
  if (!length(shared)) {
    return(list(time = numeric(), value = numeric(),
                flat = stats::median(value)))
  }
  keep <- time %in% shared
  medians <- tapply(value[keep], time[keep], stats::median)
  list(time = as.numeric(names(medians)), value = as.numeric(medians),
       flat = stats::median(value))
}

# `rule = 2` here, unlike in `.interpolate_trajectory()`. This curve is the
# whole cohort rather than one donor, so holding its first or last value beyond
# the ends is a statement about the study's own observed range, not an invented
# measurement of a patient nobody followed.
.interpolate_cohort_median <- function(median_trajectory, target_time) {
  time <- median_trajectory$time
  if (length(time) < 2L) {
    value <- if (length(time) == 1L) median_trajectory$value[[1L]]
             else median_trajectory$flat
    return(rep(value, length(target_time)))
  }
  stats::approx(time, median_trajectory$value, xout = target_time,
                ties = mean, rule = 2)$y
}

.endpoint_noise_scale <- function(data, roles, endpoint_name, transform) {
  observed <- .observation_rows(data, roles, require_present = TRUE)
  endpoint <- .endpoint(data, roles)
  values <- .transform_dv(
    data[[roles$dv]][observed & endpoint == endpoint_name], transform
  )
  if (identical(transform$method, "log_offset")) return(1)
  scale <- stats::sd(values, na.rm = TRUE)
  if (!is.finite(scale) || scale <= 0) {
    scale <- max(abs(stats::median(values, na.rm = TRUE)) * 0.1, 0.01)
  }
  scale
}

.synthesize_trajectories <- function(skeleton, data, roles, donor_indices,
                                     weights, profiles, subject_noise_sd,
                                     residual_noise_sd, residual_phi,
                                     warnings,
                                     censoring_by_endpoint = list()) {
  allowed <- .observation_rows(skeleton, roles)
  present <- allowed & !is.na(skeleton[[roles$dv]])
  endpoint <- .endpoint(skeleton, roles)
  source_observed <- .observation_rows(data, roles, require_present = TRUE)
  source_endpoint <- .endpoint(data, roles)

  for (endpoint_name in unique(endpoint[present])) {
    target_rows <- which(present & endpoint == endpoint_name)
    target_time <- .aligned_time(skeleton, roles)[target_rows]
    transform <- profiles$transforms[[endpoint_name]]
    if (is.null(transform)) {
      transform <- .choose_transform(
        data[[roles$dv]][source_observed & source_endpoint == endpoint_name]
      )
    }
    donor_matrix <- vapply(donor_indices, function(donor) {
      trajectory <- .donor_trajectory(
        data, roles, profiles$subject_rows[[donor]], endpoint_name, transform
      )
      .interpolate_trajectory(trajectory, target_time)
    }, numeric(length(target_time)))
    if (is.null(dim(donor_matrix))) {
      donor_matrix <- matrix(donor_matrix, ncol = length(donor_indices))
    }
    if (subject_noise_sd == 0 && residual_noise_sd == 0 &&
        ncol(donor_matrix) > 1L) {
      donor_spread <- apply(donor_matrix, 1L, function(values) {
        values <- values[is.finite(values)]
        if (length(values) < 2L) Inf else diff(range(values))
      })
      if (length(donor_spread) &&
          all(donor_spread <= sqrt(.Machine$double.eps))) {
        warnings$add(paste0(
          "Endpoint `", endpoint_name,
          "` had indistinguishable donor trajectories with noise disabled; ",
          "an exact source-shaped trajectory was mathematically unavoidable."
        ))
      }
    }
    blended <- apply(donor_matrix, 1L, .weighted_available, weights = weights)
    missing_blend <- !is.finite(blended)
    if (any(missing_blend)) {
      fallback <- .cohort_median_trajectory(
        data, roles, endpoint_name, transform, source_observed, source_endpoint
      )
      blended[missing_blend] <- .interpolate_cohort_median(
        fallback, target_time[missing_blend]
      )
      # Deliberately carries no count. The collector dedups on the exact string
      # and this fires per subject, so a count turns one line into one line per
      # avatar -- 14 near-identical paragraphs on `warfarin` alone.
      warnings$add(sprintf(paste0(
        "Endpoint `%s` has generated observations at times no donor was ",
        "measured at; the cohort's median trajectory was used for those rows. ",
        "Expected wherever subjects were sampled on different schedules."
      ), endpoint_name))
    }

    scale <- .endpoint_noise_scale(data, roles, endpoint_name, transform)
    shift <- stats::rnorm(1L, sd = subject_noise_sd * scale)
    residual <- .ar1_noise(
      length(blended), residual_phi, residual_noise_sd * scale
    )
    # The back-transformed blend is the latent value: what the subject would
    # have measured with no assay limit. DV, CENS, and LIMIT are then
    # reconstructed from it together, so the three always agree.
    generated <- .inverse_dv(blended + shift + residual, transform)
    # A discrete endpoint is put back on its own scale before anything else
    # reads the value, so censoring and the emitted DV see the same number. The
    # blend is a weighted mean and the noise terms are continuous, so without
    # this a 0/1 endpoint leaves the source's range entirely (`SIM-045`).
    generated <- .snap_endpoint_values(
      generated, profiles$value_types[[endpoint_name]]
    )
    skeleton[[roles$dv]][target_rows] <- generated
    censoring <- censoring_by_endpoint[[endpoint_name]]
    if (!is.null(censoring)) {
      skeleton <- .censor_latent(skeleton, target_rows, generated, roles,
                                 public = censoring)
    }
  }
  skeleton
}

.source_uses_standard_mdv <- function(data, roles) {
  if (is.null(roles$mdv)) return(FALSE)
  expected_observed <- .is_zero(data[[roles$evid]]) & !is.na(data[[roles$dv]])
  actual_observed <- .is_zero(data[[roles$mdv]])
  all(expected_observed == actual_observed)
}

.derive_standard_mdv <- function(skeleton, roles) {
  observed <- .is_zero(skeleton[[roles$evid]]) &
    !is.na(skeleton[[roles$dv]])
  template <- skeleton[[roles$mdv]]
  values <- ifelse(observed, 0, 1)
  if (is.factor(template)) {
    skeleton[[roles$mdv]] <- factor(
      as.character(values), levels = levels(template),
      ordered = is.ordered(template)
    )
  } else if (is.character(template)) {
    skeleton[[roles$mdv]] <- as.character(values)
  } else if (is.integer(template)) {
    skeleton[[roles$mdv]] <- as.integer(values)
  } else {
    skeleton[[roles$mdv]] <- as.numeric(values)
  }
  skeleton
}

# Advice that depends on what the caller already did. Three alerts used to end
# with "declare a `nominal_time` role" unconditionally, which is worse than
# useless to somebody who declared it on the first line of their script: it
# reads as though the run did not notice. When the grid is already nominal the
# remaining levers are the pool split and the floor, so say those instead.
.schedule_advice <- function(source_roles, min_pattern_share, grid) {
  nominal <- identical(grid, "nominal")
  strata <- source_roles$strata
  c(
    if (!nominal) {
      paste("declare a `nominal_time` role so coarsening puts more patients on",
            "the same visits")
    },
    if (identical(grid, "mixed")) {
      paste("fill in `nominal_time` on the rows that are missing it -- the run",
            "reports how many were snapped to an inferred grid instead")
    },
    if (length(strata)) {
      sprintf(paste("drop `strata` (%s) so the arms share one pool",
                    "of visit sets instead of one each"),
              paste(strata, collapse = ", "))
    },
    sprintf("lower `min_pattern_share` from %d", as.integer(min_pattern_share))
  )
}

# A deliberately loud, immediate alert for what a run cannot fix on its own --
# a source with fewer subjects than the donor floor, a schedule the grid could
# not collapse. Unlike the collected end-of-run warning() (which
# `suppressWarnings()` removes -- as the demo did, hiding it), the message()
# banner survives suppressWarnings and prints at once, in red on an interactive
# console.
#
# `message()` is the only channel that prints. A condition is also signalled so
# a caller can react programmatically, but it is deliberately NOT a `warning`
# subclass, and that detail is the whole point: `warning()`, and equally a
# signalled condition carrying the `warning` class, is picked up by knitr's
# calling handler and rendered *in addition to* the message, so every alert
# appeared twice in a knitted report -- which is exactly the volume this
# formatting exists to cut. A plain `synpmx_alert` condition is invisible to
# knitr, to `suppressWarnings()`, and to `options(warn = 2)`; handle it with
# `withCallingHandlers(synpmx_alert = ...)` when you want to collect them.
.loud_warn <- function(title, what, why = NULL, fix = NULL, kind = "ALERT") {
  text <- .alert_text(title, what, why, fix, kind)
  banner <- if (interactive()) paste0("\033[1;31m", text, "\033[0m") else text
  # The trailing blank line is what makes a run emitting three of these
  # skimmable rather than a wall.
  message(banner, "\n")
  signalCondition(structure(
    class = c("synpmx_alert", "condition"),
    list(message = paste0(text, "\n"), call = NULL)
  ))
  invisible(NULL)
}

# "3 of 32 patients (9%)", the phrase every schedule alert opens with.
.count_phrase <- function(k, n, noun = "patient") {
  sprintf("%d of %d %s%s (%.0f%%)", k, n, noun, if (k == 1L) "" else "s",
          if (n) 100 * k / n else 0)
}

# Choose the donors whose trajectories are blended onto the anchor's event
# skeleton, in two stages, both confined to the anchor's own administration
# route (`.route_key()`; never crossed, at any stage, for any shortfall):
#
#   1. Same-schedule donors -- identical event signature -- taken nearest-first
#      in profile space, up to `k`.
#   2. If stage 1 yields fewer than `k`, the shortfall is filled with the
#      nearest remaining route-compatible subjects in profile space, whatever
#      their dose or schedule.
#
# "Nearest" is Euclidean distance between retained PCA profile coordinates
# (`.neighbor_distances()`), ties broken by subject index so selection is
# deterministic under a fixed seed. Stage-2 donors' measurements are mapped onto
# the anchor's own observation times by interpolation, so the avatar keeps its
# anchor's regimen while its values are averaged across >= k real patients.
#
# A shortfall that survives both stages means the anchor's route arm holds fewer
# than k + 1 subjects. There is no legal donor left to borrow, so the caller
# drops such anchors before generation and alerts; reaching this function with
# no compatible donor at all is the residual case where dropping was impossible.
# Source subjects whose event structure is extreme on the high side -- a
# follow-up or dose count more than `mult` times the cohort's 90th percentile.
# These are the anchors that would give an avatar a structurally extreme
# skeleton (the long tail a reader notices); the default screen keeps them out
# of the anchor pool. Anchoring on the 90th percentile, not the median, means
# ordinary spread does not trip it: only a subject well beyond the high end of
# normal is excluded. (wbcSim follow-up: 90th percentile ~648 h, so the cut sits
# near 1300 h and drops only the 1730/4580 h subjects, not the ordinary ~650 h.)
# A stratum below this many source patients is not balanced exactly, because
# reproducing its size exactly would disclose that size. A cohort with one
# patient in a joint arm x cohort cell gets exactly one avatar there on every
# seed, which says "that cell held one patient" as clearly as printing it. Small
# strata keep the old behavior instead -- their slots are drawn with replacement
# from the pooled small-stratum members, so the count varies by seed.
#
# Three, not two: a floor of 2 leaves pairs exactly reproduced, and a pair is
# the case `min_pattern_share` already treats as too rare to reuse elsewhere.
.strata_balance_floor <- 3L

# The number of avatars to build on each stratum's anchors.
#
# Sampling anchors with replacement across the whole pool -- which is what this
# replaces -- leaves stratum sizes to the binomial: a balanced six-arm design of
# 30 patients per arm came back between 21 and 39 per arm depending only on the
# seed. Arm membership was always exact (an avatar never leaves its anchor's
# stratum), but the *balance* was not, and any downstream summary that groups by
# arm silently inherits the wobble.
#
# Targets are proportions of the SOURCE cohort, so they still mean something
# when `n_subjects` differs from the source size, and largest-remainder rounding
# makes them sum to `n_subjects` exactly. Strata below the floor above, and
# strata with no eligible anchor left after screening, fall back into a pooled
# remainder drawn as before.
# Who an avatar may be re-anchored on when its own anchor cannot be masked.
#
# Inside its own stratum, because the anchor's rows are copied into the output
# with its `strata` values attached: an avatar moved across arms leaves the run
# reporting the balance it planned while the finished table holds something
# else. Drawing from the whole pool -- which is what this replaces -- scrambled
# arm sizes in proportion to how often re-anchoring fired, and nothing said so.
#
# The whole pool is still the fallback when the stratum has nobody else to
# offer. Emitting an avatar wearing one real patient's visit set is a
# disclosure; an arm off its target is a fidelity cost, and the caller counts
# it as `strata_crossed` so the trade is visible rather than silent.
.reanchor_candidates <- function(allowed, anchor, strata_key, target_stratum) {
  candidates <- setdiff(allowed, anchor)
  if (is.null(target_stratum) || !length(candidates)) {
    return(list(candidates = candidates, crossed = FALSE))
  }
  within <- candidates[strata_key[candidates] == target_stratum]
  if (length(within)) {
    list(candidates = within, crossed = FALSE)
  } else {
    list(candidates = candidates, crossed = TRUE)
  }
}

.strata_targets <- function(strata_key, allowed, n_subjects,
                            floor = .strata_balance_floor) {
  source_sizes <- table(strata_key)
  eligible <- unique(strata_key[allowed])
  # Selected as a logical rather than by name, because a blank stratum label is
  # a name `[` cannot look up: `source_sizes[""]` is NA, not the count, and the
  # NA reaches `spare` as a missing condition rather than as a wrong number.
  keep <- source_sizes >= floor & names(source_sizes) %in% eligible
  if (!any(keep)) return(NULL)

  share <- as.numeric(source_sizes[keep]) / length(strata_key)
  exact <- share * n_subjects
  target <- floor(exact)
  # Largest remainder: hand out the seats floor() dropped, to the strata with
  # the biggest fractional parts. The balanced strata claim round(sum(exact))
  # seats in total, so whatever is left over stays with the small strata.
  spare <- round(sum(exact)) - sum(target)
  if (spare > 0) {
    order_by_remainder <- order(exact - target, decreasing = TRUE)
    take <- utils::head(order_by_remainder, spare)
    target[take] <- target[take] + 1L
  }
  target <- as.integer(target)
  names(target) <- names(source_sizes)[keep]
  target[target > 0L]
}

.structural_outlier_anchors <- function(source, roles, mult = 2) {
  subjects <- .unique_in_order(source[[roles$id]])
  key <- factor(as.character(source[[roles$id]]),
                levels = as.character(subjects))
  observed <- .observation_rows(source, roles, require_present = TRUE)
  dosed <- .dose_rows(source, roles)
  time <- suppressWarnings(as.numeric(source[[roles$time]]))
  follow_up <- vapply(split(ifelse(observed, time, NA_real_), key), function(v) {
    v <- v[is.finite(v)]
    if (length(v)) max(v) else NA_real_
  }, numeric(1))[as.character(subjects)]
  n_doses <- as.numeric(tapply(as.integer(dosed), key, sum)[
    as.character(subjects)
  ])
  high <- function(x) {
    finite <- x[is.finite(x)]
    if (!length(finite)) return(rep(FALSE, length(x)))
    reference <- stats::quantile(finite, 0.90, names = FALSE)
    if (!is.finite(reference) || reference <= 0) return(rep(FALSE, length(x)))
    is.finite(x) & x > mult * reference
  }
  which(high(follow_up) | high(n_doses))
}

.select_donors <- function(anchor, profiles, k, warnings,
                           max_donor_weight = 0.50) {
  target <- as.integer(k)
  routes <- profiles$routes %||% rep("none", length(profiles$subjects))
  compatible <- setdiff(which(routes == routes[anchor]), anchor)
  if (!length(compatible)) {
    # No subject shares the anchor's route, and route is never crossed, so the
    # avatar is a noised copy of this one patient. The caller drops such anchors
    # where it can and alerts loudly where every anchor is in this position.
    return(list(indices = anchor, distances = 0, weights = 1))
  }

  nearest <- function(pool) {
    if (!length(pool)) return(integer())
    pool[order(.neighbor_distances(profiles$coordinates, anchor, pool), pool)]
  }

  same_group <- intersect(
    which(profiles$signatures == profiles$signatures[anchor]), compatible
  )
  chosen <- nearest(same_group)
  if (length(chosen) > target) chosen <- chosen[seq_len(target)]

  if (length(chosen) < target) {
    fill <- nearest(setdiff(compatible, chosen))
    need <- target - length(chosen)
    chosen <- c(chosen, fill[seq_len(min(length(fill), need))])
    if (length(same_group) < target) {
      warnings$add(paste0(
        "Fewer than ", target, " same-schedule donors were available for at ",
        "least one subject; the nearest donors from other dose/schedule groups ",
        "on the same route were borrowed to reach the floor, so some ",
        "measurements are blended across doses."
      ))
    }
  }

  chosen_distances <- .neighbor_distances(profiles$coordinates, anchor, chosen)
  if (length(chosen_distances) &&
      all(chosen_distances <= sqrt(.Machine$double.eps))) {
    warnings$add(
      "Duplicated subject profiles produced zero neighbor distances; epsilon-stabilized randomized weights were used."
    )
  }
  list(
    indices = chosen,
    distances = chosen_distances,
    weights = .randomized_weights(chosen_distances, max_donor_weight)
  )
}

#' Synthesize a structurally faithful PMX dataset (AVATAR-style)
#'
#' Samples complete subject event templates and fills them with AVATAR-like,
#' endpoint-specific blends of compatible subjects' baseline covariates and
#' longitudinal measurements. Event-control fields such as EVID, AMT, RATE,
#' CMT, and DVID are never averaged or independently generated.
#'
#' This is an AVATAR-inspired adaptation, not an exact implementation of
#' published AVATAR software. It creates synthetic data for model-workflow
#' exploration. It does not provide formal anonymization or preserve scientific
#' parameter or covariate-response relationships.
#'
#' @details
#' Donors are selected in two stages, both confined to the anchor's own
#' administration route, which is never crossed: same-signature donors first,
#' taken nearest-first by Euclidean distance between retained PCA profile
#' coordinates, then --- if that yields fewer than `k` --- the nearest remaining
#' route-compatible subjects regardless of dose or schedule.
#'
#' For the selected donors, randomized raw weights are
#' `Exp(1) / max(distance, epsilon) * 2^(-randomized_rank)`. They are normalized
#' and then capped so that *no* donor exceeds `max_donor_weight`, the excess
#' being redistributed proportionally among the donors still below the cap until
#' none is over. A cap below `1/K` for `K` donors cannot be satisfied and
#' relaxes to `1/K`, i.e. uniform weights. The same subject weights are used for
#' covariates and all endpoints; weights are renormalized locally when a donor
#' lacks a requested endpoint/time value.
#'
#' Positive-like endpoints use an offset log scale and are constrained to be
#' nonnegative after back-transformation. Other endpoints use the identity
#' scale. Transform choices and interpolation alignment are recorded in the
#' returned `pmx_settings` attribute.
#'
#' @param data A source PMX data frame or tibble.
#' @param roles Explicit roles from [pmx_roles()]. Columns listed in
#'   `roles$exclude` are omitted from the generated output.
#' @param n_subjects Number of synthetic subjects. `NULL` retains the source count.
#' @param seed Reproducibility seed. The caller's random-number state is
#'   restored on exit.
#' @param event_method Event generation method. The prototype supports
#'   `"template"`.
#' @param dv_method Measurement method. The prototype supports
#'   `"avatar_blend"`.
#' @param k Number of real patients blended into each synthetic subject
#'   (default 5). Same-schedule donors are used first; when a subject's
#'   dose/schedule group holds fewer than `k`, the nearest subjects from other
#'   groups *on the same administration route* are borrowed to reach `k`,
#'   blending measurements across doses. Route is never crossed, so a route arm
#'   holding fewer than `k + 1` subjects cannot reach the floor at all; those
#'   subjects are dropped from the anchor pool with a loud alert, and the
#'   synthetic cohort does not represent that arm.
#' @param pca_variance Fraction of usable profile variance retained for
#'   neighborhood distances.
#' @param subject_noise_sd Nonnegative subject perturbation multiplier.
#' @param residual_noise_sd Nonnegative within-trajectory noise multiplier.
#' @param residual_phi AR(1) correlation in observation order, strictly between
#'   -1 and 1.
#' @param time_jitter Standard deviation for coherent tied-time jitter. Zero,
#'   the default, leaves the event template's times unchanged. This is a realism
#'   control, **not** a privacy control: every jittered time is clamped inside
#'   its own Voronoi cell, so no value of `time_jitter` moves a visit more than
#'   half a gap from where the source subject's visit was, and the source
#'   schedule stays recoverable. Use `coarsen_time` for that.
#' @param screen When `TRUE` (default), a source subject whose follow-up length
#'   or dose count is more than twice the cohort's 90th percentile is not used as
#'   an anchor, so no avatar inherits an extreme skeleton (the long tail a reader
#'   notices). Anchoring the cut on the 90th percentile, not the median, means
#'   ordinary spread is kept: only a subject well beyond the high end of normal
#'   is dropped. Only these structural axes are screened; dose magnitude (which
#'   weight-based dosing makes noisy) and DV (which is blended, not copied) are
#'   not. A source with no extreme subject is unaffected. Set `FALSE` to anchor
#'   on every subject. For a fuller, tunable screen of the generated output, see
#'   [flag_identifiable_subjects()] and [remediate_identifiable_subjects()].
#' @param coarsen_time When `TRUE` (default), source times are collapsed onto a
#'   shared visit grid before generation, and per-visit deviations are pooled
#'   across the cohort and resampled independently onto each avatar. The grid is
#'   the `nominal_time` role where one is declared, and K-means centres of the
#'   pooled times otherwise. This is the mechanism that stops an avatar from
#'   carrying one real subject's exact visit schedule: the event skeleton is
#'   copied verbatim from a single anchor, and under actual recorded times almost
#'   every subject is alone in its event-signature class, so the copy is
#'   identifying. Snapping is many-to-one and *destroys* the deviation rather
#'   than perturbing it, which is what distinguishes this from `time_jitter`. A
#'   source already on nominal time has no deviation to remove or restore, so its
#'   output is unchanged. Run [skeleton_uniqueness()] on the source to see how
#'   much this has to do, and what it leaves behind. The cost is timing
#'   fidelity: an avatar's deviation from nominal is drawn from the cohort, not
#'   inherited, so `TIME` no longer pairs with its `DV` as precisely as the
#'   source did. Set `FALSE` to keep exact source timing.
#' @param min_pattern_share How many source subjects must share an attendance
#'   pattern before an avatar may be given it. Default 2; `1` restores copying
#'   the anchor's own pattern.
#'
#'   Two is the smallest value that means something, and what it means is
#'   precise: **no synthetic patient carries a schedule unique to a real
#'   patient.** An attacker who links a reproduced pattern to a participant gets
#'   at least two candidates, never one. Higher values hide more and are harder
#'   to state — "shared by at least three" is more conservative but no more
#'   defensible — and they cost sharply more, because attendance patterns are
#'   distributed with one common pattern and a long tail of singletons rather
#'   than a populated middle.
#'
#'   Once `coarsen_time` has put every subject on a shared visit grid, what
#'   remains of a schedule is which of those visits each subject attended. Every
#'   time is then shared, so no single visit is identifying — but the
#'   *combination of absences* is, and a patient who missed weeks 2 and 3 can be
#'   singled out by a fingerprint made of gaps. No grid can fix that at any
#'   resolution, because the grid decides where the visits are and not which ones
#'   a subject has. So the pattern is sampled from ones at least
#'   `min_pattern_share` subjects hold, frequency-weighted, and whole patterns
#'   are drawn rather than individual visits, which keeps an ending monotone
#'   instead of producing implausible attend/miss/attend sequences.
#'
#'   Unlike dropping the exposed subjects, nobody leaves the cohort: a subject
#'   with a rare pattern still contributes measurements as a donor, only their
#'   distinctive absences stop being reproduced. Dose events are never sampled,
#'   since that could emit a regimen no protocol permits. Raising this hides more
#'   and flattens the cohort's missingness further; where no pattern is shared
#'   widely enough, anchors keep their own and the run alerts loudly. Pools are
#'   formed within each `strata` stratum and endpoint set.
#'
#'   **Patterns below the floor are lost, not approximated.** An ending or
#'   dose-interruption pattern held by too few patients simply will not appear in
#'   the synthetic data, and that loss is the mechanism working — it is what
#'   stops an avatar carrying a schedule traceable to one person. Because it is a
#'   real cost to the data's realism, every run reports it: the number of source
#'   patterns excluded and how many subjects held them, both as a loud alert and
#'   as `patterns_total`, `patterns_dropped` and `subjects_with_dropped_pattern`
#'   in the settings, alongside `pattern_generated_fraction` for how often an
#'   arrangement had to be invented. What survives is how much missingness there
#'   was and of what kind; what is lost is which specific visits. Check those
#'   figures before deciding the default suits your study.
#' @section Dose recomputed from a blended covariate:
#' When the dose is a fixed multiple of a baseline covariate within each
#' assigned stratum — mg/kg, mg/m^2 — that multiplier is a protocol property the
#' whole stratum shares, while the covariate is individual and is already
#' blended across donors. `synpmx_avatar()` detects this and recomputes each
#' avatar's `AMT` (and any `rate`, keeping the infusion duration) from the
#' avatar's *own* blended covariate.
#'
#' This fixes two things at once. The amount is no longer one real patient's
#' real dose, which under proportional dosing discloses that patient's weight
#' exactly. And the avatar stops violating the protocol it claims to follow:
#' previously `AMT` was copied from the anchor while covariates were blended, so
#' a cohort dosed at exactly 5 mg/kg produced avatars anywhere from 4.4 to 5.3.
#'
#' Several dose levels are handled by clustering the observed ratios rather than
#' averaging within a group, so a 1/2/3 mg/kg escalation is recognised without
#' the arm being declared — and so is **intra-patient** escalation, where the
#' level changes within a subject and no subject-constant grouping could see it.
#'
#' Detection is conservative and fails closed: the ratios must collapse onto a
#' handful of levels while the amounts they came from stay varied, so a study
#' whose dose is unrelated to the covariate produces about as many ratios as
#' amounts and is refused. The covariate and the levels found are recorded as
#' `dose_basis` and `dose_levels` in the settings, `NA` when none was found.
#'
#' @param max_donor_weight Largest share of one synthetic subject that any one
#'   real donor may contribute. The default 0.50 states simply that no single
#'   real patient is more than half of any synthetic patient.
#'
#'   The floor `k` sets how many patients are blended; this cap is what bounds
#'   any single patient's contribution, so it, not `k`, is the parameter that
#'   limits how closely an avatar can resemble one real person. Without a cap
#'   the randomized weights are strongly concentrated --- a median 58% of an
#'   avatar in one donor at `k = 5`, and about 2.4 effective donors --- so the
#'   cap is what makes the floor mean anything.
#'
#'   Two diagnostics in the returned `pmx_settings` say where a given value
#'   landed: `mean_effective_donors` is `1 / sum(w^2)`, the number of donors an
#'   avatar is effectively blended from, and `cap_binding_fraction` is how often
#'   the cap actually fired. A cap binding on nearly every subject is not a
#'   guardrail but the weighting scheme itself, with the inverse-distance term
#'   underneath it doing little; one that never fires is not protecting
#'   anything. At `k = 5` the default binds on roughly two thirds of subjects.
#' @param preserve_strata_balance Give each declared stratum ([pmx_roles()]
#'   `strata`) the same share of the synthetic cohort it holds in the source.
#'   Default `TRUE`.
#'
#'   An avatar never leaves its anchor's stratum, but anchors are sampled with
#'   replacement, so without this the *balance* is left to chance: a balanced
#'   six-arm design of 30 patients per arm came back between 21 and 39 per arm
#'   depending only on the seed, and any downstream summary grouped by arm
#'   inherits that. Targets are proportions of the source cohort, so they still
#'   hold when `n_subjects` differs from the source size.
#'
#'   **Strata holding fewer than three source patients are deliberately not
#'   balanced**, because reproducing such a stratum's size exactly would
#'   disclose it: a joint cell holding one real patient would get exactly one
#'   avatar on every seed. Their slots are drawn with replacement as before, so
#'   the count varies. `strata_balanced` and `strata_stochastic` in the returned
#'   `pmx_settings` report how many strata fell on each side.
#'
#'   `FALSE` restores the unstratified draw. With no `strata` declared the
#'   argument has no effect.
#' @param on_donor_shortfall What to do with a subject whose administration
#'   route holds fewer than `k + 1` subjects, so that no legal donor set exists
#'   for it. `"drop"` (default) omits those subjects from the anchor pool: no
#'   avatar is built on them and the synthetic cohort does not represent that
#'   arm. `"noise"` keeps them, blending however many same-route donors exist
#'   (possibly none) and relying on `subject_noise_sd` and `residual_noise_sd`
#'   for the rest --- **not recommended**, because such a synthetic subject can
#'   remain close to one real patient; screen the result with
#'   [flag_identifiable_subjects()] if you use it. `"error"` refuses to generate
#'   and names the choice. Every branch alerts loudly. When *every* route arm is
#'   below the floor, `"drop"` would leave nothing to generate, so generation
#'   proceeds as if `"noise"`.
#'
#' @return An ordinary data frame or tibble with retained source columns, order,
#'   and practical classes. A lightweight `pmx_settings` attribute records the
#'   generator choices and endpoint transformations.
#' @export
#'
#' @examples
#' source <- data.frame(
#'   ID = rep(1:3, each = 4),
#'   TIME = rep(c(0, 0, 1, 2), 3),
#'   DV = c(0, 0.2, 2, 1, 0, 0.3, 3, 1.5, 0, 0.4, 4, 2),
#'   AMT = rep(c(100, 0, 0, 0), 3),
#'   EVID = rep(c(1L, 0L, 0L, 0L), 3),
#'   CMT = rep(c(1L, 2L, 2L, 2L), 3),
#'   WT = rep(c(60, 70, 80), each = 4)
#' )
#' roles <- pmx_roles(
#'   id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
#'   cmt = "CMT", covariates = "WT"
#' )
#' synthetic <- synpmx_avatar(source, roles, n_subjects = 2, seed = 123)
#' validate_pmx(synthetic, roles)$valid
synpmx_avatar <- function(data, roles, n_subjects = NULL, seed = 123,
                     event_method = "template",
                     dv_method = "avatar_blend", k = 5,
                     pca_variance = 0.90, subject_noise_sd = 0.15,
                     residual_noise_sd = 0.05, residual_phi = 0.6,
                     time_jitter = 0, screen = TRUE, coarsen_time = TRUE,
                     min_pattern_share = 2L,
                     max_donor_weight = 0.50,
                     preserve_strata_balance = TRUE,
                     on_donor_shortfall = c("drop", "noise", "error")) {
  on_donor_shortfall <- match.arg(on_donor_shortfall)
  if (!is.logical(preserve_strata_balance) ||
      length(preserve_strata_balance) != 1L ||
      is.na(preserve_strata_balance)) {
    stop("`preserve_strata_balance` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.data.frame(data)) stop("`data` must be a data frame or tibble.",
                                 call. = FALSE)
  .assert_roles(data, roles)
  .reject_dp_only_roles(roles)
  .require_endpoint_key(data, roles)
  .note_single_endpoint(data, roles)
  # Allowlist, not blocklist: keep only what a role names, so a column the user
  # never mentioned -- a secondary identifier, a site, a randomization date --
  # cannot ride out of a real subject into synthetic data by being forgotten.
  # There is no `exclude` here by design: not naming a column is how you drop it.
  retained_names <- intersect(names(data), .retained_role_columns(roles))
  dropped <- setdiff(names(data), retained_names)
  if (length(dropped)) {
    message("synpmx_avatar(): dropped ", length(dropped),
            " undeclared column(s): ", paste(dropped, collapse = ", "),
            ".\n  Declare a column in `keep` to carry it through verbatim.")
  }
  source <- data[, retained_names, drop = FALSE]
  source_roles <- roles
  source_roles$exclude <- NULL
  class(source_roles) <- "pmx_roles"
  validate_pmx(source, source_roles, strict = TRUE)
  subjects <- .unique_in_order(source[[source_roles$id]])
  n_subjects <- .validate_generator_options(
    n_subjects, length(subjects), event_method, dv_method, k, pca_variance,
    subject_noise_sd, residual_noise_sd, residual_phi, time_jitter,
    max_donor_weight, coarsen_time, min_pattern_share
  )

  .with_local_seed(seed, {
    warnings <- .warning_collector()
    # Boundaries come from the reported data, before imputation overwrites the
    # substituted values they are read from.
    censoring_by_endpoint <- list()
    if (!is.null(source_roles$cens)) {
      endpoint_names <- unique(.endpoint(
        source, source_roles
      )[.observation_rows(source, source_roles, require_present = TRUE)])
      for (name in endpoint_names) {
        censoring_by_endpoint[[name]] <- .source_censoring(
          source, source_roles, name
        )
      }
      censored_rows <- sum(suppressWarnings(as.numeric(as.character(
        source[[source_roles$cens]]
      ))) != 0, na.rm = TRUE)
      if (censored_rows && !length(Filter(Negate(is.null),
                                          censoring_by_endpoint))) {
        warnings$add(paste0(
          "A `cens` role is declared and ", censored_rows, " row(s) are ",
          "flagged, but no censoring boundary could be read from them; ",
          "CENS was carried through without being reconstructed."
        ))
      }
      source <- .impute_censored(source, source_roles)
    }
    # Before profiling, so equivalence classes are formed on the coarsened
    # schedule: donor compatibility (`.select_donors()` stage 1) matches on the
    # event signature, which under actual times is unique per subject and forces
    # every avatar into the route-only fallback.
    coarsened <- if (coarsen_time) {
      .coarsen_source_time(source, source_roles)
    } else {
      list(source = source, deviations = numeric(), grid = "off",
           grid_derived_rows = 0L, grid_rows = 0L)
    }
    source <- coarsened$source
    time_deviations <- coarsened$deviations
    # Coarsening onto a derived grid is best-effort: where subjects hold samples
    # closer together than the spread between subjects, no grid can merge them
    # without collapsing a real visit, so the schedule stays unique and the
    # verbatim skeleton copy stays identifying. That is the safe failure, but it
    # is a silent one, and a caller who asked for `coarsen_time` and did not get
    # it should hear so rather than infer it from the absence of an alert.
    unique_schedule_n <- NA_integer_
    unique_obs_time_n <- NA_integer_
    if (coarsen_time) {
      exposure <- skeleton_uniqueness(source, source_roles)
      still_unique <- sum(exposure$unique_schedule)
      # Two causes, opposite remedies, so the alert has to say which. A subject
      # sampled at a moment nobody else was is a grid failure -- declaring
      # `nominal_time` fixes it. A patient whose every time is shared and whose
      # *set of attended visits* is unique has visits missing, and no grid however fine
      # or coarse touches that; it is the outlier screen's problem.
      unshared <- sum(exposure$n_share_rarest_time == 1L, na.rm = TRUE)
      pattern_only <- max(still_unique - unshared, 0L)
      unique_schedule_n <- still_unique
      unique_obs_time_n <- unshared
      if (unshared > 0L) {
        .loud_warn(
          "unique observation times",
          paste(.count_phrase(unshared, length(subjects)),
                "were sampled at a moment no other patient was, even after",
                "coarsening."),
          why = paste("an avatar copies its anchor's observation times",
                      "verbatim, so it wears a schedule that belongs to one",
                      "real patient."),
          fix = paste("declare a `nominal_time` role. Coarsening then snaps",
                      "visits onto the real protocol grid instead of a",
                      "guessed one.")
        )
      }
      if (pattern_only > 0L) {
        .loud_warn(
          "unique visit sets",
          paste(.count_phrase(pattern_only, length(subjects)),
                "share every individual observation time with somebody, but",
                "the set of visits they have observations at is theirs alone --",
                "a missed visit, a discontinuation, or follow-up that has not",
                "reached the later visits."),
          why = paste("no time grid can help here, however fine or coarse: a",
                      "grid decides where the visits are, not which ones a",
                      "patient turned up for."),
          fix = paste("`min_pattern_share` already stops these sets being",
                      "reused (see the run report). Screen the result with",
                      "`flag_identifiable_subjects()` if it still matters.")
        )
      }
    }
    # The detector always returns its reasoning; `dose_basis` stays NULL unless
    # it actually found a basis, so nothing downstream changes.
    dose_basis_detection <- .dose_basis_for(source, source_roles)
    dose_basis <- if (is.null(dose_basis_detection$covariate)) NULL else
      dose_basis_detection
    profiles <- .build_profiles(source, source_roles, pca_variance)
    attendance <- .attendance_pool(source, source_roles, profiles,
                                   as.integer(min_pattern_share))
    if (!is.null(attendance) && attendance$dropped_patterns > 0L) {
      # A NOTE, not an alert: this is the mechanism working. It is reported
      # because it is the one place the synthetic data is deliberately less
      # faithful than the source, and the caller decides whether that matters.
      .loud_warn(
        "rare visit sets not reused",
        sprintf(paste("%d of %d distinct visit set%s, held by %d patient%s,",
                      "%s shared by fewer than %d patients and %s given to no",
                      "avatar."),
                attendance$dropped_patterns, attendance$total_patterns,
                if (attendance$total_patterns == 1L) "" else "s",
                attendance$dropped_subjects,
                if (attendance$dropped_subjects == 1L) "" else "s",
                if (attendance$dropped_patterns == 1L) "is" else "are",
                as.integer(min_pattern_share),
                if (attendance$dropped_patterns == 1L) "is" else "are"),
        why = paste("an avatar carrying a visit set unique to one real patient",
                    "could be traced back to them. Kept instead: how many",
                    "visits were missed and of what kind -- all at the end",
                    "(follow-up ending), a run in the middle (an interruption),",
                    "or",
                    "scattered. Which specific visits were missed is not",
                    "preserved."),
        fix = paste("nothing, unless this study's interruptions matter.",
                    "`min_pattern_share = 1` copies exact visit sets and gives",
                    "up the guarantee."),
        kind = "NOTE"
      )
    }
    if (!is.null(attendance)) {
      unpooled <- setdiff(unique(attendance$group), names(attendance$pool))
      if (length(unpooled)) {
        stranded <- sum(attendance$group %in% unpooled)
        .loud_warn(
          "no shared visit set to draw from",
          sprintf(paste("in %d of %d schedule group%s, covering %d patient%s,",
                        "no visit set is shared by %d or more patients."),
                  length(unpooled), length(unique(attendance$group)),
                  if (length(unpooled) == 1L) "" else "s",
                  stranded, if (stranded == 1L) "" else "s",
                  as.integer(min_pattern_share)),
          why = paste("avatars anchored there keep their anchor's own set of",
                      "attended visits, which is exactly what",
                      "`min_pattern_share` exists to avoid reproducing."),
          fix = paste0(paste(.schedule_advice(source_roles, min_pattern_share,
                                             coarsened$grid),
                             collapse = "; "), ".")
        )
      }
    }
    # The dose side of the same guarantee. `min_pattern_share = 1` turns the
    # whole guarantee off, so it turns this off too.
    dose_identifying <- if (as.integer(min_pattern_share) <= 1L) {
      rep(FALSE, length(subjects))
    } else {
      unname(.dose_schedule_sharing(source, source_roles, min_pattern_share))
    }
    # Whether there is anywhere safe to move an avatar TO. Where every patient's
    # dose schedule is unique -- individualised infusions, per-patient
    # titration -- re-anchoring cannot help, so it is not attempted and the
    # caller is told once instead of the run silently doing nothing.
    # Truncation first: a patient who merely stopped early can be given a
    # different, safe stopping depth, so they are not unmaskable after all.
    dose_plan <- if (as.integer(min_pattern_share) <= 1L) {
      list(depth = rep(NA_integer_, length(subjects)),
           target = rep(NA_integer_, length(subjects)), cells = numeric())
    } else {
      .dose_truncation_plan(source, source_roles, min_pattern_share)
    }
    dose_identifying <- dose_identifying & is.na(dose_plan$target)
    dose_maskable <- any(!dose_identifying)
    if (!dose_maskable && any(dose_identifying)) {
      .loud_warn(
        "no dose schedule can be masked",
        sprintf(paste("all %d patients have a set of dose times nobody else",
                      "shares."),
                length(subjects)),
        why = paste("a dose schedule is copied from its anchor verbatim --",
                    "resampling dose events would emit regimens the protocol",
                    "never permitted -- so every avatar carries one real",
                    "patient's exact dosing. There is no patient to build on",
                    "instead, because every one of them is in this position."),
        fix = paste("nothing within this study. Declaring `nominal_time` helps",
                    "only if the dose times are protocol times recorded",
                    "loosely. Otherwise treat the output as individually",
                    "identifying on dosing, and use the observation-side",
                    "guarantee for what it does cover.")
      )
    } else if (any(dose_identifying)) {
      .loud_warn(
        "patients whose dose schedule nobody else shares",
        sprintf("%s have a set of dose times no other patient has.",
                .count_phrase(sum(dose_identifying), length(subjects))),
        why = paste("dose events are copied from the anchor verbatim, so an",
                    "avatar built on one of these would carry that patient's",
                    "exact dosing. No avatar is anchored on them; they still",
                    "contribute as donors, so their measurements still shape",
                    "the output."),
        fix = paste("nothing, unless those regimens need to appear in the",
                    "synthetic data. `min_pattern_share = 1` keeps them and",
                    "gives up the guarantee."),
        kind = "NOTE"
      )
    }

    new_ids <- .new_ids(source[[source_roles$id]], n_subjects)
    # Donor floor: each avatar should blend `k` real patients, borrowing across
    # dose/schedule groups to reach it. The only unfixable shortfall is a source
    # with fewer than k + 1 subjects, so there are not k others to borrow.
    available_donors <- length(subjects) - 1L
    if (available_donors < as.integer(k)) {
      .loud_warn(
        "source too small for the donor floor",
        sprintf(paste("the source has %d patient%s, so every avatar is blended",
                      "from at most %d real patient%s -- fewer than the floor",
                      "of k = %d."),
                length(subjects), if (length(subjects) == 1L) "" else "s",
                max(available_donors, 0L),
                if (available_donors == 1L) "" else "s", as.integer(k)),
        why = paste("blending across few patients leaves each avatar close to",
                    "an individual, which markedly raises",
                    "re-identifiability."),
        fix = paste("use a larger source, or treat the output as individually",
                    "identifying and keep it under the source's own access",
                    "controls.")
      )
    }
    # Keep the output from looking extreme by default: do not anchor an avatar
    # on a source subject whose event structure is far beyond the cohort, since
    # the anchor's skeleton is copied verbatim. See .structural_outlier_anchors.
    # Screening uses no randomness, so a source with no such outlier yields
    # byte-identical output to `screen = FALSE`. Turn it off to keep every
    # structure.
    allowed <- seq_along(subjects)
    # Route is an absolute barrier (see `.route_key()`), so a subject whose
    # route arm holds fewer than k + 1 subjects can never be blended to the
    # floor: there is no legal donor left to borrow. What to do about it is the
    # caller's call, because only the caller knows whether a sparsely populated
    # arm matters more than the re-identification risk of reproducing it.
    # `on_donor_shortfall` picks; every branch is loud, because each silently
    # changes either what the cohort covers or how identifying it is.
    route_size <- table(profiles$routes)
    short_arms <- names(route_size)[route_size < as.integer(k) + 1L]
    route_excluded <- which(profiles$routes %in% short_arms)
    if (length(route_excluded)) {
      arm_summary <- paste0(
        short_arms, " (n=", as.integer(route_size[short_arms]), ")",
        collapse = "; "
      )
      # The most donors any of these subjects could actually get: arm size minus
      # itself. Zero means noise on a lone patient, which is worth saying out
      # loud rather than leaving the reader to work out.
      best_available <- max(as.integer(route_size[short_arms])) - 1L
      one_subject <- length(route_excluded) == 1L
      shortfall_context <- sprintf(
        paste0("%d patient%s sit%s in %d route arm%s holding fewer than the ",
               "donor floor of k = %d: %s."),
        length(route_excluded), if (one_subject) "" else "s",
        if (one_subject) "s" else "",
        length(short_arms), if (length(short_arms) == 1L) "" else "s",
        as.integer(k), arm_summary
      )
      if (identical(on_donor_shortfall, "error")) {
        stop(shortfall_context, " Donors are never blended across routes, so ",
             if (one_subject) "this patient has" else "these patients have",
             " no legal donor set.\n  Choose how to proceed with ",
             "`on_donor_shortfall`:\n",
             "  \"drop\"  (default) omit these subjects from the anchor pool; ",
             "the synthetic cohort will not represent the arm.\n",
             "  \"noise\" keep them, ",
             if (best_available > 0L) {
               sprintf(paste0("blending the %d donor%s available within their ",
                              "own route and relying on subject and residual ",
                              "noise for the rest"),
                       best_available, if (best_available == 1L) "" else "s")
             } else {
               "with no donor at all, so only subject and residual noise"
             },
             ". Not recommended: such a subject can stay close to one real ",
             "patient.",
             call. = FALSE)
      }
      shortfall_why <- paste("donors are never blended across routes, so",
                             if (one_subject) "this patient has"
                             else "these patients have",
                             "no legal donor set.")
      if (identical(on_donor_shortfall, "noise")) {
        .loud_warn(
          "route arm below the donor floor -- kept",
          shortfall_context,
          why = paste(shortfall_why,
                      sprintf(paste("Under `on_donor_shortfall = \"noise\"`",
                                    "they were kept anyway, built from at most",
                                    "%d donor%s plus subject and residual",
                                    "noise. Such an avatar can stay close to",
                                    "one real patient."),
                              max(best_available, 0L),
                              if (best_available == 1L) "" else "s")),
          fix = paste("treat these avatars as individually identifying and run",
                      "`flag_identifiable_subjects()` on the result. The",
                      "default `on_donor_shortfall = \"drop\"` omits them",
                      "instead.")
        )
      } else if (length(route_excluded) < length(subjects)) {
        allowed <- setdiff(allowed, route_excluded)
        .loud_warn(
          "route arm below the donor floor -- dropped",
          shortfall_context,
          why = paste(shortfall_why, "They were dropped from the anchor pool,",
                      "so the synthetic cohort does not represent",
                      if (length(short_arms) == 1L) "this arm."
                      else "these arms."),
          fix = paste("nothing, if the arm does not matter. To keep it, set",
                      "`on_donor_shortfall = \"noise\"` -- not recommended,",
                      "since such an avatar can stay close to one real",
                      "patient.")
        )
      } else {
        .loud_warn(
          "every route arm is below the donor floor",
          shortfall_context,
          why = paste(shortfall_why, "Dropping every arm would leave nothing",
                      "to generate, so generation proceeded as if",
                      "`on_donor_shortfall = \"noise\"`."),
          fix = paste("treat the whole output as individually identifying, or",
                      "use a larger source.")
        )
      }
    }
    # Two mechanisms remove subjects from the anchor pool, for different reasons
    # and with different consequences, and neither was counted anywhere. Both
    # are recorded below so a run can say how many real patients it declined to
    # build on and why.
    anchors_route_excluded <- length(setdiff(seq_along(subjects), allowed))
    anchors_screened_out <- 0L
    if (isTRUE(screen)) {
      excluded <- .structural_outlier_anchors(source, source_roles)
      # Tested against `allowed`, not the whole cohort: the route floor above
      # may already have removed anchors, and screening the rest to nothing
      # would leave no anchor to sample.
      if (length(excluded) && length(setdiff(allowed, excluded))) {
        anchors_screened_out <- length(intersect(allowed, excluded))
        allowed <- setdiff(allowed, excluded)
      } else if (length(excluded)) {
        warning("Screening would exclude every eligible source subject as a ",
                "structural outlier; it was skipped for this call.",
                call. = FALSE)
      }
    }
    # Anchors, drawn within stratum when there are strata to preserve.
    strata_balanced <- 0L
    strata_stochastic <- 0L
    strata_key <- if (length(source_roles$strata)) {
      keys <- lapply(source_roles$strata, function(column) {
        vapply(profiles$subject_rows, function(rows) {
          as.character(.first_present(source[[column]][rows]))
        }, character(1))
      })
      do.call(paste, c(keys, list(sep = "\r")))
    } else NULL
    targets <- if (isTRUE(preserve_strata_balance) && !is.null(strata_key)) {
      .strata_targets(strata_key, allowed, n_subjects)
    } else NULL
    if (length(targets)) {
      strata_balanced <- length(targets)
      strata_stochastic <- length(unique(strata_key)) - strata_balanced
      pooled <- allowed[!strata_key[allowed] %in% names(targets)]
      anchors <- unlist(c(
        # By position for the same reason as in `.strata_targets()`: a blank
        # stratum label cannot be looked up by name.
        lapply(seq_along(targets), function(i) {
          members <- allowed[strata_key[allowed] == names(targets)[[i]]]
          members[sample.int(length(members), targets[[i]], replace = TRUE)]
        }),
        if (n_subjects > sum(targets)) {
          # Whatever the balanced strata did not claim goes to the strata too
          # small to balance; with none of those left eligible, it goes back to
          # the whole pool rather than failing.
          remainder <- if (length(pooled)) pooled else allowed
          list(remainder[sample.int(length(remainder),
                                    n_subjects - sum(targets),
                                    replace = TRUE)])
        }
      ), use.names = FALSE)
      anchors <- anchors[sample.int(length(anchors))]
    } else {
      anchors <- allowed[sample.int(length(allowed), n_subjects,
                                    replace = TRUE)]
      if (!is.null(strata_key)) {
        strata_stochastic <- length(unique(strata_key))
      }
    }
    standard_mdv <- .source_uses_standard_mdv(source, source_roles)
    generated <- vector("list", n_subjects)
    # Inverse participation ratio of the donor weights, 1 / sum(w^2): the number
    # of donors an avatar is *effectively* blended from, which is what the
    # weight cap actually controls. It equals K for uniform weights and 1 for a
    # sole donor, and -- treating donors as independent -- the blend retains
    # sum(w^2) of individual variance, so it reads as a privacy floor and as the
    # between-subject-variability cost in the same number.
    effective_donors <- numeric(n_subjects)
    # How often the cap actually binds. A cap that fires on nearly every subject
    # is not a guardrail, it *is* the weighting scheme, and the inverse-distance
    # term underneath it is doing little; a cap that never fires is not
    # protecting anything. Reporting the rate is what makes that visible on real
    # data instead of inferable only by simulation.
    cap_bound <- logical(n_subjects)
    pattern_sampled <- logical(n_subjects)
    pattern_generated <- logical(n_subjects)
    # Whether the number of missed visits had to move to find a placement that
    # was not some real patient's. Recorded because it is the one thing the
    # attendance mechanism promises to preserve exactly and sometimes cannot.
    pattern_shifted <- logical(n_subjects)
    # An avatar keeping its anchor's own visit set is only a disclosure when
    # that set is one nobody else shares. Kept apart, because conflating them
    # made an ordinary shared schedule look like a failure.
    pattern_substituted <- logical(n_subjects)
    pattern_identifying <- logical(n_subjects)
    # Avatars moved to a different anchor because their first one could not be
    # given a visit set that masks it.
    pattern_reanchored <- logical(n_subjects)
    # Re-anchoring that had to leave the avatar's own stratum because nobody
    # else in it could be masked either. Counted separately from
    # `pattern_reanchored`: an in-stratum move costs nothing a reader cares
    # about, and this one costs the arm sizes.
    strata_crossed <- logical(n_subjects)
    # The visit set each avatar actually ended up with, and whose dose schedule
    # it inherited. Recorded per avatar so the end-to-end checks can be exact.
    # Snapping the finished table back onto the grid was tried and rejected: the
    # resampled deviations are comparable to the grid spacing on an irregular
    # study, so an avatar could snap into a neighbouring cell and be reported as
    # holding a schedule it was never given -- a false alarm on `wbcSim`.
    applied_key <- rep(NA_character_, n_subjects)
    dose_truncated <- logical(n_subjects)

    for (synthetic_index in seq_len(n_subjects)) {
      anchor <- anchors[synthetic_index]
      # Decide the visit set BEFORE building anything. An anchor whose own set
      # nobody shares, and for whom nothing legal can be placed, cannot be
      # masked -- any avatar built on it wears that one real patient's exact
      # pattern of which visits have observations. Rather than emit it, or drop
      # the patient from the cohort entirely, this avatar is simply anchored
      # somewhere else. Only the avatar moves; every source patient stays a
      # donor and stays available to anchor other avatars.
      choice <- NULL
      buildable <- function(key) {
        is.na(key) || !is.null(.apply_attendance(
          source[profiles$subject_rows[[anchor]], , drop = FALSE],
          source_roles, key
        ))
      }
      unsafe <- function(index) {
        (dose_maskable && isTRUE(dose_identifying[[index]])) ||
          isTRUE(choice$identifying)
      }
      # The stratum this avatar was allocated to, captured before any move. An
      # anchor carries its own `strata` values into the output, so re-anchoring
      # across arms silently undoes `.strata_targets()`: the run still reports
      # the balance it planned while the finished table holds something else.
      # Measured on a study with 8 arms and 6 unmaskable schedule groups, a
      # 19-patient arm came back with 5.
      target_stratum <- if (is.null(strata_key)) NULL else strata_key[[anchor]]
      if (!is.null(attendance) || dose_maskable) {
        for (attempt in seq_len(12L)) {
          choice <- if (is.null(attendance)) NULL else
            .attendance_choice(anchor, attendance)
          if (!is.null(choice) && !buildable(choice$key)) {
            choice$key <- NA_character_
            choice$identifying <- isTRUE(attendance$identifying[[anchor]])
          }
          if (!unsafe(anchor)) break
          move <- .reanchor_candidates(allowed, anchor, strata_key,
                                       target_stratum)
          if (move$crossed) strata_crossed[synthetic_index] <- TRUE
          if (!length(move$candidates)) break
          anchor <- move$candidates[[
            sample.int(length(move$candidates), 1L)
          ]]
          pattern_reanchored[synthetic_index] <- TRUE
        }
        anchors[synthetic_index] <- anchor
      }
      skeleton <- source[profiles$subject_rows[[anchor]], , drop = FALSE]
      original_order <- seq_len(nrow(skeleton))
      if (!is.null(attendance)) {
        applied <- if (is.na(choice$key)) NULL else
          .apply_attendance(skeleton, source_roles, choice$key)
        applied_key[synthetic_index] <- if (is.null(applied)) {
          attendance$keys[[anchor]]      # kept its anchor's own set
        } else choice$key
        if (!is.null(applied)) {
          skeleton <- applied
          original_order <- seq_len(nrow(skeleton))
          pattern_sampled[synthetic_index] <- TRUE
          pattern_generated[synthetic_index] <- choice$generated
          pattern_shifted[synthetic_index] <- isTRUE(choice$shifted)
          pattern_substituted[synthetic_index] <- isTRUE(choice$substituted)
        } else if (isTRUE(attendance$identifying[[anchor]])) {
          # Either nothing was drawn, or the drawn set could not be built on
          # this anchor. Both leave the anchor's own set in place, which is a
          # disclosure only when that set is one nobody else shares.
          pattern_identifying[synthetic_index] <- TRUE
        }
      }
      if (!is.na(dose_plan$target[[anchor]])) {
        skeleton <- .apply_dose_truncation(skeleton, source_roles,
                                           dose_plan$cells,
                                           dose_plan$target[[anchor]])
        original_order <- seq_len(nrow(skeleton))
        dose_truncated[synthetic_index] <- TRUE
      }
      skeleton <- .jitter_skeleton_time(skeleton, source_roles, time_jitter)
      if (length(time_deviations)) {
        skeleton <- .offset_unique_times(skeleton, source_roles, function(n) {
          sample(time_deviations, n, replace = TRUE)
        })
      }
      donors <- .select_donors(anchor, profiles, k, warnings, max_donor_weight)
      effective_donors[synthetic_index] <- 1 / sum(donors$weights^2)
      cap_bound[synthetic_index] <- length(donors$weights) > 1L &&
        max(donors$weights) >= max_donor_weight - sqrt(.Machine$double.eps)
      dose_level <- .dose_levels_for(skeleton, source_roles, dose_basis)
      skeleton <- .synthesize_covariates(
        skeleton, source, source_roles, donors$indices, donors$weights, profiles,
        subject_noise_sd
      )
      skeleton <- .apply_dose_basis(skeleton, source_roles, dose_basis,
                                    dose_level)
      skeleton <- .synthesize_trajectories(
        skeleton, source, source_roles, donors$indices, donors$weights, profiles,
        subject_noise_sd, residual_noise_sd, residual_phi, warnings,
        censoring_by_endpoint
      )
      skeleton <- .recompute_tad(skeleton, source_roles)
      if (standard_mdv) skeleton <- .derive_standard_mdv(skeleton, source_roles)
      new_id <- new_ids[synthetic_index]
      if (is.factor(skeleton[[source_roles$id]])) {
        new_label <- as.character(new_id)
        levels(skeleton[[source_roles$id]]) <- unique(c(
          levels(skeleton[[source_roles$id]]), new_label
        ))
        skeleton[[source_roles$id]][] <- new_label
      } else {
        skeleton[[source_roles$id]][] <- new_id
      }
      row_order <- order(skeleton[[source_roles$time]], original_order)
      generated[[synthetic_index]] <- skeleton[row_order, , drop = FALSE]
    }

    # An avatar that drew no visit set keeps its anchor's, which is one real
    # patient's exact pattern of absences -- precisely what `min_pattern_share`
    # exists to stop being reproduced. Two things cause it and neither alerted
    # before: a schedule group with no set shared widely enough (that one warns
    # above, from `attendance$pool`), and a group where every legal placement of
    # the shape happens to already be somebody's, so `.place_attendance()`
    # rejects all of them. The second is invisible until the run is over, which
    # is why the check is here rather than next to the pool. A small residue is
    # not worth interrupting for; a large one changes what the output is.
    # An avatar keeping its anchor's own visit set is NOT in itself a problem:
    # if six patients share that set, copying it discloses nothing. What matters
    # is the subset where the anchor's set was one nobody else had and nothing
    # could be substituted, because that pattern of absences belongs to exactly
    # one real person. Alerting on the wrong one of these made an ordinary run
    # look like a failure.

    result <- do.call(rbind, generated)
    rownames(result) <- NULL

    # The end-to-end check, run on the finished table rather than inferred from
    # the mechanisms that built it.
    #
    # Every leak found in this area so far leaked because each mechanism was
    # correct on its own terms and nobody asked the whole question afterwards:
    # placements were rejected properly but the caller then copied the anchor;
    # the pool was built properly but a one-patient group never got one. This
    # asks the question directly -- does any avatar's set of (endpoint, time)
    # cells match a source patient's set that fewer than `min_pattern_share`
    # patients share? -- so a future mechanism that gets it wrong is caught here
    # instead of in a report somebody reads months later.
    # The end-to-end checks, on what each avatar actually received.
    #
    # Every leak found in this area leaked because each mechanism was correct on
    # its own terms and nobody asked the whole question afterwards: placements
    # were rejected properly but the caller then copied the anchor; the pool was
    # built properly but a one-patient group never got one; a borrowed set could
    # not be built on its anchor and the failure was silent. These ask the
    # question directly, so the next mechanism that gets it wrong fails a test
    # rather than surfacing in a report months later.
    identifying_visit_sets <- 0L
    if (!is.null(attendance)) {
      rare_source <- unique(attendance$keys[attendance$identifying])
      rare_source <- rare_source[nzchar(rare_source)]
      identifying_visit_sets <- sum(!is.na(applied_key) &
                                      applied_key %in% rare_source)
    }
    # Dose events are copied from the anchor verbatim -- resampling them would
    # emit regimens the protocol never permitted -- so an avatar's dose schedule
    # simply IS its anchor's, and the question is whether that anchor's was one
    # nobody shared.
    identifying_dose_schedules <- sum(dose_identifying[anchors])
    # What declining to build on those patients COST. Not building on a patient
    # whose dose schedule nobody shares is the only safe answer, but it removes
    # that regimen from the output entirely, and until now nothing said so: a
    # cohort of nineteen patients on three doses, one on two and one on one came
    # out as twenty-one patients on three doses, silently.
    dose_regimens_source <- length(unique(
      .dose_schedule_keys(source, source_roles)
    ))
    dose_regimens_represented <- length(unique(
      .dose_schedule_keys(source, source_roles)[
        as.character(subjects)[unique(anchors)]
      ]
    ))
    if (identifying_visit_sets > 0L) {
      .loud_warn(
        "avatars carrying a visit set nobody else shares",
        sprintf(paste("%d of %d avatars were emitted holding a set of visits",
                      "that exactly one real patient has."),
                identifying_visit_sets, n_subjects),
        why = paste("that pattern of which visits have observations belongs to",
                    "one real person, so an avatar carrying it can be traced",
                    "back to them. Nothing legal was available to put in its",
                    "place."),
        fix = paste0(paste(.schedule_advice(source_roles, min_pattern_share,
                                            coarsened$grid),
                           collapse = "; "),
                     ". Until then, treat these avatars as individually ",
                     "identifying.")
      )
    }

    result <- .restore_schema(result, source, source_roles)
    # Measured on the finished table rather than trusted from the mechanism, for
    # the same reason as the two guarantees above: every part of this was
    # verified on its own terms once and the whole question was asked of nobody.
    # A discrete endpoint holding a value the source could not have held is a
    # defect in this package, not a property of the study, so it is a warning
    # rather than an alert.
    endpoint_type_violations <- .endpoint_type_violations(
      result, source_roles, profiles$value_types
    )
    if (sum(endpoint_type_violations) > 0L) {
      offending <- names(endpoint_type_violations)[
        endpoint_type_violations > 0L
      ]
      warnings$add(paste0(
        "Endpoint(s) ", paste(offending, collapse = ", "),
        " were generated with values off their source scale (",
        sum(endpoint_type_violations), " row(s)); please report this."
      ))
    }
    settings <- list(
      seed = as.integer(seed),
      n_subjects = n_subjects,
      event_method = event_method,
      dv_method = dv_method,
      k = as.integer(k),
      pca_variance = pca_variance,
      subject_noise_sd = subject_noise_sd,
      residual_noise_sd = residual_noise_sd,
      residual_phi = residual_phi,
      time_jitter = time_jitter,
      # Who was available to build on. `source_subjects` is the cohort as given;
      # the two exclusions are real patients this run declined to anchor any
      # avatar on, and `anchors_available` is what was left to sample from.
      # Cohort *size* is unaffected -- `n_subjects` avatars are still built --
      # so an excluded arm shows up as coverage lost, not as a shorter dataset.
      # The `screen` flag itself is deliberately *not* recorded. A clean source
      # must return byte-identical output whether screening is on or off
      # (`test-avatar-screen.R`), and storing the flag would break that identity
      # on the settings attribute while the data stayed the same. The count
      # below carries the information that matters and is 0 either way when
      # nothing was excluded.
      source_subjects = length(subjects),
      anchors_screened_out = as.integer(anchors_screened_out),
      anchors_route_excluded = as.integer(anchors_route_excluded),
      anchors_available = length(allowed),
      preserve_strata_balance = preserve_strata_balance,
      strata_balanced = as.integer(strata_balanced),
      strata_stochastic = as.integer(strata_stochastic),
      # Avatars that ended up in an arm other than the one they were allocated
      # to. Any value above 0 is why C3 does not add up.
      strata_crossed = as.integer(sum(strata_crossed)),
      coarsen_time = coarsen_time,
      min_pattern_share = as.integer(min_pattern_share),
      pattern_sampled_fraction = mean(pattern_sampled),
      # Of the sampled patterns, how many were built from a shape rather than
      # reused from a real subject.
      pattern_generated_fraction = mean(pattern_generated),
      # Of the avatars given a visit set, how many had the miss count moved.
      pattern_shifted_fraction = mean(pattern_shifted),
      # Avatars whose anchor's own set was too rare to reuse and was swapped
      # for the group's most widely held one.
      pattern_substituted_fraction = mean(pattern_substituted),
      # The number that must be zero: avatars emitted carrying a visit set no
      # real patient shares. Everything above is fidelity; this is disclosure.
      pattern_identifying_fraction = mean(pattern_identifying),
      # Measured on the finished table, not inferred from the mechanisms. This
      # is the number the guarantee is stated in terms of, and it must be 0.
      identifying_visit_sets = as.integer(identifying_visit_sets),
      identifying_dose_schedules = as.integer(identifying_dose_schedules),
      dose_truncated_fraction = mean(dose_truncated),
      dose_regimens_source = as.integer(dose_regimens_source),
      dose_regimens_represented = as.integer(dose_regimens_represented),
      pattern_reanchored_fraction = mean(pattern_reanchored),
      # What the floor cost: source attendance patterns excluded from the pool,
      # and how many real subjects held them.
      # Note that at the default floor of 2 a dropped pattern is *by definition*
      # held by exactly one subject, so `patterns_dropped` and
      # `subjects_with_dropped_pattern` are necessarily equal there. They
      # diverge only at 3 or more.
      patterns_total = if (is.null(attendance)) NA_integer_ else
        as.integer(attendance$total_patterns),
      patterns_dropped = if (is.null(attendance)) 0L else
        as.integer(attendance$dropped_patterns),
      subjects_with_dropped_pattern = if (is.null(attendance)) 0L else
        as.integer(attendance$dropped_subjects),
      time_grid = coarsened$grid,
      time_grid_derived_rows = as.integer(coarsened$grid_derived_rows %||% 0L),
      time_grid_rows = as.integer(coarsened$grid_rows %||% 0L),
      # How many source subjects still hold a schedule nobody else shares, once
      # coarsening has done what it can, and why. Recorded rather than only
      # alerted so a report can tabulate it without recomputing.
      dose_basis = if (is.null(dose_basis)) NA_character_ else
        dose_basis$covariate,
      # Why the answer above is what it is. Without this, "not detected" and
      # "detection was never attempted" are the same blank in a report, and a
      # caller whose study *is* weight-based has no way to tell which happened.
      dose_basis_note = dose_basis_detection$note,
      dose_levels = if (is.null(dose_basis)) NA_real_ else dose_basis$levels,
      # TRUE when the caller named the covariate rather than the run inferring
      # it. The two behave differently -- declared holds each dose row's own
      # ratio, inferred snaps to shared protocol levels -- so a report that
      # says "yes" needs to say which.
      dose_basis_declared = !is.null(source_roles$dose_covariate),
      unique_schedule_n = unique_schedule_n,
      unique_obs_time_n = unique_obs_time_n,
      unique_visit_set_n = if (is.na(unique_schedule_n)) NA_integer_ else
        unique_schedule_n - unique_obs_time_n,
      time_deviation_sd = if (length(time_deviations)) {
        stats::sd(time_deviations)
      } else 0,
      roles = unclass(source_roles),
      endpoint_transforms = profiles$transforms,
      # Which endpoints were treated as discrete, what levels they were snapped
      # to, and whether the caller declared that or the run inferred it.
      # `pmx_endpoint_types()` prints the same decision as a table.
      endpoint_value_types = profiles$value_types,
      endpoint_type_violations = as.integer(sum(endpoint_type_violations)),
      alignment = paste(
        "time relative to first positive dose within compatible schedules;",
        "normalized observation-window fallback"
      ),
      compatible_event_groups = length(unique(profiles$signatures)),
      routes = length(unique(profiles$routes)),
      on_donor_shortfall = on_donor_shortfall,
      max_donor_weight = max_donor_weight,
      cap_binding_fraction = mean(cap_bound),
      mean_effective_donors = mean(effective_donors),
      min_effective_donors = min(effective_donors),
      warnings = warnings$messages
    )
    attr(result, "pmx_settings") <- settings
    validate_pmx(result, source_roles, strict = TRUE)
    if (length(warnings$messages)) {
      warning(
        "Synthetic generation used documented small-group/profile fallbacks:\n",
        paste(vapply(warnings$messages, .wrap_plain, character(1),
                     initial = "- ", prefix = "  "),
              collapse = "\n"),
        call. = FALSE
      )
    }
    result
  })
}
