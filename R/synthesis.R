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
  none <- list(source = source, deviations = numeric(), grid = "none")
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
  list(source = source, deviations = deviations, grid = kind)
}

# The stratum a subject was assigned to: the combination of every
# `subject_properties` column, or one shared stratum when none is declared.
# Protocol properties -- what dose level the arm received, which visits the
# schedule called for -- are constant within a stratum and differ between them,
# so it is the right grouping for both the dose basis below and the attendance
# pattern pool. It is deliberately *not* a blending barrier; only `.route_key()`
# is.
.subject_strata <- function(data, roles) {
  if (!length(roles$subject_properties)) {
    return(rep("all", nrow(data)))
  }
  do.call(paste, c(lapply(roles$subject_properties, function(column) {
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
.detect_dose_basis <- function(source, roles, tolerance = 0.02,
                               max_levels = 10L) {
  if (is.null(roles$amt) || !length(roles$covariates)) return(NULL)
  amount <- suppressWarnings(as.numeric(source[[roles$amt]]))
  usable <- .dose_rows(source, roles) & is.finite(amount) & amount > 0
  if (sum(usable) < 3L) return(NULL)
  amount <- amount[usable]
  distinct_amounts <- length(unique(signif(amount, 8)))
  # Flat dosing: nothing varies, so nothing to explain and nothing identifying
  # about the amount.
  if (distinct_amounts < 2L) return(NULL)

  for (covariate in roles$covariates) {
    if (is.factor(source[[covariate]])) next
    values <- suppressWarnings(as.numeric(source[[covariate]]))
    if (!is.numeric(values)) next
    values <- values[usable]
    if (any(!is.finite(values) | values <= 0)) next

    levels <- .dose_ratio_levels(amount / values, tolerance)
    if (is.null(levels)) next
    # Ratios must be far more concentrated than the amounts they came from.
    # A protocol has a handful of dose levels and many patients, so genuine
    # proportional dosing collapses dozens of distinct amounts onto a few
    # ratios; a study where dose is unrelated to the covariate produces about as
    # many ratios as amounts and is refused here.
    if (length(levels) > max_levels) next
    if (length(levels) * 2L > sum(usable)) next
    if (distinct_amounts < 2L * length(levels)) next
    return(list(covariate = covariate, levels = levels))
  }
  NULL
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
  level[dosed] <- basis$levels[
    max.col(-abs(outer(ratio, basis$levels, "-")), ties.method = "first")
  ]
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
# correlation that makes dropout dropout -- once a subject discontinues they are
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
  paste(sort(unique(paste0(endpoint, "@",
                           format(time, digits = 12, trim = TRUE)))),
        collapse = ";")
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
    shapes <- vapply(names(counts), function(key) {
      .attendance_shape(key, cells)
    }, character(1))
    shape_counts <- tapply(as.integer(counts), shapes, sum)
    shape_kept <- shape_counts[shape_counts >= min_pattern_share]

    # What is lost is now what neither an exact pattern nor its shape can carry.
    lost <- !kept & !(shapes %in% names(shape_kept))
    dropped_patterns <- dropped_patterns + sum(lost)
    dropped_subjects <- dropped_subjects + sum(as.integer(counts[lost]))
    if (!length(shape_kept)) next

    pool[[name]] <- list(
      cells = cells,
      shapes = names(shape_kept), shape_weight = as.numeric(shape_kept),
      # Exact patterns that clear the floor on their own, indexed by shape, so a
      # real pattern is preferred wherever reusing one is safe.
      exact = split(names(counts)[kept], shapes[kept]),
      exact_weight = split(as.numeric(counts[kept]), shapes[kept]),
      # Patterns too rare to reuse. A generated placement that lands on one of
      # these is rejected, so "nothing was copied" does not quietly weaken into
      # "something was reproduced by accident".
      rare = names(counts)[!kept]
    )
  }
  list(pool = pool, group = group, keys = keys,
       total_patterns = total_patterns, dropped_patterns = dropped_patterns,
       dropped_subjects = dropped_subjects)
}

# The union of every cell any subject in the group was observed at, in time
# order. Time first, so "trailing" means late in the study -- dropout -- rather
# than late in an alphabet.
.attendance_cells <- function(keys) {
  cells <- unique(unlist(strsplit(keys, ";", fixed = TRUE), use.names = FALSE))
  cells <- cells[nzchar(cells)]
  time <- suppressWarnings(as.numeric(sub("^.*@", "", cells)))
  cells[order(time, cells)]
}

# Describe a pattern by *how* it departs from full attendance rather than by
# which cells it holds: the number of missed visits and their arrangement.
#
#   complete   attended everything
#   trailing   every miss is at the end -- dropout or early discontinuation
#   block      the misses are contiguous but not terminal -- an interruption
#   scattered  anything else
#
# This is the whole point of the fallback. The exact-cell definition cannot see
# that "missed one visit, early" and "missed one visit, late" are the same kind
# of event, so each is a group of one and both are discarded. As a shape they are
# one group of two.
.attendance_shape <- function(key, cells) {
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

.attendance_key_from <- function(cells, missed) {
  paste(sort(cells[-missed]), collapse = ";")
}

# Place a shape's misses on the grid without copying anyone's arrangement, and
# refuse any placement that happens to reproduce a pattern too rare to have been
# reusable. Without that rejection the guarantee would weaken from "no synthetic
# patient carries a schedule unique to a real one" to merely "nothing was copied
# on purpose", and an attacker cannot tell those apart.
.place_attendance <- function(cells, shape, rare, tries = 24L) {
  parts <- strsplit(shape, "|", fixed = TRUE)[[1L]]
  kind <- parts[[1L]]
  n_miss <- as.integer(parts[[2L]])
  n <- length(cells)
  if (kind == "complete" || n_miss < 1L) return(paste(sort(cells), collapse = ";"))
  if (n_miss >= n) return(NA_character_)

  for (attempt in seq_len(tries)) {
    missed <- switch(
      kind,
      trailing = seq.int(n - n_miss + 1L, n),
      block = {
        # Contiguous but not terminal, so a block stays distinguishable from
        # dropout.
        last_start <- n - n_miss
        if (last_start < 1L) return(NA_character_)
        start <- sample.int(last_start, 1L)
        seq.int(start, start + n_miss - 1L)
      },
      scattered = sort(sample.int(n, n_miss)),
      return(NA_character_)
    )
    if (kind == "scattered") {
      # Keep the kind honest: a scattered draw that came out contiguous or
      # terminal is a different shape, so redraw.
      if (all(diff(missed) == 1L) || all(missed > n - n_miss)) next
    }
    key <- .attendance_key_from(cells, missed)
    if (!(key %in% rare)) return(key)
  }
  NA_character_
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
                generated = FALSE))
  }
  key <- .place_attendance(available$cells, shape, available$rare)
  list(key = key, generated = !is.na(key))
}

# Rebuild the skeleton's observation rows to match a sampled pattern, keeping
# every dose row exactly as the anchor had it. Rows for a visit the anchor did
# not attend are cloned from one of its own observation rows for that endpoint,
# so compartment, MDV, and every carried column stay coherent; the DV itself is
# filled by the donor blend afterwards, exactly as for a copied skeleton.
.apply_attendance <- function(skeleton, roles, key) {
  observed <- .observation_rows(skeleton, roles, require_present = TRUE)
  if (!any(observed) || !nzchar(key)) return(skeleton)
  wanted <- strsplit(key, ";", fixed = TRUE)[[1L]]
  wanted_endpoint <- sub("@.*$", "", wanted)
  wanted_time <- as.numeric(sub("^.*@", "", wanted))
  endpoint <- .endpoint(skeleton, roles)

  template <- vapply(unique(wanted_endpoint), function(name) {
    candidate <- which(observed & endpoint == name)
    if (length(candidate)) candidate[1L] else NA_integer_
  }, integer(1))
  names(template) <- unique(wanted_endpoint)
  if (anyNA(template)) return(skeleton)

  built <- skeleton[template[wanted_endpoint], , drop = FALSE]
  built[[roles$time]] <- wanted_time
  rbind(skeleton[!observed, , drop = FALSE], built)
}

# TAD is time since the most recent dose, so it goes stale the moment TIME
# moves. Recomputed from the generated skeleton rather than carried over from
# the anchor, where it would describe a schedule the avatar no longer has.
.recompute_tad <- function(skeleton, roles) {
  if (is.null(roles$tad)) return(skeleton)
  time <- suppressWarnings(as.numeric(skeleton[[roles$time]]))
  dose_time <- sort(time[.dose_rows(skeleton, roles) & is.finite(time)])
  if (!length(dose_time)) return(skeleton)
  index <- findInterval(time, dose_time)
  # Before the first dose there is no elapsed time to report; zero is the only
  # value `validate_pmx()` accepts and the only one that is not a fiction.
  tad <- ifelse(index >= 1L, time - dose_time[pmax(index, 1L)], 0)
  tad[!is.finite(tad)] <- 0
  skeleton[[roles$tad]] <- pmax(tad, 0)
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

.interpolate_trajectory <- function(trajectory, target_time) {
  time <- trajectory$time
  value <- trajectory$value
  if (!length(time)) return(rep(NA_real_, length(target_time)))
  if (length(unique(time)) == 1L) return(rep(mean(value), length(target_time)))

  absolute <- stats::approx(time, value, xout = target_time,
                            ties = mean, rule = 1)$y
  missing <- !is.finite(absolute)
  if (any(missing) && length(unique(target_time)) > 1L) {
    target_range <- range(target_time)
    donor_range <- range(time)
    target_fraction <- (target_time[missing] - target_range[1L]) /
      diff(target_range)
    mapped_time <- donor_range[1L] + target_fraction * diff(donor_range)
    absolute[missing] <- stats::approx(
      time, value, xout = mapped_time, ties = mean, rule = 2
    )$y
  }
  absolute
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
      fallback <- .transform_dv(
        data[[roles$dv]][source_observed & source_endpoint == endpoint_name],
        transform
      )
      fallback <- stats::median(fallback[is.finite(fallback)], na.rm = TRUE)
      if (!is.finite(fallback)) {
        stop("No usable DV values exist for endpoint `", endpoint_name, "`.",
             call. = FALSE)
      }
      blended[missing_blend] <- fallback
      warnings$add(paste0(
        "Endpoint `", endpoint_name,
        "` required a dataset-median interpolation fallback."
      ))
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

# A deliberately loud, immediate alert for the one case pooling cannot fix: a
# source with fewer subjects than the donor floor, so every avatar is blended
# from too few real patients. Unlike the collected end-of-run warning() (which
# `suppressWarnings()` removes -- as the demo did, hiding it), the message()
# banner survives suppressWarnings and prints at once, in red on an interactive
# console. A real warning condition is also raised so it is catchable and shows
# in non-interactive logs.
.loud_warn <- function(msg) {
  banner <- paste0("SYNPMX ALERT: ", msg)
  if (interactive()) banner <- paste0("\033[1;31m", banner, "\033[0m")
  message(banner)
  warning(msg, call. = FALSE, immediate. = TRUE)
  invisible(NULL)
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
#'   are drawn rather than individual visits, which keeps dropout monotone
#'   instead of producing implausible attend/miss/attend sequences.
#'
#'   Unlike dropping the exposed subjects, nobody leaves the cohort: a subject
#'   with a rare pattern still contributes measurements as a donor, only their
#'   distinctive absences stop being reproduced. Dose events are never sampled,
#'   since that could emit a regimen no protocol permits. Raising this hides more
#'   and flattens the cohort's missingness further; where no pattern is shared
#'   widely enough, anchors keep their own and the run alerts loudly. Pools are
#'   formed within each `subject_properties` stratum and endpoint set.
#'
#'   **Patterns below the floor are lost, not approximated.** A dropout or
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
                     on_donor_shortfall = c("drop", "noise", "error")) {
  on_donor_shortfall <- match.arg(on_donor_shortfall)
  if (!is.data.frame(data)) stop("`data` must be a data frame or tibble.",
                                 call. = FALSE)
  .assert_roles(data, roles)
  .reject_dp_only_roles(roles)
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
      list(source = source, deviations = numeric(), grid = "off")
    }
    source <- coarsened$source
    time_deviations <- coarsened$deviations
    # Coarsening onto a derived grid is best-effort: where subjects hold samples
    # closer together than the spread between subjects, no grid can merge them
    # without collapsing a real visit, so the schedule stays unique and the
    # verbatim skeleton copy stays identifying. That is the safe failure, but it
    # is a silent one, and a caller who asked for `coarsen_time` and did not get
    # it should hear so rather than infer it from the absence of an alert.
    exposure_alone <- NA_integer_
    exposure_unique_moment <- NA_integer_
    if (coarsen_time) {
      exposure <- skeleton_uniqueness(source, source_roles)
      still_alone <- sum(exposure$alone)
      # Two causes, opposite remedies, so the alert has to say which. A subject
      # observed at a moment nobody else was is a grid failure -- declaring
      # `nominal_time` fixes it. A subject whose every time is shared and whose
      # *attendance pattern* is unique is dropout, and no grid however fine or
      # coarse touches that; it is the outlier screen's problem.
      unshared <- sum(exposure$min_time_share == 1L, na.rm = TRUE)
      pattern_only <- max(still_alone - unshared, 0L)
      exposure_alone <- still_alone
      exposure_unique_moment <- unshared
      if (unshared > 0L) {
        .loud_warn(sprintf(
          paste0("`coarsen_time = TRUE` left %d of %d source subject%s observed ",
                 "at a moment no other subject shares, so the grid could not ",
                 "hide the schedule and an avatar anchored there carries it. ",
                 "Declare a `nominal_time` role to snap to the protocol grid ",
                 "exactly; `scripts/measure_skeleton_uniqueness.R` shows what ",
                 "the grid did and did not collapse."),
          unshared, length(subjects), if (unshared == 1L) "" else "s"
        ))
      }
      if (pattern_only > 0L) {
        .loud_warn(sprintf(
          paste0("%d of %d source subject%s share every observation time with ",
                 "others but hold a unique *pattern* of which visits were ",
                 "attended -- dropout, discontinuation, or missed visits. ",
                 "Coarsening cannot change this, because the times are already ",
                 "shared. Screen those subjects with ",
                 "`flag_identifiable_subjects()` and ",
                 "`remediate_identifiable_subjects()` if the pattern matters."),
          pattern_only, length(subjects), if (pattern_only == 1L) "" else "s"
        ))
      }
    }
    dose_basis <- .detect_dose_basis(source, source_roles)
    profiles <- .build_profiles(source, source_roles, pca_variance)
    attendance <- .attendance_pool(source, source_roles, profiles,
                                   as.integer(min_pattern_share))
    if (!is.null(attendance) && attendance$dropped_patterns > 0L) {
      .loud_warn(sprintf(
        paste0("%d source attendance pattern(s), held by %d subject(s), were ",
               "not shared by %d or more patients and so will not appear in ",
               "the synthetic data. These are real dropout and dose-",
               "interruption patterns, and losing them is what stops an avatar ",
               "carrying a schedule traceable to one patient. To keep more of ",
               "them, lower `min_pattern_share` (2 means no synthetic patient ",
               "carries a schedule unique to a real one); `1` disables the ",
               "mechanism and copies each anchor's own pattern."),
        attendance$dropped_patterns, attendance$dropped_subjects,
        as.integer(min_pattern_share)
      ))
    }
    if (!is.null(attendance)) {
      unpooled <- setdiff(unique(attendance$group), names(attendance$pool))
      if (length(unpooled)) {
        stranded <- sum(attendance$group %in% unpooled)
        .loud_warn(sprintf(
          paste0("no attendance pattern is shared by %d or more subjects in ",
                 "%d of %d group(s), covering %d subject(s), so avatars ",
                 "anchored there keep their anchor's own pattern of attended ",
                 "visits -- which is what `min_pattern_share` exists to avoid ",
                 "reproducing. Lower `min_pattern_share`, or declare a ",
                 "`nominal_time` role so coarsening can put more subjects on ",
                 "the same visits."),
          as.integer(min_pattern_share), length(unpooled),
          length(unique(attendance$group)), stranded
        ))
      }
    }
    new_ids <- .new_ids(source[[source_roles$id]], n_subjects)
    # Donor floor: each avatar should blend `k` real patients, borrowing across
    # dose/schedule groups to reach it. The only unfixable shortfall is a source
    # with fewer than k + 1 subjects, so there are not k others to borrow.
    available_donors <- length(subjects) - 1L
    if (available_donors < as.integer(k)) {
      .loud_warn(sprintf(
        paste0("the source has %d subject%s, so every avatar is blended from ",
               "at most %d real patient%s -- fewer than the floor of %d. This ",
               "markedly raises re-identifiability; use a larger source or ",
               "treat the output as individually identifying."),
        length(subjects), if (length(subjects) == 1L) "" else "s",
        max(available_donors, 0L), if (available_donors == 1L) "" else "s",
        as.integer(k)
      ))
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
        paste0("%d subject%s in %d route arm%s below the donor floor of %d: ",
               "%s. Donors are never blended across routes, so %s no legal ",
               "donor set."),
        length(route_excluded), if (one_subject) "" else "s",
        length(short_arms), if (length(short_arms) == 1L) "" else "s",
        as.integer(k), arm_summary,
        if (one_subject) "this subject has" else "these subjects have"
      )
      if (identical(on_donor_shortfall, "error")) {
        stop(shortfall_context, "\n  Choose how to proceed with ",
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
      if (identical(on_donor_shortfall, "noise")) {
        .loud_warn(sprintf(
          paste0("%s They were kept anyway under `on_donor_shortfall = ",
                 "\"noise\"`, generated from at most %d donor%s plus subject ",
                 "and residual noise. Such avatars can stay close to one real ",
                 "patient; treat them as individually identifying, and use ",
                 "`flag_identifiable_subjects()` on the result. The default ",
                 "`on_donor_shortfall = \"drop\"` omits them instead."),
          shortfall_context, max(best_available, 0L),
          if (best_available == 1L) "" else "s"
        ))
      } else if (length(route_excluded) < length(subjects)) {
        allowed <- setdiff(allowed, route_excluded)
        .loud_warn(sprintf(
          paste0("%s They were dropped from the anchor pool, so the synthetic ",
                 "cohort does not represent %s arm%s. To keep them anyway, set ",
                 "`on_donor_shortfall = \"noise\"` -- not recommended, since ",
                 "such a subject can stay close to one real patient."),
          shortfall_context,
          if (length(short_arms) == 1L) "this" else "these",
          if (length(short_arms) == 1L) "" else "s"
        ))
      } else {
        .loud_warn(sprintf(
          paste0("%s Dropping every arm would leave nothing to generate, so ",
                 "generation proceeded as if `on_donor_shortfall = \"noise\"`. ",
                 "Treat the output as individually identifying."),
          shortfall_context
        ))
      }
    }
    if (isTRUE(screen)) {
      excluded <- .structural_outlier_anchors(source, source_roles)
      # Tested against `allowed`, not the whole cohort: the route floor above
      # may already have removed anchors, and screening the rest to nothing
      # would leave no anchor to sample.
      if (length(excluded) && length(setdiff(allowed, excluded))) {
        allowed <- setdiff(allowed, excluded)
      } else if (length(excluded)) {
        warning("Screening would exclude every eligible source subject as a ",
                "structural outlier; it was skipped for this call.",
                call. = FALSE)
      }
    }
    anchors <- allowed[sample.int(length(allowed), n_subjects, replace = TRUE)]
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

    for (synthetic_index in seq_len(n_subjects)) {
      anchor <- anchors[synthetic_index]
      skeleton <- source[profiles$subject_rows[[anchor]], , drop = FALSE]
      original_order <- seq_len(nrow(skeleton))
      if (!is.null(attendance)) {
        available <- attendance$pool[[attendance$group[[anchor]]]]
        if (!is.null(available)) {
          drawn <- .draw_attendance(available)
          if (!is.na(drawn$key)) {
            skeleton <- .apply_attendance(skeleton, source_roles, drawn$key)
            original_order <- seq_len(nrow(skeleton))
            pattern_sampled[synthetic_index] <- TRUE
            pattern_generated[synthetic_index] <- drawn$generated
          }
        }
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

    result <- do.call(rbind, generated)
    rownames(result) <- NULL
    result <- .restore_schema(result, source, source_roles)
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
      coarsen_time = coarsen_time,
      min_pattern_share = as.integer(min_pattern_share),
      pattern_sampled_fraction = mean(pattern_sampled),
      # Of the sampled patterns, how many were built from a shape rather than
      # reused from a real subject.
      pattern_generated_fraction = mean(pattern_generated),
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
      # How many source subjects still hold a schedule nobody else shares, once
      # coarsening has done what it can, and why. Recorded rather than only
      # alerted so a report can tabulate it without recomputing.
      dose_basis = if (is.null(dose_basis)) NA_character_ else
        dose_basis$covariate,
      dose_levels = if (is.null(dose_basis)) NA_real_ else dose_basis$levels,
      exposure_alone = exposure_alone,
      exposure_unique_moment = exposure_unique_moment,
      exposure_unique_pattern = if (is.na(exposure_alone)) NA_integer_ else
        exposure_alone - exposure_unique_moment,
      time_deviation_sd = if (length(time_deviations)) {
        stats::sd(time_deviations)
      } else 0,
      roles = unclass(source_roles),
      endpoint_transforms = profiles$transforms,
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
        "Synthetic generation used documented small-group/profile fallbacks:\n- ",
        paste(warnings$messages, collapse = "\n- "),
        call. = FALSE
      )
    }
    result
  })
}
