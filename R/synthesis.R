.validate_generator_options <- function(n_subjects, source_n, event_method,
                                        dv_method, k, pca_variance,
                                        subject_noise_sd, residual_noise_sd,
                                        residual_phi, time_jitter,
                                        max_donor_weight, coarsen_time) {
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
# Times are pooled and cut into cells wherever the gap between consecutive
# distinct values exceeds a threshold -- single-linkage clustering in one
# dimension, exact and O(n log n). The threshold is searched for rather than
# assumed, because there is no defensible constant: it is the largest cut that
# still preserves every subject's own visit count (see below).
.derive_time_grid <- function(times, subject_index) {
  keep <- is.finite(times)
  times <- times[keep]
  subject_index <- subject_index[keep]
  distinct <- sort(unique(times))
  if (length(distinct) < 2L) return(distinct)

  # The coarsest grid that still preserves every subject's own visit count. Two
  # pooled times are the same nominal visit exactly when merging them never puts
  # two of one subject's visits in one cell -- that constraint is the definition
  # of "too coarse", and the largest threshold satisfying it is the most
  # collapsing this can safely do.
  #
  # A fixed fraction of the smallest within-subject gap was tried first and is
  # far too conservative: one subject anywhere with a tightly spaced pair drags
  # the threshold to near zero and nothing merges at all, which
  # `scripts/measure_skeleton_uniqueness.R` shows plainly on the public
  # datasets. Searching for the constraint boundary is immune to that, because
  # one tight pair only pins the cells around itself.
  pairs <- unique(data.frame(subject = subject_index, time = times,
                             stringsAsFactors = FALSE))
  preserves_visits <- function(threshold) {
    cluster <- cumsum(c(TRUE, diff(distinct) > threshold))
    !anyDuplicated(data.frame(
      subject = pairs$subject,
      cell = cluster[match(pairs$time, distinct)],
      stringsAsFactors = FALSE
    ))
  }
  candidates <- sort(unique(diff(distinct)))
  candidates <- candidates[candidates > 0]
  if (!length(candidates)) return(distinct)
  # Binary search the largest candidate that still preserves visit counts. Cells
  # only ever merge as the threshold grows, so validity is monotone in it.
  low <- 0L
  high <- length(candidates)
  while (low < high) {
    middle <- (low + high + 1L) %/% 2L
    if (preserves_visits(candidates[middle])) low <- middle else high <- middle - 1L
  }
  if (low == 0L) return(distinct)
  cluster <- cumsum(c(TRUE, diff(distinct) > candidates[low]))
  if (max(cluster) == length(distinct)) return(distinct)
  # Centre each visit on the mean of the times actually recorded there, so a
  # visit sampled by many patients is not pulled by a sparse neighbour.
  member <- cluster[match(times, distinct)]
  as.numeric(tapply(times, member, mean))
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
    subject_properties = "carry the column with `keep`",
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
    max_donor_weight, coarsen_time
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
    if (coarsen_time) {
      still_alone <- sum(skeleton_uniqueness(source, source_roles)$alone)
      if (still_alone > 0L) {
        .loud_warn(sprintf(
          paste0("`coarsen_time = TRUE` left %d of %d source subject%s holding ",
                 "the only copy of its observation schedule, so an avatar ",
                 "anchored there still carries an identifying visit pattern. ",
                 "Declare a `nominal_time` role to snap to the protocol grid ",
                 "exactly; `scripts/measure_skeleton_uniqueness.R` shows what ",
                 "the grid did and did not collapse."),
          still_alone, length(subjects),
          if (still_alone == 1L) "" else "s"
        ))
      }
    }
    profiles <- .build_profiles(source, source_roles, pca_variance)
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

    for (synthetic_index in seq_len(n_subjects)) {
      anchor <- anchors[synthetic_index]
      skeleton <- source[profiles$subject_rows[[anchor]], , drop = FALSE]
      original_order <- seq_len(nrow(skeleton))
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
      skeleton <- .synthesize_covariates(
        skeleton, source, source_roles, donors$indices, donors$weights, profiles,
        subject_noise_sd
      )
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
      time_grid = coarsened$grid,
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
