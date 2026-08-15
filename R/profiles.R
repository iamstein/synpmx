# Documented in `avatar-algorithm.Rmd`, Step 5. Change one, change the other.
.event_signature <- function(subject_data, roles) {
  event <- !.is_zero(subject_data[[roles$evid]])
  observed <- .observation_rows(subject_data, roles)
  endpoint <- .endpoint(subject_data, roles)

  event_token <- if (any(event)) {
    pieces <- list(as.character(subject_data[[roles$evid]][event]))
    for (role in c("cmt", "dvid")) {
      if (!is.null(roles[[role]])) {
        pieces[[length(pieces) + 1L]] <-
          as.character(subject_data[[roles[[role]][[1L]]]][event])
      }
    }
    for (role in c("amt", "rate")) {
      if (!is.null(roles[[role]])) {
        value <- suppressWarnings(as.numeric(subject_data[[roles[[role]]]][event]))
        # Keep regimen magnitude in compatibility. Otherwise an anchor dose
        # can be paired with a donor trajectory from a different dose.
        pieces[[length(pieces) + 1L]] <- ifelse(
          is.finite(value), format(signif(value, 8), trim = TRUE,
                                   scientific = FALSE), "NA"
        )
      }
    }
    paste(do.call(paste, c(pieces, sep = ":")), collapse = ";")
  } else {
    "none"
  }
  start_event <- event
  if (!is.null(roles$amt)) {
    amount <- subject_data[[roles$amt]]
    positive_event <- event & !is.na(amount) & amount > 0
    if (any(positive_event)) start_event <- positive_event
  }
  start_time <- as.numeric(subject_data[[roles$time]][start_event])
  schedule_token <- if (length(start_time) <= 1L) {
    paste0("starts=", length(start_time))
  } else {
    # Two significant digits tolerate nominal-time noise while keeping clearly
    # different repeat-dose/occasion intervals in separate donor groups.
    interval <- signif(diff(start_time), 2L)
    paste0(
      "starts=", length(start_time), ":gaps=",
      paste(format(interval, trim = TRUE, scientific = FALSE), collapse = ",")
    )
  }
  endpoint_token <- paste(sort(unique(endpoint[observed])), collapse = ",")
  paste0(
    "events=", event_token, "|schedule=", schedule_token,
    "|endpoints=", endpoint_token
  )
}

# The administration route, and the one structural axis that is an absolute
# barrier: an avatar is never blended from a donor dosed by a different route.
# Every other difference -- dose size, number of doses, interval -- makes a
# donor a worse match, and the fallback in `.select_donors()` is allowed to
# accept it. A route difference is not a worse match, it is a different
# experiment: a bolus, an infusion, and an oral dose produce categorically
# different concentration shapes, so blending across them yields a trajectory no
# protocol could have produced. Read from the dosing rows only: which EVID and
# compartment the dose enters, and whether it is delivered over time. NONMEM's
# RATE < 0 (modeled rate or duration) is still an infusion. The tokens form a
# *set*, so three oral doses and five oral doses share a route -- dose count is
# a schedule difference, handled by the signature, not a route difference.
# Documented in `avatar-algorithm.Rmd`, Step 5. Change one, change the other.
.route_key <- function(subject_data, roles) {
  dosed <- .dose_rows(subject_data, roles)
  if (!any(dosed)) return("none")
  pieces <- list(as.character(subject_data[[roles$evid]][dosed]))
  if (!is.null(roles$cmt)) {
    pieces[[length(pieces) + 1L]] <-
      as.character(subject_data[[roles$cmt[[1L]]]][dosed])
  }
  delivery <- rep("bolus", sum(dosed))
  if (!is.null(roles$rate)) {
    rate <- suppressWarnings(as.numeric(subject_data[[roles$rate]][dosed]))
    delivery[is.finite(rate) & rate != 0] <- "infusion"
  }
  pieces[[length(pieces) + 1L]] <- delivery
  paste(sort(unique(do.call(paste, c(pieces, sep = ":")))), collapse = ";")
}

# Documented in `avatar-algorithm.Rmd`, Step 3. Change one, change the other.
.choose_transform <- function(values) {
  values <- values[is.finite(values)]
  if (!length(values)) {
    return(list(method = "identity", offset = 0, positive = FALSE))
  }
  negative_fraction <- mean(values < 0)
  positive <- values[values > 0]
  positive_like <- length(positive) > 0L && stats::median(values) > 0 &&
    negative_fraction <= 0.01
  if (!positive_like) {
    return(list(method = "identity", offset = 0, positive = FALSE))
  }
  offset <- if (length(positive)) min(positive) / 2 else 1e-06
  offset <- max(offset, .Machine$double.eps^0.5)
  list(method = "log_offset", offset = offset, positive = TRUE)
}

.transform_dv <- function(values, transform) {
  if (identical(transform$method, "log_offset")) {
    log(pmax(values, 0) + transform$offset)
  } else {
    values
  }
}

.inverse_dv <- function(values, transform) {
  if (identical(transform$method, "log_offset")) {
    pmax(exp(values) - transform$offset, 0)
  } else {
    values
  }
}

# Documented in `avatar-algorithm.Rmd`, Step 4. Change one, change the other.
.common_grid <- function(times, max_points = 15L) {
  times <- sort(unique(times[is.finite(times)]))
  if (length(times) <= max_points) return(times)
  grid <- as.numeric(stats::quantile(
    times, probs = seq(0, 1, length.out = max_points),
    names = FALSE, type = 8
  ))
  sort(unique(grid))
}

.trajectory_on_grid <- function(time, value, grid) {
  okay <- is.finite(time) & is.finite(value)
  time <- time[okay]
  value <- value[okay]
  if (!length(time)) return(rep(NA_real_, length(grid)))
  if (length(unique(time)) == 1L) {
    out <- rep(NA_real_, length(grid))
    closest <- which.min(abs(grid - time[1L]))
    out[closest] <- mean(value)
    return(out)
  }
  stats::approx(time, value, xout = grid, ties = mean, rule = 1)$y
}

# Documented in `avatar-algorithm.Rmd`, Step 4. Change one, change the other.
.build_profiles <- function(data, roles, pca_variance) {
  subjects <- .unique_in_order(data[[roles$id]])
  subject_rows <- lapply(subjects, function(subject) {
    data[[roles$id]] == subject
  })
  observed <- .observation_rows(data, roles, require_present = TRUE)
  endpoint <- .endpoint(data, roles)
  aligned_time <- .aligned_time(data, roles)
  endpoints <- sort(unique(endpoint[observed]))

  transforms <- stats::setNames(lapply(endpoints, function(ep) {
    .choose_transform(data[[roles$dv]][observed & endpoint == ep])
  }), endpoints)
  # Decided once here, from the source, so that generation and the run report
  # cannot disagree about which endpoints are discrete.
  value_types <- .endpoint_value_types(data, roles)
  grids <- stats::setNames(lapply(endpoints, function(ep) {
    .common_grid(aligned_time[observed & endpoint == ep])
  }), endpoints)

  features <- list()
  for (covariate in roles$covariates) {
    template <- data[[covariate]]
    values <- lapply(subject_rows, function(rows) .first_present(template[rows]))
    if (is.numeric(template) && !is.factor(template)) {
      features[[paste0("cov_", covariate)]] <- as.numeric(unlist(values))
    } else {
      character_values <- vapply(values, function(value) {
        if (!length(value) || is.na(value)) NA_character_ else as.character(value)
      }, character(1))
      lev <- if (is.factor(template)) levels(template) else
        sort(unique(character_values[!is.na(character_values)]))
      for (level in lev) {
        feature <- as.numeric(character_values == level)
        feature[is.na(character_values)] <- NA_real_
        features[[paste0("cov_", covariate, "_", make.names(level))]] <- feature
      }
    }
  }

  for (ep in endpoints) {
    grid <- grids[[ep]]
    endpoint_matrix <- t(vapply(subject_rows, function(rows) {
      selected <- rows & observed & endpoint == ep
      transformed <- .transform_dv(data[[roles$dv]][selected], transforms[[ep]])
      .trajectory_on_grid(aligned_time[selected], transformed, grid)
    }, numeric(length(grid))))
    if (!length(grid)) next
    colnames(endpoint_matrix) <- paste0(
      "dv_", make.names(ep), "_", seq_len(ncol(endpoint_matrix))
    )
    for (column in seq_len(ncol(endpoint_matrix))) {
      features[[colnames(endpoint_matrix)[column]]] <- endpoint_matrix[, column]
    }
  }

  profile <- if (length(features)) {
    as.matrix(as.data.frame(features, check.names = FALSE))
  } else {
    matrix(numeric(), nrow = length(subjects), ncol = 0L)
  }
  rownames(profile) <- as.character(subjects)

  usable <- if (ncol(profile)) colSums(is.finite(profile)) > 0 else logical()
  profile <- profile[, usable, drop = FALSE]
  medians <- if (ncol(profile)) apply(profile, 2L, stats::median, na.rm = TRUE) else numeric()
  imputed <- profile
  if (ncol(imputed)) {
    for (column in seq_len(ncol(imputed))) {
      missing <- !is.finite(imputed[, column])
      imputed[missing, column] <- medians[column]
    }
  }
  scales <- if (ncol(imputed)) apply(imputed, 2L, stats::sd) else numeric()
  variable <- is.finite(scales) & scales > sqrt(.Machine$double.eps)
  imputed <- imputed[, variable, drop = FALSE]
  scales <- scales[variable]
  centered <- if (ncol(imputed)) {
    scale(imputed, center = TRUE, scale = scales)
  } else {
    matrix(numeric(), nrow = length(subjects), ncol = 0L)
  }

  pca <- NULL
  if (nrow(centered) >= 2L && ncol(centered) >= 2L) {
    pca <- tryCatch(
      stats::prcomp(centered, center = FALSE, scale. = FALSE,
                    rank. = min(nrow(centered) - 1L, ncol(centered))),
      error = function(error) NULL
    )
  }
  if (!is.null(pca) && length(pca$sdev) && sum(pca$sdev^2) > 0) {
    explained <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
    components <- which(explained >= pca_variance)[1L]
    coordinates <- pca$x[, seq_len(components), drop = FALSE]
  } else if (ncol(centered)) {
    coordinates <- centered
  } else {
    coordinates <- matrix(0, nrow = length(subjects), ncol = 1L,
                          dimnames = list(as.character(subjects), "fallback"))
  }

  signatures <- vapply(subject_rows, function(rows) {
    .event_signature(data[rows, , drop = FALSE], roles)
  }, character(1))
  routes <- vapply(subject_rows, function(rows) {
    .route_key(data[rows, , drop = FALSE], roles)
  }, character(1))

  list(
    subjects = subjects,
    subject_rows = subject_rows,
    coordinates = coordinates,
    transforms = transforms,
    value_types = value_types,
    grids = grids,
    signatures = signatures,
    routes = routes,
    dropped_features = setdiff(colnames(profile), colnames(imputed)),
    pca = pca
  )
}

# Documented in `avatar-algorithm.Rmd`, Step 8. Change one, change the other.
.neighbor_distances <- function(coordinates, anchor, candidates) {
  if (!length(candidates)) return(numeric())
  difference <- sweep(coordinates[candidates, , drop = FALSE], 2L,
                      coordinates[anchor, ], "-")
  sqrt(rowSums(difference^2))
}

# Hold *every* donor at or below `max_weight`, not only the largest one.
# Water-filling: donors over the cap are pinned there, the mass they give up is
# redistributed proportionally among those still under it, and the pass repeats
# until none is over. Pinned donors are never topped up again, so the loop
# closes after at most one pass per donor.
#
# Capping only `which.max()` (the behavior through the 0.80 era) was adequate
# while the cap was loose, because a second donor rarely reached it. At a cap
# that actually bites, the runner-up routinely lands above it once the leader's
# excess is redistributed, so a single pass would leave the documented maximum
# violated by the very donor the redistribution created.
#
# No weight vector of length K can satisfy a cap below 1/K, so such a cap
# relaxes to exactly 1/K: uniform weights, the flattest blend available. This is
# what a small source falls back to -- with two donors and a 0.30 cap the answer
# is 0.5/0.5, not an error.
# Documented in `avatar-algorithm.Rmd`, Step 9, M5. Change one, change the other.
.cap_weights <- function(weights, max_weight) {
  n <- length(weights)
  if (!n) return(weights)
  cap <- max(max_weight, 1 / n)
  tolerance <- sqrt(.Machine$double.eps)
  pinned <- rep(FALSE, n)
  for (pass in seq_len(n)) {
    over <- !pinned & weights > cap + tolerance
    if (!any(over)) break
    pinned <- pinned | over
    weights[pinned] <- cap
    free <- !pinned
    if (!any(free)) break
    remaining <- 1 - sum(weights[pinned])
    share <- weights[free]
    weights[free] <- if (sum(share) > 0) {
      remaining * share / sum(share)
    } else {
      remaining / sum(free)
    }
  }
  weights / sum(weights)
}

# Documented in `avatar-algorithm.Rmd`, Step 9, M5. Change one, change the other.
.randomized_weights <- function(distances, max_weight = 0.50) {
  n <- length(distances)
  if (!n) return(numeric())
  if (n == 1L) return(1)
  positive <- distances[is.finite(distances) & distances > 0]
  epsilon <- if (length(positive)) max(1e-08, stats::median(positive) * 1e-06) else 1e-08
  randomized_rank <- sample.int(n)
  raw <- stats::rexp(n) / pmax(distances, epsilon) * 2^(-randomized_rank)
  if (!all(is.finite(raw)) || sum(raw) <= 0) raw <- rep(1, n)
  .cap_weights(raw / sum(raw), max_weight)
}
