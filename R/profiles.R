.event_signature <- function(subject_data, roles) {
  event <- !.is_zero(subject_data[[roles$evid]])
  observed <- .observation_rows(subject_data, roles)
  endpoint <- .endpoint(subject_data, roles)

  event_token <- if (any(event)) {
    pieces <- list(as.character(subject_data[[roles$evid]][event]))
    for (role in c("cmt", "dvid")) {
      if (!is.null(roles[[role]])) {
        pieces[[length(pieces) + 1L]] <-
          as.character(subject_data[[roles[[role]]]][event])
      }
    }
    for (role in c("amt", "rate")) {
      if (!is.null(roles[[role]])) {
        value <- subject_data[[roles[[role]]]][event]
        pieces[[length(pieces) + 1L]] <- ifelse(
          is.na(value), "NA", ifelse(value > 0, "+", ifelse(value < 0, "-", "0"))
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

  list(
    subjects = subjects,
    subject_rows = subject_rows,
    coordinates = coordinates,
    transforms = transforms,
    grids = grids,
    signatures = signatures,
    dropped_features = setdiff(colnames(profile), colnames(imputed)),
    pca = pca
  )
}

.neighbor_distances <- function(coordinates, anchor, candidates) {
  if (!length(candidates)) return(numeric())
  difference <- sweep(coordinates[candidates, , drop = FALSE], 2L,
                      coordinates[anchor, ], "-")
  sqrt(rowSums(difference^2))
}

.randomized_weights <- function(distances, max_weight = 0.80) {
  n <- length(distances)
  if (!n) return(numeric())
  if (n == 1L) return(1)
  positive <- distances[is.finite(distances) & distances > 0]
  epsilon <- if (length(positive)) max(1e-08, stats::median(positive) * 1e-06) else 1e-08
  randomized_rank <- sample.int(n)
  raw <- stats::rexp(n) / pmax(distances, epsilon) * 2^(-randomized_rank)
  if (!all(is.finite(raw)) || sum(raw) <= 0) raw <- rep(1, n)
  weights <- raw / sum(raw)
  if (max(weights) > max_weight) {
    dominant <- which.max(weights)
    remainder <- weights[-dominant]
    if (!is.finite(sum(remainder)) || sum(remainder) <= 0) {
      remainder <- rep(1, length(remainder))
    }
    weights[dominant] <- max_weight
    weights[-dominant] <- (1 - max_weight) * remainder / sum(remainder)
  }
  weights / sum(weights)
}
