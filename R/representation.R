.nearest_grid_cell <- function(x, grid) {
  if (!length(x)) return(integer())
  vapply(x, function(value) {
    if (!is.finite(value)) return(NA_integer_)
    which.min(abs(grid - value))[1L]
  }, integer(1))
}

.subject_clock <- function(data, roles, endpoint) {
  time <- as.numeric(data[[roles$time]])
  dose_rows <- .dose_rows(data, roles)
  dose_times <- sort(unique(time[dose_rows & is.finite(time)]))
  if (!length(dose_times)) dose_times <- min(time, na.rm = TRUE)
  occasion <- findInterval(time, dose_times)
  occasion[occasion < 1L] <- 1L
  if (!is.null(roles$occasion)) {
    declared_occasion <- suppressWarnings(as.numeric(data[[roles$occasion]]))
    okay <- is.finite(declared_occasion) & declared_occasion >= 1
    occasion[okay] <- as.integer(floor(declared_occasion[okay]))
  }
  qualifying <- dose_times[pmin(occasion, length(dose_times))]
  tad <- pmax(0, time - qualifying)
  if (!is.null(roles$tad)) {
    declared <- suppressWarnings(as.numeric(data[[roles$tad]]))
    okay <- is.finite(declared)
    tad[okay] <- pmax(0, declared[okay])
  }
  study <- time
  if (!is.null(roles$nominal_time)) {
    nominal <- suppressWarnings(as.numeric(data[[roles$nominal_time]]))
    okay <- is.finite(nominal)
    study[okay] <- nominal[okay]
  }
  list(
    primary = if (endpoint$alignment %in% c("dose_relative", "occasion")) {
      tad
    } else study,
    local = tad,
    study = study,
    occasion = occasion,
    dose_times = dose_times
  )
}

.bound_subject_contributions <- function(data, roles, endpoints, bounds,
                                         contribution_limits) {
  id <- data[[roles$id]]
  subjects <- .unique_in_order(id)
  out <- vector("list", length(subjects))
  for (i in seq_along(subjects)) {
    rows <- which(!is.na(id) & id == subjects[[i]])
    rows <- rows[seq_len(min(length(rows), contribution_limits$max_rows))]
    piece <- data[rows, , drop = FALSE]
    piece[[roles$time]] <- .clip(piece[[roles$time]], bounds$time)
    if (!is.null(roles$nominal_time)) {
      piece[[roles$nominal_time]] <- .clip(
        piece[[roles$nominal_time]], bounds$time
      )
    }
    if (!is.null(roles$tad)) {
      piece[[roles$tad]] <- .clip(
        piece[[roles$tad]], c(0, diff(bounds$time))
      )
    }
    if (!is.null(roles$amt)) {
      piece[[roles$amt]] <- .clip(piece[[roles$amt]], bounds$amt)
    }
    if (!is.null(roles$rate)) {
      piece[[roles$rate]] <- .clip(piece[[roles$rate]], bounds$rate)
    }
    endpoint_name <- .endpoint_name_for_rows(piece, roles, endpoints)
    for (name in names(endpoints)) {
      selected <- !is.na(endpoint_name) & endpoint_name == name
      piece[[roles$dv]][selected] <- .clip(
        piece[[roles$dv]][selected], bounds$endpoints[[name]]
      )
      if (!is.null(roles$limit)) {
        limit_bound <- bounds$limit[[name]] %||% bounds$endpoints[[name]]
        piece[[roles$limit]][selected] <- .clip(
          piece[[roles$limit]][selected], limit_bound
        )
      }
    }
    for (name in intersect(roles$covariates, names(bounds$covariates))) {
      piece[[name]] <- .clip(piece[[name]], bounds$covariates[[name]])
    }
    observed <- .observation_rows(piece, roles, require_present = TRUE)
    unknown <- observed & is.na(endpoint_name)
    if (any(unknown)) {
      stop("An observed DVID is not covered by the public endpoint declarations.",
           call. = FALSE)
    }
    observations <- stats::setNames(vector("list", length(endpoints)),
                                    names(endpoints))
    for (name in names(endpoints)) {
      candidate <- which(observed & endpoint_name == name)
      limit <- endpoints[[name]]$observation_limit
      if (length(candidate)) {
        occasion <- .subject_clock(piece, roles, endpoints[[name]])$occasion
        candidate <- candidate[occasion[candidate] <=
                                 contribution_limits$max_occasions]
        groups <- split(candidate, occasion[candidate])
        candidate <- unlist(lapply(groups, function(index) {
          index[seq_len(min(length(index), limit))]
        }), use.names = FALSE)
      }
      observations[[name]] <- candidate
    }
    dose <- which(.dose_rows(piece, roles))
    dose <- dose[seq_len(min(length(dose), contribution_limits$max_doses))]
    out[[i]] <- list(data = piece, observations = observations, doses = dose)
  }
  out
}

.base_event_features <- function(subjects, roles, bounds,
                                 contribution_limits) {
  names <- c(
    "has_event", "dose_count", "mean_interval", "mean_amount",
    "mean_rate", "has_infusion", "mean_duration", "occasion_count",
    "observation_count"
  )
  matrix <- matrix(0, nrow = length(subjects), ncol = length(names),
                   dimnames = list(NULL, names))
  time_span <- diff(bounds$time)
  for (i in seq_along(subjects)) {
    subject <- subjects[[i]]
    data <- subject$data
    time <- .clip(data[[roles$time]], bounds$time)
    dose <- subject$doses
    event <- .event_rows(data, roles)
    matrix[i, "has_event"] <- as.numeric(any(event))
    matrix[i, "dose_count"] <- length(dose) / contribution_limits$max_doses
    matrix[i, "occasion_count"] <- min(length(dose),
      contribution_limits$max_occasions) / contribution_limits$max_occasions
    matrix[i, "observation_count"] <-
      sum(vapply(subject$observations, length, integer(1))) /
      contribution_limits$max_rows
    if (length(dose) > 1L) {
      matrix[i, "mean_interval"] <-
        min(max(mean(diff(time[dose])) / time_span, 0), 1)
    }
    if (!is.null(roles$amt) && length(dose)) {
      amount <- .to_unit(data[[roles$amt]][dose], bounds$amt)
      amount <- amount[is.finite(amount)]
      if (length(amount)) matrix[i, "mean_amount"] <- mean(amount)
    }
    if (!is.null(roles$rate) && length(dose)) {
      rate <- suppressWarnings(as.numeric(data[[roles$rate]][dose]))
      rate_unit <- .to_unit(rate, bounds$rate)
      rate_unit <- rate_unit[is.finite(rate_unit)]
      if (length(rate_unit)) matrix[i, "mean_rate"] <- mean(rate_unit)
      matrix[i, "has_infusion"] <- as.numeric(any(is.finite(rate) & rate > 0))
      durations <- numeric()
      event_rows <- which(event)
      for (start in dose) {
        later <- event_rows[event_rows > start]
        stop_rows <- later[
          suppressWarnings(as.numeric(data[[roles$rate]][later])) < 0 |
            (!is.null(roles$amt) &
             suppressWarnings(as.numeric(data[[roles$amt]][later])) < 0)
        ]
        if (length(stop_rows)) {
          durations <- c(durations, time[stop_rows[1L]] - time[start])
        } else if (!is.null(roles$amt)) {
          amount_value <- suppressWarnings(as.numeric(data[[roles$amt]][start]))
          rate_value <- suppressWarnings(as.numeric(data[[roles$rate]][start]))
          if (is.finite(amount_value) && amount_value > 0 &&
              is.finite(rate_value) && rate_value > 0) {
            durations <- c(durations, amount_value / rate_value)
          }
        }
      }
      durations <- durations[is.finite(durations) & durations >= 0]
      if (length(durations)) {
        matrix[i, "mean_duration"] <-
          min(max(mean(durations) / time_span, 0), 1)
      }
    }
  }
  pmin(pmax(matrix, 0), 1)
}

.subject_property_spec <- function(roles, public_design) {
  property_names <- roles$subject_properties
  if (!length(property_names)) {
    return(list(names = character(), levels = list(), strata = list()))
  }
  levels <- stats::setNames(lapply(property_names, function(name) {
    column <- .schema_column(public_design$schema, name)
    values <- if (!is.null(column) && "factor" %in% column$class) {
      column$levels
    } else {
      public_design$category_levels[[name]]
    }
    unique(as.character(values))
  }), property_names)
  grid <- expand.grid(
    levels, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  strata <- lapply(seq_len(nrow(grid)), function(index) {
    stats::setNames(as.list(as.character(grid[index, , drop = TRUE])),
                    property_names)
  })
  list(names = property_names, levels = levels, strata = strata)
}

.subject_property_index <- function(data, spec) {
  if (!length(spec$names)) return(NA_integer_)
  indices <- vapply(spec$names, function(name) {
    value <- data[[name]][which(!is.na(data[[name]]))[1L]]
    match(as.character(value), spec$levels[[name]])
  }, integer(1))
  if (anyNA(indices)) {
    stop(
      "A subject-property value is outside the declared public category levels.",
      call. = FALSE
    )
  }
  multiplier <- cumprod(c(1L, utils::head(lengths(spec$levels), -1L)))
  as.integer(1L + sum((indices - 1L) * multiplier))
}

.event_features <- function(subjects, roles, bounds, contribution_limits,
                            public_design) {
  base <- .base_event_features(
    subjects, roles, bounds, contribution_limits
  )
  property_spec <- .subject_property_spec(roles, public_design)
  if (!length(property_spec$names)) {
    return(list(
      matrix = base,
      sensitivity = ncol(base),
      map = list(base_fields = colnames(base), properties = NULL)
    ))
  }

  strata <- lapply(seq_along(property_spec$strata), function(index) {
    prefix <- sprintf("property_stratum_%03d", index)
    list(
      values = property_spec$strata[[index]],
      count_field = paste0(prefix, "__count"),
      event_fields = stats::setNames(
        paste(prefix, colnames(base), sep = "__"), colnames(base)
      )
    )
  })
  columns <- unlist(lapply(strata, function(stratum) {
    c(stratum$count_field, unname(stratum$event_fields))
  }), use.names = FALSE)
  matrix <- matrix(
    0, nrow = length(subjects), ncol = length(columns),
    dimnames = list(NULL, columns)
  )
  for (index in seq_along(subjects)) {
    stratum_index <- .subject_property_index(
      subjects[[index]]$data, property_spec
    )
    stratum <- strata[[stratum_index]]
    matrix[index, stratum$count_field] <- 1
    matrix[index, unname(stratum$event_fields)] <- base[index, ]
  }
  list(
    matrix = matrix,
    # One subject contributes to exactly one stratum: one count indicator plus
    # the bounded base event vector. Dimension grows with declared strata, but
    # per-subject L1 sensitivity does not.
    sensitivity = 1 + ncol(base),
    map = list(
      base_fields = colnames(base),
      properties = list(names = property_spec$names, strata = strata)
    )
  )
}

.timing_features <- function(subjects, roles, endpoints,
                             contribution_limits) {
  map <- list()
  columns <- character()
  for (name in names(endpoints)) {
    cells <- seq_along(endpoints[[name]]$grid)
    grid_presence <- paste(name, "presence", cells, sep = "__")
    occasion_presence <- character()
    occasion_count <- character()
    if (endpoints[[name]]$alignment %in% c("dose_relative", "occasion")) {
      occasion_presence <- paste(
        name, "occasion_presence",
        seq_len(contribution_limits$max_occasions), sep = "__"
      )
      occasion_count <- paste(
        name, "occasion_count",
        seq_len(contribution_limits$max_occasions), sep = "__"
      )
    }
    map[[name]] <- list(
      grid_presence = grid_presence,
      occasion_presence = occasion_presence,
      occasion_count = occasion_count,
      observation_limit = endpoints[[name]]$observation_limit
    )
    columns <- c(
      columns, grid_presence, occasion_presence, occasion_count
    )
  }
  matrix <- matrix(0, nrow = length(subjects), ncol = length(columns),
                   dimnames = list(NULL, columns))
  for (i in seq_along(subjects)) {
    subject <- subjects[[i]]
    for (name in names(endpoints)) {
      rows <- subject$observations[[name]]
      if (!length(rows)) next
      endpoint <- endpoints[[name]]
      clock <- .subject_clock(subject$data, roles, endpoint)
      cells <- unique(.nearest_grid_cell(clock$primary[rows], endpoint$grid))
      cells <- cells[!is.na(cells)]
      matrix[i, map[[name]]$grid_presence[cells]] <- 1
      if (length(map[[name]]$occasion_count)) {
        counts <- table(factor(
          clock$occasion[rows],
          levels = seq_len(contribution_limits$max_occasions)
        ))
        matrix[i, map[[name]]$occasion_presence] <- as.numeric(counts > 0)
        matrix[i, map[[name]]$occasion_count] <- pmin(
          as.numeric(counts) / endpoints[[name]]$observation_limit, 1
        )
      }
    }
  }
  list(matrix = matrix, map = map)
}

.trajectory_features <- function(subjects, roles, endpoints) {
  map <- list()
  columns <- character()
  for (name in names(endpoints)) {
    cells <- seq_along(endpoints[[name]]$grid)
    presence <- paste(name, "presence", cells, sep = "__")
    value <- paste(name, "value", cells, sep = "__")
    local_presence <- local_value <- character()
    if (endpoints[[name]]$alignment == "hybrid") {
      local_presence <- paste(name, "local_presence", cells, sep = "__")
      local_value <- paste(name, "local_value", cells, sep = "__")
    }
    map[[name]] <- list(
      presence = presence, value = value,
      local_presence = local_presence, local_value = local_value
    )
    columns <- c(columns, presence, value, local_presence, local_value)
  }
  matrix <- matrix(0, nrow = length(subjects), ncol = length(columns),
                   dimnames = list(NULL, columns))
  for (i in seq_along(subjects)) {
    subject <- subjects[[i]]
    for (name in names(endpoints)) {
      rows <- subject$observations[[name]]
      if (!length(rows)) next
      endpoint <- endpoints[[name]]
      clock <- .subject_clock(subject$data, roles, endpoint)
      cell <- .nearest_grid_cell(clock$primary[rows], endpoint$grid)
      working <- .transform_endpoint(subject$data[[roles$dv]][rows], endpoint)
      working_bounds <- .endpoint_working_bounds(endpoint)
      value <- .to_unit(working, working_bounds)
      for (j in seq_along(endpoint$grid)) {
        selected <- which(cell == j & is.finite(value))
        if (length(selected)) {
          matrix[i, map[[name]]$presence[j]] <- 1
          matrix[i, map[[name]]$value[j]] <- mean(value[selected])
        }
      }
      if (endpoint$alignment == "hybrid") {
        local_cell <- .nearest_grid_cell(clock$local[rows],
                                         endpoint$local_grid)
        overall <- mean(value[is.finite(value)])
        for (j in seq_along(endpoint$local_grid)) {
          selected <- which(local_cell == j & is.finite(value))
          if (length(selected)) {
            matrix[i, map[[name]]$local_presence[j]] <- 1
            centered <- mean(value[selected]) - overall + 0.5
            matrix[i, map[[name]]$local_value[j]] <-
              min(max(centered, 0), 1)
          }
        }
      }
    }
  }
  list(matrix = matrix, map = map)
}

.covariate_features <- function(subjects, roles, bounds, public_design) {
  if (!length(roles$covariates)) {
    return(list(matrix = matrix(numeric(), nrow = length(subjects), ncol = 0),
                map = list()))
  }
  schema <- public_design$schema
  map <- list()
  columns <- character()
  for (name in roles$covariates) {
    column <- .schema_column(schema, name)
    forced_categorical <- length(public_design$category_levels[[name]]) > 0L
    numeric <- !forced_categorical && !is.null(column) &&
      ("numeric" %in% column$class || "integer" %in% column$class) &&
      !("factor" %in% column$class)
    if (numeric) {
      if (is.null(bounds$covariates[[name]])) {
        stop("A public bound is required for numeric covariate `", name, "`.",
             call. = FALSE)
      }
      fields <- paste(name, c("presence", "value", "square"), sep = "__")
      map[[name]] <- list(type = "numeric", fields = fields)
    } else {
      levels <- column$levels %||% public_design$category_levels[[name]]
      if (is.null(levels) || !length(levels)) {
        stop("Public category levels are required for covariate `", name, "`.",
             call. = FALSE)
      }
      fields <- paste(name, "category", seq_along(levels), sep = "__")
      map[[name]] <- list(type = "categorical", fields = fields,
                          levels = as.character(levels))
    }
    columns <- c(columns, fields)
  }
  matrix <- matrix(0, nrow = length(subjects), ncol = length(columns),
                   dimnames = list(NULL, columns))
  for (i in seq_along(subjects)) {
    data <- subjects[[i]]$data
    for (name in roles$covariates) {
      value <- data[[name]][which(!is.na(data[[name]]))[1L]]
      if (!length(value) || is.na(value)) next
      spec <- map[[name]]
      if (spec$type == "numeric") {
        unit <- .to_unit(value, bounds$covariates[[name]])
        matrix[i, spec$fields] <- c(1, unit, unit^2)
      } else {
        index <- match(as.character(value), spec$levels)
        if (!is.na(index)) matrix[i, spec$fields[index]] <- 1
      }
    }
  }
  list(matrix = matrix, map = map)
}

.censoring_features <- function(subjects, roles, endpoints, bounds) {
  if (is.null(roles$cens)) {
    return(list(matrix = matrix(numeric(), nrow = length(subjects), ncol = 0),
                map = list()))
  }
  map <- list()
  columns <- character()
  for (name in names(endpoints)) {
    fields <- paste(name,
                    c("observed", "left", "right", "interval",
                      "boundary", "other_boundary"), sep = "__")
    map[[name]] <- fields
    columns <- c(columns, fields)
  }
  matrix <- matrix(0, nrow = length(subjects), ncol = length(columns),
                   dimnames = list(NULL, columns))
  for (i in seq_along(subjects)) {
    subject <- subjects[[i]]
    for (name in names(endpoints)) {
      rows <- subject$observations[[name]]
      if (!length(rows)) next
      data <- subject$data
      cens <- suppressWarnings(as.numeric(as.character(data[[roles$cens]][rows])))
      okay <- cens %in% c(-1, 0, 1)
      if (!any(okay)) next
      cens <- cens[okay]
      selected_rows <- rows[okay]
      other <- if (is.null(roles$limit)) rep(NA_real_, length(cens)) else
        suppressWarnings(as.numeric(data[[roles$limit]][selected_rows]))
      dv <- suppressWarnings(as.numeric(data[[roles$dv]][selected_rows]))
      interval <- is.finite(other) & cens != 0
      left <- cens == 1 & !interval
      right <- cens == -1 & !interval
      fields <- map[[name]]
      matrix[i, fields[1L]] <- 1
      matrix[i, fields[2L]] <- mean(left)
      matrix[i, fields[3L]] <- mean(right)
      matrix[i, fields[4L]] <- mean(interval)
      boundary_bounds <- bounds$limit[[name]] %||% bounds$endpoints[[name]]
      if (any(cens != 0 & is.finite(dv))) {
        matrix[i, fields[5L]] <- mean(.to_unit(
          dv[cens != 0 & is.finite(dv)], boundary_bounds
        ))
      }
      if (any(interval & is.finite(other))) {
        matrix[i, fields[6L]] <- mean(.to_unit(
          other[interval & is.finite(other)], boundary_bounds
        ))
      }
    }
  }
  list(matrix = matrix, map = map)
}

.release_matrix_sum <- function(accountant, query, matrix, epsilon,
                                sensitivity = ncol(matrix)) {
  if (!ncol(matrix)) return(numeric())
  matrix <- pmin(pmax(matrix, 0), 1)
  .private_release(
    accountant, query, colSums(matrix), sensitivity = sensitivity,
    epsilon = epsilon
  ) |>
    stats::setNames(colnames(matrix))
}

.decode_event <- function(released, count, bounds, contribution_limits) {
  mean <- pmin(pmax(released / max(count, 1), 0), 1)
  list(
    event_probability = mean[["has_event"]],
    n_doses = max(0L, as.integer(round(
      mean[["dose_count"]] * contribution_limits$max_doses
    ))),
    dose_interval = mean[["mean_interval"]] * diff(bounds$time),
    dose_amount = if (is.null(bounds$amt)) 0 else
      .from_unit(mean[["mean_amount"]], bounds$amt),
    dose_rate = if (is.null(bounds$rate)) 0 else
      .from_unit(mean[["mean_rate"]], bounds$rate),
    infusion_probability = mean[["has_infusion"]],
    infusion_duration = mean[["mean_duration"]] * diff(bounds$time),
    occasions = max(1L, as.integer(round(
      mean[["occasion_count"]] * contribution_limits$max_occasions
    ))),
    observation_count = mean[["observation_count"]] *
      contribution_limits$max_rows
  )
}

.decode_property_events <- function(released, event_map, count, bounds,
                                    contribution_limits) {
  properties <- event_map$properties
  if (is.null(properties)) {
    return(list(
      event = .decode_event(
        released[event_map$base_fields], count, bounds,
        contribution_limits
      ),
      subject_properties = NULL
    ))
  }

  global <- stats::setNames(rep(0, length(event_map$base_fields)),
                            event_map$base_fields)
  strata <- lapply(properties$strata, function(stratum) {
    released_count <- max(as.numeric(released[[stratum$count_field]]), 0)
    event_values <- released[unname(stratum$event_fields)]
    names(event_values) <- names(stratum$event_fields)
    global <<- global + event_values
    list(
      values = stratum$values,
      released_count = released_count,
      probability = released_count,
      event = .decode_event(
        event_values, released_count, bounds, contribution_limits
      )
    )
  })
  probability <- vapply(strata, `[[`, numeric(1), "probability")
  if (!sum(probability) > 0) probability[] <- 1
  probability <- probability / sum(probability)
  for (index in seq_along(strata)) {
    strata[[index]]$probability <- unname(probability[index])
  }
  list(
    event = .decode_event(global, count, bounds, contribution_limits),
    subject_properties = list(
      names = properties$names,
      strata = strata
    )
  )
}

.decode_timing <- function(released, count, timing_map) {
  lapply(timing_map, function(spec) {
    grid_probability <- pmin(pmax(
      released[spec$grid_presence] / max(count, 1), 0
    ), 1)
    released_occasion_presence <- if (length(spec$occasion_presence)) {
      pmax(released[spec$occasion_presence], 0)
    } else numeric()
    occasion_presence_probability <- if (
      length(released_occasion_presence)
    ) {
      pmin(released_occasion_presence / max(count, 1), 1)
    } else numeric()
    occasion_observation_count <- if (length(spec$occasion_count)) {
      pmin(pmax(
        released[spec$occasion_count] /
          pmax(released_occasion_presence, 0.25) *
          spec$observation_limit,
        0
      ), spec$observation_limit)
    } else numeric()
    list(
      grid_probability = unname(grid_probability),
      occasion_presence_probability = unname(occasion_presence_probability),
      occasion_observation_count = unname(occasion_observation_count)
    )
  })
}

.fill_unoccupied_curve <- function(value_unit, presence, grid,
                                   threshold = 0.25) {
  occupied <- which(presence > threshold & is.finite(value_unit))
  if (!length(occupied)) return(rep(0.5, length(value_unit)))
  if (length(occupied) == 1L) return(rep(value_unit[occupied], length(value_unit)))
  # Empty cells mean "no released support here", not a measurement at the
  # midpoint of the endpoint domain. Interpolation is source-free
  # post-processing of the released occupied cells and avoids artificial
  # troughs on wide log-transformed domains.
  stats::approx(
    x = grid[occupied], y = value_unit[occupied], xout = grid,
    rule = 2, ties = mean
  )$y
}

.decode_trajectories <- function(released, trajectory_map, endpoints) {
  out <- list()
  for (name in names(endpoints)) {
    spec <- trajectory_map[[name]]
    presence <- pmax(released[spec$presence], 0)
    value_sum <- pmax(released[spec$value], 0)
    value_unit <- ifelse(presence > 0.25,
                         value_sum / pmax(presence, 1e-8), NA_real_)
    value_unit <- pmin(pmax(value_unit, 0), 1)
    endpoint <- endpoints[[name]]
    value_unit <- .fill_unoccupied_curve(
      value_unit, presence, endpoint$grid
    )
    working <- .from_unit(value_unit, .endpoint_working_bounds(endpoint))
    item <- list(
      grid = endpoint$grid,
      mean_working = working,
      released_presence = presence
    )
    if (endpoint$alignment == "hybrid") {
      local_presence <- pmax(released[spec$local_presence], 0)
      local_sum <- pmax(released[spec$local_value], 0)
      local_unit <- ifelse(local_presence > 0.25,
                           local_sum / pmax(local_presence, 1e-8), 0.5)
      item$local_grid <- endpoint$local_grid
      item$local_excursion_unit <- pmin(pmax(local_unit, 0), 1) - 0.5
    }
    out[[name]] <- item
  }
  out
}

.decode_covariates <- function(released, count, covariate_map, bounds) {
  out <- list()
  for (name in names(covariate_map)) {
    spec <- covariate_map[[name]]
    if (spec$type == "numeric") {
      fields <- spec$fields
      presence <- max(released[[fields[1L]]], 0)
      center <- if (presence > 0.25) released[[fields[2L]]] / presence else 0.5
      second <- if (presence > 0.25) released[[fields[3L]]] / presence else
        center^2 + 0.05^2
      center <- min(max(center, 0), 1)
      variance <- min(max(second - center^2, 0.01^2), 0.25)
      out[[name]] <- list(
        type = "numeric",
        center = .from_unit(center, bounds$covariates[[name]]),
        scale = sqrt(variance) * diff(bounds$covariates[[name]]),
        bounds = bounds$covariates[[name]],
        presence_probability = min(max(presence / max(count, 1), 0), 1)
      )
    } else {
      probability <- pmax(released[spec$fields], 0)
      if (!sum(probability) > 0) probability[] <- 1
      probability <- probability / sum(probability)
      out[[name]] <- list(type = "categorical", levels = spec$levels,
                          probability = as.numeric(probability))
    }
  }
  out
}

.decode_censoring <- function(released, count, censor_map, endpoints, bounds) {
  out <- list()
  for (name in names(censor_map)) {
    fields <- censor_map[[name]]
    observed <- max(released[[fields[1L]]], 0)
    denominator <- max(observed, 0.25)
    rates <- pmin(pmax(released[fields[2:4]] / denominator, 0), 1)
    total <- sum(rates)
    if (total > 0.95) rates <- rates * 0.95 / total
    boundary_bounds <- bounds$limit[[name]] %||% bounds$endpoints[[name]]
    boundary <- .from_unit(
      min(max(released[[fields[5L]]] / denominator, 0), 1),
      boundary_bounds
    )
    other <- .from_unit(
      min(max(released[[fields[6L]]] / denominator, 0), 1),
      boundary_bounds
    )
    out[[name]] <- list(
      left = unname(rates[1L]), right = unname(rates[2L]),
      interval = unname(rates[3L]), boundary = boundary,
      other_boundary = other
    )
  }
  out
}

.fit_population_summaries <- function(data, roles, endpoints, bounds,
                                      contribution_limits, budget_allocation,
                                      epsilon, accountant, public_design) {
  subjects <- .bound_subject_contributions(
    data, roles, endpoints, bounds, contribution_limits
  )
  fractions <- unlist(budget_allocation, use.names = TRUE)
  group_epsilon <- epsilon * fractions
  overshoot <- sum(group_epsilon) - epsilon
  if (overshoot > 0) {
    # IEEE arithmetic can round a mathematically valid allocation slightly
    # upward. Charge that rounding to the largest group so the mechanisms,
    # not merely the report, remain within the exact requested budget.
    index <- which.max(group_epsilon)
    group_epsilon[index] <- group_epsilon[index] - overshoot
  }
  active <- c(
    subject_count = TRUE, event = TRUE, timing = TRUE,
    covariates = length(roles$covariates) > 0L,
    endpoints = TRUE, censoring = !is.null(roles$cens)
  )
  missing_budget <- names(active)[active & fractions[names(active)] <= 0]
  if (length(missing_budget)) {
    stop("Active private summary groups need positive budget: ",
         paste(missing_budget, collapse = ", "), ".", call. = FALSE)
  }
  count_release <- .private_release(
    accountant, "subject_count", length(subjects), sensitivity = 1,
    epsilon = group_epsilon[["subject_count"]]
  )
  private_count <- max(as.numeric(count_release), 1)

  event_matrix <- .event_features(
    subjects, roles, bounds, contribution_limits, public_design
  )
  event_release <- .release_matrix_sum(
    accountant, "event_and_regimen", event_matrix$matrix,
    group_epsilon[["event"]], sensitivity = event_matrix$sensitivity
  )
  timing_features <- .timing_features(
    subjects, roles, endpoints, contribution_limits
  )
  timing_release <- .release_matrix_sum(
    accountant, "endpoint_timing", timing_features$matrix,
    group_epsilon[["timing"]]
  )
  trajectory_features <- .trajectory_features(subjects, roles, endpoints)
  trajectory_release <- .release_matrix_sum(
    accountant, "endpoint_trajectories", trajectory_features$matrix,
    group_epsilon[["endpoints"]]
  )
  covariate_features <- .covariate_features(
    subjects, roles, bounds, public_design
  )
  covariate_release <- if (ncol(covariate_features$matrix)) {
    .release_matrix_sum(
      accountant, "baseline_covariates", covariate_features$matrix,
      group_epsilon[["covariates"]]
    )
  } else numeric()
  censor_features <- .censoring_features(subjects, roles, endpoints, bounds)
  censor_release <- if (ncol(censor_features$matrix)) {
    .release_matrix_sum(
      accountant, "censoring", censor_features$matrix,
      group_epsilon[["censoring"]]
    )
  } else numeric()

  decoded_event <- .decode_property_events(
    event_release, event_matrix$map, private_count, bounds,
    contribution_limits
  )
  list(
    private_subject_count = private_count,
    event = decoded_event$event,
    subject_properties = decoded_event$subject_properties,
    timing = .decode_timing(timing_release, private_count,
                            timing_features$map),
    trajectories = .decode_trajectories(
      trajectory_release, trajectory_features$map, endpoints
    ),
    covariates = .decode_covariates(
      covariate_release, private_count, covariate_features$map, bounds
    ),
    censoring = if (length(censor_release)) .decode_censoring(
      censor_release, private_count, censor_features$map, endpoints, bounds
    ) else list()
  )
}
