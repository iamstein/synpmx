.resolved_regimen <- function(model, event = NULL) {
  design <- model$public$design
  property_conditioned <- !is.null(event)
  event <- event %||% model$population$event
  bounds <- model$public$bounds
  n_doses <- design$n_doses %||% event$n_doses
  if (!length(n_doses) || !is.finite(n_doses)) n_doses <- 0L
  n_doses <- max(0L, as.integer(round(n_doses)))
  interval <- design$dose_interval %||% event$dose_interval
  if (!is.finite(interval) || interval <= 0) {
    interval <- max(diff(bounds$time) / max(n_doses, 1L),
                    sqrt(.Machine$double.eps))
  }
  dose_times <- design$dose_times
  if (is.null(dose_times) && n_doses > 0L) {
    dose_times <- bounds$time[1L] + seq.int(0L, n_doses - 1L) * interval
  }
  if (!is.null(dose_times)) {
    dose_times <- dose_times[dose_times >= bounds$time[1L] &
                               dose_times <= bounds$time[2L]]
    n_doses <- length(dose_times)
  }
  amount <- design$dose_amount %||% event$dose_amount
  if (!is.finite(amount)) amount <- 0
  if (!is.null(bounds$amt)) amount <- .clip(amount, bounds$amt)
  rate <- design$dose_rate %||% event$dose_rate
  if (!is.finite(rate)) rate <- 0
  if (!is.null(bounds$rate)) rate <- .clip(rate, bounds$rate)
  duration <- design$infusion_duration %||% event$infusion_duration
  if (!is.finite(duration) || duration < 0) duration <- 0
  infusion <- n_doses > 0L && rate > 0 && duration > 0 &&
    (!is.null(design$dose_rate) || event$infusion_probability >= 0.5)
  list(
    n_doses = n_doses, dose_times = dose_times %||% numeric(),
    interval = interval, amount = amount, rate = rate,
    duration = duration, infusion = infusion,
    observation_count = event$observation_count,
    amount_jitter_sd = if (
      property_conditioned || !is.null(model$public$roles$assigned_dose)
    ) 0 else 0.08
  )
}

.new_public_row <- function(model) {
  schema <- model$public$design$schema
  defaults <- model$public$design$defaults %||% list()
  out <- stats::setNames(vector("list", length(.schema_names(schema))),
                         .schema_names(schema))
  for (name in names(out)) {
    out[[name]] <- if (!is.null(defaults[[name]])) defaults[[name]][1L] else NA
  }
  out
}

.set_row_role <- function(row, roles, role, value) {
  column <- roles[[role]]
  if (!is.null(column)) row[[column]] <- value
  row
}

.event_row <- function(model, time, amount, rate, cmt, evid, interval,
                       endpoint, occasion, tie_order, stop = FALSE) {
  roles <- model$public$roles
  row <- .new_public_row(model)
  row <- .set_row_role(row, roles, "time", time)
  row <- .set_row_role(row, roles, "nominal_time", time)
  row <- .set_row_role(row, roles, "tad", 0)
  row <- .set_row_role(row, roles, "occasion", occasion)
  event_dv <- model$public$design$defaults[[roles$dv]] %||% 0
  row <- .set_row_role(row, roles, "dv", event_dv)
  row <- .set_row_role(row, roles, "amt", if (stop) -abs(amount) else amount)
  row <- .set_row_role(row, roles, "rate", if (stop) -abs(rate) else rate)
  row <- .set_row_role(row, roles, "evid", evid)
  row <- .set_row_role(row, roles, "cmt", cmt)
  row <- .set_row_role(row, roles, "dvid", endpoint$dvid)
  row <- .set_row_role(row, roles, "mdv", 1)
  row <- .set_row_role(row, roles, "cens", 0)
  row <- .set_row_role(row, roles, "limit", NA_real_)
  row <- .set_row_role(row, roles, "addl", 0)
  row <- .set_row_role(row, roles, "ii", interval)
  row$.endpoint_name <- NA_character_
  row$.nominal_internal <- time
  row$.tad_internal <- 0
  row$.occasion_internal <- occasion
  row$.tie_order <- tie_order
  row
}

.observation_row <- function(model, endpoint_name, nominal, actual, tad,
                             occasion, tie_order) {
  roles <- model$public$roles
  endpoint <- model$public$endpoints[[endpoint_name]]
  row <- .new_public_row(model)
  row <- .set_row_role(row, roles, "time", actual)
  row <- .set_row_role(row, roles, "nominal_time", nominal)
  row <- .set_row_role(row, roles, "tad", tad)
  row <- .set_row_role(row, roles, "occasion", occasion)
  row <- .set_row_role(row, roles, "dv", NA_real_)
  row <- .set_row_role(row, roles, "amt", 0)
  row <- .set_row_role(row, roles, "rate", 0)
  row <- .set_row_role(row, roles, "evid", 0)
  cmt <- endpoint$cmt %||%
    model$public$design$endpoint_cmt[[endpoint_name]] %||% NA
  row <- .set_row_role(row, roles, "cmt", cmt)
  row <- .set_row_role(row, roles, "dvid", endpoint$dvid)
  row <- .set_row_role(row, roles, "mdv", 0)
  row <- .set_row_role(row, roles, "cens", 0)
  row <- .set_row_role(row, roles, "limit", NA_real_)
  row <- .set_row_role(row, roles, "addl", 0)
  row <- .set_row_role(row, roles, "ii", 0)
  row$.endpoint_name <- endpoint_name
  row$.nominal_internal <- nominal
  row$.tad_internal <- tad
  row$.occasion_internal <- occasion
  row$.tie_order <- tie_order
  row
}

.ensure_grid_presence <- function(probability, minimum = 3L) {
  selected <- which(stats::runif(length(probability)) < probability)
  target <- min(length(probability), minimum)
  if (length(selected) < target) {
    ranked <- order(probability, decreasing = TRUE)
    selected <- sort(unique(c(selected, ranked[seq_len(target)])))
  }
  selected
}

.timing_grid_probability <- function(timing, cells) {
  probability <- if (is.list(timing)) timing$grid_probability else timing
  if (!is.numeric(probability) || length(probability) != cells ||
      any(!is.finite(probability))) {
    return(rep(1, cells))
  }
  pmin(pmax(probability, 0), 1)
}

.weighted_cells <- function(cells, size, probability) {
  if (size <= 0L || !length(cells)) return(integer())
  if (size >= length(cells)) return(cells)
  probability <- pmax(probability, sqrt(.Machine$double.eps))
  cells[sample.int(
    length(cells), size = size, replace = FALSE, prob = probability
  )]
}

.select_timing_cells <- function(probability, count) {
  count <- min(length(probability), max(0L, as.integer(count)))
  if (!count) return(integer())
  probability <- pmin(pmax(probability, 0), 1)
  selected <- which(stats::runif(length(probability)) < probability)

  if (length(selected) > count) {
    certain <- selected[
      probability[selected] >= 1 - sqrt(.Machine$double.eps)
    ]
    if (length(certain) >= count) {
      selected <- .weighted_cells(
        certain, count, rep(1, length(certain))
      )
    } else {
      candidates <- setdiff(selected, certain)
      selected <- c(
        certain,
        .weighted_cells(
          candidates, count - length(certain), probability[candidates]
        )
      )
    }
  }
  if (length(selected) < count) {
    candidates <- setdiff(seq_along(probability), selected)
    selected <- c(
      selected,
      .weighted_cells(
        candidates, count - length(selected), probability[candidates]
      )
    )
  }
  sort(selected)
}

.inferred_occasion_counts <- function(timing, occasions) {
  count <- if (is.list(timing)) timing$occasion_observation_count else NULL
  if (!is.numeric(count) || !length(count) || any(!is.finite(count))) {
    return(NULL)
  }
  c(count, rep(0, max(0L, occasions - length(count))))[seq_len(occasions)]
}

.inferred_occasion_presence <- function(timing, occasions) {
  probability <- if (is.list(timing)) {
    timing$occasion_presence_probability
  } else NULL
  if (!is.numeric(probability) || !length(probability) ||
      any(!is.finite(probability))) {
    count <- .inferred_occasion_counts(timing, occasions)
    if (is.null(count)) return(NULL)
    probability <- pmin(count, 1)
  }
  c(
    probability,
    rep(0, max(0L, occasions - length(probability)))
  )[seq_len(occasions)] |>
    pmin(1) |>
    pmax(0)
}

.generated_local_grid <- function(endpoint, selected, regimen,
                                  occasion = NULL, timing = NULL) {
  grid <- endpoint$grid[selected]
  if (!isTRUE(endpoint$grid_automatic) || regimen$n_doses <= 1L ||
      !is.finite(regimen$interval) || regimen$interval <= 0) {
    return(grid)
  }
  horizon <- endpoint$grid_horizon %||% max(endpoint$grid)
  if (!is.finite(horizon) || horizon <= 0) return(grid)
  # Nonterminal local profiles must remain before the following dose. The last
  # occasion has no such boundary, so retain the released occupied horizon;
  # this supports a terminal washout profile longer than one dose interval
  # without reading or copying a source visit schedule.
  target <- regimen$interval * (1 - 1e-6)
  if (!is.null(occasion) && occasion == regimen$n_doses) {
    probability <- .timing_grid_probability(timing, length(endpoint$grid))
    occupied <- which(probability > sqrt(.Machine$double.eps))
    if (length(occupied)) target <- max(endpoint$grid[occupied])
  }
  grid * target / horizon
}

.observation_group_weight <- function(rows, data, model) {
  endpoint_name <- data$.endpoint_name[rows[1L]]
  endpoint <- model$public$endpoints[[endpoint_name]]
  if (!endpoint$alignment %in% c("dose_relative", "occasion")) {
    timing <- model$population$timing[[endpoint_name]]
    probability <- .timing_grid_probability(timing, length(endpoint$grid))
    return(max(sum(probability), 1))
  }
  occasion <- as.integer(data$.occasion_internal[rows[1L]])
  public_schedule <- model$public$design$endpoint_occasion_grids[[endpoint_name]]
  if (!is.null(public_schedule)) {
    return(max(length(public_schedule[[as.character(occasion)]]), 1L))
  }
  inferred <- .inferred_occasion_counts(
    model$population$timing[[endpoint_name]], occasion
  )
  if (is.null(inferred)) return(length(rows))
  max(inferred[occasion], sqrt(.Machine$double.eps))
}

.trim_observations_to_private_count <- function(data, model, target = NULL) {
  observation <- which(!is.na(data$.endpoint_name))
  if (!length(observation)) return(data)
  target <- target %||% model$population$event$observation_count
  if (!is.numeric(target) || length(target) != 1L || !is.finite(target)) {
    return(data)
  }

  endpoint <- data$.endpoint_name[observation]
  local <- vapply(endpoint, function(name) {
    model$public$endpoints[[name]]$alignment %in%
      c("dose_relative", "occasion")
  }, logical(1))
  group <- endpoint
  group[local] <- paste(
    endpoint[local], data$.occasion_internal[observation][local], sep = "__"
  )
  groups <- split(observation, group)
  sizes <- vapply(groups, length, integer(1))

  # Preserve at least one observation for every endpoint/occasion represented
  # by the generated skeleton, then honor the privately released total as
  # closely as structural validity permits.
  target <- min(length(observation), max(length(groups), as.integer(round(target))))
  allocation <- rep(1L, length(groups))
  names(allocation) <- names(groups)
  remaining <- target - sum(allocation)
  weights <- vapply(groups, .observation_group_weight, numeric(1),
                    data = data, model = model)
  while (remaining > 0L && any(allocation < sizes)) {
    eligible <- which(allocation < sizes)
    deficit <- weights[eligible] - allocation[eligible]
    if (any(deficit > 0)) {
      chosen <- eligible[[which.max(deficit)]]
    } else {
      capacity_score <- (sizes[eligible] - allocation[eligible]) *
        pmax(weights[eligible], sqrt(.Machine$double.eps))
      chosen <- eligible[[which.max(capacity_score)]]
    }
    allocation[chosen] <- allocation[chosen] + 1L
    remaining <- remaining - 1L
  }

  keep <- which(is.na(data$.endpoint_name))
  for (name in names(groups)) {
    rows <- groups[[name]]
    rows <- rows[order(data$.nominal_internal[rows])]
    count <- allocation[[name]]
    endpoint <- model$public$endpoints[[data$.endpoint_name[rows[1L]]]]
    if (endpoint$alignment %in% c("dose_relative", "occasion")) {
      # A sparse local occasion uses its first post-dose cell. Denser profiles
      # are selected from the privacy-accounted grid-presence probabilities;
      # filling rows from the beginning would systematically delete late PK
      # samples whenever the released per-subject count is shorter than the
      # public grid.
      selected <- if (count == 1L && length(rows) > 1L) {
        postdose <- which(data$.tad_internal[rows] >
                            sqrt(.Machine$double.eps))
        if (length(postdose)) postdose[1L] else 1L
      } else {
        timing <- model$population$timing[[data$.endpoint_name[rows[1L]]]]
        grid_probability <- .timing_grid_probability(
          timing, length(endpoint$grid)
        )
        cell <- .nearest_grid_cell(
          data$.tad_internal[rows], endpoint$grid
        )
        .select_timing_cells(grid_probability[cell], count)
      }
    } else {
      anchors <- unique(as.integer(round(seq(
        1L, length(rows), length.out = min(count, 3L)
      ))))
      selected <- anchors
      if (length(selected) < count) {
        available <- setdiff(seq_along(rows), selected)
        extra <- available[sample.int(
          length(available), count - length(selected), replace = FALSE
        )]
        selected <- sort(c(selected, extra))
      }
    }
    keep <- c(keep, rows[selected])
  }
  data[sort(keep), , drop = FALSE]
}

.build_subject_skeleton <- function(model, regimen) {
  endpoints <- model$public$endpoints
  design <- model$public$design
  bounds <- model$public$bounds
  rows <- list()
  tie <- 0L
  first_endpoint <- endpoints[[1L]]
  if (regimen$n_doses > 0L) {
    for (occasion in seq_along(regimen$dose_times)) {
      tie <- tie + 1L
      dose_time <- regimen$dose_times[occasion]
      amount <- regimen$amount
      if (!is.null(bounds$amt) && amount > 0) {
        amount <- .clip(
          amount * exp(stats::rnorm(1L, sd = regimen$amount_jitter_sd)),
          bounds$amt
        )
      }
      generated_rate <- if (regimen$infusion && regimen$amount != 0) {
        regimen$rate * amount / regimen$amount
      } else if (regimen$infusion) regimen$rate else 0
      if (!is.null(bounds$rate)) generated_rate <- .clip(generated_rate, bounds$rate)
      rows[[length(rows) + 1L]] <- .event_row(
        model, dose_time, amount,
        generated_rate,
        design$dose_cmt, design$dose_evid, regimen$interval,
        first_endpoint, occasion, tie
      )
      if (regimen$infusion) {
        tie <- tie + 1L
        rows[[length(rows) + 1L]] <- .event_row(
          model, dose_time + regimen$duration, amount, generated_rate,
          design$dose_cmt, design$dose_evid, regimen$interval,
          first_endpoint, occasion, tie, stop = TRUE
        )
      }
    }
  }
  for (name in names(endpoints)) {
    endpoint <- endpoints[[name]]
    timing <- model$population$timing[[name]]
    probability <- .timing_grid_probability(timing, length(endpoint$grid))
    occasion_schedule <- design$endpoint_occasion_grids[[name]]
    inferred_counts <- if (
      is.null(occasion_schedule) &&
      endpoint$alignment %in% c("dose_relative", "occasion")
    ) {
      .inferred_occasion_counts(timing, length(regimen$dose_times))
    } else NULL
    selected <- if (is.null(occasion_schedule)) {
      if (!is.null(inferred_counts)) seq_along(endpoint$grid) else
        .ensure_grid_presence(probability)
    } else integer()
    if (endpoint$alignment %in% c("dose_relative", "occasion")) {
      origins <- regimen$dose_times
      if (!length(origins)) origins <- bounds$time[1L]
      active_occasions <- seq_along(origins)
      if (is.null(occasion_schedule)) {
        if (!is.null(inferred_counts)) {
          inferred_presence <- .inferred_occasion_presence(
            timing, length(origins)
          )
          active_occasions <- which(
            stats::runif(length(origins)) < inferred_presence
          )
          if (!length(active_occasions)) {
            active_occasions <- which.max(
              inferred_presence * pmax(inferred_counts, 1)
            )
          }
        }
      }
      for (occasion in seq_along(origins)) {
        if (!occasion %in% active_occasions) next
        local_grid <- if (is.null(occasion_schedule)) {
          .generated_local_grid(
            endpoint, selected, regimen, occasion = occasion,
            timing = timing
          )
        } else {
          occasion_schedule[[as.character(occasion)]]
        }
        if (is.null(local_grid) || !length(local_grid)) next
        nominal <- origins[occasion] + local_grid
        keep <- nominal >= bounds$time[1L] & nominal <= bounds$time[2L]
        if (occasion < length(origins)) {
          keep <- keep & nominal < origins[occasion + 1L]
        }
        nominal <- nominal[keep]
        local <- local_grid[keep]
        actual <- nominal
        for (j in seq_along(nominal)) {
          tie <- tie + 1L
          rows[[length(rows) + 1L]] <- .observation_row(
            model, name, nominal[j], actual[j],
            max(0, actual[j] - origins[occasion]), occasion, tie
          )
        }
      }
    } else {
      nominal <- endpoint$grid[selected]
      nominal <- nominal[nominal >= bounds$time[1L] & nominal <= bounds$time[2L]]
      actual <- nominal
      occasions <- if (length(regimen$dose_times)) {
        pmax(1L, findInterval(actual, regimen$dose_times))
      } else rep(1L, length(actual))
      tad <- if (length(regimen$dose_times)) {
        actual - regimen$dose_times[pmin(occasions, length(regimen$dose_times))]
      } else actual - bounds$time[1L]
      for (j in seq_along(nominal)) {
        tie <- tie + 1L
        rows[[length(rows) + 1L]] <- .observation_row(
          model, name, nominal[j], actual[j], max(0, tad[j]),
          occasions[j], tie
        )
      }
    }
  }
  if (!length(rows)) stop("The private model generated no event or observation rows.", call. = FALSE)
  skeleton <- as.data.frame(do.call(rbind, lapply(rows, function(x) {
    as.data.frame(x, stringsAsFactors = FALSE, optional = TRUE)
  })), stringsAsFactors = FALSE, optional = TRUE)
  skeleton <- .fill_assigned_dose(skeleton, model$public$roles)
  skeleton <- .trim_observations_to_private_count(
    skeleton, model, regimen$observation_count
  )
  .coherent_actual_times(skeleton, model, regimen)
}

.coherent_actual_times <- function(data, model, regimen) {
  roles <- model$public$roles
  bounds <- model$public$bounds$time
  observation <- !is.na(data$.endpoint_name)
  nominal <- as.numeric(data$.nominal_internal[observation])
  actual <- .jitter_times(nominal, model$public$design$time_jitter_sd)
  if (length(regimen$dose_times)) {
    occasion <- as.integer(data$.occasion_internal[observation])
    origin <- regimen$dose_times[pmin(occasion, length(regimen$dose_times))]
    next_dose <- vapply(occasion, function(x) {
      if (x < length(regimen$dose_times)) regimen$dose_times[x + 1L] else Inf
    }, numeric(1))
    tied_to_dose <- nominal %in% regimen$dose_times
    actual[tied_to_dose] <- nominal[tied_to_dose]
    local_alignment <- vapply(data$.endpoint_name[observation], function(name) {
      model$public$endpoints[[name]]$alignment %in% c("dose_relative", "occasion")
    }, logical(1))
    actual[local_alignment] <- pmax(actual[local_alignment], origin[local_alignment])
    # Keep a dose-relative observation inside the occasion that generated it.
    # Equality with the next dose would cause standard interval assignment to
    # relabel it as an observation from the following occasion.
    upper <- next_dose
    finite_next <- is.finite(next_dose)
    upper[finite_next] <- next_dose[finite_next] -
      pmax(abs(next_dose[finite_next] - origin[finite_next]), 1) * 1e-8
    actual[local_alignment] <- pmin(actual[local_alignment], upper[local_alignment])
  }
  actual <- pmin(pmax(actual, bounds[1L]), bounds[2L])
  data[[roles$time]][observation] <- actual
  data$.tad_internal[observation] <- if (length(regimen$dose_times)) {
    occasion <- as.integer(data$.occasion_internal[observation])
    pmax(0, actual - regimen$dose_times[pmin(occasion,
                                            length(regimen$dose_times))])
  } else pmax(0, actual - bounds[1L])
  if (!is.null(roles$tad)) data[[roles$tad]][observation] <-
    data$.tad_internal[observation]
  data
}

.smooth_private_curve <- function(x) {
  if (length(x) < 3L) return(x)
  out <- x
  out[2:(length(x) - 1L)] <-
    (x[1:(length(x) - 2L)] + 2 * x[2:(length(x) - 1L)] +
       x[3:length(x)]) / 4
  out
}

.project_unimodal_at <- function(x, peak) {
  n <- length(x)
  if (n < 3L) return(x)
  peak <- min(max(as.integer(peak), 1L), n)
  left_index <- seq_len(peak)
  right_index <- peak:n
  left <- stats::isoreg(left_index, x[left_index])$yf
  right <- -stats::isoreg(seq_along(right_index), -x[right_index])$yf
  peak_value <- max(left[length(left)], right[1L])
  left[length(left)] <- peak_value
  right[1L] <- peak_value
  c(left[-length(left)], right)
}

.best_unimodal_projection <- function(x) {
  if (length(x) < 3L || any(!is.finite(x))) {
    return(list(values = x, peak = which.max(x), relative_error = 0))
  }
  candidates <- lapply(seq_along(x), function(peak) {
    .project_unimodal_at(x, peak)
  })
  error <- vapply(candidates, function(candidate) {
    sum((candidate - x)^2)
  }, numeric(1))
  peak <- which.min(error)
  values <- candidates[[peak]]
  scale <- max(diff(range(x)), sqrt(.Machine$double.eps))
  list(
    values = values, peak = peak,
    relative_error = max(abs(values - x)) / scale
  )
}

.occasion_ar1_noise <- function(data, rows, endpoint, sd) {
  if (!length(rows) || sd == 0) return(numeric(length(rows)))
  if (!endpoint$alignment %in% c("dose_relative", "occasion")) {
    return(.ar1_noise(length(rows), phi = 0.65, sd = sd))
  }
  occasion <- as.integer(data$.occasion_internal[rows])
  out <- numeric(length(rows))
  for (value in unique(occasion)) {
    index <- which(occasion == value)
    out[index] <- .ar1_noise(length(index), phi = 0.65, sd = sd)
  }
  out
}

.stabilize_unimodal_profiles <- function(latent, mean_working, clock,
                                          occasion) {
  out <- latent
  for (value in unique(occasion)) {
    index <- which(occasion == value)
    if (length(index) < 3L) next
    ordered <- index[order(clock[index])]
    peak <- which.max(mean_working[ordered])
    out[ordered] <- .project_unimodal_at(out[ordered], peak)
  }
  out
}

.generate_endpoint_values <- function(data, model) {
  roles <- model$public$roles
  for (name in names(model$public$endpoints)) {
    rows <- which(data$.endpoint_name == name)
    if (!length(rows)) next
    endpoint <- model$public$endpoints[[name]]
    trajectory <- model$population$trajectories[[name]]
    curve <- .smooth_private_curve(trajectory$mean_working)
    curve_shape <- .best_unimodal_projection(curve)
    stabilize_shape <-
      endpoint$alignment %in% c("dose_relative", "occasion") &&
      curve_shape$relative_error <= 0.10
    if (stabilize_shape) curve <- curve_shape$values
    if (endpoint$alignment %in% c("dose_relative", "occasion")) {
      clock <- as.numeric(data$.tad_internal[rows])
      mean_working <- stats::approx(
        trajectory$grid, curve, xout = clock, rule = 2, ties = mean
      )$y
    } else {
      clock <- if (!is.null(roles$nominal_time)) {
        as.numeric(data[[roles$nominal_time]][rows])
      } else as.numeric(data[[roles$time]][rows])
      mean_working <- stats::approx(
        trajectory$grid, curve, xout = clock, rule = 2, ties = mean
      )$y
      if (endpoint$alignment == "hybrid") {
        tad <- as.numeric(data$.tad_internal[rows])
        excursion <- stats::approx(
          trajectory$local_grid, trajectory$local_excursion_unit,
          xout = tad, rule = 2, ties = mean
        )$y
        mean_working <- mean_working + excursion *
          diff(.endpoint_working_bounds(endpoint)) * 0.35
      }
    }
    working_range <- diff(.endpoint_working_bounds(endpoint))
    subject_shift <- stats::rnorm(
      1L, sd = endpoint$subject_sd * working_range * 0.15
    )
    occasion_shift <- numeric(length(rows))
    if (endpoint$alignment == "occasion") {
      occasion <- as.integer(data$.occasion_internal[rows])
      shifts <- stats::rnorm(length(unique(occasion)),
                             sd = endpoint$subject_sd * working_range * 0.08)
      occasion_shift <- shifts[match(occasion, unique(occasion))]
    }
    residual <- .occasion_ar1_noise(
      data, rows, endpoint,
      sd = endpoint$residual_sd * working_range * 0.10
    )
    latent_working <- mean_working + subject_shift + occasion_shift + residual
    if (stabilize_shape) {
      latent_working <- .stabilize_unimodal_profiles(
        latent_working, mean_working, clock,
        as.integer(data$.occasion_internal[rows])
      )
    }
    latent <- .inverse_endpoint(
      latent_working, endpoint
    )
    data[[roles$dv]][rows] <- latent
    data <- .apply_censoring(data, rows, name, latent, model)
  }
  data
}

.apply_censoring <- function(data, rows, endpoint_name, latent, model) {
  roles <- model$public$roles
  if (is.null(roles$cens)) return(data)
  endpoint <- model$public$endpoints[[endpoint_name]]
  private <- model$population$censoring[[endpoint_name]] %||%
    list(left = 0, right = 0, interval = 0,
         boundary = NA_real_, other_boundary = NA_real_)
  public <- endpoint$censoring %||% list()
  cens <- rep(0, length(rows))
  limit <- rep(NA_real_, length(rows))
  reported <- latent

  if (!is.null(public$interval)) {
    lower <- public$interval[1L]
    upper <- public$interval[2L]
    selected <- latent >= lower & latent <= upper
    cens[selected] <- 1
    reported[selected] <- upper
    limit[selected] <- lower
  } else if (private$interval > 0 && is.finite(private$boundary) &&
             is.finite(private$other_boundary)) {
    lower <- min(private$boundary, private$other_boundary)
    upper <- max(private$boundary, private$other_boundary)
    selected <- latent >= lower & latent <= upper &
      stats::runif(length(rows)) < private$interval
    cens[selected] <- 1
    reported[selected] <- upper
    limit[selected] <- lower
  }
  left_boundary <- public$left %||% private$boundary
  if (is.finite(left_boundary)) {
    selected <- cens == 0 & latent <= left_boundary &
      (if (!is.null(public$left)) TRUE else
         stats::runif(length(rows)) < private$left)
    cens[selected] <- 1
    reported[selected] <- left_boundary
  }
  right_boundary <- public$right %||% private$boundary
  if (is.finite(right_boundary)) {
    selected <- cens == 0 & latent >= right_boundary &
      (if (!is.null(public$right)) TRUE else
         stats::runif(length(rows)) < private$right)
    cens[selected] <- -1
    reported[selected] <- right_boundary
  }
  data[[roles$dv]][rows] <- reported
  data[[roles$cens]][rows] <- cens
  if (!is.null(roles$limit)) data[[roles$limit]][rows] <- limit
  data
}

.generate_covariates <- function(data, model) {
  roles <- model$public$roles
  for (name in roles$covariates) {
    spec <- model$population$covariates[[name]]
    if (is.null(spec)) next
    if (spec$type == "numeric") {
      value <- .clip(stats::rnorm(1L, spec$center, max(spec$scale, 1e-8)),
                     spec$bounds)
    } else {
      value <- sample(spec$levels, 1L, prob = spec$probability)
    }
    data[[name]][] <- value
  }
  data
}

.sample_subject_property <- function(model) {
  properties <- model$population$subject_properties
  if (is.null(properties) || !length(properties$strata)) {
    return(list(values = NULL, event = NULL))
  }
  probability <- vapply(
    properties$strata, `[[`, numeric(1), "probability"
  )
  probability <- pmax(probability, 0)
  if (!sum(probability) > 0) probability[] <- 1
  stratum <- properties$strata[[sample.int(
    length(properties$strata), 1L, prob = probability
  )]]
  list(values = stratum$values, event = stratum$event)
}

.generate_subject_properties <- function(data, property) {
  if (is.null(property$values)) return(data)
  for (name in names(property$values)) {
    data[[name]][] <- property$values[[name]]
  }
  data
}

.fill_assigned_dose <- function(data, roles) {
  if (is.null(roles$assigned_dose)) return(data)
  occasion <- as.integer(data$.occasion_internal)
  amount <- suppressWarnings(as.numeric(data[[roles$amt]]))
  evid <- as.character(data[[roles$evid]])
  event <- !is.na(evid) & !evid %in% c("0", "0.0") &
    is.finite(amount) & amount > 0
  assigned <- rep(NA_real_, nrow(data))
  for (value in unique(occasion[event])) {
    nominal <- amount[which(event & occasion == value)[1L]]
    assigned[occasion == value] <- nominal
  }
  data[[roles$assigned_dose]] <- assigned
  data
}

.repair_derived_time <- function(data, roles, regimen) {
  if (!length(regimen$dose_times)) return(data)
  time <- as.numeric(data[[roles$time]])
  occasion <- pmax(1L, findInterval(time, regimen$dose_times))
  occasion <- pmin(occasion, length(regimen$dose_times))
  tad <- pmax(0, time - regimen$dose_times[occasion])
  if (!is.null(roles$occasion)) data[[roles$occasion]] <- occasion
  if (!is.null(roles$tad)) data[[roles$tad]] <- tad
  data
}

#' Generate a new PMX event dataset from a fitted private model
#'
#' This function never reads source data. Repeated calls are post-processing and
#' do not alter or consume the fitted model's privacy accounting.
#'
#' @param private_model A fitted model from [fit_private_pmx()].
#' @param n_subjects Optional positive number of new synthetic subjects. By default,
#'   generation uses the privacy-accounted subject-count release stored in the
#'   fitted model (the exact source count for the explicitly public fixture
#'   backend).
#' @param seed Ordinary reproducibility seed for post-processing generation.
#'
#' @return An ordinary data frame or tibble with the declared public schema and
#' a lightweight `pmx_privacy` attribute.
#' @export
generate_pmx <- function(private_model, n_subjects = NULL, seed = 123) {
  validate_private_model(private_model, strict = TRUE)
  if (is.null(n_subjects)) {
    n_subjects <- max(
      1L, as.integer(round(private_model$population$private_subject_count))
    )
  }
  n_subjects <- .positive_integer(n_subjects, "n_subjects")
  .with_local_seed(seed, {
    model <- private_model
    roles <- model$public$roles
    schema <- model$public$design$schema
    ids <- .new_public_ids(schema, roles$id, n_subjects)
    generated <- vector("list", n_subjects)
    for (i in seq_len(n_subjects)) {
      property <- .sample_subject_property(model)
      regimen <- .resolved_regimen(model, property$event)
      subject <- .build_subject_skeleton(model, regimen)
      subject <- .generate_endpoint_values(subject, model)
      subject <- .generate_covariates(subject, model)
      subject <- .generate_subject_properties(subject, property)
      subject[[roles$id]] <- rep(ids[i], nrow(subject))
      subject <- .repair_derived_time(subject, roles, regimen)
      order <- order(as.numeric(subject[[roles$time]]), subject$.tie_order)
      subject <- subject[order, , drop = FALSE]
      subject$.endpoint_name <- NULL
      subject$.nominal_internal <- NULL
      subject$.tad_internal <- NULL
      subject$.occasion_internal <- NULL
      subject$.tie_order <- NULL
      generated[[i]] <- subject
    }
    result <- do.call(rbind, generated)
    rownames(result) <- NULL
    result <- .restore_public_schema(result, schema, roles)
    validate_pmx(result, roles, endpoints = model$public$endpoints,
                 strict = TRUE)
    attr(result, "pmx_privacy") <- list(
      release_id = model$ledger$release_id,
      formal_dp = model$privacy$formal_dp,
      epsilon = model$privacy$epsilon, delta = model$privacy$delta,
      generation_seed = as.integer(seed),
      post_processing = TRUE,
      statement = if (isTRUE(model$privacy$formal_dp)) {
        paste0("Generated from a subject-level (", model$privacy$epsilon,
               ", ", model$privacy$delta,
               ")-differentially private model.")
      } else "Generated from an explicitly public-source fixture model."
    )
    result
  })
}
