.assert_pair <- function(x, name, finite = TRUE) {
  okay <- is.numeric(x) && length(x) == 2L && !anyNA(x) &&
    (!finite || all(is.finite(x))) && x[1L] < x[2L]
  if (!okay) {
    stop("`", name, "` must be two increasing finite numeric bounds.",
         call. = FALSE)
  }
  as.numeric(x)
}

.assert_named_pairs <- function(x, name, allow_null = TRUE) {
  if (is.null(x) && allow_null) return(list())
  if (!is.list(x) || is.null(names(x)) || any(!nzchar(names(x))) ||
      anyDuplicated(names(x))) {
    stop("`", name, "` must be a uniquely named list of bound pairs.",
         call. = FALSE)
  }
  lapply(seq_along(x), function(i) {
    .assert_pair(x[[i]], paste0(name, "$", names(x)[i]))
  }) |>
    stats::setNames(names(x))
}

#' Declare public numeric domains for private PMX fitting
#'
#' Bounds must be chosen without inspecting confidential patient values. They
#' define clipping domains and therefore enter the sensitivity argument.
#'
#' @param time Bounds for actual study time.
#' @param endpoints Named list of DV bounds, one pair per endpoint declaration.
#' @param amt,rate Optional amount and rate bounds.
#' @param covariates Named list of numeric covariate bounds.
#' @param limit Named list of censoring-limit bounds by endpoint.
#'
#' @return A public `pmx_bounds` configuration object.
#' @export
pmx_bounds <- function(time, endpoints, amt = NULL, rate = NULL,
                       covariates = NULL, limit = NULL) {
  out <- list(
    time = .assert_pair(time, "time"),
    endpoints = .assert_named_pairs(endpoints, "endpoints", FALSE),
    amt = if (is.null(amt)) NULL else .assert_pair(amt, "amt"),
    rate = if (is.null(rate)) NULL else .assert_pair(rate, "rate"),
    covariates = .assert_named_pairs(covariates, "covariates"),
    limit = .assert_named_pairs(limit, "limit")
  )
  structure(out, class = "pmx_bounds")
}

.positive_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < 1 || x != floor(x) || x > .Machine$integer.max) {
    stop("`", name, "` must be one positive integer.", call. = FALSE)
  }
  as.integer(x)
}

#' Declare subject contribution limits
#'
#' @param max_rows,max_doses,max_occasions,max_timing_cells Positive public
#'   limits applied independently to every source subject.
#' @param max_observations_per_endpoint A positive scalar or a named vector
#'   giving per-endpoint limits.
#'
#' @return A `pmx_contribution_limits` object.
#' @export
pmx_contribution_limits <- function(max_rows, max_doses, max_occasions,
                                    max_observations_per_endpoint,
                                    max_timing_cells = 12L) {
  observations <- max_observations_per_endpoint
  if (!is.numeric(observations) || !length(observations) ||
      anyNA(observations) || any(!is.finite(observations)) ||
      any(observations < 1) || any(observations != floor(observations))) {
    stop("`max_observations_per_endpoint` must contain positive integers.",
         call. = FALSE)
  }
  observations <- as.integer(observations)
  names(observations) <- names(max_observations_per_endpoint)
  structure(list(
    max_rows = .positive_integer(max_rows, "max_rows"),
    max_doses = .positive_integer(max_doses, "max_doses"),
    max_occasions = .positive_integer(max_occasions, "max_occasions"),
    max_observations_per_endpoint = observations,
    max_timing_cells = {
      value <- .positive_integer(max_timing_cells, "max_timing_cells")
      if (value < 2L) stop("`max_timing_cells` must be at least two.",
                           call. = FALSE)
      value
    }
  ), class = "pmx_contribution_limits")
}

#' Allocate an epsilon budget across private summary groups
#'
#' Values are fractions of the requested epsilon. Active groups must receive a
#' positive fraction and the total may not exceed one.
#'
#' @param subject_count,event,timing,covariates,endpoints,censoring Nonnegative
#'   budget fractions.
#'
#' @return A `pmx_budget_allocation` object.
#' @export
pmx_budget_allocation <- function(subject_count, event, timing, covariates,
                                  endpoints, censoring) {
  values <- c(
    subject_count = subject_count, event = event, timing = timing,
    covariates = covariates, endpoints = endpoints, censoring = censoring
  )
  if (!is.numeric(values) || anyNA(values) || any(!is.finite(values)) ||
      any(values < 0)) {
    stop("Budget fractions must be finite and nonnegative.", call. = FALSE)
  }
  if (sum(values) > 1 || sum(values) <= 0) {
    stop("Budget fractions must sum to a value in (0, 1].", call. = FALSE)
  }
  structure(as.list(values), class = "pmx_budget_allocation")
}

#' Capture a schema asserted to be public
#'
#' This helper records names and practical classes only. Calling it is an
#' assertion that those metadata and factor/category levels are public; do not
#' derive them from confidential data unless their release has been approved.
#'
#' @param data A public schema template.
#' @param exclude Columns omitted from generated output.
#'
#' @return A `pmx_schema` object.
#' @export
pmx_schema <- function(data, exclude = NULL) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  keep <- setdiff(names(data), exclude)
  columns <- lapply(keep, function(name) {
    x <- data[[name]]
    list(
      name = name,
      class = class(x),
      typeof = typeof(x),
      levels = if (is.factor(x)) levels(x) else NULL,
      ordered = is.ordered(x)
    )
  })
  structure(list(columns = columns, data_class = class(data)),
            class = "pmx_schema")
}

#' Declare public event-design information
#'
#' Every supplied value is treated as public protocol or schema metadata and
#' is not charged to the privacy budget. Omitted regimen and timing quantities
#' are estimated through budgeted private summaries.
#'
#' @param schema A schema from [pmx_schema()]. It is required.
#' @param dose_times,dose_interval,n_doses,dose_amount,dose_rate Optional public
#'   regimen overrides. Normally omit these so the budgeted event summaries
#'   infer the regimen from the input data. Supply them only when they are
#'   independently public protocol facts, never by inspecting confidential
#'   records outside the private fit.
#' @param infusion_duration Optional public infusion duration.
#' @param dose_evid,dose_cmt Public event and dose-compartment values.
#' @param endpoint_grids Optional named list of fixed discretization bases on
#'   each endpoint's declared scientific clock. These are not sampling
#'   schedules. When omitted, the package constructs generic bases from public
#'   bounds and contribution limits and privately learns their occupancy.
#' @param endpoint_occasion_grids Optional named list of endpoint-specific,
#'   public dose-occasion sampling schedules. This exceptional override should
#'   be used only for an independently public protocol. Each endpoint entry is
#'   a list named by positive occasion number; omitted occasions generate no
#'   observations.
#' @param endpoint_cmt Named list or vector of public observation compartments.
#' @param category_levels Named lists of allowed values for categorical
#'   covariates and subject properties. Supplying levels forces even a numeric
#'   covariate to be treated categorically. Factor levels may instead come from
#'   the public schema. These levels are public domains, not values discovered
#'   by inspecting confidential records outside the private fit.
#' @param defaults Named values for otherwise unmodeled public columns.
#' @param time_jitter_sd Nonnegative generation-time jitter scale, expressed as
#'   a fraction of the closest public grid spacing.
#' @param subject_count Optional independently public source subject count.
#'   When omitted, generation uses the privacy-accounted fitted count.
#'
#' @return A `pmx_public_design` object.
#' @export
pmx_public_design <- function(schema, dose_times = NULL,
                              dose_interval = NULL, n_doses = NULL,
                              dose_amount = NULL, dose_rate = NULL,
                              infusion_duration = NULL, dose_evid = 1,
                              dose_cmt = 1, endpoint_grids = NULL,
                              endpoint_occasion_grids = NULL,
                              endpoint_cmt = NULL, category_levels = NULL,
                              defaults = NULL,
                              time_jitter_sd = 0.02,
                              subject_count = NULL) {
  if (!inherits(schema, "pmx_schema")) {
    stop("`schema` must be created by `pmx_schema()`.", call. = FALSE)
  }
  if (!is.null(dose_times)) {
    if (!is.numeric(dose_times) || anyNA(dose_times) ||
        any(!is.finite(dose_times)) || any(diff(dose_times) < 0)) {
      stop("`dose_times` must be finite and nondecreasing.", call. = FALSE)
    }
    dose_times <- as.numeric(dose_times)
  }
  scalar_numeric <- list(
    dose_interval = dose_interval, dose_amount = dose_amount,
    dose_rate = dose_rate, infusion_duration = infusion_duration
  )
  for (name in names(scalar_numeric)) {
    value <- scalar_numeric[[name]]
    if (!is.null(value) && (!is.numeric(value) || length(value) != 1L ||
                            is.na(value) || !is.finite(value))) {
      stop("`", name, "` must be one finite number or NULL.", call. = FALSE)
    }
  }
  if (!is.null(dose_interval) && dose_interval <= 0) {
    stop("`dose_interval` must be positive when supplied.", call. = FALSE)
  }
  if (!is.null(infusion_duration) && infusion_duration < 0) {
    stop("`infusion_duration` must be nonnegative when supplied.",
         call. = FALSE)
  }
  if (!is.null(n_doses)) n_doses <- .positive_integer(n_doses, "n_doses")
  if (!is.null(subject_count)) {
    subject_count <- .positive_integer(subject_count, "subject_count")
  }
  if (!is.numeric(time_jitter_sd) || length(time_jitter_sd) != 1L ||
      is.na(time_jitter_sd) || !is.finite(time_jitter_sd) ||
      time_jitter_sd < 0) {
    stop("`time_jitter_sd` must be one finite nonnegative number.",
         call. = FALSE)
  }
  if (is.null(endpoint_grids)) endpoint_grids <- list()
  if (!is.list(endpoint_grids)) {
    stop("`endpoint_grids` must be a named list.", call. = FALSE)
  }
  if (!is.null(endpoint_cmt)) {
    if (is.null(names(endpoint_cmt)) || any(!nzchar(names(endpoint_cmt)))) {
      stop("`endpoint_cmt` must be named by endpoint.", call. = FALSE)
    }
    endpoint_cmt <- as.list(endpoint_cmt)
  } else endpoint_cmt <- list()
  for (name in names(endpoint_grids)) {
    grid <- endpoint_grids[[name]]
    if (!is.numeric(grid) || !length(grid) || anyNA(grid) ||
        any(!is.finite(grid)) || any(diff(grid) <= 0)) {
      stop("Each endpoint grid must be finite and strictly increasing.",
           call. = FALSE)
    }
    endpoint_grids[[name]] <- as.numeric(grid)
  }
  if (is.null(endpoint_occasion_grids)) endpoint_occasion_grids <- list()
  if (!is.list(endpoint_occasion_grids) ||
      (length(endpoint_occasion_grids) &&
       (is.null(names(endpoint_occasion_grids)) ||
        any(!nzchar(names(endpoint_occasion_grids))) ||
        anyDuplicated(names(endpoint_occasion_grids))))) {
    stop("`endpoint_occasion_grids` must be a uniquely named list.",
         call. = FALSE)
  }
  for (endpoint_name in names(endpoint_occasion_grids)) {
    schedule <- endpoint_occasion_grids[[endpoint_name]]
    occasions <- suppressWarnings(as.integer(names(schedule)))
    if (!is.list(schedule) || !length(schedule) || is.null(names(schedule)) ||
        any(!nzchar(names(schedule))) || anyNA(occasions) ||
        any(occasions < 1L) || anyDuplicated(occasions) ||
        any(names(schedule) != as.character(occasions))) {
      stop(
        "Each endpoint occasion schedule must be a list named by positive occasion number.",
        call. = FALSE
      )
    }
    normalized <- vector("list", length(schedule))
    names(normalized) <- as.character(occasions)
    for (i in seq_along(schedule)) {
      grid <- schedule[[i]]
      if (!is.numeric(grid) || !length(grid) || anyNA(grid) ||
          any(!is.finite(grid)) || any(grid < 0) || any(diff(grid) <= 0)) {
        stop(
          "Endpoint occasion grids must be nonnegative, finite, and strictly increasing.",
          call. = FALSE
        )
      }
      normalized[[i]] <- as.numeric(grid)
    }
    endpoint_occasion_grids[[endpoint_name]] <- normalized
  }
  if (is.null(category_levels)) category_levels <- list()
  if (!is.list(category_levels)) {
    stop("`category_levels` must be a named list.", call. = FALSE)
  }
  for (name in names(category_levels)) {
    values <- category_levels[[name]]
    if (!length(values) || anyNA(values)) {
      stop("Public category levels must be nonempty and nonmissing.",
           call. = FALSE)
    }
    category_levels[[name]] <- unique(as.character(values))
  }
  structure(list(
    schema = schema, dose_times = dose_times,
    dose_interval = dose_interval, n_doses = n_doses,
    dose_amount = dose_amount, dose_rate = dose_rate,
    infusion_duration = infusion_duration, dose_evid = dose_evid,
    dose_cmt = dose_cmt, endpoint_grids = endpoint_grids,
    endpoint_occasion_grids = endpoint_occasion_grids,
    endpoint_cmt = endpoint_cmt, category_levels = category_levels,
    defaults = defaults,
    time_jitter_sd = as.numeric(time_jitter_sd),
    subject_count = subject_count
  ), class = "pmx_public_design")
}
