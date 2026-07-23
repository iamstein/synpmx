#' Declare endpoint scientific-clock behavior
#'
#' Endpoint alignment is public scientific metadata. It controls how a small
#' fixed-dimensional trajectory summary is built and how new trajectories are
#' generated; it is not a PK or PD model.
#'
#' @param dvid Public DVID value for this endpoint, or `NULL` when no DVID
#'   column is used.
#' @param alignment One of `"dose_relative"`, `"study_time"`, `"occasion"`,
#'   or `"hybrid"`.
#' @param transform One of `"log"`, `"identity"`, or `"auto"`. `"auto"` uses
#'   only the public DV bounds: a nonnegative domain uses an offset log scale.
#' @param shape One of `"occasion"` or `"global"`; this is a broad public shape
#'   expectation, not a fitted structural model.
#' @param units Optional public unit label.
#' @param grid Optional strictly increasing public grid on the declared clock.
#'   This is a discretization basis, not a sampling schedule. When omitted, a
#'   generic basis is constructed from public bounds and contribution limits;
#'   sampling-cell occupancy is then learned from the fitted data.
#' @param cmt Optional public observation compartment value.
#' @param subject_sd,residual_sd Public generation variability multipliers.
#' @param censoring Optional public list with `left`, `right`, or a two-value
#'   `interval`. Source-dependent censoring frequencies are separately private.
#'
#' @return A `pmx_endpoint` declaration.
#' @export
pmx_endpoint <- function(dvid = NULL,
                         alignment,
                         transform = c("auto", "log", "identity"),
                         shape, units = NULL,
                         grid = NULL, cmt = NULL, subject_sd = 0.20,
                         residual_sd = 0.08, censoring = NULL) {
  alignment <- match.arg(
    alignment, c("dose_relative", "study_time", "occasion", "hybrid")
  )
  transform <- match.arg(transform)
  shape <- match.arg(shape, c("occasion", "global"))
  if (!is.null(dvid) && (length(dvid) != 1L || is.na(dvid))) {
    stop("`dvid` must be one public endpoint value or NULL.", call. = FALSE)
  }
  if (!is.null(units) && (!is.character(units) || length(units) != 1L ||
                          is.na(units) || !nzchar(units))) {
    stop("`units` must be one non-empty string or NULL.", call. = FALSE)
  }
  if (!is.null(grid)) {
    if (!is.numeric(grid) || length(grid) < 2L || anyNA(grid) ||
        any(!is.finite(grid)) || any(diff(grid) <= 0)) {
      stop("`grid` must be finite and strictly increasing.", call. = FALSE)
    }
    grid <- as.numeric(grid)
  }
  for (name in c("subject_sd", "residual_sd")) {
    value <- get(name)
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        !is.finite(value) || value < 0) {
      stop("`", name, "` must be one finite nonnegative number.",
           call. = FALSE)
    }
  }
  if (!is.null(censoring)) {
    if (!is.list(censoring) ||
        !all(names(censoring) %in% c("left", "right", "interval"))) {
      stop("`censoring` must be a list containing left, right, or interval.",
           call. = FALSE)
    }
    if (!is.null(censoring$left) &&
        (!is.numeric(censoring$left) || length(censoring$left) != 1L ||
         !is.finite(censoring$left))) {
      stop("`censoring$left` must be one finite boundary.", call. = FALSE)
    }
    if (!is.null(censoring$right) &&
        (!is.numeric(censoring$right) || length(censoring$right) != 1L ||
         !is.finite(censoring$right))) {
      stop("`censoring$right` must be one finite boundary.", call. = FALSE)
    }
    if (!is.null(censoring$interval)) {
      censoring$interval <- .assert_pair(censoring$interval,
                                         "censoring$interval")
    }
  }
  structure(list(
    dvid = dvid, alignment = alignment, transform = transform,
    shape = shape, units = units, grid = grid, cmt = cmt,
    subject_sd = as.numeric(subject_sd), residual_sd = as.numeric(residual_sd),
    censoring = censoring
  ), class = "pmx_endpoint")
}

#' @export
print.pmx_endpoint <- function(x, ...) {
  cat("PMX endpoint:\n",
      "  dvid: ", if (is.null(x$dvid)) "<implicit>" else x$dvid, "\n",
      "  alignment: ", x$alignment, "\n",
      "  transform: ", x$transform, "\n",
      "  shape: ", x$shape, "\n", sep = "")
  invisible(x)
}

.normalize_endpoints <- function(endpoints, roles, bounds, public_design,
                                 contribution_limits) {
  if (!is.list(endpoints) || !length(endpoints) || is.null(names(endpoints)) ||
      any(!nzchar(names(endpoints))) || anyDuplicated(names(endpoints)) ||
      !all(vapply(endpoints, inherits, logical(1), "pmx_endpoint"))) {
    stop("`endpoints` must be a uniquely named list of `pmx_endpoint()` objects.",
         call. = FALSE)
  }
  if (is.null(roles$dvid) && length(endpoints) != 1L) {
    stop("Exactly one endpoint is required when no DVID role is declared.",
         call. = FALSE)
  }
  if (is.null(roles$dvid) && !is.null(endpoints[[1L]]$dvid)) {
    stop("Endpoint `dvid` must be NULL when no DVID role is declared.",
         call. = FALSE)
  }
  if (!is.null(roles$dvid)) {
    values <- vapply(endpoints, function(x) {
      if (is.null(x$dvid)) NA_character_ else as.character(x$dvid)
    }, character(1))
    if (anyNA(values) || anyDuplicated(values)) {
      stop("Every endpoint needs a unique public `dvid` value.", call. = FALSE)
    }
  }
  if (!setequal(names(endpoints), names(bounds$endpoints))) {
    stop("Endpoint names and `bounds$endpoints` names must match exactly.",
         call. = FALSE)
  }
  unknown_schedules <- setdiff(
    names(public_design$endpoint_occasion_grids), names(endpoints)
  )
  if (length(unknown_schedules)) {
    stop("Public occasion grids refer to undeclared endpoints: ",
         paste(unknown_schedules, collapse = ", "), ".", call. = FALSE)
  }
  incompatible_schedules <- names(public_design$endpoint_occasion_grids)[
    !vapply(endpoints[names(public_design$endpoint_occasion_grids)], function(x) {
      x$alignment %in% c("dose_relative", "occasion")
    }, logical(1))
  ]
  if (length(incompatible_schedules)) {
    stop("Public occasion grids require dose-relative or occasion alignment: ",
         paste(incompatible_schedules, collapse = ", "), ".", call. = FALSE)
  }
  observation_limits <- contribution_limits$max_observations_per_endpoint
  if (length(observation_limits) > 1L &&
      (is.null(names(observation_limits)) ||
       !all(names(endpoints) %in% names(observation_limits)))) {
    stop("Named observation contribution limits are required for all endpoints.",
         call. = FALSE)
  }

  for (name in names(endpoints)) {
    endpoint <- endpoints[[name]]
    grid <- endpoint$grid
    if (is.null(grid)) grid <- public_design$endpoint_grids[[name]]
    automatic_grid <- is.null(grid)
    if (is.null(grid)) {
      cells <- contribution_limits$max_timing_cells
      if (endpoint$alignment %in% c("dose_relative", "occasion")) {
        # This is a generic, source-independent basis rather than a disclosed
        # visit schedule. With no public interval, the public maximum number of
        # occasions supplies a coarse local horizon. Log spacing gives PK-like
        # clocks useful early-time resolution without inspecting source times.
        generic_interval <- diff(bounds$time) /
          max(contribution_limits$max_occasions - 1L, 1L)
        upper <- min(
          diff(bounds$time),
          public_design$dose_interval %||% generic_interval
        )
        grid <- expm1(seq(
          0, log1p(upper), length.out = cells
        ))
      } else {
        grid <- seq(bounds$time[1L], bounds$time[2L], length.out = cells)
      }
    }
    endpoint$grid <- unique(as.numeric(grid))
    endpoint$grid_automatic <- automatic_grid
    endpoint$grid_horizon <- if (
      automatic_grid &&
      endpoint$alignment %in% c("dose_relative", "occasion")
    ) max(endpoint$grid) else NULL
    if (length(endpoint$grid) > contribution_limits$max_timing_cells) {
      index <- unique(round(seq(1, length(endpoint$grid),
                                length.out = contribution_limits$max_timing_cells)))
      endpoint$grid <- endpoint$grid[index]
    }
    if (length(endpoint$grid) < 2L || any(diff(endpoint$grid) <= 0)) {
      stop("Every endpoint requires at least two increasing public grid cells.",
           call. = FALSE)
    }
    endpoint$bound <- bounds$endpoints[[name]]
    endpoint$transform_resolved <- if (endpoint$transform == "auto") {
      if (endpoint$bound[1L] >= 0) "log" else "identity"
    } else endpoint$transform
    endpoint$offset <- if (endpoint$transform_resolved == "log") {
      max(diff(endpoint$bound) * 1e-6,
          if (endpoint$bound[1L] > 0) endpoint$bound[1L] * 0.5 else 0,
          sqrt(.Machine$double.eps))
    } else 0
    endpoint$observation_limit <- if (length(observation_limits) == 1L) {
      unname(observation_limits)
    } else unname(observation_limits[[name]])
    if (endpoint$alignment == "hybrid") {
      upper <- min(diff(bounds$time),
                   public_design$dose_interval %||% diff(bounds$time))
      endpoint$local_grid <- seq(0, upper,
                                 length.out = length(endpoint$grid))
    }
    endpoints[[name]] <- endpoint
  }
  endpoints
}

.endpoint_name_for_rows <- function(data, roles, endpoints) {
  if (is.null(roles$dvid)) return(rep(names(endpoints)[1L], nrow(data)))
  values <- as.character(data[[.dvid_primary(roles)]])
  declared <- vapply(endpoints, function(x) as.character(x$dvid), character(1))
  names(declared) <- names(endpoints)
  names(declared)[match(values, declared)]
}

.clip <- function(x, bounds) pmin(pmax(as.numeric(x), bounds[1L]), bounds[2L])

.to_unit <- function(x, bounds) (.clip(x, bounds) - bounds[1L]) / diff(bounds)

.from_unit <- function(x, bounds) bounds[1L] + pmin(pmax(x, 0), 1) * diff(bounds)

.transform_endpoint <- function(x, endpoint) {
  x <- .clip(x, endpoint$bound)
  if (endpoint$transform_resolved == "log") log(x + endpoint$offset) else x
}

.inverse_endpoint <- function(x, endpoint) {
  out <- if (endpoint$transform_resolved == "log") {
    exp(x) - endpoint$offset
  } else x
  .clip(out, endpoint$bound)
}

.endpoint_working_bounds <- function(endpoint) {
  .transform_endpoint(endpoint$bound, endpoint)
}
