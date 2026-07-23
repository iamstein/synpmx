# Baseline covariates -------------------------------------------------------
#
# Covariates exist mainly so that covariate-handling pipeline code (joins,
# filters, covariate models) has columns to run against. Fidelity is secondary.
#
# Each covariate costs exactly one budget slice, regardless of its number of
# levels: a continuous covariate releases one clipped mean, and a categorical
# covariate releases one level-count vector whose L1 sensitivity is one, because
# adding or removing a subject changes exactly one level's count by one.

#' Declare one public baseline covariate
#'
#' A covariate is either continuous, with a public plausible `range` used for
#' clipping, or categorical, with public `levels`. The range or level set must
#' be chosen without inspecting the confidential data.
#'
#' @param range Two increasing numbers bracketing a continuous covariate.
#' @param levels Character levels of a categorical covariate.
#' @param source Required provenance string.
#'
#' @return A `pmx_covariate`.
#' @export
pmx_covariate <- function(range = NULL, levels = NULL, source) {
  if (missing(source) || !is.character(source) || length(source) != 1L ||
      !nzchar(trimws(source))) {
    stop("`source` is required for a covariate.", call. = FALSE)
  }
  if (is.null(range) == is.null(levels)) {
    stop("Supply exactly one of `range` (continuous) or `levels` ",
         "(categorical).", call. = FALSE)
  }
  if (!is.null(range)) {
    if (!is.numeric(range) || length(range) != 2L || anyNA(range) ||
        any(!is.finite(range)) || range[1L] >= range[2L]) {
      stop("`range` must be two increasing finite numbers.", call. = FALSE)
    }
    return(structure(list(type = "continuous", range = as.numeric(range),
                          source = source), class = "pmx_covariate"))
  }
  levels <- as.character(levels)
  if (!length(levels) || anyNA(levels) || anyDuplicated(levels)) {
    stop("`levels` must be unique non-missing labels.", call. = FALSE)
  }
  structure(list(type = "categorical", levels = levels, source = source),
            class = "pmx_covariate")
}

#' Collect public covariate declarations
#'
#' @param ... Named [pmx_covariate()] objects. Each name is a column in the
#'   source data and in the generated output.
#'
#' @return A `pmx_covariates` object, or `NULL` if nothing is supplied.
#' @export
pmx_covariates <- function(...) {
  covariates <- list(...)
  if (!length(covariates)) return(NULL)
  if (is.null(names(covariates)) || any(!nzchar(names(covariates))) ||
      anyDuplicated(names(covariates))) {
    stop("`pmx_covariates()` needs uniquely named `pmx_covariate()` objects.",
         call. = FALSE)
  }
  if (!all(vapply(covariates, inherits, logical(1), "pmx_covariate"))) {
    stop("Every element must come from `pmx_covariate()`.", call. = FALSE)
  }
  structure(covariates, class = "pmx_covariates")
}

#' @export
print.pmx_covariates <- function(x, ...) {
  cat("Public covariates\n")
  for (name in names(x)) {
    cov <- x[[name]]
    if (cov$type == "continuous") {
      cat(sprintf("  %s: continuous [%g, %g]\n", name, cov$range[1L],
                  cov$range[2L]))
    } else {
      cat(sprintf("  %s: categorical {%s}\n", name,
                  paste(cov$levels, collapse = ", ")))
    }
  }
  invisible(x)
}

# One value per subject: the first non-missing entry in that subject's rows.
.subject_covariate <- function(values, id) {
  vapply(split(values, id), function(v) {
    v <- v[!is.na(v)]
    if (length(v)) v[[1L]] else NA
  }, values[[1L]][NA][1L])
}

# Release a differentially private summary of each declared covariate. Returns a
# list of per-covariate summaries used later during generation.
.covariate_summaries <- function(data, id, covariates, accountant, per_query) {
  if (is.null(covariates)) return(NULL)
  summaries <- list()
  for (name in names(covariates)) {
    cov <- covariates[[name]]
    if (is.null(data[[name]])) {
      stop("Covariate column `", name, "` is not in the data.", call. = FALSE)
    }
    if (cov$type == "continuous") {
      per_subject <- .subject_covariate(
        suppressWarnings(as.numeric(data[[name]])), id
      )
      unit <- .to_unit(per_subject[is.finite(per_subject)], cov$range)
      total <- .private_release(accountant, paste0("covariate_", name),
                                sum(unit), sensitivity = 1, epsilon = per_query)
      mean_unit <- min(max(as.numeric(total) / max(length(unit), 1), 0), 1)
      summaries[[name]] <- list(type = "continuous", range = cov$range,
                                mean = .from_unit(mean_unit, cov$range))
    } else {
      per_subject <- as.character(.subject_covariate(
        as.character(data[[name]]), id
      ))
      counts <- as.numeric(table(factor(per_subject, levels = cov$levels)))
      # One subject occupies one level, so the whole count vector has L1
      # sensitivity one.
      released <- .private_release(accountant, paste0("covariate_", name),
                                   counts, sensitivity = 1, epsilon = per_query)
      released <- pmax(as.numeric(released), 0)
      if (!sum(released) > 0) released[] <- 1
      summaries[[name]] <- list(type = "categorical", levels = cov$levels,
                                prob = released / sum(released))
    }
  }
  summaries
}

# Draw a covariate table, one row per generated subject. Uses released summaries
# when present (calibrated mode) and the public declaration otherwise (prior
# mode). Continuous spread is a public assumption; only the centre is calibrated.
.draw_covariate_table <- function(covariates, summaries, n) {
  if (is.null(covariates)) return(NULL)
  out <- list()
  for (name in names(covariates)) {
    cov <- covariates[[name]]
    summary <- summaries[[name]]
    if (cov$type == "continuous") {
      centre <- summary$mean %||% mean(cov$range)
      spread <- diff(cov$range) / 6            # public: range spans ~6 SD
      values <- stats::rnorm(n, centre, spread)
      out[[name]] <- .clip(values, cov$range)
    } else {
      prob <- summary$prob %||% rep(1 / length(cov$levels), length(cov$levels))
      out[[name]] <- sample(cov$levels, n, replace = TRUE, prob = prob)
    }
  }
  as.data.frame(out, stringsAsFactors = FALSE)
}
