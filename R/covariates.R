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

#' Declare bootstrap-resampled covariates by column name
#'
#' A low-ceremony alternative to [pmx_covariates()] for a long list of
#' covariates whose fidelity does not matter. Instead of a public range or level
#' set per column, the columns are named and their values are drawn directly
#' from the source data: a uniform draw over the (clipped) observed range for
#' continuous columns, and a proportional resample for categorical ones. Column
#' type is detected from the data at fit time.
#'
#' This is the approach used by Novartis's `synadam`, and it is **not**
#' differentially private: it exposes the data-derived support of each column. A
#' model that uses it is marked as having non-private covariates, and its
#' privacy report says so. Use it only inside a trusted environment, and never
#' when the covariate columns may cross a trust boundary.
#'
#' @param names Character vector of covariate column names.
#' @param clip Two probabilities giving the quantiles a continuous column is
#'   clipped to before its range is taken, so the exact minimum and maximum are
#'   not exposed. Defaults to the 1st and 99th percentiles. Pass `NULL` to use
#'   the raw observed minimum and maximum, matching `synadam` exactly.
#'
#' @return A `pmx_covariates` object of bootstrap covariates.
#' @export
pmx_covariates_auto <- function(names, clip = c(0.01, 0.99)) {
  if (!is.character(names) || !length(names) || anyNA(names) ||
      any(!nzchar(names)) || anyDuplicated(names)) {
    stop("`names` must be unique non-empty covariate column names.",
         call. = FALSE)
  }
  if (!is.null(clip)) {
    if (!is.numeric(clip) || length(clip) != 2L || anyNA(clip) ||
        any(clip < 0) || any(clip > 1) || clip[1L] >= clip[2L]) {
      stop("`clip` must be two increasing probabilities in [0, 1], or NULL.",
           call. = FALSE)
    }
    clip <- as.numeric(clip)
  }
  covariates <- stats::setNames(lapply(names, function(nm) {
    structure(list(type = "bootstrap", clip = clip), class = "pmx_covariate")
  }), names)
  structure(covariates, class = "pmx_covariates")
}

# TRUE if any covariate is bootstrap-resampled, which makes the covariate block
# non-private and must be surfaced in the privacy report.
.covariates_have_bootstrap <- function(covariates) {
  !is.null(covariates) &&
    any(vapply(covariates, function(c) identical(c$type, "bootstrap"),
               logical(1)))
}

#' @export
print.pmx_covariates <- function(x, ...) {
  cat("Covariates\n")
  for (name in names(x)) {
    cov <- x[[name]]
    if (cov$type == "continuous") {
      cat(sprintf("  %s: continuous [%g, %g] (public, DP)\n", name,
                  cov$range[1L], cov$range[2L]))
    } else if (cov$type == "categorical") {
      cat(sprintf("  %s: categorical {%s} (public, DP)\n", name,
                  paste(cov$levels, collapse = ", ")))
    } else {
      clip <- if (is.null(cov$clip)) "min/max" else
        sprintf("%g-%g quantile", cov$clip[1L], cov$clip[2L])
      cat(sprintf("  %s: bootstrap, %s (NOT DP)\n", name, clip))
    }
  }
  if (.covariates_have_bootstrap(x)) {
    cat("  note: bootstrap covariates are resampled from the data and are ",
        "not differentially private.\n", sep = "")
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

# Summarize each declared covariate. DP-declared covariates (continuous or
# categorical, with a public range or levels) go through the accountant with
# sensitivity one. Bootstrap covariates are summarized directly from the data,
# consume no budget, and are not differentially private.
.covariate_summaries <- function(data, id, covariates, accountant, per_query,
                                 denominator = NULL) {
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
      denominator <- denominator %||% length(unit)
      mean_unit <- min(max(as.numeric(total) / max(denominator, 1), 0), 1)
      summaries[[name]] <- list(type = "continuous", range = cov$range,
                                mean = .from_unit(mean_unit, cov$range))
    } else if (cov$type == "categorical") {
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
    } else {
      summaries[[name]] <- .bootstrap_summary(
        .subject_covariate(data[[name]], id), cov$clip
      )
    }
  }
  summaries
}

# synadam-style summary of one column, computed directly from the data with no
# privacy accounting. Continuous columns keep a (clipped) range for a uniform
# draw; other columns keep their observed values for a proportional resample.
.bootstrap_summary <- function(values, clip) {
  numeric <- suppressWarnings(as.numeric(values))
  is_continuous <- is.numeric(values) ||
    (mean(is.finite(numeric)) > 0.9 &&
       length(unique(numeric[is.finite(numeric)])) > 10L)
  if (is_continuous) {
    finite <- numeric[is.finite(numeric)]
    if (!length(finite)) {
      return(list(type = "bootstrap_continuous", range = c(0, 1),
                  integer = FALSE))
    }
    bounds <- if (is.null(clip)) {
      range(finite)
    } else {
      unname(stats::quantile(finite, clip, names = FALSE, type = 7))
    }
    if (bounds[1L] >= bounds[2L]) bounds <- range(finite)
    list(type = "bootstrap_continuous", range = as.numeric(bounds),
         integer = all(finite == round(finite)))
  } else {
    observed <- as.character(values)
    observed <- observed[!is.na(observed)]
    if (!length(observed)) observed <- NA_character_
    list(type = "bootstrap_categorical", values = observed)
  }
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
    if (cov$type == "bootstrap") {
      if (is.null(summary)) {
        stop("Bootstrap covariate `", name, "` needs the data. Declare it in ",
             "`synpmx_calibrated()`, not in prior-mode `synpmx_prior()`.",
             call. = FALSE)
      }
      if (summary$type == "bootstrap_continuous") {
        # synadam: a uniform draw over the (clipped) observed range.
        values <- stats::runif(n, summary$range[1L], summary$range[2L])
        out[[name]] <- if (summary$integer) round(values) else values
      } else {
        # Proportional resample, so the observed level frequencies are kept.
        out[[name]] <- sample(summary$values, n, replace = TRUE)
      }
    } else if (cov$type == "continuous") {
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
