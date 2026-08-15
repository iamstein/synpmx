# What kind of values an endpoint takes ---------------------------------------
#
# Blending is a weighted mean, and a weighted mean of five donors' zeros and
# ones is a number between them; the subject and residual noise terms then carry
# it outside the source range entirely. On `xgxr::mad` that turned a 0/1
# endpoint into 600 distinct values spanning -0.13 to 1.08, an ordinal 1/2/3
# endpoint into values up to 4.69, and an integer count into fractions
# (`SIM-045`). Nothing caught it: class restoration works on the DV *column*,
# and a PMX table keeps every endpoint in the same numeric column, so there is
# nothing per-endpoint left to restore.
#
# Deciding this from the data is safe here in a way that dose-basis inference is
# not. That inference asks whether a ratio *clusters*, which a study can fail to
# do for many reasons, and a wrong answer rewrites every dose amount. This one
# asks only whether every recorded value of an endpoint is a whole number, which
# the data answers outright, and the response to a false positive is to
# reproduce the granularity the source already had. It still fails closed on too
# little evidence, and `pmx_roles(endpoint_types = )` overrides it in either
# direction.
#
# The cost is worth stating where a reader will meet it: snapping a binary or
# ordinal endpoint back onto its levels means generated values ARE source
# values. Nothing else is available -- a 0/1 endpoint has no third value to emit
# -- so what protects a patient on a discrete endpoint is the schedule and
# visit-set machinery, never the distinctness of the number.

# The reported vocabulary. `"integer"` rather than `"count"` because the
# evidence is "every value is a whole number" and nothing more: warfarin's `pca`
# is a percentage recorded without decimals, not a count of anything. `"count"`
# is accepted when declaring, since that is what a pharmacometrician calls the
# endpoint they have in mind, and is reported as `"integer"`.
.endpoint_type_choices <- c("continuous", "binary", "ordinal", "integer")
.endpoint_type_synonyms <- c(count = "integer")

.canonical_endpoint_type <- function(type) {
  if (type %in% names(.endpoint_type_synonyms)) {
    return(unname(.endpoint_type_synonyms[[type]]))
  }
  type
}

# The most levels an endpoint may have and still be treated as a scale to snap
# onto. Past this the values are rounded instead, which is the same repair
# without pinning output to the exact set of values the source happened to hold.
.endpoint_ordinal_max_levels <- 12L

# Below this many observed values the answer is "not enough evidence", not
# "continuous by measurement". Six whole numbers are as likely to be a coarse
# assay as an ordinal scale.
.endpoint_type_min_observations <- 10L

.whole_numbers <- function(values) {
  all(abs(values - round(values)) <= 1e-8)
}

# Documented in `avatar-algorithm.Rmd`, Step 3. Change one, change the other.
.infer_endpoint_type <- function(values) {
  values <- values[is.finite(values)]
  if (length(values) < .endpoint_type_min_observations) {
    return(list(
      type = "continuous", levels = NULL, nonnegative = FALSE,
      reason = sprintf(
        "only %d observed value(s), too few to call; left continuous",
        length(values)
      )
    ))
  }
  if (!.whole_numbers(values)) {
    return(list(type = "continuous", levels = NULL, nonnegative = FALSE,
                reason = "not every observed value is a whole number"))
  }
  levels <- sort(unique(round(values)))
  nonnegative <- min(levels) >= 0
  if (length(levels) < 2L) {
    return(list(type = "continuous", levels = NULL, nonnegative = nonnegative,
                reason = "only one distinct observed value"))
  }
  if (identical(as.numeric(levels), c(0, 1))) {
    return(list(type = "binary", levels = levels, nonnegative = TRUE,
                reason = "every observed value is 0 or 1"))
  }
  if (length(levels) <= .endpoint_ordinal_max_levels) {
    return(list(
      type = "ordinal", levels = levels, nonnegative = nonnegative,
      reason = sprintf("%d whole-number levels: %s", length(levels),
                       paste(levels, collapse = ", "))
    ))
  }
  list(
    type = "integer", levels = NULL, nonnegative = nonnegative,
    reason = sprintf("%d whole-number levels, from %s to %s",
                     length(levels), min(levels), max(levels))
  )
}

# A declared type still reads the source for its level set, because the levels
# are what "ordinal" means operationally and asking the caller to type them out
# would only be a second chance to get them wrong.
.declared_endpoint_type <- function(type, values, endpoint_name) {
  values <- values[is.finite(values)]
  levels <- sort(unique(values))
  nonnegative <- length(levels) > 0L && min(levels) >= 0
  if (identical(type, "continuous")) {
    return(list(type = type, levels = NULL, nonnegative = nonnegative,
                reason = "declared in `pmx_roles(endpoint_types = )`"))
  }
  if (identical(type, "integer")) {
    return(list(type = type, levels = NULL, nonnegative = nonnegative,
                reason = "declared in `pmx_roles(endpoint_types = )`"))
  }
  if (identical(type, "binary")) {
    if (!length(levels)) {
      stop("Endpoint `", endpoint_name, "` is declared binary but has no ",
           "observed values.", call. = FALSE)
    }
    if (length(levels) > 2L || !all(levels %in% c(0, 1))) {
      stop("Endpoint `", endpoint_name, "` is declared binary but its ",
           "observed values are ", paste(utils::head(levels, 5L),
                                         collapse = ", "),
           if (length(levels) > 5L) ", ..." else "",
           ". Declare it as `ordinal` if it is a scale.", call. = FALSE)
    }
    return(list(type = type, levels = c(0, 1), nonnegative = TRUE,
                reason = "declared in `pmx_roles(endpoint_types = )`"))
  }
  # Ordinal. Snapping onto a large set of source values would emit source
  # measurements one at a time rather than a scale, so the cap is a refusal
  # rather than a silent widening.
  if (length(levels) > .endpoint_ordinal_max_levels) {
    stop("Endpoint `", endpoint_name, "` is declared ordinal but holds ",
         length(levels), " distinct observed values, more than the ",
         .endpoint_ordinal_max_levels, " a scale is allowed. Declare it as ",
         "`count` to round to whole numbers, or `continuous` to leave it ",
         "alone.", call. = FALSE)
  }
  if (length(levels) < 2L) {
    stop("Endpoint `", endpoint_name, "` is declared ordinal but holds fewer ",
         "than two distinct observed values.", call. = FALSE)
  }
  list(type = type, levels = levels, nonnegative = nonnegative,
       reason = "declared in `pmx_roles(endpoint_types = )`")
}

# One specification per endpoint present in `data`, declared where the caller
# said so and inferred everywhere else.
.endpoint_value_types <- function(data, roles) {
  observed <- .observation_rows(data, roles, require_present = TRUE)
  endpoint <- .endpoint(data, roles)
  names_present <- sort(unique(endpoint[observed]))
  declared <- roles$endpoint_types
  if (!is.null(declared)) {
    unknown <- setdiff(names(declared), names_present)
    if (length(unknown)) {
      stop("`endpoint_types` names endpoint(s) not present in the data: ",
           paste(unknown, collapse = ", "), ". Endpoints here are: ",
           paste(names_present, collapse = ", "), ".", call. = FALSE)
    }
  }
  values_of <- function(name) {
    suppressWarnings(as.numeric(
      data[[roles$dv]][observed & endpoint == name]
    ))
  }
  stats::setNames(lapply(names_present, function(name) {
    values <- values_of(name)
    spec <- if (!is.null(declared) && name %in% names(declared)) {
      .declared_endpoint_type(
        .canonical_endpoint_type(unname(declared[[name]])), values, name
      )
    } else {
      .infer_endpoint_type(values)
    }
    spec$declared <- !is.null(declared) && name %in% names(declared)
    spec
  }), names_present)
}

# Nearest level, or nearest whole number. `findInterval()` on the midpoints
# keeps this linear in the number of values and independent of the number of
# levels; a value exactly on a midpoint goes to the higher level.
# Documented in `avatar-algorithm.Rmd`, Step 11. Change one, change the other.
.snap_endpoint_values <- function(values, spec) {
  if (is.null(spec) || identical(spec$type, "continuous")) return(values)
  finite <- which(is.finite(values))
  if (!length(finite)) return(values)
  if (!is.null(spec$levels) && length(spec$levels) > 1L) {
    levels <- spec$levels
    midpoints <- (levels[-1L] + levels[-length(levels)]) / 2
    values[finite] <- levels[findInterval(values[finite], midpoints) + 1L]
    return(values)
  }
  values[finite] <- round(values[finite])
  if (isTRUE(spec$nonnegative)) {
    values[finite] <- pmax(values[finite], 0)
  }
  values
}

# Observation rows of a discrete endpoint whose value is not one the source
# could have held. Censored rows are excluded: censoring reports a boundary
# rather than a measurement, so a boundary that is not a level is the censoring
# working, not this mechanism failing.
.endpoint_type_violations <- function(data, roles, specs) {
  observed <- .observation_rows(data, roles, require_present = TRUE)
  if (!is.null(roles$cens) && !is.null(data[[roles$cens]])) {
    observed <- observed & .is_zero(data[[roles$cens]])
  }
  endpoint <- .endpoint(data, roles)
  counts <- vapply(names(specs), function(name) {
    spec <- specs[[name]]
    if (identical(spec$type, "continuous")) return(0L)
    values <- suppressWarnings(as.numeric(
      data[[roles$dv]][observed & endpoint == name]
    ))
    values <- values[is.finite(values)]
    if (!length(values)) return(0L)
    bad <- if (!is.null(spec$levels)) {
      !(values %in% spec$levels)
    } else {
      abs(values - round(values)) > 1e-8 |
        (isTRUE(spec$nonnegative) & values < 0)
    }
    sum(bad)
  }, integer(1))
  stats::setNames(counts, names(specs))
}

#' What kind of values each endpoint takes
#'
#' Reports, per endpoint, whether its values are continuous or discrete, and
#' where that answer came from. [synpmx_avatar()] blends real trajectories, and
#' a weighted mean of several patients' zeros and ones is a number between them,
#' so a discrete endpoint would come back continuous unless the generated values
#' are snapped back onto the levels the source used. This is the function that
#' decides which endpoints that applies to.
#'
#' An endpoint is called:
#'
#' - `"binary"` when every observed value is 0 or 1;
#' - `"ordinal"` when every observed value is a whole number and there are at
#'   most 12 distinct ones, which are then the scale;
#' - `"integer"` when every observed value is a whole number and there are more
#'   levels than that, so generated values are rounded rather than snapped onto
#'   a scale. Counts are the usual case, and `"count"` may be used to declare
#'   one, but the evidence is only that the values are whole numbers;
#' - `"continuous"` otherwise, including when there are fewer than 10 observed
#'   values, which is too few to call.
#'
#' Declare `endpoint_types` in [pmx_roles()] to override any of it. The `reason`
#' column says what the data showed, so a study whose endpoint was called wrongly
#' can be corrected without guessing at the rule.
#'
#' Snapping a binary or ordinal endpoint means generated values are source
#' values: a 0/1 endpoint has no third value to emit. What protects a patient on
#' a discrete endpoint is the visit-set and dose-schedule machinery, not the
#' distinctness of any one number.
#'
#' @param data PMX data.
#' @param roles Explicit roles from [pmx_roles()].
#'
#' @return A `pmx_endpoint_types` data frame with one row per endpoint and
#'   columns `endpoint`, `type`, `levels`, `decided_by`, and `reason`. Marked
#'   `"restricted_not_releasable"`: the level set is read from real data.
#' @seealso [pmx_roles()], [synpmx_avatar()], [synpmx_scorecard()].
#' @export
#' @examples
#' data <- pmx_simulated_fixture(30)
#' roles <- pmx_roles(
#'   id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
#'   cmt = "CMT", dvid = "DVID", covariates = "WT"
#' )
#' pmx_endpoint_types(data, roles)
pmx_endpoint_types <- function(data, roles) {
  .assert_roles(data, roles)
  specs <- .endpoint_value_types(data, roles)
  out <- data.frame(
    endpoint = names(specs),
    type = vapply(specs, function(s) s$type, character(1)),
    levels = vapply(specs, function(s) {
      if (is.null(s$levels)) "--" else paste(s$levels, collapse = ", ")
    }, character(1)),
    decided_by = vapply(specs, function(s) {
      if (isTRUE(s$declared)) "declared" else "inferred"
    }, character(1)),
    reason = vapply(specs, function(s) s$reason, character(1)),
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out <- structure(out, class = c("pmx_endpoint_types", "data.frame"))
  .mark_release(out, "restricted_not_releasable")
}

#' @export
print.pmx_endpoint_types <- function(x, ...) {
  cat("Endpoint value types\n\n")
  print(as.data.frame(x), row.names = FALSE)
  cat("\n", .wrap_plain(paste(
    "Generated values on a `binary` or `ordinal` endpoint are snapped to the",
    "levels above, and on an `integer` endpoint rounded to whole numbers.",
    "Override with `pmx_roles(endpoint_types = )`. Source-derived; not",
    "releasable unless separately public or privately budgeted."
  )), "\n", sep = "")
  invisible(x)
}

#' @exportS3Method knitr::knit_print
knit_print.pmx_endpoint_types <- function(x, ...) {
  knitr::knit_print(knitr::kable(as.data.frame(x), row.names = FALSE))
}
