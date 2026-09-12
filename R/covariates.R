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
#' A covariate is either continuous, with a public plausible `range`, or
#' categorical, with public `levels`. The range or level set must be chosen
#' without inspecting the confidential data.
#'
#' @section The range is a bound, and the distribution is separate:
#' `range` is the clipping bound that makes a differentially private release of
#' this covariate possible, so it should be generous: a value outside it is
#' pulled to the edge. `median` and `cv` say where the population actually sits
#' inside that bound. Supplying neither leaves the generator with only the
#' range to work from, and it falls back to a normal draw centred on the
#' midpoint with the range spanning six standard deviations. That couples two
#' unrelated things -- widening the bound for safety also widens the
#' distribution -- and for a right-skewed covariate such as body weight the
#' shape is wrong as well. State `median` and `cv` for anything whose
#' distribution matters, and see [pmx_covariates_reference()] for defensible
#' starting values.
#'
#' For a categorical covariate, `prob` is the same point: without it every
#' level is equally likely, which is rarely what a trial looked like.
#'
#' @param range Two increasing numbers bracketing a continuous covariate.
#' @param levels Character levels of a categorical covariate.
#' @param source Required provenance string.
#' @param median Typical value of a continuous covariate. For the default
#'   lognormal this is the median, equivalently the geometric mean, following
#'   the same convention as a population parameter in
#'   [pmx_structural_model()].
#' @param cv Coefficient of variation of a continuous covariate, as a
#'   proportion: `0.18` for 18 percent. Requires `median`.
#' @param distribution `"lognormal"` (the default when `median` and `cv` are
#'   given) or `"normal"`. Lognormal is right-skewed and cannot go negative,
#'   which suits body size and clearance; normal suits age and a laboratory
#'   value that is roughly symmetric.
#' @param prob Probability per level of a categorical covariate, in the order
#'   of `levels`. Normalized if it does not sum to one.
#' @param integer Round draws to whole numbers. For a covariate recorded as a
#'   count of years or a score.
#'
#' @return A `pmx_covariate`.
#' @seealso [pmx_covariates()] to collect them,
#'   [pmx_covariates_reference()] for reference-population values,
#'   [pmx_covariates_auto()] to resample from the data instead.
#' @export
#' @examples
#' # Body weight: a generous bound, and the distribution stated separately.
#' pmx_covariate(range = c(35, 160), median = 75, cv = 0.18,
#'               source = "protocol inclusion criteria")
#'
#' # A trial that enrolled seven men for every three women.
#' pmx_covariate(levels = c("M", "F"), prob = c(0.7, 0.3),
#'               source = "protocol enrollment targets")
pmx_covariate <- function(range = NULL, levels = NULL, source, median = NULL,
                          cv = NULL, distribution = NULL, prob = NULL,
                          integer = FALSE) {
  if (missing(source) || !is.character(source) || length(source) != 1L ||
      !nzchar(trimws(source))) {
    stop("`source` is required for a covariate.", call. = FALSE)
  }
  if (is.null(range) == is.null(levels)) {
    stop("Supply exactly one of `range` (continuous) or `levels` ",
         "(categorical).", call. = FALSE)
  }
  if (!is.logical(integer) || length(integer) != 1L || is.na(integer)) {
    stop("`integer` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(range)) {
    if (!is.null(prob)) {
      stop("`prob` describes a categorical covariate; a continuous one takes ",
           "`median` and `cv`.", call. = FALSE)
    }
    if (!is.numeric(range) || length(range) != 2L || anyNA(range) ||
        any(!is.finite(range)) || range[1L] >= range[2L]) {
      stop("`range` must be two increasing finite numbers.", call. = FALSE)
    }
    range <- as.numeric(range)
    if (is.null(median) && !is.null(cv)) {
      stop("`cv` needs `median`: a coefficient of variation is a proportion ",
           "of the typical value, so there is nothing for it to scale.",
           call. = FALSE)
    }
    if (!is.null(median)) {
      if (!is.numeric(median) || length(median) != 1L || !is.finite(median)) {
        stop("`median` must be one finite number.", call. = FALSE)
      }
      median <- as.numeric(median)
      if (median < range[1L] || median > range[2L]) {
        stop("`median` (", median, ") is outside `range` [", range[1L], ", ",
             range[2L], "]. Every draw would be clipped to one edge.",
             call. = FALSE)
      }
    }
    if (!is.null(cv)) {
      if (!is.numeric(cv) || length(cv) != 1L || !is.finite(cv) || cv <= 0) {
        stop("`cv` must be one positive number, as a proportion.",
             call. = FALSE)
      }
      cv <- as.numeric(cv)
    }
    distribution <- .covariate_distribution(distribution, median, cv, range)
    return(structure(list(type = "continuous", range = range, median = median,
                          cv = cv, distribution = distribution,
                          integer = integer, source = source),
                     class = "pmx_covariate"))
  }
  if (!is.null(median) || !is.null(cv)) {
    stop("`median` and `cv` describe a continuous covariate; a categorical ",
         "one takes `prob`.", call. = FALSE)
  }
  levels <- as.character(levels)
  if (!length(levels) || anyNA(levels) || anyDuplicated(levels)) {
    stop("`levels` must be unique non-missing labels.", call. = FALSE)
  }
  if (!is.null(prob)) {
    if (!is.numeric(prob) || length(prob) != length(levels) || anyNA(prob) ||
        any(!is.finite(prob)) || any(prob < 0) || !sum(prob) > 0) {
      stop("`prob` must be one non-negative number per level, not all zero.",
           call. = FALSE)
    }
    prob <- as.numeric(prob) / sum(prob)
  }
  structure(list(type = "categorical", levels = levels, prob = prob,
                 source = source), class = "pmx_covariate")
}

# A lognormal median is the geometric mean, which is what `.draw_subject_params`
# already assumes for a population parameter, so it is the default wherever a
# median and a spread are stated. `normal` has to be asked for, and only a
# stated spread makes either choice meaningful.
.covariate_distribution <- function(distribution, median, cv, range) {
  if (is.null(distribution)) {
    return(if (is.null(median) || is.null(cv)) "range" else "lognormal")
  }
  distribution <- match.arg(distribution, c("lognormal", "normal"))
  if (is.null(median) || is.null(cv)) {
    stop("`distribution = \"", distribution, "\"` needs `median` and `cv`. ",
         "Without them only the range is available and the draw is normal ",
         "over it.", call. = FALSE)
  }
  if (identical(distribution, "lognormal") && range[1L] < 0) {
    stop("A lognormal covariate cannot take a negative `range` lower bound.",
         call. = FALSE)
  }
  distribution
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

# Reference covariate distributions ------------------------------------------
#
# Round numbers for a broad adult population, chosen so that covariate-handling
# code has plausible columns to run against. They are the same kind of input as
# the built-in structural models: illustrative on purpose, and not a rendering
# of any particular trial. Two things are worth knowing about how far they
# travel.
#
# The CV is the durable half. Measured across the adult studies in this
# package's own surveys, body weight runs at a CV of 14 to 18 percent whatever
# the median is, over medians from 70 to 117 kg. So the shape and the spread
# transfer between populations and the median does not, which is why `medians`
# is the argument this function offers. `test-covariates.R` holds the reference
# median against those studies rather than leaving the claim in a comment.
#
# Age, sex and race are not population constants at all. They are set by the
# protocol's inclusion criteria and by where it enrolled: the three studies
# here with an age run at medians of 31, 46 and 61, and the three with a sex
# split at 87/13, 50/50 and 70/30. The entries below exist so a column is
# there, and anyone who cares what is in it should state it.
.covariate_reference <- list(
  WT = list(
    quantity = "body weight", unit = "kg", range = c(35, 160),
    median = 75, cv = 0.18, distribution = "lognormal"
  ),
  HT = list(
    quantity = "height", unit = "cm", range = c(130, 210),
    median = 170, cv = 0.06, distribution = "lognormal"
  ),
  BMI = list(
    quantity = "body mass index", unit = "kg/m2", range = c(15, 55),
    median = 26, cv = 0.18, distribution = "lognormal"
  ),
  BSA = list(
    quantity = "body surface area", unit = "m2", range = c(1.1, 2.9),
    median = 1.9, cv = 0.12, distribution = "lognormal"
  ),
  AGE = list(
    quantity = "age", unit = "years", range = c(18, 90),
    median = 45, cv = 0.30, distribution = "normal", integer = TRUE
  ),
  CRCL = list(
    quantity = "creatinine clearance", unit = "mL/min", range = c(15, 200),
    median = 100, cv = 0.25, distribution = "lognormal"
  ),
  EGFR = list(
    quantity = "estimated glomerular filtration rate",
    unit = "mL/min/1.73m2", range = c(15, 180),
    median = 95, cv = 0.25, distribution = "lognormal"
  ),
  ALB = list(
    quantity = "serum albumin", unit = "g/L", range = c(25, 55),
    median = 42, cv = 0.12, distribution = "normal"
  ),
  SEX = list(
    quantity = "sex", levels = c("M", "F"), prob = c(0.5, 0.5)
  ),
  RACE = list(
    quantity = "race",
    levels = c("White", "Black or African American", "Asian", "Other"),
    prob = c(0.65, 0.12, 0.15, 0.08)
  )
)

# Column names that occur in the wild for a quantity the catalogue holds. Only
# the ones actually seen: `WEIGHTB` in `case1_pkpd` and `mad`, `HGT` in
# `nimoData`. Anything else is named explicitly, which is what the named-vector
# form of `names` is for.
.covariate_reference_aliases <- c(
  WEIGHTB = "WT", WEIGHT = "WT", WGT = "WT", BW = "WT",
  HGT = "HT", HEIGHT = "HT",
  GENDER = "SEX", AGEY = "AGE", CLCR = "CRCL", CRCLN = "CRCL"
)

#' Reference distributions for the usual baseline covariates
#'
#' Builds [pmx_covariate()] declarations for body size, age, sex and the common
#' renal and hepatic laboratory values, so that a public-model generator has
#' plausible covariate columns without each distribution being elicited by
#' hand. The values are round numbers for a broad adult population.
#'
#' @section What these numbers are, and are not:
#' They are the same kind of input as the built-in structural models in
#' [pmx_structural_model()]: illustrative on purpose, so that covariate-handling
#' code has something to run against. They are not a rendering of any
#' particular trial and were not measured from one.
#'
#' The spread is the half that travels. Body weight holds much the same
#' coefficient of variation across adult populations whatever its median is, so
#' the shape and the CV are reusable and the median is not: an ordinary adult
#' cohort and a cohort selected for obesity differ in the median alone.
#' `medians` is how to move it, and a regression test holds the reference
#' median against the adult studies in this package's own surveys.
#'
#' Age, sex and race are not population constants. The protocol's inclusion
#' criteria and where it enrolled decide them, and a phase 1 healthy-volunteer
#' cohort, a renal-impairment study and an oncology trial have nothing in
#' common here. The entries put a column in the data. State your own if
#' anything downstream reads it. No dataset in this package carries a race
#' column at all.
#'
#' Nothing here fits a neonatal or paediatric cohort, where body size differs
#' from an adult by more than an order of magnitude. Declare those with
#' [pmx_covariate()].
#'
#' Before generated data crosses a trust boundary, replace these with the
#' protocol's own criteria or a published description of the population, and
#' record that in each covariate's `source`.
#'
#' @param names Covariate column names. A plain vector takes each name as the
#'   quantity, so `c("WT", "AGE")` is body weight and age. A named vector maps
#'   a column to a quantity, so `c(BWT = "WT")` puts the weight distribution in
#'   a column called `BWT`. A handful of common spellings resolve on their own:
#'   `WEIGHTB`, `WEIGHT`, `WGT` and `BW` to `WT`, `HGT` and `HEIGHT` to `HT`,
#'   `GENDER` to `SEX`, `CLCR` to `CRCL`.
#' @param medians Optional named numeric replacing the reference median for
#'   those columns, in the units listed by [pmx_covariate_reference_table()].
#'   Named by output column, so `c(BWT = 82)` goes with `names = c(BWT = "WT")`.
#' @param source Provenance recorded on every covariate built here. The default
#'   says these are package reference values rather than a study's own.
#'
#' @return A `pmx_covariates` object.
#' @seealso [pmx_covariate_reference_table()] for what is available,
#'   [pmx_covariate()] to state one by hand,
#'   [pmx_covariates_auto()] to resample from the data instead.
#' @export
#' @examples
#' pmx_covariates_reference(c("WT", "AGE", "SEX"))
#'
#' # A heavier cohort, in a column named the way the study names it.
#' pmx_covariates_reference(c(WEIGHTB = "WT"), medians = c(WEIGHTB = 117))
#'
#' # `WEIGHTB` resolves on its own, so the mapping is only needed to rename.
#' pmx_covariates_reference("WEIGHTB")
pmx_covariates_reference <- function(names, medians = NULL,
                                     source = paste(
                                       "synpmx reference distribution for a",
                                       "broad adult population; a package",
                                       "default, not measured from any study"
                                     )) {
  if (!is.character(names) || !length(names) || anyNA(names) ||
      any(!nzchar(names))) {
    stop("`names` must be non-empty covariate column names.", call. = FALSE)
  }
  columns <- base::names(names) %||% names
  columns[!nzchar(columns)] <- names[!nzchar(columns)]
  if (anyDuplicated(columns)) {
    stop("`names` must give each output column once.", call. = FALSE)
  }
  quantities <- unname(names)
  unknown <- !quantities %in% base::names(.covariate_reference)
  aliased <- unknown & quantities %in% base::names(.covariate_reference_aliases)
  quantities[aliased] <- .covariate_reference_aliases[quantities[aliased]]
  missing <- setdiff(quantities, base::names(.covariate_reference))
  if (length(missing)) {
    stop("No reference distribution for ", paste(missing, collapse = ", "),
         ". Available: ",
         paste(base::names(.covariate_reference), collapse = ", "),
         ". Declare anything else with `pmx_covariate()`, or map a column to ",
         "one of these with a named vector such as `c(BWT = \"WT\")`.",
         call. = FALSE)
  }
  if (!is.null(medians)) {
    if (!is.numeric(medians) || is.null(base::names(medians)) ||
        anyNA(medians) || any(!is.finite(medians))) {
      stop("`medians` must be a named numeric vector.", call. = FALSE)
    }
    unplaced <- setdiff(base::names(medians), columns)
    if (length(unplaced)) {
      stop("`medians` names a column that is not being built: ",
           paste(unplaced, collapse = ", "), ".", call. = FALSE)
    }
  }
  built <- list()
  for (i in seq_along(columns)) {
    column <- columns[[i]]
    spec <- .covariate_reference[[quantities[[i]]]]
    provenance <- paste0(source, " (", spec$quantity,
                         if (!is.null(spec$unit)) paste0(", ", spec$unit),
                         ")")
    if (is.null(spec$levels)) {
      median <- spec$median
      if (!is.null(medians) && column %in% base::names(medians)) {
        median <- unname(medians[[column]])
        # A replaced median may sit outside a range chosen around the default,
        # and `pmx_covariate()` would refuse it. Widen the bound to keep it,
        # since the bound is a safety margin rather than a claim.
        spec$range <- c(min(spec$range[1L], median / 2),
                        max(spec$range[2L], median * 2))
      }
      built[[column]] <- pmx_covariate(
        range = spec$range, source = provenance, median = median,
        cv = spec$cv, distribution = spec$distribution,
        integer = isTRUE(spec$integer)
      )
    } else {
      built[[column]] <- pmx_covariate(
        levels = spec$levels, prob = spec$prob, source = provenance
      )
    }
  }
  structure(built, class = "pmx_covariates")
}

#' What [pmx_covariates_reference()] holds
#'
#' One row per available covariate, with its units and the distribution it
#' builds. Read it before relying on a reference value: the units are the ones
#' this package writes, and a study recording height in metres rather than
#' centimetres needs its own declaration.
#'
#' @return A data frame.
#' @seealso [pmx_covariates_reference()]
#' @export
#' @examples
#' pmx_covariate_reference_table()
pmx_covariate_reference_table <- function() {
  rows <- lapply(base::names(.covariate_reference), function(name) {
    spec <- .covariate_reference[[name]]
    data.frame(
      covariate = name,
      quantity = spec$quantity,
      unit = spec$unit %||% NA_character_,
      distribution = if (is.null(spec$levels)) spec$distribution else
        "categorical",
      typical = if (is.null(spec$levels)) as.character(spec$median) else
        paste(sprintf("%s %.0f%%", spec$levels, 100 * spec$prob),
              collapse = ", "),
      cv = if (is.null(spec$levels)) sprintf("%.0f%%", 100 * spec$cv) else
        NA_character_,
      range = if (is.null(spec$levels))
        sprintf("%g to %g", spec$range[1L], spec$range[2L]) else
          NA_character_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
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
#' privacy report says so. Use it only where the source data's own access
#' controls and confidentiality obligations still apply, and never when the
#' covariate columns may reach anyone the source data could not.
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
      shape <- if (identical(cov$distribution, "range")) {
        sprintf("normal about the midpoint %g, range as 6 SD",
                mean(cov$range))
      } else {
        sprintf("%s, median %g, CV %.0f%%", cov$distribution, cov$median,
                100 * cov$cv)
      }
      cat(sprintf("  %s: continuous [%g, %g]; %s (public, DP)\n", name,
                  cov$range[1L], cov$range[2L], shape))
    } else if (cov$type == "categorical") {
      levels <- if (is.null(cov$prob)) {
        paste(cov$levels, collapse = ", ")
      } else {
        paste(sprintf("%s %.0f%%", cov$levels, 100 * cov$prob),
              collapse = ", ")
      }
      cat(sprintf("  %s: categorical {%s} (public, DP)\n", name, levels))
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
      # Calibrated mode has released a centre; prior mode has whatever the
      # declaration stated, and the range midpoint only if it stated nothing.
      centre <- summary$mean %||% cov$median %||% mean(cov$range)
      values <- if (identical(cov$distribution, "lognormal")) {
        # What a release carries is an arithmetic mean of clipped values, and
        # what a lognormal is centred on is its median. Dividing by
        # sqrt(1 + cv^2) converts the one to the other, so the generated
        # column's mean is the quantity that was actually measured. A declared
        # median needs no conversion because it is already one.
        median <- if (is.null(summary$mean)) centre else
          centre / sqrt(1 + cov$cv^2)
        median * stats::rlnorm(n, 0, sqrt(log(1 + cov$cv^2)))
      } else if (identical(cov$distribution, "normal")) {
        stats::rnorm(n, centre, centre * cov$cv)
      } else {
        # Nothing stated but the bound, so the bound doubles as the spread.
        stats::rnorm(n, centre, diff(cov$range) / 6)
      }
      values <- .clip(values, cov$range)
      out[[name]] <- if (isTRUE(cov$integer)) round(values) else values
    } else {
      prob <- summary$prob %||% cov$prob %||%
        rep(1 / length(cov$levels), length(cov$levels))
      out[[name]] <- sample(cov$levels, n, replace = TRUE, prob = prob)
    }
  }
  as.data.frame(out, stringsAsFactors = FALSE)
}
