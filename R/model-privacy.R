# What leaves the study, and how much any one patient moves it --------------
#
# `vignettes/articles/pmxmodel-privacy.Rmd` states every protection and check
# below for a reviewer, and `pmxmodel-algorithm.Rmd` Step 5 says where each
# sits in the algorithm. No formal privacy guarantee is made by anything here.

# The fingerprint ------------------------------------------------------------
#
# pmxmodel-algorithm.Rmd, Step 5 (REV-057). What leaves the environment that
# holds the study is built from an allowlist: the fields `.model_generate()`
# reads, and inside each of them only the parts it reads. Everything else on a
# fitted model is a diagnostic for the person who ran the fit -- the AIC of
# every candidate, the non-compartmental starting values, the covariate
# correlations with the individual random effects, the design notes, the
# timings -- and stays on the fit `synpmx_model_estimate()` returned.
#
# An allowlist rather than a list of things to remove, because a field added to
# the fitted model later is then left behind by default. Removing things is how
# the individual random effects reached synthetic data (REV-056): they were
# removed from one copy and travelled in another.

.fingerprint_top <- c("structural", "parameters", "pk_models", "pd", "arms",
                  "dosing", "visits", "cells", "discrete", "covariates",
                  "covariate_effects", "schema", "roles", "endpoints",
                  "quantification_floor", "n_source", "settings", "privacy")
.fingerprint_pk <- c("endpoint", "structural", "parameters", "effects",
                 "movement")
.fingerprint_parameters <- c("fixed", "omega", "residual")
.fingerprint_movement <- c("moved", "omega_moved", "tolerance")
.fingerprint_pd <- c("pd", "typical", "baseline_cv", "residual", "arms")
.fingerprint_dosing <- c("planned", "levels", "reduction", "interruption",
                     "discontinuation", "patients")
.fingerprint_endpoints <- c("pk", "pd", "discrete", "decided_by")
.fingerprint_settings <- c("min_arm_patients", "min_category_patients")

#' What leaves the study: the fingerprint of a fitted model
#'
#' The part of a fitted model that synthetic data is generated from, and
#' nothing else. [synpmx_model_estimate()] returns a fit that also carries
#' diagnostics for the person who ran it: the candidate comparison and its AIC
#' values, the starting values, the correlations between covariates and the
#' individual random effects, the design notes and timings. None of those is
#' needed to generate, so none of them is in the fingerprint.
#' [synpmx_model_generate()] reduces whatever it is given to its fingerprint
#' before generating, and attaches the fingerprint, never the fit, to the data
#' it returns.
#'
#' Save the fingerprint, not the fit, when it has to be carried out of the
#' environment that holds the study: [synpmx_model_generate()] accepts either
#' and produces the same data from both. `model_release()` is another name for
#' the same function.
#'
#' @param fitted_model A `pmx_fitted_model` from [synpmx_model_estimate()], or
#'   a fingerprint, which is returned unchanged.
#' @return A `pmx_model_fingerprint`, which is also a `pmx_fitted_model`.
#' @seealso [model_privacy_checks()], [synpmx_model_generate()],
#'   [model_report()].
#' @export
model_fingerprint <- function(fitted_model) {
  if (!inherits(fitted_model, "pmx_fitted_model")) {
    stop("`fitted_model` must come from `synpmx_model_estimate()`.",
         call. = FALSE)
  }
  if (inherits(fitted_model, "pmx_model_fingerprint")) return(fitted_model)
  x <- unclass(fitted_model)
  keep <- function(object, fields) object[intersect(fields, names(object))]
  out <- keep(x, .fingerprint_top)
  if (!is.null(out$parameters)) {
    out$parameters <- keep(out$parameters, .fingerprint_parameters)
  }
  out$pk_models <- lapply(out$pk_models, function(model) {
    model <- keep(model, .fingerprint_pk)
    model$parameters <- keep(model$parameters, .fingerprint_parameters)
    if (!is.null(model$movement)) {
      model$movement <- keep(model$movement, .fingerprint_movement)
    }
    model
  })
  shape <- function(one) {
    one <- keep(one, .fingerprint_pd)
    if (length(one$arms)) one$arms <- lapply(one$arms, shape)
    one
  }
  out$pd <- lapply(out$pd, shape)
  out$dosing <- lapply(out$dosing, function(entry) {
    if (!is.null(entry$planned)) return(keep(entry, .fingerprint_dosing))
    lapply(entry, keep, fields = .fingerprint_dosing)
  })
  out$endpoints <- keep(out$endpoints, .fingerprint_endpoints)
  out$settings <- keep(out$settings, .fingerprint_settings)
  out$privacy <- .fingerprint_privacy(x$privacy)
  structure(out, class = c("pmx_model_fingerprint", "pmx_fitted_model"))
}

#' @rdname model_fingerprint
#' @export
model_release <- model_fingerprint

#' @export
print.pmx_model_fingerprint <- function(x, ...) {
  cat("A model fingerprint, from model_fingerprint()\n")
  cat("Every number synthetic data is generated from, and nothing else.\n\n")
  print(model_report(x))
  invisible(x)
}

# Paths in a fingerprint that the allowlist does not name, which should be none.
# Read by the first privacy check, so that a field added to the fingerprint by
# hand, or a fingerprint built by something other than `model_fingerprint()`, is
# caught rather than trusted.
.fingerprint_strays <- function(fingerprint) {
  strays <- setdiff(names(fingerprint), .fingerprint_top)
  for (name in names(fingerprint$pk_models)) {
    model <- fingerprint$pk_models[[name]]
    strays <- c(strays,
                paste0("pk_models$", name, "$", setdiff(names(model),
                                                        .fingerprint_pk)),
                paste0("pk_models$", name, "$parameters$",
                       setdiff(names(model$parameters), .fingerprint_parameters)))
  }
  for (name in names(fingerprint$pd)) {
    strays <- c(strays, paste0("pd$", name, "$",
                               setdiff(names(fingerprint$pd[[name]]),
                                       .fingerprint_pd)))
  }
  strays <- strays[!grepl("\\$$", strays)]
  # The individual random effects anywhere at all, under any name a fit has
  # used for them.
  walk <- function(object, path) {
    if (!is.list(object) || is.data.frame(object)) return(character())
    found <- character()
    for (name in names(object)) {
      here <- paste0(path, "$", name)
      if (grepl("^eta", name)) found <- c(found, here)
      found <- c(found, walk(object[[name]], here))
    }
    found
  }
  unique(c(strays, walk(unclass(fingerprint), "fingerprint")))
}

# Single-patient influence -----------------------------------------------------
#
# pmxmodel-algorithm.Rmd, Step 5 (SIM-056). How far each released estimate
# would move if one patient were left out of the study: the largest such move,
# over every patient. Differential privacy asks the same question over every
# possible dataset and calls the answer sensitivity; asked of the one study at
# hand it is local sensitivity, and in regression the same idea is Cook's
# distance. It is a reading, not a guarantee: it describes this study only.
#
# Nothing is refitted. A population fit's typical value on the log scale sits
# close to the mean of the patients' individual log parameters, so leaving out
# patient i moves it by about `eta_i / (n - 1)`; and the between-subject
# variance sits close to the mean of the squared random effects, so leaving out
# patient i scales it by `((S - eta_i^2) / (n - 1)) / (S / n)`, where `S` is
# the sum of squares. Shrinkage pulls every random effect toward zero, which
# understates both, so the shrinkage is reported beside them. The PD baselines
# and the covariate summaries need no approximation: each is recomputed with the
# patient left out. A PD shape's other parameters are refitted with the patient
# left out, and judged on standard errors as well as percent
# (`.pd_shape_influence()`).
#
# The thresholds are a move of 15 (a review) and 30 (a failure), in percent for
# a typical value and in points of the between-subject SD on the log scale,
# which is close to a CV in percent, for a spread. On a parameter whose
# between-subject SD is 0.27, a patient 3 SD out in a cohort of 32 moves the
# typical value by about 3% and the spread by about 4 points, and on the public
# studies the package ships fits for, the largest move is under 15 everywhere
# except one PD endpoint's spread. A patient given a thousand times the recorded
# dose moves spreads past 50. The verdicts are recomputed whenever the reading
# is read, so a stored fit never carries a threshold the code has moved past.
.influence_thresholds <- c(review = 15, fail = 30)
# A PD shape parameter's move in standard errors, where that is the scale that
# governs it (`.pd_shape_influence()`).
.influence_se_thresholds <- c(review = 1, fail = 2)

# A move as a multiple of its unit's review threshold, so rows in different
# units can be ranked and judged together: 1 is a review and 2 a failure on
# either scale.
.influence_severity <- function(change, unit = "%") {
  review <- ifelse(unit == "SE", .influence_se_thresholds[["review"]],
                   .influence_thresholds[["review"]])
  change / review
}

.influence_verdict <- function(change, unit = "%") {
  severity <- .influence_severity(change, unit)
  ifelse(!is.finite(severity), "not computable",
         ifelse(severity >= 2, "FAIL",
                ifelse(severity >= 1, "review", "pass")))
}

# The thresholds in words, for messages and the check's criterion.
.influence_threshold_text <- function(level) {
  sprintf("%g (%g standard error%s for a PD shape parameter)",
          .influence_thresholds[[level]], .influence_se_thresholds[[level]],
          if (.influence_se_thresholds[[level]] == 1) "" else "s")
}

.influence_row <- function(group, name, quantity, released, change, unit,
                           n, shrinkage = NA_real_, who = NA_character_) {
  out <- data.frame(group = group, name = name, quantity = quantity,
                    released = released, change = change, unit = unit,
                    patients = as.integer(n), shrinkage = shrinkage,
                    stringsAsFactors = FALSE)
  attr(out, "who") <- who
  out
}

# `rbind()` drops attributes, so the patient behind each row is carried beside
# the rows and put back on the bound table. It is read once, by the warning at
# estimation, and removed before the table is stored.
.bind_influence <- function(rows) {
  who <- unlist(lapply(rows, function(r) attr(r, "who") %||%
                         rep(NA_character_, nrow(r))), use.names = FALSE)
  out <- do.call(rbind, rows)
  attr(out, "who") <- who
  out
}

# One PK model, from its individual random effects.
.pk_influence <- function(endpoint, etas, omega, fixed) {
  if (is.null(etas) || !nrow(etas)) return(NULL)
  rows <- list()
  for (parameter in rownames(omega)) {
    column <- paste0("eta.", parameter)
    if (!column %in% names(etas)) next
    eta <- etas[[column]]
    eta <- eta[is.finite(eta)]
    n <- length(eta)
    if (n < 3L) next
    variance <- omega[parameter, parameter]
    ids <- rownames(etas)[is.finite(etas[[column]])] %||%
      rep(NA_character_, n)
    pull <- abs(eta - mean(eta))
    typical <- 100 * (exp(max(pull) / (n - 1)) - 1)
    total <- sum(eta^2)
    moves <- if (total > 0 && is.finite(variance) && variance > 0) {
      ratio <- ((total - eta^2) / (n - 1)) / (total / n)
      100 * sqrt(variance) * abs(1 - sqrt(pmax(ratio, 0)))
    } else rep(0, n)
    spread <- max(moves)
    shrinkage <- if (is.finite(variance) && variance > 0) {
      1 - stats::sd(eta) / sqrt(variance)
    } else NA_real_
    label <- paste0(endpoint, ": ", parameter)
    rows[[length(rows) + 1L]] <- .influence_row(
      "PK", label, "typical value", unname(fixed[[parameter]]), typical, "%",
      n, shrinkage, who = ids[[which.max(pull)]])
    rows[[length(rows) + 1L]] <- .influence_row(
      "PK", label, "between-subject SD", sqrt(variance), spread, "points", n,
      shrinkage, who = ids[[which.max(moves)]])
  }
  if (!length(rows)) NULL else .bind_influence(rows)
}

# One PD shape, from each subject's own baseline: the typical baseline moves
# by about the patient's distance from the mean on the log scale over n - 1,
# and the between-subject spread is the SD of the log baselines, recomputed
# without the patient.
#
# The subjects behind the two largest moves ride along as an attribute, for the
# estimation to read and remove; the stored reading is numbers only.
.pd_influence <- function(baselines) {
  keep <- is.finite(baselines) & baselines > 0
  b <- log(baselines[keep])
  ids <- names(baselines)[keep] %||% rep(NA_character_, length(b))
  n <- length(b)
  if (n < 3L) return(NULL)
  spread <- stats::sd(b)
  without <- vapply(seq_len(n), function(i) stats::sd(b[-i]), numeric(1))
  pull <- abs(b - mean(b))
  moves <- abs(spread - without)
  out <- c(typical = 100 * (exp(max(pull) / (n - 1)) - 1),
           spread = 100 * max(moves), patients = n)
  attr(out, "who") <- c(typical = ids[[which.max(pull)]],
                        spread = ids[[which.max(moves)]])
  out
}

# The PD shape's parameters other than the baseline: a linear shape's slope, an
# exponential's plateau and rate. The shapes are least-squares curves fitted in
# milliseconds, so each subject is left out and the curve refitted, with no
# approximation.
#
# A percentage is the wrong unit for some of them: a slope near zero turns a
# negligible move into hundreds of percent. So the move is also measured in
# standard errors, and only a move large on both scales counts: at least 15% of
# the estimate and at least one standard error is a review, 30% and two a
# failure. The row reports whichever scale is the smaller multiple of its
# threshold, so it reads in percent for a well-determined estimate and in
# standard errors for one near zero.
#
# The standard error is a jackknife over subjects from the same refits. The
# fit's own standard error treats every observation as independent and a
# subject contributes many, so it is too small. And for each subject it is
# taken over the others, since an outlying subject inflates a jackknife that
# includes them and would set the yardstick they are measured against.
.pd_shape_influence <- function(shape, time, value, subject) {
  parameters <- setdiff(names(shape$typical), "baseline")
  ids <- .unique_in_order(as.character(subject))
  n <- length(ids)
  if (!length(parameters) || n < 4L) return(NULL)
  refit <- function(keep) {
    t <- time[keep]
    v <- value[keep]
    model <- switch(shape$pd,
      linear = try(stats::lm(v ~ t), silent = TRUE),
      exponential = try(stats::nls(
        v ~ plateau + (baseline - plateau) * exp(-rate * pmax(t, 0)),
        start = as.list(shape$typical)), silent = TRUE),
      NULL)
    if (is.null(model) || inherits(model, "try-error")) {
      return(stats::setNames(rep(NA_real_, length(parameters)), parameters))
    }
    coefficients <- stats::coef(model)
    if (identical(shape$pd, "linear")) names(coefficients) <- c("baseline", "slope")
    coefficients[parameters]
  }
  without <- vapply(ids, function(id) refit(as.character(subject) != id),
                    numeric(length(parameters)))
  without <- matrix(without, nrow = length(parameters),
                    dimnames = list(parameters, ids))
  rows <- lapply(parameters, function(parameter) {
    estimate <- unname(shape$typical[[parameter]])
    loo <- without[parameter, ]
    ok <- is.finite(loo)
    if (sum(ok) < 4L) return(NULL)
    loo <- loo[ok]
    move <- abs(loo - estimate)
    se <- vapply(seq_along(loo), function(i) {
      others <- loo[-i]
      m <- length(others)
      sqrt((m - 1) / m * sum((others - mean(others))^2))
    }, numeric(1))
    percent <- 100 * move / abs(estimate)
    standard <- move / se
    # The smaller multiple of its own review threshold is the one that counts.
    by_percent <- percent / .influence_thresholds[["review"]] <=
      standard / .influence_se_thresholds[["review"]]
    by_percent[is.na(by_percent)] <- TRUE
    change <- ifelse(by_percent, percent, standard)
    unit <- ifelse(by_percent, "%", "SE")
    severity <- .influence_severity(change, unit)
    worst <- which.max(ifelse(is.finite(severity), severity, -Inf))
    data.frame(parameter = parameter, released = estimate,
               change = change[[worst]], unit = unit[[worst]],
               patients = length(loo), who = names(loo)[[worst]],
               stringsAsFactors = FALSE)
  })
  rows <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(rows)) return(NULL)
  who <- rows$who
  rows$who <- NULL
  attr(rows, "who") <- who
  rows
}

# Every continuous covariate, recomputed without each patient in turn. The
# summary is the trimmed one the covariate model stores, so a patient at either
# extreme moves it little, which is the point of the trimming.
.covariate_influence <- function(source, roles, covariates) {
  baseline <- .baseline_covariates(source, roles)
  rows <- list()
  for (column in names(covariates)) {
    spec <- covariates[[column]]
    if (!spec$kind %in% c("lognormal", "normal")) next
    x <- suppressWarnings(as.numeric(baseline[[column]]))
    finite <- is.finite(x)
    ids <- names(baseline[[column]])[finite] %||% rep(NA_character_, sum(finite))
    x <- x[finite]
    n <- length(x)
    if (n < 3L) next
    scale <- if (identical(spec$kind, "lognormal")) log(x) else x
    whole <- .trimmed_moments(scale)
    loo <- vapply(seq_len(n), function(i) .trimmed_moments(scale[-i]),
                  numeric(2))
    location <- abs(loo["mean", ] - whole[["mean"]])
    spread <- abs(loo["sd", ] - whole[["sd"]])
    if (identical(spec$kind, "lognormal")) {
      rows[[length(rows) + 1L]] <- .influence_row(
        "covariate", column, "geometric mean", exp(spec$meanlog),
        100 * (exp(max(location)) - 1), "%", n,
        who = ids[[which.max(location)]])
      rows[[length(rows) + 1L]] <- .influence_row(
        "covariate", column, "SD on the log scale", spec$sdlog,
        100 * max(spread), "points", n, who = ids[[which.max(spread)]])
    } else {
      sd <- whole[["sd"]]
      rows[[length(rows) + 1L]] <- .influence_row(
        "covariate", column, "mean", spec$mean,
        if (sd > 0) 100 * max(location) / sd else 0, "% of the SD", n,
        who = ids[[which.max(location)]])
      rows[[length(rows) + 1L]] <- .influence_row(
        "covariate", column, "SD", spec$sd,
        if (sd > 0) 100 * max(spread) / sd else 0, "% of the SD", n,
        who = ids[[which.max(spread)]])
    }
  }
  if (!length(rows)) NULL else .bind_influence(rows)
}

# The whole table, one row per released estimate, with its verdict.
.model_influence <- function(pk_models, etas, pd, source, roles, covariates) {
  rows <- list()
  for (endpoint in names(pk_models)) {
    model <- pk_models[[endpoint]]
    rows[[length(rows) + 1L]] <- .pk_influence(
      endpoint, etas[[endpoint]], model$parameters$omega,
      model$parameters$fixed)
  }
  for (endpoint in names(pd)) {
    shapes <- if (length(pd[[endpoint]]$arms)) {
      stats::setNames(pd[[endpoint]]$arms,
                      paste0(endpoint, " (", .arm_label(names(pd[[endpoint]]$arms)),
                             ")"))
    } else stats::setNames(list(pd[[endpoint]]), endpoint)
    for (label in names(shapes)) {
      shape <- shapes[[label]]
      reading <- shape$influence
      if (is.null(reading)) next
      who <- attr(reading, "who") %||% c(typical = NA, spread = NA)
      rows[[length(rows) + 1L]] <- .bind_influence(list(
        .influence_row("PD", label, "typical baseline",
                       unname(shape$typical[["baseline"]]),
                       reading[["typical"]], "%", reading[["patients"]],
                       who = unname(who[["typical"]])),
        .influence_row("PD", label, "between-subject SD",
                       shape$baseline_cv %||% 0, reading[["spread"]],
                       "points", reading[["patients"]],
                       who = unname(who[["spread"]]))))
      extra <- shape$shape_influence
      if (is.data.frame(extra) && nrow(extra)) {
        who <- attr(extra, "who") %||% rep(NA_character_, nrow(extra))
        rows[[length(rows) + 1L]] <- .bind_influence(lapply(
          seq_len(nrow(extra)), function(i) {
            # The released value is the rounded one on the shape.
            .influence_row("PD", label, extra$parameter[[i]],
                           unname(shape$typical[[extra$parameter[[i]]]]),
                           extra$change[[i]],
                           extra$unit[[i]], extra$patients[[i]],
                           who = who[[i]])
          }))
      }
    }
  }
  rows[[length(rows) + 1L]] <- .covariate_influence(source, roles, covariates)
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(NULL)
  out <- .bind_influence(rows)
  who <- attr(out, "who")
  out$verdict <- .influence_verdict(out$change, out$unit)
  rownames(out) <- NULL
  attr(out, "who") <- who
  out
}

# The worst row, in a sentence.
.influence_summary <- function(influence) {
  if (!is.data.frame(influence) || !nrow(influence)) {
    return(list(verdict = "not applicable",
                result = "no released estimate to read"))
  }
  finite <- influence[is.finite(influence$change), , drop = FALSE]
  if (!nrow(finite)) {
    return(list(verdict = "not applicable",
                result = "no released estimate to read"))
  }
  worst <- finite[which.max(.influence_severity(finite$change,
                                                 finite$unit)), ]
  list(verdict = .influence_verdict(worst$change, worst$unit),
       result = sprintf("largest: %s %s, %.3g %s", worst$name, worst$quantity,
                        worst$change, worst$unit))
}

# What the fingerprint keeps of the influence reading: the verdict, and a number
# only where the verdict is not a pass. A passing fingerprint says that every
# estimate moved less than the threshold and nothing more, so the one
# data-dependent thing it carries is the verdict itself. A fingerprint that does
# not pass names its worst estimate, because it is not meant to leave as it is.
.fingerprint_privacy <- function(privacy) {
  if (is.null(privacy)) return(NULL)
  summary <- .influence_summary(privacy$influence)
  frequencies <- .frequency_summary(privacy$frequencies)
  list(influence = list(
    verdict = summary$verdict,
    result = if (identical(summary$verdict, "pass")) {
      paste("every released estimate moves less than",
            .influence_threshold_text("review"))
    } else summary$result),
    # The frequency verdict travels the same way: the counts stay on the fit.
    frequencies = list(
      verdict = frequencies$verdict,
      result = if (identical(frequencies$verdict, "pass")) {
        "every released frequency rests on at least the floor of patients"
      } else frequencies$result),
    thresholds = .influence_thresholds)
}

# The checks -----------------------------------------------------------------

#' Privacy checks on a fitted model
#'
#' Five checks on what a fitted model releases, each with its pass criterion.
#' The first three ask whether the fingerprint holds anything about one patient
#' that it should not: a per-patient table, an identifier, or a single
#' patient's value. The fourth recounts the smallest group of patients behind
#' each kind of released frequency -- arm sizes, attendance, dose-change rates,
#' categorical levels, carried values, and the levels the schema keeps for
#' factor columns and binary or ordinal endpoints -- on both sides of it, since
#' "one patient missed this visit" discloses as much as "one patient came". The
#' fifth asks how far any one patient moves a released estimate, from the
#' individual random effects for the population model, from each subject's
#' baseline for the PD shapes, and by leaving each patient out for the
#' covariate summaries; by default estimation has already left out any patient
#' who moved one by 15 or more.
#'
#' None of this is a formal privacy guarantee. The fifth check describes this
#' study only, and an adversary who knows every other patient in it can detect
#' a move smaller than any threshold. The article
#' [Privacy protections in the PMX model generator](https://iamstein.github.io/synpmx/articles/pmxmodel-privacy.html)
#' states what each check establishes and what it does not.
#'
#' @param fitted_model A `pmx_fitted_model` from [synpmx_model_estimate()], or
#'   a fingerprint from [model_fingerprint()]. A fingerprint carries only the verdict of
#'   the fifth check; the full fit carries the table behind it.
#' @return A `pmx_privacy_checks` data frame with columns `check`, `question`,
#'   `result`, `criterion` and `verdict`. On a full fit, the reading behind P5
#'   is the `influence` attribute, one row per released estimate with the
#'   released value, the largest move one patient causes, its unit and its
#'   verdict; and the recount behind P4 is the `frequencies` attribute, one row
#'   per kind of released frequency with the smallest group behind it, its floor
#'   and how many values were changed to meet it. No row is about a patient.
#' @seealso [model_fingerprint()], [synpmx_model_estimate()], [synpmx_scorecard()].
#' @export
model_privacy_checks <- function(fitted_model) {
  if (!inherits(fitted_model, "pmx_fitted_model")) {
    stop("`fitted_model` must come from `synpmx_model_estimate()`.",
         call. = FALSE)
  }
  fingerprint <- model_fingerprint(fitted_model)
  k_arm <- fingerprint$settings$min_arm_patients %||% 3L
  k_level <- fingerprint$settings$min_category_patients %||% 3L

  strays <- .fingerprint_strays(fingerprint)
  p1 <- list(
    result = if (!length(strays)) "none" else paste(strays, collapse = ", "),
    ok = !length(strays))

  schema <- fingerprint$schema
  id_column <- fingerprint$roles$id
  id_levels <- length(levels(schema$prototypes[[id_column]]))
  offset <- schema$id_offset %||% 0
  round_offset <- offset == 0 || abs(log10(offset) - round(log10(offset))) < 1e-9
  identifiers <- c(
    if (!is.null(schema$id_levels)) "source ID levels in the schema",
    if (id_levels) paste(id_levels, "source ID labels on the ID column"),
    if (!round_offset) "the largest source ID as the synthetic ID offset")
  p2 <- list(result = if (!length(identifiers)) "none" else
    paste(identifiers, collapse = "; "), ok = !length(identifiers))

  floors <- unlist(fingerprint$quantification_floor)
  on_series <- vapply(floors, function(value) {
    isTRUE(abs(value - .round_down_125(value)) <= 1e-9 * value)
  }, logical(1))
  medians <- names(fingerprint$covariates)[vapply(fingerprint$covariates, function(s) {
    !is.null(s$median)
  }, logical(1))]
  references <- as.numeric(unlist(lapply(
    c(fingerprint$covariate_effects,
      unlist(lapply(fingerprint$pk_models, function(m) m$effects),
             recursive = FALSE)),
    function(effect) effect$reference)))
  unrounded <- references[abs(references - signif(references, 1L)) >
                            1e-9 * abs(references)]
  # A factor level of a carried column that no arm carries (REV-067): the label
  # of an arm dropped for being too small, or a value one patient held.
  uncarried <- Filter(function(column) {
    prototype <- schema$prototypes[[column]]
    carried_values <- unlist(lapply(schema$arm_values, function(values) {
      as.character(values[[column]])
    }), use.names = FALSE)
    is.factor(prototype) && length(setdiff(levels(prototype), carried_values))
  }, schema$carried %||% character())
  values <- c(
    if (!all(on_series)) paste("floor not on the 1-2-5 series:",
                               paste(names(floors)[!on_series], collapse = ", ")),
    if (length(medians)) paste("covariate median stored:",
                               paste(medians, collapse = ", ")),
    if (length(unrounded)) "allometric reference not rounded",
    if (length(uncarried)) paste("factor level no arm carries:",
                                 paste(uncarried, collapse = ", ")))
  p3 <- list(result = if (!length(values)) "none" else
    paste(values, collapse = "; "), ok = !length(values))

  sizes <- as.integer(fingerprint$arms$sizes)
  frequencies <- if (inherits(fitted_model, "pmx_model_fingerprint")) {
    fingerprint$privacy$frequencies
  } else if (!is.null(fitted_model$privacy$frequencies)) {
    .frequency_summary(fitted_model$privacy$frequencies)
  } else NULL
  p4 <- if (!is.null(frequencies)) {
    list(result = frequencies$result,
         ok = if (identical(frequencies$verdict, "pass")) TRUE else
           if (identical(frequencies$verdict, "FAIL")) FALSE else NA)
  } else {
    list(result = sprintf(
      "smallest arm %d; no frequency record: re-estimate with this version",
      min(sizes)), ok = min(sizes) >= k_arm)
  }

  influence <- fitted_model$privacy$influence
  if (is.data.frame(influence)) {
    influence$verdict <- .influence_verdict(influence$change, influence$unit)
  }
  summary <- if (inherits(fitted_model, "pmx_model_fingerprint")) {
    fingerprint$privacy$influence %||%
      list(verdict = "not applicable", result = "no influence record")
  } else if (is.null(fitted_model$privacy)) {
    list(verdict = "not applicable",
         result = "no influence record: re-estimate with this version")
  } else .influence_summary(influence)
  left_out <- fitted_model$privacy$left_out$patients %||% 0L
  if (left_out > 0L) {
    summary$result <- paste0(summary$result, "; ", left_out, " left out")
  }

  verdict <- function(ok) {
    if (isTRUE(ok)) "pass" else if (isFALSE(ok)) "FAIL" else "not applicable"
  }
  out <- data.frame(
    check = c("P1", "P2", "P3", "P4", "P5"),
    question = c(
      "No per-patient table is released",
      "No source identifier is released",
      "No single patient's value is released",
      "Every released frequency rests on several patients, on both sides",
      "No single patient moves a released estimate far"),
    result = c(p1$result, p2$result, p3$result, p4$result, summary$result),
    criterion = c(
      "the fingerprint holds only the fields generation reads, none of them per patient",
      "no source ID label, and no source ID value as the synthetic ID offset",
      paste("no order statistic: floors coarsened, no median, minimum or",
            "maximum; no factor level of a carried column that no arm carries"),
      sprintf(paste("each side of every released frequency none or at",
                    "least %d patients (%d for categorical levels)"),
              k_arm, k_level),
      sprintf(paste("pass under %g, review under %g, FAIL from %g (%% or",
                    "points); %g and %g standard errors for a PD shape",
                    "parameter"),
              .influence_thresholds[["review"]],
              .influence_thresholds[["fail"]],
              .influence_thresholds[["fail"]],
              .influence_se_thresholds[["review"]],
              .influence_se_thresholds[["fail"]])),
    verdict = c(verdict(p1$ok), verdict(p2$ok), verdict(p3$ok),
                verdict(p4$ok), summary$verdict),
    stringsAsFactors = FALSE)
  attr(out, "influence") <- if (is.data.frame(influence)) influence else NULL
  attr(out, "frequencies") <- if (is.data.frame(fitted_model$privacy$frequencies))
    fitted_model$privacy$frequencies else NULL
  class(out) <- c("pmx_privacy_checks", "data.frame")
  out
}

#' @export
print.pmx_privacy_checks <- function(x, ...) {
  plain <- as.data.frame(x)
  plain <- plain[, intersect(c("check", "question", "result", "verdict"),
                             names(plain)), drop = FALSE]
  attr(plain, "influence") <- NULL
  attr(plain, "frequencies") <- NULL
  print(plain, row.names = FALSE, right = FALSE)
  influence <- attr(x, "influence")
  if (!is.null(influence) && nrow(influence)) {
    cat("\nHow far one patient moves each released estimate (P5):\n")
    shown <- influence[, c("group", "name", "quantity", "released", "change",
                           "unit", "verdict")]
    # Three significant figures, written out: one column holds a clearance
    # beside a variance a thousand times smaller, and a shared exponent would
    # make every row of it unreadable.
    figures <- function(x) vapply(x, function(v) {
      format(signif(v, 3L), scientific = FALSE, drop0trailing = TRUE,
             trim = TRUE)
    }, character(1))
    shown$released <- figures(shown$released)
    shown$change <- figures(shown$change)
    print(shown, row.names = FALSE, right = FALSE)
  }
  frequencies <- attr(x, "frequencies")
  if (!is.null(frequencies) && nrow(frequencies)) {
    cat("\nThe smallest group of patients behind each kind of released",
        "frequency (P4):\n")
    print(frequencies[, c("quantity", "smallest", "threshold", "adjusted")],
          row.names = FALSE, right = FALSE)
  }
  invisible(x)
}

# Rounding -------------------------------------------------------------------
#
# pmxmodel-algorithm.Rmd, Step 5 (REV-063, REV-074). Every released PK and PD
# estimate is carried to one significant figure, and every covariate summary to
# two. Two figures moved a parameter by 1% on average and changed nothing an
# attacker could measure; one moves it by 9% on average and 33% at worst, and
# takes a quarter to half of a membership-inference attack's edge over chance,
# because its step is as wide as an estimate's sampling variation. Covariates
# keep two, because one would put a mean height of 172 cm at 200. Applied after
# the acceptance checks, which run on the fit as estimated, and the estimates
# at full precision stay on the fit, outside the fingerprint.
.fingerprint_digits <- 1L

.round_parameters <- function(parameters) {
  if (is.null(parameters)) return(parameters)
  parameters$fixed <- signif(parameters$fixed, .fingerprint_digits)
  parameters$omega[] <- signif(parameters$omega, .fingerprint_digits)
  for (field in c("cv", "sd")) {
    if (!is.null(parameters$residual[[field]])) {
      parameters$residual[[field]] <- signif(parameters$residual[[field]],
                                             .fingerprint_digits)
    }
  }
  parameters
}

.round_shape <- function(shape) {
  if (is.null(shape)) return(shape)
  shape$typical <- signif(shape$typical, .fingerprint_digits)
  if (!is.null(shape$baseline_cv)) {
    shape$baseline_cv <- signif(shape$baseline_cv, .fingerprint_digits)
  }
  if (!is.null(shape$residual$sd)) {
    shape$residual$sd <- signif(shape$residual$sd, .fingerprint_digits)
  }
  if (length(shape$arms)) shape$arms <- lapply(shape$arms, .round_shape)
  shape
}

# A lognormal covariate is rounded on its own scale -- the geometric mean to two
# figures, 69 kg rather than 68.7 -- rather than on the log scale, where two
# figures of 4.23 would move the geometric mean by 3%.
.round_covariate <- function(spec) {
  if (identical(spec$kind, "lognormal")) {
    spec$meanlog <- log(signif(exp(spec$meanlog), 2L))
    spec$sdlog <- signif(spec$sdlog, 2L)
  } else if (identical(spec$kind, "normal")) {
    spec$mean <- signif(spec$mean, 2L)
    spec$sd <- signif(spec$sd, 2L)
  }
  spec
}

# Saying so at estimation. A failure is a warning, because a fingerprint whose
# estimates rest on one patient should not leave as it is, and a patient given
# far more drug than the record says is the commonest way to get one. A review
# is a message: a small cohort moves its estimates by more than a large one, and
# that is worth knowing rather than worth stopping for.
#
# The warning names the patient where the random effects say who it is. It is
# printed in the environment that holds the study and is stored nowhere: the
# fit keeps the reading per estimate, never per patient.
.model_announce_influence <- function(influence, quiet = FALSE) {
  summary <- .influence_summary(influence)
  who <- attr(influence, "who")
  culprit <- if (!is.null(who) && nrow(influence)) {
    severity <- .influence_severity(influence$change, influence$unit)
    worst <- which.max(ifelse(is.finite(severity), severity, -Inf))
    who[[worst]]
  } else NA_character_
  if (identical(summary$verdict, "FAIL")) {
    warning(.condition_text(
      "One patient moves a released estimate by ",
      .influence_threshold_text("fail"), " or more (", summary$result,
      if (!is.na(culprit)) paste0("; patient `", culprit, "`") else "", ").",
      why = paste("The fitted model then describes that patient as much as the",
                  "cohort. A patient who received far more or less drug than",
                  "the dosing record says produces exactly this."),
      fix = paste("Check that patient's records, or leave them out, before",
                  "sharing data generated from this fit.",
                  "`model_privacy_checks()` lists every estimate.")),
      call. = FALSE)
  } else if (identical(summary$verdict, "review") && !quiet) {
    message("One patient moves a released estimate by ",
            .influence_threshold_text("review"), " or more (", summary$result,
            "); see `model_privacy_checks()`.")
  }
  invisible(summary)
}

# The subjects behind a PD shape's influence reading leave with the reading's
# attribute; the shape the fit stores keeps the numbers only (SIM-092).
.strip_influence_who <- function(shape) {
  if (is.null(shape)) return(shape)
  if (!is.null(shape$influence)) attr(shape$influence, "who") <- NULL
  if (!is.null(shape$shape_influence)) {
    attr(shape$shape_influence, "who") <- NULL
  }
  if (length(shape$arms)) shape$arms <- lapply(shape$arms, .strip_influence_who)
  shape
}

# Released frequencies ---------------------------------------------------------
#
# pmxmodel-algorithm.Rmd, Step 5 (REV-065, REV-066). Every frequency the fingerprint
# carries is a share of some group of patients, and the threshold rule of
# statistical disclosure control asks that the patients on each side of it be
# none or at least `k`: "one patient missed this visit" says as much about a
# patient as "one patient came". This recounts, from the source, the smallest
# group behind each kind of released frequency after the rules that enforce it
# have run, and how many values those rules changed. One row per kind; nothing
# in it is about a patient.
.frequency_audit <- function(source, roles, subject_group, arm_models, cells,
                             covariates, estimate_source, discrete, schema,
                             k_arm, k_level, discrete_by_arm = FALSE) {
  smallest_of <- function(x) if (length(x)) min(x) else NA_real_
  row <- function(quantity, smallest, adjusted, threshold, rule) {
    data.frame(quantity = quantity, smallest = smallest,
               adjusted = as.integer(adjusted), threshold = threshold,
               rule = rule, stringsAsFactors = FALSE)
  }
  rows <- list(row("patients in an arm", smallest_of(arm_models$sizes), 0L,
                   k_arm, "a smaller arm is dropped before estimation"))

  attendance <- arm_models$audit$attendance
  if (!is.null(attendance)) {
    open <- attendance[!attendance$rounded & attendance$attenders > 0 &
                         attendance$misses > 0, , drop = FALSE]
    masked <- attendance$masked %||% rep(0, nrow(attendance))
    rows[[length(rows) + 1L]] <- row(
      "patients on either side of an endpoint's attendance rate",
      smallest_of(pmin(open$attenders, open$misses)),
      sum(attendance$rounded) + sum(masked), k_arm,
      paste("one rate per endpoint; fewer patients missing any of its visits",
            "leaves it at 1, and an arm with fewer at a visit does not have it"))
  }
  rates <- arm_models$audit$rates
  if (!is.null(rates)) {
    moving <- rates[rates$patients_after > 0, , drop = FALSE]
    rows[[length(rows) + 1L]] <- row(
      "patients with the dose change behind a rate",
      smallest_of(moving$patients_after), sum(rates$action == "set to zero"),
      k_arm, "pooled over arms; a rate resting on fewer is set to zero")
  }

  baseline <- .baseline_covariates(estimate_source, roles)
  level_holders <- numeric()
  excluded <- 0L
  for (column in names(covariates)) {
    values <- as.character(baseline[[column]])
    values <- values[!is.na(values)]
    if (!identical(covariates[[column]]$kind, "categorical")) {
      if (identical(covariates[[column]]$kind, "missing") && length(values) &&
          !is.numeric(baseline[[column]])) {
        excluded <- excluded + length(unique(values))
      }
      next
    }
    kept <- covariates[[column]]$levels
    level_holders <- c(level_holders, vapply(kept, function(level) {
      sum(values == level)
    }, numeric(1)))
    excluded <- excluded + length(setdiff(unique(values), kept))
  }
  rows[[length(rows) + 1L]] <- row(
    "patients holding a categorical covariate level",
    smallest_of(level_holders), excluded, k_level,
    "a rarer level is excluded")

  ids <- as.character(source[[roles$id]])
  subjects <- .unique_in_order(source[[roles$id]])
  if (length(discrete) && any(!vapply(unlist(discrete, recursive = FALSE),
                                      is.null, logical(1)))) {
    nominal <- suppressWarnings(as.numeric(source[[roles$nominal_time]]))
    planned <- source
    planned[[roles$time]] <- nominal
    aligned <- .aligned_time(planned, roles)
    observed <- .observation_rows(source, roles, require_present = TRUE)
    endpoint <- .endpoint(source, roles)
    arm_of <- stats::setNames(subject_group, as.character(subjects))
    row_arm <- arm_of[ids]
    discrete_holders <- numeric()
    folded <- 0L
    # One marginal per endpoint, read over every visit and every arm that
    # shares it, or one per arm under `pd_by_arm` (REV-071, REV-072).
    groups <- if (discrete_by_arm) as.list(names(discrete)) else
      list(names(discrete))
    for (group in groups) {
      for (name in unique(cells$endpoint)) {
        index <- which(cells$endpoint == name)
        at <- rep(FALSE, nrow(source))
        marginal <- NULL
        for (i in index) {
          sharing <- group[!vapply(group, function(arm) {
            is.null(discrete[[arm]][[i]])
          }, logical(1))]
          if (!length(sharing)) next
          marginal <- discrete[[sharing[[1L]]]][[i]]
          at <- at | (observed & row_arm %in% sharing & endpoint == name &
                        abs(aligned - cells$time[i]) < sqrt(.Machine$double.eps))
        }
        if (is.null(marginal)) next
        values <- suppressWarnings(as.numeric(source[[roles$dv]][at]))
        holders <- ids[at]
        present <- !is.na(values)
        discrete_holders <- c(discrete_holders, vapply(marginal$levels,
          function(level) {
            length(unique(holders[present & values == level]))
          }, numeric(1)))
        folded <- folded + length(setdiff(unique(values[present]),
                                          marginal$levels))
      }
    }
    rows[[length(rows) + 1L]] <- row(
      "patients holding a discrete endpoint's level",
      smallest_of(discrete_holders), folded, k_level,
      "a rarer level is folded into the endpoint's most common one")
  }

  carried <- intersect(setdiff(roles$keep, roles$strata), names(source))
  if (length(carried)) {
    first_row <- vapply(subjects, function(subject) {
      which(!is.na(source[[roles$id]]) & source[[roles$id]] == subject)[1L]
    }, integer(1))
    carried_holders <- numeric()
    removed <- 0L
    for (arm in names(schema$arm_values)) {
      members <- first_row[subject_group == arm]
      for (column in carried) {
        stored <- schema$arm_values[[arm]][[column]]
        own <- as.character(source[[column]][members])
        if (is.null(stored) || is.na(stored)) {
          if (!is.na(own[[1L]])) removed <- removed + 1L
          next
        }
        carried_holders <- c(carried_holders,
                             sum(own == as.character(stored), na.rm = TRUE))
      }
    }
    rows[[length(rows) + 1L]] <- row(
      "patients holding an arm's keep value",
      smallest_of(carried_holders), removed, k_level,
      "a value fewer hold is written as missing")
  }

  # The levels the schema keeps for a factor column (REV-067), other than the
  # subject identifier, which keeps none, a covariate, counted above, and a
  # stratum, whose levels are the arms counted first.
  factor_holders <- numeric()
  dropped <- 0L
  for (column in names(schema$prototypes)) {
    prototype <- schema$prototypes[[column]]
    if (!is.factor(prototype) || identical(column, roles$id) ||
        column %in% c(roles$covariates, roles$strata) ||
        !column %in% names(source)) next
    values <- as.character(source[[column]])
    factor_holders <- c(factor_holders, vapply(levels(prototype),
      function(level) {
        length(unique(ids[!is.na(values) & values == level]))
      }, numeric(1)))
    if (is.factor(source[[column]])) {
      dropped <- dropped + length(setdiff(levels(source[[column]]),
                                          levels(prototype)))
    }
  }
  if (length(factor_holders) || dropped) {
    rows[[length(rows) + 1L]] <- row(
      "patients holding a level of a factor column",
      smallest_of(factor_holders), dropped, k_level,
      "a rarer level, or one no arm carries, is not stored")
  }

  # The levels the schema keeps for each binary or ordinal endpoint's value
  # type (REV-068), against every level the study recorded.
  recorded <- .endpoint_value_types(source, roles)
  observed <- .observation_rows(source, roles, require_present = TRUE)
  endpoint <- .endpoint(source, roles)
  dv <- suppressWarnings(as.numeric(source[[roles$dv]]))
  endpoint_holders <- numeric()
  dropped <- 0L
  for (name in names(recorded)) {
    if (!length(recorded[[name]]$levels)) next
    stored <- schema$endpoint_specs[[name]]$levels
    at <- which(observed & endpoint == name & is.finite(dv))
    endpoint_holders <- c(endpoint_holders, vapply(stored, function(level) {
      length(unique(ids[at][abs(dv[at] - level) <= 1e-8]))
    }, numeric(1)))
    dropped <- dropped + length(recorded[[name]]$levels) - length(stored)
  }
  if (length(endpoint_holders) || dropped) {
    rows[[length(rows) + 1L]] <- row(
      "patients holding a discrete endpoint's level in the schema",
      smallest_of(endpoint_holders), dropped, k_level,
      "a rarer level is not stored")
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# The audit in a sentence, and its verdict.
.frequency_summary <- function(frequencies) {
  if (is.null(frequencies) || !nrow(frequencies)) {
    return(list(verdict = "not applicable", result = "no frequency record"))
  }
  known <- frequencies[is.finite(frequencies$smallest), , drop = FALSE]
  short <- known[known$smallest < known$threshold, , drop = FALSE]
  adjusted <- frequencies[!is.na(frequencies$adjusted) &
                            frequencies$adjusted > 0, , drop = FALSE]
  smallest <- if (nrow(known)) known[which.min(known$smallest), ] else NULL
  result <- paste0(
    if (nrow(short)) {
      paste0("below the floor: ", paste(sprintf("%s (%g)", short$quantity,
                                                short$smallest),
                                        collapse = "; "))
    } else if (!is.null(smallest)) {
      sprintf("smallest group %g", smallest$smallest)
    } else "nothing to count",
    if (nrow(adjusted)) {
      sprintf("; %d adjusted", sum(adjusted$adjusted))
    } else "")
  list(verdict = if (nrow(short)) "FAIL" else "pass", result = result)
}
