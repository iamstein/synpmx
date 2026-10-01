# What leaves the study, and how much any one patient moves it --------------
#
# `vignettes/articles/pmxmodel-privacy.Rmd` states every protection and check
# below for a reviewer, and `pmxmodel-algorithm.Rmd` Step 5 says where each
# sits in the algorithm. No formal privacy guarantee is made by anything here.

# The release ----------------------------------------------------------------
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

.release_top <- c("structural", "parameters", "pk_models", "pd", "arms",
                  "dosing", "visits", "cells", "discrete", "covariates",
                  "covariate_effects", "schema", "roles", "endpoints",
                  "quantification_floor", "n_source", "settings", "privacy")
.release_pk <- c("endpoint", "structural", "parameters", "effects",
                 "movement")
.release_parameters <- c("fixed", "omega", "residual")
.release_movement <- c("moved", "omega_moved", "tolerance")
.release_pd <- c("pd", "typical", "baseline_cv", "residual", "arms")
.release_dosing <- c("planned", "levels", "reduction", "interruption",
                     "discontinuation", "patients")
.release_endpoints <- c("pk", "pd", "discrete", "decided_by")
.release_settings <- c("min_arm_patients", "min_category_patients")

#' What leaves the study: the release of a fitted model
#'
#' The part of a fitted model that synthetic data is generated from, and
#' nothing else. [synpmx_model_estimate()] returns a fit that also carries
#' diagnostics for the person who ran it: the candidate comparison and its AIC
#' values, the starting values, the correlations between covariates and the
#' individual random effects, the design notes and timings. None of those is
#' needed to generate, so none of them is in the release.
#' [synpmx_model_generate()] attaches the release, never the fit, to the data
#' it returns.
#'
#' Save the release, not the fit, when a fingerprint has to be carried out of
#' the environment that holds the study: [synpmx_model_generate()] accepts
#' either and produces the same data from both.
#'
#' @param fitted_model A `pmx_fitted_model` from [synpmx_model_estimate()], or
#'   a release, which is returned unchanged.
#' @return A `pmx_model_release`, which is also a `pmx_fitted_model`.
#' @seealso [model_privacy_checks()], [synpmx_model_generate()],
#'   [model_report()].
#' @export
model_release <- function(fitted_model) {
  if (!inherits(fitted_model, "pmx_fitted_model")) {
    stop("`fitted_model` must come from `synpmx_model_estimate()`.",
         call. = FALSE)
  }
  if (inherits(fitted_model, "pmx_model_release")) return(fitted_model)
  x <- unclass(fitted_model)
  keep <- function(object, fields) object[intersect(fields, names(object))]
  out <- keep(x, .release_top)
  if (!is.null(out$parameters)) {
    out$parameters <- keep(out$parameters, .release_parameters)
  }
  out$pk_models <- lapply(out$pk_models, function(model) {
    model <- keep(model, .release_pk)
    model$parameters <- keep(model$parameters, .release_parameters)
    if (!is.null(model$movement)) {
      model$movement <- keep(model$movement, .release_movement)
    }
    model
  })
  shape <- function(one) {
    one <- keep(one, .release_pd)
    if (length(one$arms)) one$arms <- lapply(one$arms, shape)
    one
  }
  out$pd <- lapply(out$pd, shape)
  out$dosing <- lapply(out$dosing, function(entry) {
    if (!is.null(entry$planned)) return(keep(entry, .release_dosing))
    lapply(entry, keep, fields = .release_dosing)
  })
  out$endpoints <- keep(out$endpoints, .release_endpoints)
  out$settings <- keep(out$settings, .release_settings)
  out$privacy <- .release_privacy(x$privacy)
  structure(out, class = c("pmx_model_release", "pmx_fitted_model"))
}

#' @export
print.pmx_model_release <- function(x, ...) {
  cat("A model release, from model_release()\n")
  cat("Every number synthetic data is generated from, and nothing else.\n\n")
  print(model_report(x))
  invisible(x)
}

# Paths in a release that the allowlist does not name, which should be none.
# Read by the first privacy check, so that a field added to the release by
# hand, or a release built by something other than `model_release()`, is
# caught rather than trusted.
.release_strays <- function(release) {
  strays <- setdiff(names(release), .release_top)
  for (name in names(release$pk_models)) {
    model <- release$pk_models[[name]]
    strays <- c(strays,
                paste0("pk_models$", name, "$", setdiff(names(model),
                                                        .release_pk)),
                paste0("pk_models$", name, "$parameters$",
                       setdiff(names(model$parameters), .release_parameters)))
  }
  for (name in names(release$pd)) {
    strays <- c(strays, paste0("pd$", name, "$",
                               setdiff(names(release$pd[[name]]),
                                       .release_pd)))
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
  unique(c(strays, walk(unclass(release), "release")))
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
# patient left out.
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

.influence_verdict <- function(change) {
  ifelse(!is.finite(change), "not computable",
         ifelse(change >= .influence_thresholds[["fail"]], "FAIL",
                ifelse(change >= .influence_thresholds[["review"]], "review",
                       "pass")))
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
.pd_influence <- function(baselines) {
  b <- log(baselines[is.finite(baselines) & baselines > 0])
  n <- length(b)
  if (n < 3L) return(NULL)
  spread <- stats::sd(b)
  without <- vapply(seq_len(n), function(i) stats::sd(b[-i]), numeric(1))
  c(typical = 100 * (exp(max(abs(b - mean(b))) / (n - 1)) - 1),
    spread = 100 * max(abs(spread - without)), patients = n)
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
    x <- x[is.finite(x)]
    n <- length(x)
    if (n < 3L) next
    scale <- if (identical(spec$kind, "lognormal")) log(x) else x
    whole <- .trimmed_moments(scale)
    loo <- vapply(seq_len(n), function(i) .trimmed_moments(scale[-i]),
                  numeric(2))
    if (identical(spec$kind, "lognormal")) {
      rows[[length(rows) + 1L]] <- .influence_row(
        "covariate", column, "geometric mean", exp(spec$meanlog),
        100 * (exp(max(abs(loo["mean", ] - whole[["mean"]]))) - 1), "%", n)
      rows[[length(rows) + 1L]] <- .influence_row(
        "covariate", column, "SD on the log scale", spec$sdlog,
        100 * max(abs(loo["sd", ] - whole[["sd"]])), "points", n)
    } else {
      sd <- whole[["sd"]]
      rows[[length(rows) + 1L]] <- .influence_row(
        "covariate", column, "mean", spec$mean,
        if (sd > 0) 100 * max(abs(loo["mean", ] - whole[["mean"]])) / sd else 0,
        "% of the SD", n)
      rows[[length(rows) + 1L]] <- .influence_row(
        "covariate", column, "SD", spec$sd,
        if (sd > 0) 100 * max(abs(loo["sd", ] - sd)) / sd else 0,
        "% of the SD", n)
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
      rows[[length(rows) + 1L]] <- .bind_influence(list(
        .influence_row("PD", label, "typical baseline",
                       unname(shape$typical[["baseline"]]),
                       reading[["typical"]], "%", reading[["patients"]]),
        .influence_row("PD", label, "between-subject SD",
                       shape$baseline_cv %||% 0, reading[["spread"]],
                       "points", reading[["patients"]])))
    }
  }
  rows[[length(rows) + 1L]] <- .covariate_influence(source, roles, covariates)
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(NULL)
  out <- .bind_influence(rows)
  who <- attr(out, "who")
  out$verdict <- .influence_verdict(out$change)
  rownames(out) <- NULL
  attr(out, "who") <- who
  out
}

# The worst row, in a sentence.
.influence_summary <- function(influence) {
  if (is.null(influence) || !nrow(influence)) {
    return(list(verdict = "not applicable",
                result = "no released estimate to read"))
  }
  finite <- influence[is.finite(influence$change), , drop = FALSE]
  if (!nrow(finite)) {
    return(list(verdict = "not applicable",
                result = "no released estimate to read"))
  }
  worst <- finite[which.max(finite$change), ]
  list(verdict = .influence_verdict(worst$change),
       result = sprintf("largest: %s %s, %.3g %s", worst$name, worst$quantity,
                        worst$change, worst$unit))
}

# What the release keeps of the influence reading: the verdict, and a number
# only where the verdict is not a pass. A passing release says that every
# estimate moved less than the threshold and nothing more, so the one
# data-dependent thing it carries is the verdict itself. A release that does
# not pass names its worst estimate, because it is not meant to leave as it is.
.release_privacy <- function(privacy) {
  if (is.null(privacy)) return(NULL)
  summary <- .influence_summary(privacy$influence)
  list(influence = list(
    verdict = summary$verdict,
    result = if (identical(summary$verdict, "pass")) {
      sprintf("every released estimate moves less than %g",
              .influence_thresholds[["review"]])
    } else summary$result),
    thresholds = .influence_thresholds)
}

# The checks -----------------------------------------------------------------

#' Privacy checks on a fitted model
#'
#' Five checks on what a fitted model releases, each with its pass criterion.
#' The first four ask whether the release holds anything about one patient
#' that it should not: a per-patient table, an identifier, a single patient's
#' value, or a frequency resting on too few patients. The fifth asks how far
#' any one patient moves a released estimate, from the individual random
#' effects for the population model, from each subject's baseline for the PD
#' shapes, and by leaving each patient out for the covariate summaries.
#'
#' None of this is a formal privacy guarantee. The fifth check describes this
#' study only, and an adversary who knows every other patient in it can detect
#' a move smaller than any threshold. The article
#' [Privacy protections in the PMX model generator](https://iamstein.github.io/synpmx/articles/pmxmodel-privacy.html)
#' states what each check establishes and what it does not.
#'
#' @param fitted_model A `pmx_fitted_model` from [synpmx_model_estimate()], or
#'   a release from [model_release()]. A release carries only the verdict of
#'   the fifth check; the full fit carries the table behind it.
#' @return A `pmx_privacy_checks` data frame with columns `check`, `question`,
#'   `result`, `criterion` and `verdict`. On a full fit, the per-estimate
#'   reading behind check P5 is the `influence` attribute: one row per released
#'   estimate, with the released value, the largest move one patient causes,
#'   its unit and its verdict. No row is about a patient.
#' @seealso [model_release()], [synpmx_model_estimate()], [synpmx_scorecard()].
#' @export
model_privacy_checks <- function(fitted_model) {
  if (!inherits(fitted_model, "pmx_fitted_model")) {
    stop("`fitted_model` must come from `synpmx_model_estimate()`.",
         call. = FALSE)
  }
  release <- model_release(fitted_model)
  k_arm <- release$settings$min_arm_patients %||% 3L
  k_level <- release$settings$min_category_patients %||% 3L

  strays <- .release_strays(release)
  p1 <- list(
    result = if (!length(strays)) "none" else paste(strays, collapse = ", "),
    ok = !length(strays))

  schema <- release$schema
  id_column <- release$roles$id
  id_levels <- length(levels(schema$prototypes[[id_column]]))
  offset <- schema$id_offset %||% 0
  round_offset <- offset == 0 || abs(log10(offset) - round(log10(offset))) < 1e-9
  identifiers <- c(
    if (!is.null(schema$id_levels)) "source ID levels in the schema",
    if (id_levels) paste(id_levels, "source ID labels on the ID column"),
    if (!round_offset) "the largest source ID as the synthetic ID offset")
  p2 <- list(result = if (!length(identifiers)) "none" else
    paste(identifiers, collapse = "; "), ok = !length(identifiers))

  floors <- unlist(release$quantification_floor)
  on_series <- vapply(floors, function(value) {
    isTRUE(abs(value - .round_down_125(value)) <= 1e-9 * value)
  }, logical(1))
  medians <- names(release$covariates)[vapply(release$covariates, function(s) {
    !is.null(s$median)
  }, logical(1))]
  references <- as.numeric(unlist(lapply(
    c(release$covariate_effects,
      unlist(lapply(release$pk_models, function(m) m$effects),
             recursive = FALSE)),
    function(effect) effect$reference)))
  unrounded <- references[abs(references - signif(references, 1L)) >
                            1e-9 * abs(references)]
  values <- c(
    if (!all(on_series)) paste("floor not on the 1-2-5 series:",
                               paste(names(floors)[!on_series], collapse = ", ")),
    if (length(medians)) paste("covariate median stored:",
                               paste(medians, collapse = ", ")),
    if (length(unrounded)) "allometric reference not rounded")
  p3 <- list(result = if (!length(values)) {
    paste0("none: ", if (length(floors)) "floors on the 1-2-5 series, " else "",
           "covariates summarized without their extremes, ",
           "no median, minimum or maximum stored")
  } else paste(values, collapse = "; "), ok = !length(values))

  sizes <- as.integer(release$arms$sizes)
  p4 <- list(
    result = sprintf("smallest arm %d; visit slots, dose levels and %s rest on at least %d",
                     min(sizes), "categorical levels", min(k_arm, k_level)),
    ok = min(sizes) >= k_arm)

  influence <- fitted_model$privacy$influence
  if (is.data.frame(influence)) {
    influence$verdict <- .influence_verdict(influence$change)
  }
  summary <- if (inherits(fitted_model, "pmx_model_release")) {
    release$privacy$influence %||%
      list(verdict = "not applicable", result = "no influence record")
  } else if (is.null(fitted_model$privacy)) {
    list(verdict = "not applicable",
         result = "no influence record: re-estimate with this version")
  } else .influence_summary(influence)

  verdict <- function(ok) if (isTRUE(ok)) "pass" else "FAIL"
  out <- data.frame(
    check = c("P1", "P2", "P3", "P4", "P5"),
    question = c(
      "No per-patient table is released",
      "No source identifier is released",
      "No single patient's value is released",
      "Every released frequency rests on several patients",
      "No single patient moves a released estimate far"),
    result = c(p1$result, p2$result, p3$result, p4$result, summary$result),
    criterion = c(
      "the release holds only the fields generation reads, none of them per patient",
      "no source ID label, and no source ID value as the synthetic ID offset",
      "no order statistic: floors coarsened, no median, minimum or maximum",
      sprintf("every arm at least %d patients", k_arm),
      sprintf("pass under %g, review under %g, FAIL from %g (%% or points)",
              .influence_thresholds[["review"]],
              .influence_thresholds[["fail"]],
              .influence_thresholds[["fail"]])),
    verdict = c(verdict(p1$ok), verdict(p2$ok), verdict(p3$ok),
                verdict(p4$ok), summary$verdict),
    stringsAsFactors = FALSE)
  attr(out, "influence") <- if (is.data.frame(influence)) influence else NULL
  class(out) <- c("pmx_privacy_checks", "data.frame")
  out
}

#' @export
print.pmx_privacy_checks <- function(x, ...) {
  plain <- as.data.frame(x)
  plain <- plain[, intersect(c("check", "question", "result", "verdict"),
                             names(plain)), drop = FALSE]
  attr(plain, "influence") <- NULL
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
  invisible(x)
}

# Rounding -------------------------------------------------------------------
#
# pmxmodel-algorithm.Rmd, Step 5 (REV-063). Every released estimate is carried
# to two significant figures. Leaving out a typical patient moves a typical
# value by less than one rounding step, so for most patients the released
# number is the same with or without them; an outlier moves a variance by more
# than a step, which is what the influence reading above is for. Applied after
# the acceptance checks, which run on the fit as estimated.
.round_parameters <- function(parameters) {
  if (is.null(parameters)) return(parameters)
  parameters$fixed <- signif(parameters$fixed, 2L)
  parameters$omega[] <- signif(parameters$omega, 2L)
  for (field in c("cv", "sd")) {
    if (!is.null(parameters$residual[[field]])) {
      parameters$residual[[field]] <- signif(parameters$residual[[field]], 2L)
    }
  }
  parameters
}

.round_shape <- function(shape) {
  if (is.null(shape)) return(shape)
  shape$typical <- signif(shape$typical, 2L)
  if (!is.null(shape$baseline_cv)) {
    shape$baseline_cv <- signif(shape$baseline_cv, 2L)
  }
  if (!is.null(shape$residual$sd)) {
    shape$residual$sd <- signif(shape$residual$sd, 2L)
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

# Saying so at estimation. A failure is a warning, because a release whose
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
    worst <- which.max(ifelse(is.finite(influence$change), influence$change,
                              -Inf))
    who[[worst]]
  } else NA_character_
  if (identical(summary$verdict, "FAIL")) {
    warning(.condition_text(
      "One patient moves a released estimate by ",
      .influence_thresholds[["fail"]], " or more (", summary$result,
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
            .influence_thresholds[["review"]], " or more (", summary$result,
            "); see `model_privacy_checks()`.")
  }
  invisible(summary)
}
