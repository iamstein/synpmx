# The PMX model generator: the fitted-model object, and the gates that guard it.
#
# The generator estimates a small set of linear PK models, picks one on AIC, and
# draws new subjects by simulating from it. Everything that is not the
# concentration-time curve -- dose reductions, skipped cycles, discontinuation,
# visit attendance, arms, covariates, censoring -- comes from the shared models
# in `R/dose-visit-models.R`, unchanged from what `synpmx_pca_summarize()` uses.
#
# Nothing here fits anything: this file is the object and the input gates, built
# before any estimation exists so that the half needing no fitter is under test
# on its own. `synpmx_model_estimate()` is deliberately absent rather than
# exported and empty; it arrives with a fit behind it, and until then these
# gates are called by name from the tests. `vignettes/model-algorithm.Rmd` will
# document the steps, and these comments will name its sections once it does.

# The gates.
#
# Two of the seven cannot be decided from the inputs alone and arrive with the
# steps that decide them: whether any endpoint is the drug concentration, and
# whether any candidate converged. The five here read the inputs.

# A covariance matrix fitted to a handful of subjects describes those subjects.
# PCA's floor is 10 and this one is higher because a parameter estimate
# concentrates on its cohort faster than a score does: the score is one
# subject's coordinates in a basis everybody shares, and the fixed effect is a
# statement about the population that has only these people in it. A threshold
# rather than an accounting; what would replace it with a number is `SIM-056`.
.model_require_subjects <- function(n_source, minimum) {
  minimum <- .positive_integer(minimum, "min_subjects")
  if (n_source < minimum) {
    stop("`synpmx_model_estimate()` needs at least ", minimum,
         " subjects to fit a population model; this source has ", n_source,
         ". A covariance matrix fitted to fewer describes those subjects ",
         "rather than a population. `synpmx_avatar()` and `synpmx_pca()` need ",
         "no identifiable structure and have lower floors.", call. = FALSE)
  }
  invisible(TRUE)
}

# The same fork `synpmx_pca()` answers the same way, and for a second reason
# here: the dosing and visit models sit on the nominal grid, and the estimation
# step reads recorded times against recorded dosing histories. Both axes are
# needed, so both have to be declared. Inferring the grid would be a statement
# about the protocol that only the caller is in a position to make.
#
# Unlike PCA this does not replace `roles$time` with the nominal column. A
# population PK fit is a statement about time after the dose that was actually
# given, so the source is returned untouched and the two axes stay separate.
.model_require_nominal_time <- function(source, roles) {
  if (is.null(roles$nominal_time)) {
    stop("`synpmx_model_estimate()` requires `nominal_time` in `pmx_roles()`. ",
         "The dosing and visit models sit on the nominal grid, and inferring ",
         "that grid from recorded times is a statement about the protocol ",
         "that only you can make. Add the protocol's planned times as a ",
         "column and declare it.", call. = FALSE)
  }
  nominal <- suppressWarnings(as.numeric(source[[roles$nominal_time]]))
  relevant <- .observation_rows(source, roles, require_present = TRUE) |
    .event_rows(source, roles)
  missing <- relevant & !is.finite(nominal)
  if (any(missing)) {
    stop("`nominal_time` is missing on ", sum(missing), " of ", sum(relevant),
         " dose and observation rows. Every row the dosing and visit models ",
         "read needs a nominal time; fill them in or drop those rows before ",
         "calling.", call. = FALSE)
  }
  invisible(TRUE)
}

# Below this the fit reports parameters that came from the starting values.
# Counted on the nominal grid rather than on recorded times, so that "bin" needs
# no width chosen here: the grid is the protocol's own sampling slots, and a
# study sampled at six distinct times after dose has six of them however
# precisely the clock recorded each visit. Recorded times would make this a
# count of how noisy the clock was.
.model_time_coverage <- function(source, roles) {
  nominal <- suppressWarnings(as.numeric(source[[roles$nominal_time]]))
  planned <- source
  planned[[roles$time]] <- nominal
  tad <- .derived_tad(planned, roles)
  observed <- .observation_rows(source, roles, require_present = TRUE) &
    is.finite(tad)
  length(unique(tad[observed]))
}

.model_require_time_coverage <- function(source, roles, minimum) {
  minimum <- .positive_integer(minimum, "min_time_bins")
  bins <- .model_time_coverage(source, roles)
  if (bins < minimum) {
    stop("`synpmx_model_estimate()` needs observations at ", minimum,
         " distinct nominal times after a dose to identify a linear model; ",
         "this source has ", bins, ". Below that the fit reports the ",
         "starting values. `synpmx_avatar()` and `synpmx_pca()` need no ",
         "identifiable structure and will run on this study.", call. = FALSE)
  }
  invisible(TRUE)
}

# `pmx_structural_model()` demands a `source` string because a structural model
# entering the differentially private path is treated as data-independent: no
# budget is charged for it. A `pmx_fitted_model` holds fixed effects estimated
# from the confidential study, so accepting one there would spend nothing for
# information taken from the data. The two classes are disjoint, which is what
# makes the existing `inherits(model, "pmx_structural_model")` checks refuse it;
# this names the reason rather than letting a generic type error stand for it.
#
# What no check reaches is a caller reading numbers off a fitted model and
# typing them into `pmx_structural_model(typical = )` by hand. That is `REV-042`
# and is disclosed rather than gated.
.reject_fitted_model <- function(x, argument, what) {
  if (inherits(x, "pmx_fitted_model")) {
    stop("`", argument, "` is a `pmx_fitted_model`, which `", what,
         "` cannot accept. Its parameters were estimated from the ",
         "confidential study, so treating them as a public input would ",
         "charge no privacy budget for information taken from the data. ",
         "Supply a `pmx_structural_model()` whose values come from a source ",
         "outside this dataset.", call. = FALSE)
  }
  invisible(TRUE)
}

# The object -----------------------------------------------------------------
#
# Two halves: what `nlmixr2` estimated, and the apparatus the shared dosing and
# visit models build, which is the same apparatus `synpmx_pca_summarize()`
# returns and means the same thing there.
#
# No empirical Bayes estimates. They are per-subject quantities, and an object
# carrying them would be a description of each real patient in the study.
# Generation draws random effects from `parameters$omega` instead.
.pmx_fitted_model <- function(structural, candidates, parameters, endpoints,
                              arms, dosing, visits, schema, roles, settings,
                              n_source, cells = NULL, pd = list(),
                              covariate_effects = list(), covariates = list(),
                              discrete = list(), design = NULL,
                              correlations = NULL, censoring = NULL) {
  if (!structural %in% .pk_models) {
    stop("`structural` must be one of: ", paste(.pk_models, collapse = ", "),
         ".", call. = FALSE)
  }
  if (!is.data.frame(candidates) ||
      !all(c("model", "converged", "aic", "note") %in% names(candidates))) {
    stop("`candidates` must be a data frame with columns model, converged, ",
         "aic and note.", call. = FALSE)
  }
  if (!structural %in% candidates$model[which(candidates$converged)]) {
    stop("`structural` names a model that is not among the converged ",
         "candidates.", call. = FALSE)
  }
  needed <- c("fixed", "omega", "residual")
  if (!is.list(parameters) || !all(needed %in% names(parameters))) {
    stop("`parameters` must hold ", paste(needed, collapse = ", "), ".",
         call. = FALSE)
  }
  required <- .required_pk_params[[structural]]
  missing_params <- setdiff(required, names(parameters$fixed))
  if (length(missing_params)) {
    stop("`parameters$fixed` is missing: ",
         paste(missing_params, collapse = ", "), ".", call. = FALSE)
  }
  if (!is.matrix(parameters$omega) ||
      nrow(parameters$omega) != ncol(parameters$omega) ||
      is.null(rownames(parameters$omega))) {
    stop("`parameters$omega` must be a named square matrix.", call. = FALSE)
  }
  if (!all(rownames(parameters$omega) %in% names(parameters$fixed))) {
    stop("Every random effect in `parameters$omega` needs a fixed effect of ",
         "the same name.", call. = FALSE)
  }
  if (!inherits(roles, "pmx_roles")) {
    stop("`roles` must come from `pmx_roles()`.", call. = FALSE)
  }
  structure(list(
    structural = structural,
    candidates = candidates,
    parameters = parameters,
    endpoints = endpoints,
    arms = arms,
    dosing = dosing,
    visits = visits,
    schema = schema,
    roles = roles,
    settings = settings,
    n_source = n_source,
    cells = cells,
    pd = pd,
    covariate_effects = covariate_effects,
    covariates = covariates,
    discrete = discrete,
    design = design,
    correlations = correlations,
    censoring = censoring
  ), class = "pmx_fitted_model")
}

#' @export
print.pmx_fitted_model <- function(x, ...) {
  cat("A fitted PMX model, from synpmx_model_estimate()\n\n")
  cat("  fitted on   ", x$n_source, "patients,",
      length(x$arms$arms), "arm(s)\n")
  cat("  structural  ", x$structural,
      sprintf("(chosen from %d candidate(s) on AIC)", nrow(x$candidates)),
      "\n")
  cat("  fixed       ",
      paste(sprintf("%s %.4g", names(x$parameters$fixed),
                    as.numeric(x$parameters$fixed)), collapse = ", "), "\n")
  cat("  random on   ", paste(rownames(x$parameters$omega), collapse = ", "),
      "\n")
  cat("  pk endpoint ", x$endpoints$pk %||% "none", "\n")
  # The out-of-scope statement prints with the object rather than living only
  # in the manual, because this object's contents look exactly like the output
  # of a real population analysis and will be read as one otherwise.
  cat("\n", .wrap_plain(paste(
    "These parameters are not estimates to report. They exist to make",
    "simulated profiles resemble the source study; the candidate set is too",
    "small and the covariate model too thin for any of them to answer a",
    "scientific question."
  ), "  ", "  "), "\n", sep = "")
  invisible(x)
}

#' What a fitted model carries
#'
#' An inventory of everything in a `pmx_fitted_model`, in two halves: what
#' `nlmixr2` estimated, and the dosing, visit and covariate models that are
#' summaries of the source rather than estimates. Nothing here is per-subject.
#'
#' @param fitted_model A `pmx_fitted_model` from [synpmx_model_estimate()].
#'
#' @return A `pmx_model_report` list, printed as sections.
#' @seealso [model_candidates()], [model_parameters()],
#'   [synpmx_model_estimate()].
#' @export
model_report <- function(fitted_model) {
  stopifnot(inherits(fitted_model, "pmx_fitted_model"))
  structure(list(
    structural = fitted_model$structural,
    n_source = fitted_model$n_source,
    endpoints = fitted_model$endpoints,
    design = fitted_model$design,
    parameters = fitted_model$parameters,
    covariate_effects = fitted_model$covariate_effects,
    correlations = fitted_model$correlations,
    censoring = fitted_model$censoring,
    pd = fitted_model$pd,
    arms = fitted_model$arms,
    dosing = fitted_model$dosing,
    cells = fitted_model$cells,
    settings = fitted_model$settings
  ), class = "pmx_model_report")
}

#' @export
print.pmx_model_report <- function(x, ...) {
  cat("What this fitted model carries\n\n")
  cat("Estimated by nlmixr2\n")
  cat("  structural model  ", x$structural, "\n")
  cat("  fixed effects     ",
      paste(sprintf("%s %.4g", names(x$parameters$fixed),
                    as.numeric(x$parameters$fixed)), collapse = ", "), "\n")
  cat("  between-subject   ",
      paste(sprintf("%s %.3g", rownames(x$parameters$omega),
                    sqrt(diag(x$parameters$omega))), collapse = ", "),
      "(as SD on the log scale)\n")
  cat("  residual error    ", x$parameters$residual$kind,
      sprintf("%.3g", x$parameters$residual$cv %||% x$parameters$residual$sd),
      "\n")
  cat("  covariate effects ",
      if (length(x$covariate_effects)) {
        paste(vapply(names(x$covariate_effects), function(parameter) {
          effect <- x$covariate_effects[[parameter]]
          sprintf("%s ~ (%s/%.4g)^%.2f", parameter, effect$covariate,
                  effect$reference, effect$exponent)
        }, character(1)), collapse = ", ")
      } else "none", "\n")
  if (length(x$pd)) {
    cat("  pd shapes         ",
        paste(sprintf("%s: %s", names(x$pd),
                      vapply(x$pd, function(s) s$pd, character(1))),
              collapse = ", "), "\n")
  }

  # The assay limit is reported with the fit rather than left in the schema,
  # because how much of an endpoint sits below it is how much of the fit is a
  # statement about the imputation rather than about measurements.
  if (!is.null(x$censoring) && any(x$censoring$imputed > 0)) {
    cat("  below the limit  ",
        paste(sprintf("%s %d of %d (%.0f%%) imputed below %.4g",
                      x$censoring$endpoint, x$censoring$imputed,
                      x$censoring$observations, 100 * x$censoring$fraction,
                      x$censoring$limit)[x$censoring$imputed > 0],
              collapse = "; "), "\n")
  }

  cat("\nSummarized from the source, not estimated\n")
  cat("  cohort            ", x$n_source, "patients in",
      length(x$arms$arms), "arm(s)\n")
  cat("  visit model       ", nrow(x$cells), "grid cells over",
      length(unique(x$cells$endpoint)), "endpoint(s)\n")
  varies <- vapply(x$dosing, function(d) {
    d$discontinuation > 0 || d$interruption > 0 || d$reduction > 0
  }, logical(1))
  cat("  dosing model      ",
      sprintf("%d planned cycle(s) per arm", stats::median(
        vapply(x$dosing, function(d) nrow(d$planned), integer(1)))),
      if (any(varies)) sprintf("| %d of %d arm(s) reduce, skip or stop early",
                               sum(varies), length(varies)) else
        "| no reductions, skips or early stops", "\n")

  cat("\nHow the concentration endpoint was decided\n")
  cat("  endpoint          ", x$endpoints$pk, sprintf("(%s)",
                                                      x$endpoints$decided_by), "\n")
  if (!is.null(x$endpoints$signals)) {
    print(x$endpoints$signals, row.names = FALSE)
  }
  if (!is.null(x$design)) {
    cat("  design            ", x$design$reason, "\n")
    # The candidate set is one-compartment. Where the sampling would support a
    # distribution phase, say so, because asking for it is the caller's move.
    if (isTRUE(x$design$richness$rich) && !grepl("^2cmt", x$structural)) {
      cat("  also available    ",
          sprintf("the sampling would support a two-compartment model (%s): ask for it with `pk = \"2cmt_%s\"`",
                  sprintf("median %g distinct times after a dose, %g after the peak",
                          x$design$richness$per_subject,
                          x$design$richness$after_peak),
                  if (grepl("oral", x$structural)) "oral" else "iv"), "\n")
    }
  }

  # The correlations an unmodelled covariate relationship shows up in. A
  # covariate that influences the real profiles and is not in the model is
  # generated independently of them, and this is the only place that says so.
  if (!is.null(x$correlations) && nrow(x$correlations)) {
    strongest <- x$correlations[order(-abs(x$correlations$correlation)), ,
                                drop = FALSE]
    cat("\nCovariate against the individual random effects\n")
    print(utils::head(strongest, 5L), row.names = FALSE, digits = 2)
    cat("\n", .wrap_plain(paste(
      "A covariate that moves with a random effect and is not in the model",
      "above is generated independently of the profiles, so the synthetic",
      "data carries no relationship between them. `synpmx_avatar()` keeps",
      "those relationships without modelling them."
    ), "  ", "  "), "\n", sep = "")
  }
  invisible(x)
}

#' The candidate models the selection was made from
#'
#' Every candidate the design admitted, whether or not it converged, with the
#' AIC it was compared on. A candidate that failed keeps its reason, so a search
#' that came down to one survivor does not look like a search that had one
#' candidate.
#'
#' @param fitted_model A `pmx_fitted_model` from [synpmx_model_estimate()].
#' @return A data frame with columns `model`, `converged`, `aic` and `note`.
#' @seealso [model_report()], [model_parameters()].
#' @export
model_candidates <- function(fitted_model) {
  stopifnot(inherits(fitted_model, "pmx_fitted_model"))
  fitted_model$candidates
}

#' The estimated parameters
#'
#' Fixed effects, the between-subject covariance matrix and the residual error.
#' Not estimates to report: see [synpmx_model_estimate()].
#'
#' @param fitted_model A `pmx_fitted_model` from [synpmx_model_estimate()].
#' @return A list with `fixed`, `omega` and `residual`.
#' @seealso [model_report()], [model_candidates()].
#' @export
model_parameters <- function(fitted_model) {
  stopifnot(inherits(fitted_model, "pmx_fitted_model"))
  fitted_model$parameters[c("fixed", "omega", "residual")]
}
