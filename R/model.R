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
                              n_source) {
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
    n_source = n_source
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
