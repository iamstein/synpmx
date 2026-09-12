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
#
# It warns rather than refuses. Whether a covariance matrix estimated from
# eighteen people is fit for the purpose at hand is a judgement about the
# purpose, which the count cannot make on the caller's behalf -- and the study
# small enough to trip this is often the one most in need of a synthetic copy.
# What the warning has to carry is the part nothing downstream will say: a
# scorecard reads whether the output copies anybody or changed the study's
# shape, and a small-cohort fit does neither.
# The loudest thing this package can say at estimation time, because it is the
# one failure that otherwise looks exactly like success: an AIC, a parameter
# table, a generated dataset, and a clean scorecard, all resting on numbers no
# optimizer ever touched. See `.model_fit_movement()`.
.model_warn_unmoved <- function(structural, movement) {
  changes <- movement$changes
  warning(.condition_text(
    "THE FIT DID NOT MOVE. `", structural, "` returned every parameter within ",
    sprintf("%.0f%%", 100 * movement$tolerance),
    " of its starting value, so the numbers below are starting values rather ",
    "than estimates:",
    items = sprintf("%-10s %-12s -> %-12s (%.2f%%)", changes$parameter,
                    sprintf("%.4g", changes$start),
                    sprintf("%.4g", changes$estimate),
                    100 * changes$relative_change),
    why = paste("`nlmixr2` still reports an objective and an AIC, and nothing",
                "downstream will contradict them: the generator simulates from",
                "these values and the scorecard passes, because it asks",
                "whether the output copies anybody or changed the study's",
                "shape and a fit that never moved does neither."),
    fix = paste("Do not generate from this fit. Check the starting values in",
                "`model_report()`, and prefer `synpmx_avatar()` or",
                "`synpmx_pca()`, which fit no structural model.")),
    call. = FALSE)
  invisible(TRUE)
}

.model_note_subjects <- function(n_source, minimum) {
  minimum <- .positive_integer(minimum, "min_subjects")
  if (n_source < minimum) {
    warning(.condition_text(
      "`synpmx_model_estimate()` is fitting a population model to ", n_source,
      " subjects, below `min_subjects` = ", minimum, ".",
      why = paste0("The fixed effects and the covariance matrix describe those ",
                   n_source, " subjects rather than a population, and nothing ",
                   "downstream will tell you so: the scorecard asks whether ",
                   "the output copies anybody or changed the study's shape, ",
                   "and a small-cohort fit does neither."),
      fix = paste("`synpmx_avatar()` and `synpmx_pca()` need no identifiable",
                  "structure.")), call. = FALSE)
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
    stop(.condition_text(
      "`synpmx_model_estimate()` requires `nominal_time` in `pmx_roles()`.",
      why = paste("The dosing and visit models sit on the nominal grid, and",
                  "inferring that grid from recorded times is a statement",
                  "about the protocol that only you can make."),
      fix = "Add the protocol's planned times as a column and declare it."),
      call. = FALSE)
  }
  nominal <- suppressWarnings(as.numeric(source[[roles$nominal_time]]))
  relevant <- .observation_rows(source, roles, require_present = TRUE) |
    .event_rows(source, roles)
  missing <- relevant & !is.finite(nominal)
  if (any(missing)) {
    stop(.missing_nominal_time_message(source, roles, missing, relevant,
                                      "the dosing and visit models read"),
         call. = FALSE)
  }
  invisible(TRUE)
}

# How much of a curve the sampling shows. Counted on the nominal grid rather
# than on recorded times, so that "bin" needs no width chosen here: the grid is
# the protocol's own sampling slots, and a study sampled at six distinct times
# after dose has six of them however precisely the clock recorded each visit.
# Recorded times would make this a count of how noisy the clock was.
.model_time_coverage <- function(source, roles) {
  nominal <- suppressWarnings(as.numeric(source[[roles$nominal_time]]))
  planned <- source
  planned[[roles$time]] <- nominal
  tad <- .derived_tad(planned, roles)
  observed <- .observation_rows(source, roles, require_present = TRUE) &
    is.finite(tad)
  length(unique(tad[observed]))
}

# What a count of zero actually means. "0 distinct nominal times after a dose"
# reads as a statement about the sampling schedule, and it almost never is one:
# a study with no post-dose sample would not have been run. It is a statement
# about the role columns -- nothing the roles select is an observation, or
# nothing is a dose -- so the message names which of the two, with the column it
# read, rather than leaving the caller to guess what was being counted.
#
# Each branch writes its own opening clause rather than sharing one. "Found no
# observation recorded after a dose" is exactly right where the two exist and
# never meet, and misleading in front of a study that has no dose records at
# all: it sends the reader to look at the sampling times of samples that are
# fine. The finding a caller acts on is the first thing the sentence says.
.model_time_coverage_shortfall <- function(source, roles) {
  observed <- .observation_rows(source, roles, require_present = TRUE)
  dosed <- .dose_rows(source, roles)
  if (!any(observed)) {
    return(paste0("found nothing to fit: none of the ", nrow(source),
                  " rows is an observation. An observation is `", roles$evid,
                  "` 0",
                  if (!is.null(roles$mdv)) paste0(" and `", roles$mdv, "` 0"),
                  " with `", roles$dv, "` recorded. Check the columns named ",
                  "in `pmx_roles()` against the data."))
  }
  if (!any(dosed)) {
    # What the columns actually hold, because "no row is a dose" reads as a
    # claim about the study and the fix is in the data: every value the event
    # column takes, and how many rows carry an amount anyway. A dataset that
    # marks its doses with the amount alone, or that lost its dose rows to a
    # filter upstream, is the difference between those two numbers.
    return(paste0("found no dose records, so there is no concentration-time ",
                  "curve to fit: ", sum(observed), " observations are present ",
                  "and no row is a dose. A dose is a non-zero `", roles$evid,
                  "`",
                  if (!is.null(roles$amt))
                    paste0(" with a positive `", roles$amt, "`"),
                  ". `", roles$evid, "` holds ",
                  .value_counts(source[[roles$evid]]),
                  if (!is.null(roles$amt)) {
                    amount <- suppressWarnings(as.numeric(source[[roles$amt]]))
                    paste0(", and ", sum(is.finite(amount) & amount > 0),
                           " of ", nrow(source), " rows have a positive `",
                           roles$amt, "`")
                  },
                  ". Check the columns named in `pmx_roles()` against the ",
                  "data, and whether the dose records were filtered out ",
                  "before the call. A study that genuinely has no dosing ",
                  "cannot have a population model fitted to it; ",
                  "`synpmx_avatar()` and `synpmx_pca()` need no structure ",
                  "and will run on it."))
  }
  paste0("found no observation recorded after a dose, so there is no ",
         "concentration-time curve to fit: ", sum(observed),
         " observations and ", sum(dosed), " dose rows are present, but no ",
         "subject holds an observation after a dose of its own. Check that ",
         "both carry the same `", roles$id, "`.")
}

# Sparse sampling is a weak fit rather than a refusal: below `min_time_bins` the
# fixed effects sit close to their starting values and the generated profiles
# are a plausible shape rather than this study's. That is worth saying loudly
# and is the caller's call to make, so it warns. There is no simpler structural
# model to fall back to -- the default candidate set is already the
# one-compartment model, and below it there is no concentration-time curve at
# all. Nothing to fit at all is still an error, because no fit follows.
.model_require_time_coverage <- function(source, roles, minimum) {
  minimum <- .positive_integer(minimum, "min_time_bins")
  bins <- .model_time_coverage(source, roles)
  if (bins == 0L) {
    stop("`synpmx_model_estimate()` ",
         .model_time_coverage_shortfall(source, roles), call. = FALSE)
  }
  if (bins < minimum) {
    warning(.condition_text(
      "`synpmx_model_estimate()` has observations at ", bins,
      " distinct nominal time(s) after a dose, below `min_time_bins` = ",
      minimum, ".",
      why = paste("A one-compartment model is not identifiable from that",
                  "sampling: expect parameters close to their starting values",
                  "and simulated profiles that are a plausible shape rather",
                  "than an estimate of this study's."),
      fix = paste("`synpmx_avatar()` and `synpmx_pca()` carry sparse sampling",
                  "without fitting a structure to it.")), call. = FALSE)
  }
  invisible(TRUE)
}

# How long something took, written the way a person would say it. Seconds to
# one decimal below a minute, then minutes and seconds, because "just over four
# minutes" is what a caller comparing runs actually wants and "254.7 s" is not.
# `.subject_strata()` joins the strata columns with a control character, which
# is right for a key and wrong in anything a person reads: a terminal treats it
# as a carriage return and prints the arm on top of itself, so "Placebo\r0"
# comes out as "0lacebo". The arm keeps its key internally and is labelled with
# the columns joined readably. Every generator's output goes through this.
.arm_label <- function(arm) gsub("\r", " / ", as.character(arm), fixed = TRUE)

# The distinct values of one column with how often each occurs, for a message
# that has to say what a column contained rather than what it should have.
.value_counts <- function(x, limit = 8L) {
  counts <- table(as.character(x), useNA = "ifany")
  labels <- names(counts)
  labels[is.na(labels)] <- "NA"
  shown <- seq_len(min(length(counts), limit))
  paste0(paste(sprintf("%s (%d rows)", labels[shown],
                       as.integer(counts[shown])), collapse = ", "),
         if (length(counts) > limit) paste0(", and ", length(counts) - limit,
                                            " more"))
}

# How the shape was arrived at. A search only happened where more than one
# shape had a degree of freedom left to be judged on; where one did, saying
# "chosen on AIC" would describe a comparison that never ran, and where the
# endpoint is a level the note on the row is the honest account.
# Why an endpoint was read as the concentration, from the signals that decided
# it. `post_dose` and `proportional` are the two that decide, so they are always
# stated; `compartment` and `shape` only break ties and are mentioned only where
# they passed.
# One signal to a line, in the words of what was read rather than of the field
# that holds it, and only for the signals that were readable. Semicolons ran the
# four together into a sentence that had to be parsed to be counted, and two of
# them said less than they knew: "measured where the doses go" is a statement
# about the `cmt` column, and the shape is read inside one dose interval and not
# across the study.
.model_endpoint_reason <- function(signals, endpoint) {
  if (is.null(signals)) return("the signals did not separate the endpoints")
  row <- signals[signals$endpoint == endpoint, , drop = FALSE]
  if (!nrow(row)) return("the signals did not separate the endpoints")
  c(
    if (isTRUE(row$post_dose[[1L]])) "absent before the first dose" else
      "present before the first dose",
    # Silent where proportionality could not be read. A study with one dose
    # level, or an amount per patient, has nothing to compare across doses;
    # `.model_classify_endpoints()` waives the requirement there and the other
    # signals decide, so a line about a signal that did not speak is one more
    # thing to read past on the way to the ones that did.
    if (is.na(row$proportional[[1L]])) NULL else
      if (isTRUE(row$proportional[[1L]]))
        "dose-proportional: the peak scales with the dose" else
          "not dose-proportional",
    if (isTRUE(row$compartment[[1L]]))
      "recorded in the compartment the doses go into",
    if (isTRUE(row$shape[[1L]]))
      "rises to one peak and comes back down within one dose interval"
  )
}

.model_pd_selection <- function(shape) {
  table <- shape$candidates
  if (is.null(table)) return(NULL)
  converged <- table$shape[table$converged]
  note <- table$note[table$converged & nzchar(table$note)]
  if (length(converged) > 1L) {
    return(paste0("; chosen on AIC from ", paste(converged, collapse = ", ")))
  }
  if (length(note)) return(paste0("; ", note[[1L]]))
  paste0("; the only shape these observations admit")
}

.model_duration <- function(seconds) {
  if (!is.finite(seconds)) return("unknown")
  if (seconds < 60) return(sprintf("%.1f s", seconds))
  minutes <- floor(seconds / 60)
  sprintf("%d min %.0f s", minutes, seconds - 60 * minutes)
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
    stop(.condition_text(
      "`", argument, "` is a `pmx_fitted_model`, which `", what,
      "` cannot accept.",
      why = paste("Its parameters were estimated from the confidential study,",
                  "so treating them as a public input would charge no privacy",
                  "budget for information taken from the data."),
      fix = paste("Supply a `pmx_structural_model()` whose values come from a",
                  "source outside this dataset.")), call. = FALSE)
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
                              correlations = NULL, censoring = NULL,
                              quantification_floor = NULL, timing = NULL,
                              movement = NULL, fit_subjects = NULL,
                              start_param = NULL, pk_models = NULL,
                              dose_records = NULL) {
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
    censoring = censoring,
    quantification_floor = quantification_floor,
    timing = timing,
    movement = movement,
    fit_subjects = fit_subjects,
    start_param = start_param,
    # How many dose records drove each concentration endpoint's fit. What the
    # report says with it is the point: a study with two concentrations either
    # shares one administration between them or has to tell two drugs' doses
    # apart, and the count is where a reader sees which happened.
    dose_records = dose_records,
    # One entry per concentration endpoint. `structural` and `parameters` above
    # are the first of them, which is every study that fits one concentration.
    # A model assembled without the list -- a hand-built fixture, or a fit
    # stored before there could be more than one -- gets a one-entry list built
    # from those fields, so every reader can index by endpoint unconditionally.
    pk_models = pk_models %||% stats::setNames(list(list(
      endpoint = endpoints$pk[[1L]], structural = structural,
      parameters = parameters, candidates = candidates, movement = movement,
      design = design, start_param = start_param,
      effects = covariate_effects,
      seconds = timing$fit %||% NA_real_)), endpoints$pk[[1L]]),
    endpoints = endpoints
  ), class = "pmx_fitted_model")
}

# Printing the object prints the whole inventory, because everything on it is
# an input to the simulation and a reader who has to call a second function to
# see half of them will read half of them. `model_report()` is the same content
# as a list, for a caller reading a number out rather than looking at it.
#' @export
print.pmx_fitted_model <- function(x, ...) {
  cat("A fitted PMX model, from synpmx_model_estimate()\n")
  cat("Everything below is an input to `synpmx_model_generate()`.\n\n")
  cat("  candidates fitted ", nrow(x$candidates),
      sprintf("(%s selected on AIC)", x$structural), "\n")
  cat("\n")
  print(model_report(x))
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
    movement = fitted_model$movement,
    fit_subjects = fitted_model$fit_subjects,
    start_param = fitted_model$start_param,
    pk_models = fitted_model$pk_models,
    covariate_effects = fitted_model$covariate_effects,
    correlations = fitted_model$correlations,
    censoring = fitted_model$censoring,
    quantification_floor = fitted_model$quantification_floor,
    timing = fitted_model$timing,
    pd = fitted_model$pd,
    arms = fitted_model$arms,
    dosing = fitted_model$dosing,
    visits = fitted_model$visits,
    cells = fitted_model$cells,
    dose_records = fitted_model$dose_records,
    roles = fitted_model$roles,
    covariates = fitted_model$covariates,
    discrete = fitted_model$discrete,
    schema = fitted_model$schema,
    settings = fitted_model$settings
  ), class = "pmx_model_report")
}

#' @export
print.pmx_model_report <- function(x, ...) {
  # One label column for the whole report, wrapped under itself, because these
  # lines carry sentences rather than single numbers.
  field <- function(label, ...) {
    cat(.wrap_plain(paste0(...), sprintf("  %-18s ", label),
                    strrep(" ", 21L)), "\n", sep = "")
  }
  # A shape, one parameter to a line. The generator draws every one of these
  # numbers, and as a wrapped paragraph a reader looking for one of them had to
  # find it mid-sentence, three lines down, beside a different endpoint's.
  shape_block <- function(shape, indent) {
    lines <- sprintf("%s%-16s %.4g", indent, names(shape$typical),
                     as.numeric(shape$typical))
    lines <- c(lines,
               sprintf("%s%-16s %.3g (SD on the log baseline)", indent,
                       "between-subject", shape$baseline_cv %||% 0),
               sprintf("%s%-16s %s %.3g", indent, "residual",
                       shape$residual$kind, shape$residual$sd))
    selection <- .model_pd_selection(shape)
    if (!is.null(selection)) {
      lines <- c(lines, paste0(indent, sub("^; ", "", selection)))
    }
    cat(paste(lines, collapse = "\n"), "\n", sep = "")
  }
  # Simplest first. The trial and what was done to the patients in it, then
  # the shapes fitted to the endpoints that need least, then the population
  # model -- so a reader meets a summary of the source before a structural
  # model of it. The `did not move` banner keeps the top whatever the order,
  # because it is the one line that says do not use this fit at all.
  # Before the numbers, not after: a reader who stops at the parameter block
  # has to have been told already that it is not a parameter block.
  if (!is.null(x$movement) && !isTRUE(x$movement$moved)) {
    changes <- x$movement$changes
    cat("!! THE FIT DID NOT MOVE !!\n")
    cat("  ", .wrap_plain(paste0(
      "Every parameter came back within ",
      sprintf("%.0f%%", 100 * x$movement$tolerance), " of its starting value, ",
      "so everything under `Estimated by nlmixr2` below is a starting value ",
      "rather than an estimate. Do not generate from this fit."
    ), "", "  "), "\n", sep = "")
    cat(sprintf("    %-10s %-12s -> %-12s (%.2f%%)\n", changes$parameter,
                sprintf("%.4g", changes$start),
                sprintf("%.4g", changes$estimate),
                100 * changes$relative_change), sep = "")
    cat("\n")
  }
  cat("Summarized from the source, not estimated\n")
  # One arm to a line. A study whose arms are named for their regimen -- a
  # priming dose, a maintenance dose, a route and a population, all in the
  # label -- runs to sixty characters an arm, and seven of those wrapped into
  # a paragraph is not a list of arms any reader can count.
  field("cohort", x$n_source, " patients in ", length(x$arms$arms), " arm(s)")
  cat(sprintf("%s%s (%d)\n", strrep(" ", 21L),
              .arm_label(names(x$arms$sizes)), as.integer(x$arms$sizes)),
      sep = "")
  # How the doses are given, where the study says. A schedule of the same
  # amounts means something different given over an hour than given as a bolus,
  # and different again given subcutaneously, so the report says which.
  # One entry per arm, or one per arm per drug where `dose_endpoints` split the
  # doses, labelled so that a reader sees which schedule belongs to which drug.
  dosing <- list()
  for (arm in names(x$dosing)) {
    models <- .dose_group_models(x$dosing[[arm]])
    for (k in seq_along(models)) {
      label <- if (is.null(names(models))) arm else
        paste0(arm, " / ", names(models)[[k]])
      dosing[[label]] <- models[[k]]
    }
  }
  routes <- .unique_in_order(unlist(lapply(dosing, function(d) d$planned$route)))
  routes <- routes[!is.na(routes)]
  durations <- unlist(lapply(dosing, function(d) {
    with(d$planned, ifelse(rate > 0, amt / rate, NA_real_))
  }))
  durations <- durations[is.finite(durations) & durations > 0]
  if (length(routes) || length(durations)) {
    field("dose routes",
          if (length(routes)) paste(routes, collapse = " and ") else
            "one route, undeclared",
          if (length(durations)) {
            sprintf("; infused over %s h", if (length(unique(round(durations, 6))) == 1L)
              signif(durations[[1L]], 4) else
                paste(signif(range(durations), 4), collapse = " to "))
          } else if (length(routes)) "; given as a bolus" else NULL)
  }

  # The three rates are the whole model of missed doses and reductions, and a
  # reader looking for "what happens to the dosing" has to be able to find them
  # by name. Reported per arm, because they are per arm.
  rates <- data.frame(
    arm = names(dosing),
    reduce = vapply(dosing, function(d) d$reduction, numeric(1)),
    skip = vapply(dosing, function(d) d$interruption, numeric(1)),
    stop_early = vapply(dosing, function(d) d$discontinuation, numeric(1)),
    levels = vapply(dosing, function(d) length(d$levels), integer(1)),
    stringsAsFactors = FALSE
  )
  if (any(rates$reduce > 0 | rates$skip > 0 | rates$stop_early > 0)) {
    field("dose changes", "per planned cycle:")
    cat(sprintf("      %-18s reduce %.0f%%, skip %.0f%%, stop early %.0f%% (%d dose level(s))\n",
                .arm_label(rates$arm), 100 * rates$reduce, 100 * rates$skip,
                100 * rates$stop_early, rates$levels), sep = "")
  } else {
    field("dose changes", "none")
  }

  # Attendance is the model of a missed observation, and it is one probability
  # per endpoint per nominal time per arm -- the fraction of that arm with an
  # observation there. Said in those words: the first readers of this report
  # took "cell(s) of the visit grid" for a count of nominal times, which it is
  # not, because endpoints are not all measured at the same times. The spread is
  # what a reader needs beside the median, because a slot nobody misses and a
  # slot half the arm misses are the same line otherwise.
  attendance <- unlist(lapply(x$visits, function(v) as.numeric(v$probability)))
  field("visit grid", length(unique(x$cells$endpoint)), " endpoint(s) at ",
        length(unique(x$cells$time)), " nominal time(s), ", nrow(x$cells),
        " slot(s) in all")
  field("visit attendance", if (length(attendance)) sprintf(
    "median %.0f%% of an arm attends a slot (%.0f%% to %.0f%%)",
    100 * stats::median(attendance), 100 * min(attendance),
    100 * max(attendance)) else "no attendance model")

  if (length(x$covariates)) {
    field("covariates", paste0(
      paste(sprintf("%s %s", names(x$covariates),
                    vapply(x$covariates, function(spec) spec$kind,
                           character(1))),
            collapse = ", "),
      ", drawn once for the whole study, independently of the profiles"))
  }
  drawn <- .unique_in_order(unlist(lapply(x$discrete %||% list(),
    function(arm) x$cells$endpoint[!vapply(arm, is.null, logical(1))])))
  if (length(drawn)) {
    field("discrete endpoints", paste(drawn, collapse = ", "),
          ": drawn from each arm's recorded frequencies at each visit, ",
          "not simulated")
  }
  if (!is.null(x$schema)) {
    field("columns emitted", paste(x$schema$columns %||% character(0),
                                   collapse = ", "))
  }

  # The assay limit and the emission floor are what confused the first readers
  # of this report, because both are one terse number about values near zero and
  # they mean opposite things: one is what the source reported and the fit was
  # given instead, the other is what the generator refuses to write. Each is
  # spelled out in a sentence rather than compressed into a label.
  censored <- if (is.null(x$censoring)) NULL else
    x$censoring[x$censoring$censored > 0, , drop = FALSE]
  if ((!is.null(censored) && nrow(censored)) || length(x$quantification_floor)) {
    cat("\nValues at the lower limit of what was observed\n")
  }
  if (!is.null(censored) && nrow(censored)) {
    cat("  Reported below the assay limit:\n")
    cat(sprintf("    %-18s %d of %d (%.0f%%) below %.4g (the limit)\n",
                censored$endpoint, censored$censored, censored$observations,
                100 * censored$fraction, censored$limit), sep = "")
  }
  # One line per endpoint, the same shape as the censored lines above it, and
  # carrying its own sentence rather than a paragraph over the group: an
  # endpoint with a declared limit and an endpoint with a floor are the same
  # kind of fact about the same kind of number, and reading the second one
  # meant reading a wrapped preamble first to find out which it was.
  floors <- unlist(x$quantification_floor)
  cat(sprintf("    %-18s %.4g, half the smallest value seen, no assay limit\n",
              names(floors), floors), sep = "")

  # The shape name alone is not the fit: the generator draws every one of these
  # numbers, so a report that is an inventory of its inputs has to show them.
  # `constant` and `linear` come from `lm()` and `exponential` from `nls()`,
  # all on study time from the first dose.
  if (length(x$pd)) {
    cat("\nEach non-PK continuous endpoint, fitted as constant, linear,",
        "or exponential\n")
    for (name in names(x$pd)) {
      shape <- x$pd[[name]]
      # One line per arm where the shape was fitted per arm, because that is
      # then six fits rather than one and the arms are the whole point of
      # asking for it.
      if (length(shape$arms)) {
        field(name, "one shape per arm (`pd_by_arm = TRUE`):")
        for (arm in names(shape$arms)) {
          own <- shape$arms[[arm]]
          field("", "  ", .arm_label(arm), ": ", own$pd)
          shape_block(own, strrep(" ", 25L))
        }
        next
      }
      field(name, shape$pd)
      shape_block(shape, strrep(" ", 23L))
    }
  }

  # Which endpoint carries the structural model, and on what grounds. The four
  # signals behind an inferred answer are on the object at
  # `fit$endpoints$signals`; as a grid of bare logicals they read as a puzzle,
  # so what prints here is the decision in words. Before the model rather than
  # after it: the reader is told which endpoint was fitted, and why that one,
  # before they are shown the fit to it.
  cat("\nPK endpoint for the PopPK model\n")
  for (endpoint in x$endpoints$pk) {
    if (identical(x$endpoints$decided_by, "declared")) {
      field(endpoint, "declared through `endpoint_roles`")
    } else {
      field(endpoint, "inferred from the following data characteristics:")
      for (reason in .model_endpoint_reason(x$endpoints$signals, endpoint)) {
        cat(.wrap_plain(reason, strrep(" ", 23L), strrep(" ", 25L)), "\n",
            sep = "")
      }
    }
  }
  # Only where the doses were NOT separated, because only then is there
  # something to say: every dose record drives every concentration, which is
  # right for a parent and its metabolite, which share one administration, and
  # wrong for two drugs given together, where each model is fitted against the
  # other drug's doses as well as its own -- and silent until it was printed.
  # Declared, the split is the unsurprising case and its record count is a
  # number nobody reads.
  if (length(x$endpoints$pk) > 1L && !is.null(x$dose_records) &&
      is.null(x$roles$dose_endpoints)) {
    field("doses", "all ", x$dose_records[[1L]], " dose record(s) drive every ",
          "endpoint above: right for a parent and its metabolite, wrong for ",
          "two drugs given together, which `dose_endpoints` in `pmx_roles()` ",
          "separates")
  }
  if (!is.null(x$design)) {
    field("route", x$design$route, ": ", x$design$reason)
    if (isTRUE(x$design$richness$rich) && !grepl("^2cmt", x$structural)) {
      field("also available",
            sprintf("2cmt_%s, which the sampling would support: median %g distinct times after a dose, %g after the peak",
                    if (grepl("oral", x$structural)) "oral" else "iv",
                    x$design$richness$per_subject,
                    x$design$richness$after_peak))
    }
  }

  cat("\nThe PopPK model\n\n")
  models <- x$pk_models
  if (!length(models)) models <- list(list(structural = x$structural,
                                           parameters = x$parameters,
                                           movement = x$movement,
                                           start_param = x$start_param))
  for (k in seq_along(models)) {
    own <- models[[k]]
    cat(if (length(models) > 1L)
      sprintf("Estimated by nlmixr2: `%s`\n", own$endpoint) else
        "Estimated by nlmixr2\n")
    cat("  structural model  ", own$structural, "\n")
    if (k == 1L && !is.null(x$fit_subjects)) {
      fs <- x$fit_subjects
      field("fitted on", if (fs$fitted < fs$of) paste0(
        fs$fitted, " of ", fs$of, " patients with a concentration, drawn in ",
        "proportion to the arms under `max_fit_subjects` = ", fs$cap, "; the ",
        "dosing, visit and covariate models below read the whole study") else
          paste0("all ", fs$fitted, " patients with a concentration"))
    }
    if (!is.null(own$movement) && !isTRUE(own$movement$moved)) {
      changes <- own$movement$changes
      cat("  !! THE FIT DID NOT MOVE !!\n")
      cat("  ", .wrap_plain(paste0(
        "Every parameter came back within ",
        sprintf("%.0f%%", 100 * own$movement$tolerance),
        " of its starting value, so the numbers below are starting values ",
        "rather than estimates. Do not generate from this fit."
      ), "", "  "), "\n", sep = "")
      cat(sprintf("    %-10s %-12s -> %-12s (%.2f%%)\n", changes$parameter,
                  sprintf("%.4g", changes$start),
                  sprintf("%.4g", changes$estimate),
                  100 * changes$relative_change), sep = "")
    }
    cat("  fixed effects     ",
        paste(sprintf("%s %.4g", names(own$parameters$fixed),
                      as.numeric(own$parameters$fixed)), collapse = ", "), "\n")
    cat("  between-subject   ",
        paste(sprintf("%s %.3g", rownames(own$parameters$omega),
                      sqrt(diag(own$parameters$omega))), collapse = ", "),
        "(as SD on the log scale)\n")
    if (!is.null(own$movement) && isTRUE(own$movement$moved) &&
        isFALSE(own$movement$omega_moved)) {
      field("", .wrap_plain(paste0(
        "!! Every between-subject term is within ",
        sprintf("%.0f%%", 100 * own$movement$tolerance),
        " of its starting value of ", .model_eta_init, " while the fixed ",
        "effects moved: the fit found a population mean and did not estimate ",
        "its spread. Synthetic subjects will be spread by the starting value, ",
        "not by this study.")))
    }
    if (length(own$start_param)) {
      field("starting values", paste(sprintf("%s %.4g", names(own$start_param),
                                             own$start_param), collapse = ", "),
            " declared through `start_param`; the rest were read off the ",
            "cohort's median profile")
    }
    cat("  residual error    ", own$parameters$residual$kind,
        sprintf("%.3g", own$parameters$residual$cv %||%
                  own$parameters$residual$sd), "\n")
  }
  # What the wait was. A fit is the slow thing this package does, and a caller
  # deciding whether to change a setting and run it again asks this first.
  if (!is.null(x$timing)) {
    field("time to fit", .model_duration(x$timing$fit),
          sprintf(" (%s for the whole call)",
                  .model_duration(x$timing$total)))
    # Which fit the wait was. Every candidate is a separate compiled
    # population fit, so a search that ran four of them waited four times --
    # and one to a line, because the number a reader is after is the one
    # candidate that took the minutes, which a wrapped list buries.
    if (!is.null(x$timing$candidates) && nrow(x$timing$candidates)) {
      rows <- x$timing$candidates
      field("", "nlmixr2")
      cat(sprintf("%s%-26s %s%s\n", strrep(" ", 23L), rows$model,
                  vapply(rows$seconds, .model_duration, character(1)),
                  ifelse(rows$converged, "", " (did not converge)")),
          sep = "")
    }
  }
  # Only where something was fitted. `covariate_effects` is `"none"` by
  # default, and a line reading "none" on every report says nothing.
  if (length(x$covariate_effects)) {
    field("covariate effects",
          paste(vapply(names(x$covariate_effects), function(parameter) {
            effect <- x$covariate_effects[[parameter]]
            sprintf("%s ~ (%s/%.4g)^%.2f", parameter, effect$covariate,
                    effect$reference, effect$exponent)
          }, character(1)), collapse = ", "))
  }
  # Which endpoint carries the structural model, and on what grounds. The four
  # signals behind an inferred answer are on the object at
  # `fit$endpoints$signals`; as a grid of bare logicals they read as a puzzle,
  # so what prints here is the decision in words.
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
#' @return A data frame with columns `model`, `converged`, `aic`, `seconds` --
#'   how long that candidate took to fit -- and `note`.
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
