# The four generation modes ---------------------------------------------------
#
# Every mode answers the same question -- "give me a synthetic dataset" -- and
# every one of these functions returns a data frame. They differ only in what
# the output is built from, and therefore in what may be claimed about it.
#
#   synpmx_avatar()      real subject templates and blended real trajectories
#   synpmx_prior()       a public model and protocol only, no data read
#   synpmx_calibrated()  a public model, magnitude corrected by a small release
#   synpmx_empirical()   a dense set of noised population summaries
#
# The last two read confidential data and spend privacy budget. They attach
# their release to the returned dataset so that further datasets can be drawn
# from the same release as post-processing, at no additional cost:
#
#   syn  <- synpmx_calibrated(...)     # spends epsilon once
#   syn2 <- synpmx_generate(syn)       # free
#
# `.release_of()` is what makes that work, and it is also why the accounting
# helpers accept a generated dataset as readily as a release object.

#' The release behind a generated dataset
#'
#' Accepts a dataset produced by [synpmx_calibrated()] or [synpmx_empirical()]
#' and returns the release attached to it, or passes a release object through
#' unchanged. This is what lets `privacy_report(syn)` work on a dataset.
#'
#' @param x A generated dataset, or a release object.
#' @param what Human-readable description of the caller, used in errors.
#'
#' @return A `pmx_calibrated_model` or `private_pmx_model`.
#' @keywords internal
.release_of <- function(x, what = "x") {
  if (inherits(x, c("pmx_calibrated_model", "private_pmx_model"))) {
    return(x)
  }
  release <- attr(x, "synpmx_release", exact = TRUE)
  if (!is.null(release)) {
    return(release)
  }
  if (is.data.frame(x)) {
    stop("`", what, "` carries no privacy release. Only datasets from ",
         "`synpmx_calibrated()` or `synpmx_empirical()` do; ",
         "`synpmx_avatar()` and `synpmx_prior()` make no formal claim.",
         call. = FALSE)
  }
  stop("`", what, "` must be a generated dataset or a release object.",
       call. = FALSE)
}

# One dataset, or a list of them, from an already-paid-for release. Splitting
# this out keeps `n_datasets` behaving identically across both private modes.
.draw_datasets <- function(release, n_subjects, seed, n_datasets, generator) {
  n_datasets <- .positive_integer(n_datasets, "n_datasets")
  seeds <- seed + seq_len(n_datasets) - 1L
  out <- lapply(seeds, function(s) {
    dataset <- generator(release, n_subjects, s)
    attr(dataset, "synpmx_release") <- release
    dataset
  })
  if (n_datasets == 1L) out[[1L]] else out
}

# Spending budget twice on the same data is nearly always a mistake rather than
# an intent, and nothing else in the package would notice. Track fits within the
# session and say so once per repeat.
.fit_registry <- new.env(parent = emptyenv())

.warn_on_repeat_fit <- function(data, epsilon, mode) {
  key <- paste(mode, epsilon, nrow(data), ncol(data),
               sum(vapply(data, function(col) {
                 sum(as.numeric(suppressWarnings(as.numeric(col))), na.rm = TRUE)
               }, numeric(1))), sep = "|")
  seen <- get0(key, envir = .fit_registry, ifnotfound = 0L)
  assign(key, seen + 1L, envir = .fit_registry)
  if (seen > 0L) {
    warning(
      "This looks like the ", seen + 1L, .ordinal_suffix(seen + 1L),
      " budget-spending fit against the same data at epsilon = ", epsilon,
      " in this session, so roughly ", (seen + 1L) * epsilon,
      " has now been spent in total.\n",
      "Drawing more datasets from one release costs nothing: use ",
      "`synpmx_generate(syn, seed = ...)`, or ask for several at once with ",
      "`n_datasets =`.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

.ordinal_suffix <- function(n) {
  if (n %% 100L %in% 11:13) return("th")
  switch(as.character(n %% 10L), "1" = "st", "2" = "nd", "3" = "rd", "th")
}

# Session-only opt-in gate for the DP engines --------------------------------
#
# The DP engines' unaudited status was documented but not enforced: calling
# synpmx_calibrated()/synpmx_empirical() was exactly as frictionless as
# calling the audited synpmx_avatar(). See REV-023 in
# design/REVIEW_BACKLOG.md. A per-call argument defaulting to off is too easy
# to flip once and forget; a warning() is too easy to suppress in a script.
# An in-session flag makes the acknowledgment a deliberate, visible action
# that must be repeated in every fresh session, script run, or CI job.
.dp_engines_enabled <- new.env(parent = emptyenv())
.dp_engines_enabled$on <- FALSE

.dp_audit_status <- paste(
  "the differentially private engines are complete and tested, but not",
  "under active development, carry known open findings (see",
  "design/REVIEW_BACKLOG.md), and have not been independently",
  "privacy-audited. See vignette(\"synpmx-privacy\") for the trust-boundary",
  "decision rule and what a production release additionally needs."
)

#' Acknowledge the DP engines' unaudited status for this session
#'
#' [synpmx_calibrated()] and [synpmx_empirical()] refuse to run until this has
#' been called at least once in the current session. This is a deliberate
#' speed bump, not a technical safeguard: the differentially private engines
#' are complete and tested, but not under active development, carry known
#' open findings (see `design/REVIEW_BACKLOG.md`), and have not been
#' independently privacy-audited. See `vignette("synpmx-privacy")` for the
#' trust-boundary decision rule and what a production release additionally
#' needs. [synpmx_avatar()] and [synpmx_prior()] make no differential-privacy
#' claim and are unaffected.
#'
#' The acknowledgment does not persist. It applies only to the R session it is
#' called in, so a fresh session, script run, or CI job must call it again.
#' `backend = "public"` calls make no DP claim either and are exempt from this
#' gate.
#'
#' @return `TRUE`, invisibly.
#' @seealso [synpmx_disable_dp_engines()]
#' @export
#' @examples
#' synpmx_enable_dp_engines()
synpmx_enable_dp_engines <- function() {
  assign("on", TRUE, envir = .dp_engines_enabled)
  message("DP engines enabled for this session: ", .dp_audit_status)
  invisible(TRUE)
}

#' Withdraw the acknowledgment from [synpmx_enable_dp_engines()]
#'
#' Mainly useful for tests that exercise the unacknowledged error path without
#' starting a fresh session.
#'
#' @return `TRUE`, invisibly.
#' @seealso [synpmx_enable_dp_engines()]
#' @export
synpmx_disable_dp_engines <- function() {
  assign("on", FALSE, envir = .dp_engines_enabled)
  invisible(TRUE)
}

.require_dp_engines_enabled <- function(fn) {
  if (isTRUE(get0("on", envir = .dp_engines_enabled, ifnotfound = FALSE))) {
    return(invisible(NULL))
  }
  stop(
    "`", fn, "()` is one of the differentially private engines: ",
    .dp_audit_status, " Call `synpmx_enable_dp_engines()` once per session ",
    "to acknowledge this and proceed.",
    call. = FALSE
  )
}


#' Generate a dataset from public inputs only
#'
#' Simulates from a public structural model and a public protocol. No
#' confidential data is read, so there is nothing to protect and no budget to
#' spend: this is `epsilon = 0`, the strongest possible guarantee.
#'
#' The typical parameter values must come from somewhere that is not the data --
#' allometric scaling from preclinical work, a published model for the compound
#' class, or the reasoning that set the starting dose. The output is exactly as
#' good as that prior.
#'
#' @param model A [pmx_structural_model()].
#' @param design A [pmx_trial_design()].
#' @param n_subjects Number of subjects. Defaults to the planned cohort total.
#' @param seed Ordinary generation seed. Unrelated to privacy noise.
#' @param dropout Fraction of subjects who discontinue early. A public
#'   assumption from the protocol.
#' @param lloq Lower limit of quantification. Observations below it are flagged
#'   `CENS = 1` with `DV` at the limit, following the Monolix convention.
#' @param covariates Optional [pmx_covariates()].
#'
#' @return A data frame in the generated event-table schema; see
#'   [pmx_generated_roles()].
#' @seealso [synpmx_avatar()], [synpmx_calibrated()], [synpmx_empirical()]
#' @export
#' @examples
#' model <- pmx_structural_model(
#'   pk = "1cmt_oral", typical = c(cl = 6, v = 35, ka = 1.5),
#'   source = "illustrative allometric scaling"
#' )
#' design <- pmx_trial_design(
#'   dose_levels = 320, cohort_sizes = 12, sampling = c(0, 1, 2, 4, 9, 24),
#'   source = "illustrative protocol"
#' )
#' syn <- synpmx_prior(model, design, n_subjects = 12, seed = 202)
#' head(syn, 3)
synpmx_prior <- function(model, design, n_subjects = NULL, seed = NULL,
                         dropout = 0, lloq = NULL, covariates = NULL) {
  if (!inherits(model, "pmx_structural_model")) {
    stop("`model` must come from `pmx_structural_model()`. To generate from a ",
         "calibrated release, use `synpmx_generate()`.", call. = FALSE)
  }
  .generate_structural(model, design = design, n_subjects = n_subjects,
                       seed = seed, dropout = dropout, lloq = lloq,
                       covariates = covariates)
}


#' Generate a dataset from a privately calibrated structural model
#'
#' Keeps a public structural model's *shape* and spends a small privacy budget
#' correcting only its *magnitude*. Each subject is reduced to a bounded
#' multiplicative correction of the model's own prediction, clipped to a public
#' prior range, and released with calibrated noise. Only a handful of numbers
#' leave the data, which is why this mode stays usable at 20 to 60 subjects.
#'
#' This is the recommended differentially private path for early-phase cohorts.
#' The tradeoff is that everything not calibrated is *asserted*: curve shape,
#' variability, residual error, and covariate relationships come from the public
#' model, so the output is only as realistic as that model.
#'
#' **This function spends privacy budget.** Calling it again spends the budget
#' again. To draw further datasets from the release you have already paid for,
#' use [synpmx_generate()] or ask for several at once with `n_datasets`.
#'
#' @section Maintenance status:
#' A secondary, provided-as-is path. [synpmx_avatar()] is the primary, maintained
#' method. The differentially private modes are complete and tested but not under
#' active development, carry known open findings, and have not been independently
#' privacy-audited. Use them to demonstrate the privacy/utility tradeoff, not as
#' a production release mechanism; a real regulated release needs specialist
#' review and the external OpenDP backend. Requires a session-level
#' acknowledgment: call [synpmx_enable_dp_engines()] once before this function
#' will run (unless `backend = "public"`, which makes no DP claim).
#'
#' @param data The confidential dataset.
#' @param roles A [pmx_roles()] declaration for `data`.
#' @param model A public [pmx_structural_model()].
#' @param design A public [pmx_trial_design()].
#' @param priors A [pmx_priors()] giving a public range per released correction.
#' @param epsilon The privacy budget. A governance decision, not a default;
#'   [pmx_preflight()] reports the expected fold-error before any is spent.
#' @param n_subjects Number of subjects to generate. Defaults to the released
#'   noisy count.
#' @param seed Ordinary generation seed. Unrelated to privacy noise, which is
#'   controlled by the backend and is never user-seeded.
#' @param n_datasets Number of datasets to draw from the single release. One
#'   dataset is returned directly; several are returned as a list.
#' @param covariates Optional [pmx_covariates()].
#' @param backend Privacy backend. Defaults to the validated OpenDP adapter and
#'   fails closed if it is unavailable.
#' @param public_source Assert that `data` is genuinely public. Required by, and
#'   only meaningful for, `backend = "public"`, which makes no DP claim.
#'
#' @return A data frame in the generated event-table schema, carrying its
#'   release so that [privacy_report()] and [synpmx_generate()] can read it.
#'   A list of such data frames when `n_datasets > 1`.
#' @seealso [synpmx_generate()] to draw more datasets for free,
#'   [privacy_report()] for the realized accounting.
#' @export
synpmx_calibrated <- function(data, roles, model, design, priors, epsilon,
                              n_subjects = NULL, seed = 123, n_datasets = 1L,
                              covariates = NULL, backend = "opendp",
                              public_source = FALSE) {
  if (!identical(backend, "public")) {
    .require_dp_engines_enabled("synpmx_calibrated")
  }
  .warn_on_repeat_fit(data, epsilon, "calibrated")
  release <- .fit_calibrated(
    data = data, roles = roles, model = model, design = design,
    priors = priors, epsilon = epsilon, covariates = covariates,
    backend = backend, public_source = public_source
  )
  .draw_datasets(release, n_subjects, seed, n_datasets,
                 function(r, n, s) .generate_structural(r, n_subjects = n,
                                                        seed = s))
}


#' Generate a dataset from a dense differentially private release
#'
#' Rather than asserting the curve shape, this mode measures it: it releases
#' noised summaries for the subject count, event and regimen structure,
#' observation timing, endpoint trajectories, baseline covariates, and
#' censoring, then rebuilds subjects from those summaries alone.
#'
#' It asserts less than [synpmx_calibrated()] but releases far more numbers, so
#' one epsilon is split many ways. Utility therefore collapses below a few
#' hundred subjects; this mode earns its keep on large pooled corpora.
#'
#' **This function spends privacy budget.** Calling it again spends the budget
#' again. To draw further datasets from the release you have already paid for,
#' use [synpmx_generate()] or ask for several at once with `n_datasets`.
#'
#' @inheritSection synpmx_calibrated Maintenance status
#'
#' @param data The confidential dataset.
#' @param roles A [pmx_roles()] declaration for `data`.
#' @param endpoints Named list of [pmx_endpoint()] declarations.
#' @param epsilon The privacy budget. A governance decision, not a default.
#' @param delta Additive slack in the probability bound. The implemented Laplace
#'   releases spend none, so realized accounting reports `delta = 0`.
#' @param bounds Public clipping domains from [pmx_bounds()].
#' @param public_design A [pmx_public_design()].
#' @param contribution_limits A [pmx_contribution_limits()].
#' @param budget_allocation A [pmx_budget_allocation()] splitting `epsilon`
#'   across release groups.
#' @param n_subjects Number of subjects to generate. Defaults to the released
#'   noisy count.
#' @param seed Ordinary generation seed. Unrelated to privacy noise.
#' @param n_datasets Number of datasets to draw from the single release. One
#'   dataset is returned directly; several are returned as a list.
#' @param delta_justification Required when `delta > 0`.
#' @param backend Privacy backend. Defaults to the validated OpenDP adapter and
#'   fails closed if it is unavailable.
#' @param public_source Assert that `data` is genuinely public. Required by, and
#'   only meaningful for, `backend = "public"`, which makes no DP claim.
#'
#' @return A data frame in the source event-table schema, carrying its release
#'   so that [privacy_report()] and [synpmx_generate()] can read it. A list of
#'   such data frames when `n_datasets > 1`.
#' @seealso [synpmx_generate()] to draw more datasets for free,
#'   [privacy_report()] for the realized accounting.
#' @export
synpmx_empirical <- function(data, roles, endpoints, epsilon, delta, bounds,
                             public_design, contribution_limits,
                             budget_allocation, n_subjects = NULL, seed = 123,
                             n_datasets = 1L, delta_justification = NULL,
                             backend = "opendp", public_source = FALSE) {
  if (!identical(backend, "public")) {
    .require_dp_engines_enabled("synpmx_empirical")
  }
  .warn_on_repeat_fit(data, epsilon, "empirical")
  release <- .fit_private(
    data = data, roles = roles, endpoints = endpoints, epsilon = epsilon,
    delta = delta, bounds = bounds, public_design = public_design,
    contribution_limits = contribution_limits,
    budget_allocation = budget_allocation,
    delta_justification = delta_justification, backend = backend,
    public_source = public_source
  )
  .draw_datasets(release, n_subjects, seed, n_datasets,
                 function(r, n, s) .generate_private(r, n_subjects = n,
                                                     seed = s))
}


#' Draw another dataset from a release already paid for
#'
#' Generation from an existing release is post-processing: it reads no
#' confidential data and consumes no additional privacy budget, so any number of
#' datasets may be drawn from one [synpmx_calibrated()] or [synpmx_empirical()]
#' call.
#'
#' @param x A dataset returned by [synpmx_calibrated()] or
#'   [synpmx_empirical()], or the release itself.
#' @param n_subjects Number of subjects. Defaults to the released noisy count.
#' @param seed Ordinary generation seed.
#' @param n_datasets Number of datasets to draw. One is returned directly;
#'   several are returned as a list.
#'
#' @return A data frame carrying the same release, or a list of them.
#' @seealso [synpmx_calibrated()], [synpmx_empirical()]
#' @export
synpmx_generate <- function(x, n_subjects = NULL, seed = 123, n_datasets = 1L) {
  release <- .release_of(x, "x")
  generator <- if (inherits(release, "pmx_calibrated_model")) {
    function(r, n, s) .generate_structural(r, n_subjects = n, seed = s)
  } else {
    function(r, n, s) .generate_private(r, n_subjects = n, seed = s)
  }
  .draw_datasets(release, n_subjects, seed, n_datasets, generator)
}
