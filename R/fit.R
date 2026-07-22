.validate_fit_configuration <- function(data, roles, endpoints, epsilon, delta,
                                        delta_justification, bounds,
                                        public_design, contribution_limits,
                                        budget_allocation) {
  if (!is.data.frame(data) || !nrow(data)) {
    stop("`data` must be a nonempty data frame or tibble.", call. = FALSE)
  }
  .assert_roles(data, roles)
  if (!is.numeric(epsilon) || length(epsilon) != 1L || is.na(epsilon) ||
      !is.finite(epsilon) || epsilon <= 0) {
    stop("`epsilon` must be supplied explicitly as one finite positive value.",
         call. = FALSE)
  }
  if (!is.numeric(delta) || length(delta) != 1L || is.na(delta) ||
      !is.finite(delta) || delta < 0 || delta >= 1) {
    stop("`delta` must be supplied explicitly in [0, 1).", call. = FALSE)
  }
  if (delta > 0 &&
      (is.null(delta_justification) || !is.character(delta_justification) ||
       length(delta_justification) != 1L || is.na(delta_justification) ||
       !nzchar(trimws(delta_justification)))) {
    stop("A nonempty `delta_justification` is required when delta is positive.",
         call. = FALSE)
  }
  if (!inherits(bounds, "pmx_bounds")) {
    stop("`bounds` must be created by `pmx_bounds()`.", call. = FALSE)
  }
  if (!inherits(public_design, "pmx_public_design")) {
    stop("A public schema/design from `pmx_public_design()` is required.",
         call. = FALSE)
  }
  if (!inherits(contribution_limits, "pmx_contribution_limits")) {
    stop("`contribution_limits` must come from `pmx_contribution_limits()`.",
         call. = FALSE)
  }
  if (!inherits(budget_allocation, "pmx_budget_allocation")) {
    stop("`budget_allocation` must come from `pmx_budget_allocation()`.",
         call. = FALSE)
  }
  budget_values <- unlist(budget_allocation, use.names = TRUE)
  budget_names <- c(
    "subject_count", "event", "timing", "covariates", "endpoints",
    "censoring"
  )
  if (!is.numeric(budget_values) ||
      !identical(names(budget_values), budget_names) ||
      anyNA(budget_values) || any(!is.finite(budget_values)) ||
      any(budget_values < 0) || sum(budget_values) <= 0 ||
      sum(budget_values) > 1) {
    stop("`budget_allocation` is invalid or has been modified after construction.",
         call. = FALSE)
  }
  if (!is.null(roles$amt) && is.null(bounds$amt)) {
    stop("Public `amt` bounds are required when an AMT role is declared.",
         call. = FALSE)
  }
  if (!is.null(roles$rate) && is.null(bounds$rate)) {
    stop("Public `rate` bounds are required when a RATE role is declared.",
         call. = FALSE)
  }
  if (!is.null(roles$assigned_dose) &&
      (is.null(roles$amt) || is.null(roles$occasion))) {
    stop(
      "An `assigned_dose` role requires both AMT and occasion roles.",
      call. = FALSE
    )
  }
  for (name in roles$subject_properties) {
    column <- .schema_column(public_design$schema, name)
    levels <- if (!is.null(column) && "factor" %in% column$class) {
      column$levels
    } else {
      public_design$category_levels[[name]]
    }
    if (is.null(levels) || !length(levels)) {
      stop(
        "Public category levels are required for subject property `", name,
        "`.", call. = FALSE
      )
    }
  }
  if (!is.null(roles$cens) && is.null(roles$limit) &&
      !length(bounds$limit) &&
      !any(vapply(endpoints, function(x) !is.null(x$censoring), logical(1)))) {
    warning(
      "CENS is declared without LIMIT or public censoring boundaries; only a coarse privately learned boundary will be used.",
      call. = FALSE
    )
  }
  if (!is.null(roles$cmt)) {
    missing_cmt <- names(endpoints)[vapply(names(endpoints), function(name) {
      is.null(endpoints[[name]]$cmt) &&
        is.null(public_design$endpoint_cmt[[name]])
    }, logical(1))]
    if (is.null(public_design$dose_cmt) || length(missing_cmt)) {
      stop("Public dose and endpoint CMT values are required when CMT is declared.",
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

.public_input_manifest <- function(roles, endpoints, bounds, public_design,
                                   contribution_limits, budget_allocation) {
  data.frame(
    component = c(
      "roles", "endpoint definitions", "numeric bounds", "schema",
      "protocol/design values", "contribution limits", "budget allocation",
      "generator variability"
    ),
    status = rep("public_user_asserted", 8L),
    note = c(
      "Column semantics and exclusions",
      "DVID, scientific clock, transform rule, units, broad shape",
      "Clipping domains selected without confidential extrema",
      "Names, classes, and factor/category levels",
      "Any supplied regimen and timing grid values",
      "Per-subject rows, doses, occasions, observations, and cells",
      "Fractions of requested epsilon",
      "Subject and residual variability settings"
    ),
    stringsAsFactors = FALSE
  )
}

.sanitize_design_for_release <- function(public_design, id_name) {
  design <- public_design
  index <- match(id_name, .schema_names(design$schema))
  if (!is.na(index) && "factor" %in% design$schema$columns[[index]]$class) {
    # Factor levels on an identifier are identifier values, not schema. Keep
    # the factor class but generate a fresh level set during post-processing.
    design$schema$columns[[index]]$levels <- character()
  }
  if (!is.null(design$defaults)) design$defaults[[id_name]] <- NULL
  design
}

#' Fit a subject-level differentially private PMX population generator
#'
#' This is the only stage that reads source patient data. It deterministically
#' bounds each complete subject contribution, invokes a validated DP backend for
#' every released source-dependent summary, and returns no raw patient records.
#' There is deliberately no fitting seed.
#'
#' @param data Confidential PMX event data (or a public fixture when
#'   `public_source = TRUE`).
#' @param roles Explicit column roles from [pmx_roles()].
#' @param endpoints Named endpoint declarations from [pmx_endpoint()].
#' @param epsilon,delta Explicit requested subject-level privacy parameters.
#' @param bounds Public clipping domains from [pmx_bounds()].
#' @param public_design Public schema and protocol metadata from
#'   [pmx_public_design()].
#' @param contribution_limits Public contribution limits from
#'   [pmx_contribution_limits()].
#' @param budget_allocation Explicit fractions from [pmx_budget_allocation()].
#' @param delta_justification Required justification when `delta > 0`.
#' @param backend Production fitting supports `"opendp"`. `"public"` is a
#'   noiseless structural backend allowed only with `public_source = TRUE` for
#'   public examples; it makes no privacy claim.
#' @param public_source Logical assertion that the complete input is already
#'   public. Never set this for confidential or patient data.
#'
#' @return A `private_pmx_model`. It contains public configuration, noisy
#'   population summaries, privacy accounting, and one release-ledger entry;
#'   it contains no raw IDs, rows, profiles, templates, or residuals.
#' @export
fit_private_pmx <- function(data, roles, endpoints, epsilon, delta, bounds,
                            public_design, contribution_limits,
                            budget_allocation, delta_justification = NULL,
                            backend = "opendp", public_source = FALSE) {
  .validate_fit_configuration(
    data, roles, endpoints, epsilon, delta, delta_justification, bounds,
    public_design, contribution_limits, budget_allocation
  )
  resolved_backend <- .resolve_backend(backend, public_source)

  retained_names <- setdiff(names(data), roles$exclude)
  retained <- data[, retained_names, drop = FALSE]
  retained_roles <- roles
  retained_roles$exclude <- NULL
  class(retained_roles) <- "pmx_roles"

  date_columns <- retained_names[vapply(retained, function(x) {
    inherits(x, "Date") || inherits(x, "POSIXct") || inherits(x, "POSIXlt")
  }, logical(1))]
  if (length(date_columns)) {
    stop(
      "Unmodeled Date/POSIX columns are forbidden: ",
      paste(date_columns, collapse = ", "),
      ". Exclude them explicitly before fitting.", call. = FALSE
    )
  }
  identifiers <- setdiff(.direct_identifier_names(retained_names),
                         retained_roles$id)
  if (length(identifiers)) {
    stop("Possible direct identifier columns must be explicitly excluded: ",
         paste(identifiers, collapse = ", "), ".", call. = FALSE)
  }
  .schema_matches(retained, public_design$schema)
  public_design <- .sanitize_design_for_release(
    public_design, retained_roles$id
  )
  normalized_endpoints <- .normalize_endpoints(
    endpoints, retained_roles, bounds, public_design, contribution_limits
  )
  validate_pmx(retained, retained_roles, endpoints = normalized_endpoints,
               strict = TRUE)

  accountant <- .new_accountant(epsilon, delta, resolved_backend)
  population <- .fit_population_summaries(
    retained, retained_roles, normalized_endpoints, bounds,
    contribution_limits, budget_allocation, epsilon, accountant,
    public_design
  )
  accounting <- .finalize_accounting(accountant)
  backend_record <- list(
    name = resolved_backend$name, version = resolved_backend$version,
    mechanism = resolved_backend$mechanism,
    validated = isTRUE(resolved_backend$validated),
    production = isTRUE(resolved_backend$production)
  )
  formal_dp <- isTRUE(resolved_backend$validated) &&
    isTRUE(resolved_backend$production)
  recorded_delta_justification <- if (formal_dp) {
    delta_justification %||%
      "delta = 0; the OpenDP Laplace releases are pure DP"
  } else {
    "Not applicable: the source was asserted public and no DP claim is made."
  }
  release_id <- .release_id()
  ledger <- list(
    release_id = release_id,
    created_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"),
    privacy_unit = "one subject's complete bounded longitudinal record",
    adjacency = "add-or-remove one complete subject",
    requested_epsilon = epsilon, requested_delta = delta,
    realized_epsilon = accounting$realized_epsilon,
    realized_delta = accounting$realized_delta,
    backend = backend_record$name, backend_version = backend_record$version,
    refit_notice = paste(
      "A new call to fit_private_pmx() against the same confidential data is",
      "a new release and composes with this ledger entry."
    )
  )
  warnings <- character()
  if (population$private_subject_count < 6) {
    warnings <- c(warnings, paste(
      "The private subject-count release is below six; only very broad",
      "workflow utility should be expected."
    ))
  }
  released_dimensions <- sum(accounting$entries$dimensions)
  if (released_dimensions > 6 * population$private_subject_count) {
    warnings <- c(warnings, paste(
      "The requested summary dimension is high relative to the private",
      "subject-count release; reduce grids or strengthen public assumptions."
    ))
  }

  model <- structure(list(
    version = 2L,
    engine = "subject_level_dp_population_generator",
    privacy = list(
      formal_dp = formal_dp,
      unit = "one subject's complete bounded longitudinal contribution",
      adjacency = "add-or-remove one complete subject",
      epsilon = as.numeric(epsilon), delta = as.numeric(delta),
      delta_justification = recorded_delta_justification,
      backend = backend_record,
      accounting = accounting,
      proof_assumptions = c(
        "All schema, roles, endpoint metadata, domains, levels, and protocol values recorded as public were established independently of confidential values.",
        "Every source-dependent aggregate is computed from the deterministically contribution-bounded subject representation.",
        "OpenDP's Laplace measurement and privacy map are correct for the stated L1 sensitivities and finite inputs.",
        "Basic sequential composition covers every released source-dependent computation in this fit.",
        "Generation and validation after fitting consult only this released model and public inputs.",
        "Runtime, serialization, floating-point, side-channel, governance, and public-input assertions require independent review."
      )
    ),
    public = list(
      roles = retained_roles,
      endpoints = normalized_endpoints,
      bounds = bounds,
      design = public_design,
      contribution_limits = contribution_limits,
      budget_allocation = budget_allocation,
      input_manifest = .public_input_manifest(
        retained_roles, normalized_endpoints, bounds, public_design,
        contribution_limits, budget_allocation
      )
    ),
    population = population,
    ledger = ledger,
    warnings = unique(warnings)
  ), class = "private_pmx_model")
  validate_private_model(model, strict = TRUE)
  if (length(warnings)) warning(paste(warnings, collapse = "\n"), call. = FALSE)
  model
}
