.opendp_backend <- function() {
  if (!requireNamespace("opendp", quietly = TRUE)) return(NULL)
  ns <- asNamespace("opendp")
  get <- function(name) getExportedValue("opendp", name)
  get("enable_features")("contrib")

  release <- function(value, sensitivity, epsilon) {
    if (!is.numeric(value) || any(!is.finite(value)) ||
        !is.numeric(sensitivity) || length(sensitivity) != 1L ||
        !is.finite(sensitivity) || sensitivity < 0 ||
        !is.numeric(epsilon) || length(epsilon) != 1L ||
        !is.finite(epsilon) || epsilon <= 0) {
      stop("Invalid values supplied to the OpenDP adapter.", call. = FALSE)
    }
    scale <- sensitivity / epsilon
    if (scale > 0) {
      # OpenDP checks the privacy relation with exact rational arithmetic.
      # A binary floating-point quotient can round infinitesimally below the
      # required scale and make a mathematically valid fractional allocation
      # fail closed. Inflate by a negligible amount so rounding can only make
      # the mechanism more private, never less private.
      scale <- scale * (1 + sqrt(.Machine$double.eps))
    }
    if (length(value) == 1L) {
      domain <- get("atom_domain")(nan = FALSE, .T = "f64")
      metric <- get("absolute_distance")(.T = "f64")
    } else {
      atom <- get("atom_domain")(nan = FALSE, .T = "f64")
      domain <- get("vector_domain")(atom)
      metric <- get("l1_distance")(.T = "f64")
    }
    measurement <- get("make_laplace")(
      input_domain = domain, input_metric = metric, scale = scale
    )
    if (!isTRUE(get("measurement_check")(
      measurement, as.numeric(sensitivity), as.numeric(epsilon)
    ))) {
      stop("OpenDP rejected the requested sensitivity/privacy relation.",
           call. = FALSE)
    }
    output <- get("measurement_invoke")(measurement, as.numeric(value))
    output <- as.numeric(unlist(output, use.names = FALSE))
    if (length(output) != length(value) || any(!is.finite(output))) {
      stop("OpenDP returned an invalid mechanism output.", call. = FALSE)
    }
    output
  }
  structure(list(
    name = "OpenDP", version = as.character(utils::packageVersion("opendp")),
    mechanism = paste(
      "OpenDP Laplace measurement over finite f64 values",
      "(internally exact-rational/discrete sampling)"
    ),
    production = TRUE, validated = TRUE, release = release
  ), class = "pmx_dp_backend")
}

.public_fixture_backend <- function() {
  structure(list(
    name = "public-fixture", version = as.character(utils::packageVersion(
      "pmxSynthData"
    )), mechanism = "no noise; source explicitly asserted public",
    production = FALSE, validated = FALSE,
    release = function(value, sensitivity, epsilon) as.numeric(value)
  ), class = "pmx_dp_backend")
}

.resolve_backend <- function(backend, public_source) {
  if (is.null(backend)) backend <- "opendp"
  if (!is.character(backend) || length(backend) != 1L || is.na(backend)) {
    stop("`backend` must be `\"opendp\"` or `\"public\"`.", call. = FALSE)
  }
  if (identical(backend, "public")) {
    if (!isTRUE(public_source)) {
      stop(
        "The nonprivate public-fixture backend is allowed only when ",
        "`public_source = TRUE`; it must never process confidential data.",
        call. = FALSE
      )
    }
    return(.public_fixture_backend())
  }
  if (!identical(backend, "opendp")) {
    stop("Unknown DP backend `", backend, "`.", call. = FALSE)
  }
  resolved <- .opendp_backend()
  if (is.null(resolved)) {
    stop(
      "A validated differential-privacy backend is unavailable. Install ",
      "the OpenDP R package; fitting confidential data fails closed and does ",
      "not fall back to ordinary R random noise.", call. = FALSE
    )
  }
  resolved
}

#' Inspect the differential-privacy backend
#'
#' @return A one-row data frame describing whether the production OpenDP
#'   adapter is available.
#' @export
dp_backend_status <- function() {
  backend <- tryCatch(.opendp_backend(), error = function(e) NULL)
  data.frame(
    backend = "OpenDP",
    available = !is.null(backend),
    version = if (is.null(backend)) NA_character_ else backend$version,
    production = !is.null(backend) && isTRUE(backend$production),
    stringsAsFactors = FALSE
  )
}

#' Run canonical checks against the configured DP backend
#'
#' The checks verify construction, invocation, finite output, and OpenDP's own
#' privacy-map relation. They are implementation checks, not an empirical proof
#' of differential privacy.
#'
#' @return A structured check result. The function fails closed if OpenDP is
#'   unavailable.
#' @export
run_dp_backend_tests <- function() {
  backend <- .resolve_backend("opendp", public_source = FALSE)
  scalar <- backend$release(0, sensitivity = 1, epsilon = 1)
  vector <- backend$release(c(0, 1), sensitivity = 2, epsilon = 1)
  fractional <- backend$release(
    c(0, 1, 2), sensitivity = 3, epsilon = 0.15
  )
  structure(list(
    passed = length(scalar) == 1L && length(vector) == 2L &&
      length(fractional) == 3L &&
      all(is.finite(c(scalar, vector, fractional))),
    backend = backend$name,
    version = backend$version,
    tests = c(
      "scalar_laplace", "vector_l1_laplace", "fractional_budget",
      "privacy_map"
    )
  ), class = "pmx_backend_tests")
}

.new_accountant <- function(requested_epsilon, requested_delta, backend) {
  environment <- new.env(parent = emptyenv())
  environment$requested_epsilon <- requested_epsilon
  environment$requested_delta <- requested_delta
  environment$backend <- backend
  environment$entries <- list()
  environment
}

.private_release <- function(accountant, query, value, sensitivity, epsilon) {
  if (epsilon <= 0) {
    stop("Every executed private query needs positive epsilon.", call. = FALSE)
  }
  output <- accountant$backend$release(value, sensitivity, epsilon)
  accountant$entries[[length(accountant$entries) + 1L]] <- data.frame(
    query = query, mechanism = accountant$backend$mechanism,
    epsilon = epsilon, delta = 0, sensitivity = sensitivity,
    dimensions = length(value), stringsAsFactors = FALSE
  )
  output
}

.accounting_table <- function(accountant) {
  if (!length(accountant$entries)) {
    return(data.frame(
      query = character(), mechanism = character(), epsilon = numeric(),
      delta = numeric(), sensitivity = numeric(), dimensions = integer()
    ))
  }
  do.call(rbind, accountant$entries)
}

.finalize_accounting <- function(accountant) {
  table <- .accounting_table(accountant)
  realized_epsilon <- sum(table$epsilon)
  realized_delta <- sum(table$delta)
  if (realized_epsilon > accountant$requested_epsilon ||
      realized_delta > accountant$requested_delta) {
    stop("Internal privacy accounting exceeded the requested budget.",
         call. = FALSE)
  }
  list(
    composition = "basic sequential composition of pure-DP Laplace releases",
    entries = table, realized_epsilon = realized_epsilon,
    realized_delta = realized_delta,
    unspent_epsilon = accountant$requested_epsilon - realized_epsilon,
    unspent_delta = accountant$requested_delta - realized_delta
  )
}

#' Validate a fitted private PMX population model
#'
#' Checks structural integrity, accounting, prohibited payload names, backend
#' status, and the absence of direct patient records. It cannot independently
#' prove that runtime configuration or the external backend was honest.
#'
#' @param private_model An object returned by [fit_private_pmx()].
#' @param strict Stop on any failed check when `TRUE`.
#'
#' @return A `pmx_private_validation` report.
#' @export
validate_private_model <- function(private_model, strict = FALSE) {
  checks <- data.frame(check = character(), status = character(),
                       message = character(), stringsAsFactors = FALSE)
  add <- function(check, status, message) {
    checks <<- rbind(checks, data.frame(
      check = check, status = status, message = message,
      stringsAsFactors = FALSE
    ))
  }
  if (!inherits(private_model, "private_pmx_model")) {
    add("class", "error", "Object is not a fitted private PMX model.")
  } else {
    add("class", "pass", "Private-model class and version marker found.")
    privacy <- private_model$privacy
    entries <- privacy$accounting$entries
    required_accounting <- c(
      "query", "mechanism", "epsilon", "delta", "sensitivity", "dimensions"
    )
    accounting_ok <- is.data.frame(entries) && nrow(entries) > 0L &&
      all(required_accounting %in% names(entries)) &&
      all(vapply(
        entries[c("epsilon", "delta", "sensitivity", "dimensions")],
        is.numeric, logical(1)
      )) &&
      all(is.finite(entries$epsilon)) && all(entries$epsilon > 0) &&
      all(is.finite(entries$delta)) && all(entries$delta >= 0) &&
      all(is.finite(entries$sensitivity)) && all(entries$sensitivity > 0) &&
      all(is.finite(entries$dimensions)) && all(entries$dimensions >= 1) &&
      all(entries$dimensions == floor(entries$dimensions)) &&
      identical(as.numeric(sum(entries$epsilon)),
                as.numeric(privacy$accounting$realized_epsilon)) &&
      identical(as.numeric(sum(entries$delta)),
                as.numeric(privacy$accounting$realized_delta)) &&
      is.numeric(privacy$epsilon) && length(privacy$epsilon) == 1L &&
      is.finite(privacy$epsilon) && privacy$epsilon > 0 &&
      is.numeric(privacy$delta) && length(privacy$delta) == 1L &&
      is.finite(privacy$delta) && privacy$delta >= 0 && privacy$delta < 1 &&
      privacy$accounting$realized_epsilon <= privacy$epsilon &&
      privacy$accounting$realized_delta <= privacy$delta
    if (!accounting_ok) {
      add("accounting", "error", "Privacy accounting is missing or exceeds budget.")
    } else {
      add("accounting", "pass", "Composed accounting is within budget.")
    }
    banned <- c(
      "raw_rows", "raw_data", "source_ids", "subject_profiles",
      "event_templates", "raw_residuals", "unnoised_aggregates",
      "anchors", "donors"
    )
    found <- intersect(tolower(.recursive_names(private_model)), banned)
    if (length(found)) {
      add("payload", "error", paste("Prohibited payload names:",
                                     paste(found, collapse = ", ")))
    } else {
      add("payload", "pass", "No prohibited raw-data payload names found.")
    }
    if (!is.list(private_model$population) ||
        !is.list(private_model$public)) {
      add("model", "error", "Public metadata or private population summaries are missing.")
    } else {
      add("model", "pass", "Only public metadata and released summaries are retained.")
    }
    if (!isTRUE(privacy$formal_dp)) {
      add("formal_dp", "warning",
          "This model used the explicitly public-source fixture backend; no DP claim is made or needed for that public source.")
    } else if (!isTRUE(privacy$backend$validated) ||
               !isTRUE(privacy$backend$production) ||
               !identical(privacy$backend$name, "OpenDP")) {
      add("formal_dp", "error",
          "A formal claim requires the validated production OpenDP backend.")
    } else {
      add("formal_dp", "pass", "A validated production backend was recorded.")
    }
  }
  valid <- !any(checks$status == "error")
  report <- structure(list(valid = valid, checks = checks),
                      class = "pmx_private_validation")
  if (strict && !valid) {
    stop("Private-model validation failed: ",
         paste(checks$message[checks$status == "error"], collapse = " "),
         call. = FALSE)
  }
  report
}

#' Summarize a fitted model's privacy contract
#'
#' @param private_model A fitted model from [fit_private_pmx()].
#'
#' @return A machine-readable `pmx_privacy_report` list.
#' @export
privacy_report <- function(private_model) {
  validate_private_model(private_model, strict = TRUE)
  privacy <- private_model$privacy
  structure(list(
    release_id = private_model$ledger$release_id,
    guarantee = if (isTRUE(privacy$formal_dp)) {
      paste0("Generated from a subject-level (", privacy$epsilon, ", ",
             privacy$delta, ")-differentially private model.")
    } else {
      "No DP claim: the input was explicitly asserted to be a public fixture."
    },
    formal_dp = privacy$formal_dp,
    privacy_unit = privacy$unit,
    adjacency = privacy$adjacency,
    epsilon = privacy$epsilon,
    delta = privacy$delta,
    delta_justification = privacy$delta_justification,
    contribution_limits = private_model$public$contribution_limits,
    backend = privacy$backend,
    accounting = privacy$accounting,
    public_inputs = private_model$public$input_manifest,
    proof_assumptions = privacy$proof_assumptions,
    post_processing = paste(
      "Additional generate_pmx() calls from this fitted model consume no",
      "additional privacy budget. Refitting against source data does."
    ),
    qualification = if (isTRUE(privacy$formal_dp)) {
      paste(
        "The protection is mathematically bounded rather than absolute; it is",
        "not a legal anonymity or release-authorization determination."
      )
    } else {
      paste(
        "No privacy guarantee is asserted for this public-source fixture",
        "model."
      )
    }
  ), class = "pmx_privacy_report")
}

#' @export
print.pmx_privacy_report <- function(x, ...) {
  accounting_label <- if (isTRUE(x$formal_dp)) {
    "Realized DP accounting"
  } else {
    "Illustrative query allocation (not a DP accounting claim)"
  }
  cat(x$guarantee, "\n",
      "Privacy unit: ", x$privacy_unit, "\n",
      "Adjacency: ", x$adjacency, "\n",
      "Backend: ", x$backend$name, " ", x$backend$version, "\n",
      accounting_label, ": epsilon = ",
      x$accounting$realized_epsilon, ", delta = ",
      x$accounting$realized_delta, "\n",
      x$qualification, "\n", sep = "")
  invisible(x)
}
