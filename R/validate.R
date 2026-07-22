.finish_validation <- function(checks, summary) {
  check_table <- if (length(checks)) do.call(rbind, checks) else data.frame()
  structure(list(
    valid = !nrow(check_table) || !any(check_table$status == "error"),
    checks = check_table, summary = summary
  ), class = "pmx_validation")
}

#' Validate a pharmacometric event dataset
#'
#' Checks schema usability, chronological event logic, explicit endpoint
#' semantics, derived timing fields, censoring conventions, baseline
#' constancy, subject properties, and occasion-assigned dose coherence. It does
#' not assess scientific or inferential validity.
#'
#' @param data A PMX event data frame or tibble.
#' @param roles Explicit roles from [pmx_roles()].
#' @param endpoints Optional named endpoint declarations from [pmx_endpoint()].
#' @param strict Stop when any error-level check fails.
#'
#' @return A `pmx_validation` report with `valid`, `checks`, and `summary`.
#' @export
validate_pmx <- function(data, roles, endpoints = NULL, strict = FALSE) {
  checks <- list()
  add <- function(check, status, message) {
    checks[[length(checks) + 1L]] <<- data.frame(
      check = check, status = status, message = message,
      stringsAsFactors = FALSE
    )
  }
  if (!is.data.frame(data)) {
    add("data_frame", "error", "`data` is not a data frame or tibble.")
    report <- .finish_validation(checks, list())
    if (strict) stop(report$checks$message[1L], call. = FALSE)
    return(report)
  }
  add("rows", if (nrow(data)) "pass" else "error",
      if (nrow(data)) paste(nrow(data), "rows found.") else "No rows found.")
  add("column_names", if (anyDuplicated(names(data))) "error" else "pass",
      if (anyDuplicated(names(data))) "Column names are not unique." else
        "Column names are unique.")
  role_error <- tryCatch({ .assert_roles(data, roles); NULL },
                         error = function(e) conditionMessage(e))
  if (!is.null(role_error)) {
    add("roles", "error", role_error)
    report <- .finish_validation(checks, list(rows = nrow(data)))
    if (strict) stop(role_error, call. = FALSE)
    return(report)
  }
  add("roles", "pass", "All explicit role columns are available.")

  id <- data[[roles$id]]
  time <- data[[roles$time]]
  dv <- data[[roles$dv]]
  evid <- data[[roles$evid]]
  add("id_missing", if (anyNA(id)) "error" else "pass",
      if (anyNA(id)) "ID contains missing values." else
        "ID has no missing values.")
  if (!is.numeric(time) || any(!is.finite(time))) {
    add("time", "error", "Actual time must be numeric and finite.")
  } else {
    add("time", "pass", "Actual time is numeric and finite.")
  }
  if (!is.numeric(dv)) add("dv_type", "error", "DV must be numeric.") else
    add("dv_type", "pass", "DV is numeric.")
  add("evid_missing", if (anyNA(evid)) "error" else "pass",
      if (anyNA(evid)) "EVID contains missing values." else
        "EVID has no missing values.")

  subjects <- .unique_in_order(id)
  decreasing <- vapply(subjects, function(subject) {
    rows <- which(!is.na(id) & id == subject)
    if (!is.numeric(time)) return(TRUE)
    if (is.null(roles$occasion)) {
      return(any(diff(time[rows]) < -1e-10, na.rm = TRUE))
    }
    occasion <- data[[roles$occasion]][rows]
    groups <- split(rows, occasion, drop = TRUE)
    any(vapply(groups, function(index) {
      any(diff(time[index]) < -1e-10, na.rm = TRUE)
    }, logical(1)))
  }, logical(1))
  add("row_order", if (any(decreasing)) "error" else "pass",
      if (any(decreasing)) paste(
        "TIME decreases within", sum(decreasing),
        if (is.null(roles$occasion)) "subject(s)." else
          "subject/occasion profile(s)."
      ) else if (is.null(roles$occasion)) {
        "Rows are nondecreasing in actual time within subject."
      } else {
        "Rows are nondecreasing within each subject and occasion."
      })

  allowed <- .observation_rows(data, roles)
  present <- allowed & !is.na(dv)
  event <- .event_rows(data, roles)
  if (!any(allowed)) add("observations", "error", "No observation rows found.") else
    add("observations", "pass", paste(sum(present), "observed DVs found."))
  if (is.numeric(dv) && any(!is.finite(dv[present]))) {
    add("dv_finite", "error", "Observed DV contains non-finite values.")
  } else add("dv_finite", "pass", "Observed DVs are finite.")
  add("events", if (any(event)) "pass" else "warning",
      if (any(event)) paste(sum(event), "event rows found.") else
        "No nonzero-EVID event rows were found.")

  if (!is.null(roles$dvid)) {
    missing <- allowed & is.na(data[[roles$dvid]])
    add("endpoint", if (any(missing)) "error" else "pass",
        if (any(missing)) paste(sum(missing),
                               "observation rows have missing DVID.") else
          "Every observation row has a DVID.")
  }
  if (!is.null(endpoints)) {
    endpoint_name <- .endpoint_name_for_rows(data, roles, endpoints)
    unknown <- allowed & is.na(endpoint_name)
    add("endpoint_declaration", if (any(unknown)) "error" else "pass",
        if (any(unknown)) paste(sum(unknown),
                               "observations use undeclared DVID values.") else
          "All observation endpoints are explicitly declared.")
  }

  if (!is.null(roles$nominal_time)) {
    nominal <- data[[roles$nominal_time]]
    if (!is.numeric(nominal) || any(!is.finite(nominal))) {
      add("nominal_time", "error", "Nominal time must be numeric and finite.")
    } else add("nominal_time", "pass", "Nominal time is numeric and finite.")
  }
  if (!is.null(roles$tad)) {
    tad <- suppressWarnings(as.numeric(data[[roles$tad]]))
    bad <- allowed & (!is.finite(tad) | tad < -1e-8)
    add("tad", if (any(bad)) "error" else "pass",
        if (any(bad)) "Observation TAD is missing, non-finite, or negative." else
          "Observation TAD is finite and nonnegative.")
  }
  if (!is.null(roles$occasion)) {
    occasion <- suppressWarnings(as.numeric(data[[roles$occasion]]))
    bad <- allowed & (!is.finite(occasion) | occasion < 1)
    add("occasion", if (any(bad)) "error" else "pass",
        if (any(bad)) "Observation occasion must be a positive value." else
          "Observation occasion values are positive.")
  }

  if (!is.null(roles$cens)) {
    cens <- suppressWarnings(as.numeric(as.character(data[[roles$cens]])))
    bad_value <- is.na(cens) | !(cens %in% c(-1, 0, 1))
    add("cens_values", if (any(bad_value)) "error" else "pass",
        if (any(bad_value)) "CENS must contain only -1, 0, or 1." else
          "CENS uses only the supported Monolix-style states.")
    bad_event <- event & !is.na(cens) & cens != 0
    add("cens_events", if (any(bad_event)) "error" else "pass",
        if (any(bad_event)) "Event rows cannot be censored." else
          "Event rows are uncensored.")
    censored <- allowed & cens != 0
    if (any(censored & !is.finite(dv))) {
      add("cens_boundary", "error", "Censored DV boundaries must be finite.")
    } else add("cens_boundary", "pass", "Censored DV boundaries are finite.")
    if (!is.null(roles$limit)) {
      limit <- suppressWarnings(as.numeric(data[[roles$limit]]))
      left_interval <- allowed & cens == 1 & is.finite(limit)
      right_interval <- allowed & cens == -1 & is.finite(limit)
      bad_direction <- (left_interval & limit > dv) |
        (right_interval & dv > limit)
      add("cens_limit_direction", if (any(bad_direction)) "error" else "pass",
          if (any(bad_direction))
            "Censoring interval boundaries have an invalid direction." else
            "Censoring interval boundaries are ordered consistently.")
    }
  }

  if (!is.null(roles$amt)) {
    amount <- suppressWarnings(as.numeric(data[[roles$amt]]))
    nonzero_observation <- allowed & is.finite(amount) & amount != 0
    add("observation_amount", if (any(nonzero_observation)) "warning" else "pass",
        if (any(nonzero_observation))
          "Some observation rows carry nonzero AMT; verify the convention." else
          "Observation rows carry zero or missing AMT.")
  }
  if (!is.null(roles$assigned_dose)) {
    assigned <- suppressWarnings(as.numeric(data[[roles$assigned_dose]]))
    occasion <- suppressWarnings(as.numeric(data[[roles$occasion]]))
    complete <- is.finite(assigned)
    add(
      "assigned_dose_complete",
      if (all(complete)) "pass" else "error",
      if (all(complete)) "Assigned dose is present on every row." else
        paste(sum(!complete), "row(s) have missing assigned dose.")
    )
    group <- interaction(id, occasion, drop = TRUE)
    constant <- vapply(split(assigned, group), function(value) {
      length(unique(value[is.finite(value)])) == 1L
    }, logical(1))
    add(
      "assigned_dose_constant",
      if (all(constant)) "pass" else "error",
      if (all(constant)) {
        "Assigned dose is constant within subject and occasion."
      } else {
        paste("Assigned dose varies within", sum(!constant),
              "subject/occasion profile(s).")
      }
    )
    positive_event <- event & is.finite(amount) & amount > 0
    coherent <- !positive_event | (
      is.finite(assigned) &
        abs(assigned - amount) <= 1e-8 * pmax(1, abs(amount))
    )
    add(
      "assigned_dose_event",
      if (all(coherent)) "pass" else "error",
      if (all(coherent)) "Assigned dose agrees with positive event AMT." else
        paste(sum(!coherent), "event row(s) disagree with assigned dose.")
    )
  }
  if (!is.null(roles$mdv)) {
    mdv <- data[[roles$mdv]]
    inconsistent <- (.is_zero(mdv) != allowed)
    add("mdv", if (any(inconsistent)) "error" else "pass",
        if (any(inconsistent)) "MDV is inconsistent with observation eligibility." else
          "MDV is consistent with observation eligibility.")
  }

  for (covariate in roles$covariates) {
    constant <- vapply(subjects, function(subject) {
      value <- data[[covariate]][!is.na(id) & id == subject]
      length(unique(value[!is.na(value)])) <= 1L
    }, logical(1))
    add(paste0("covariate_", covariate),
        if (all(constant)) "pass" else "error",
        if (all(constant)) paste(covariate, "is constant within subject.") else
          paste(covariate, "varies within", sum(!constant), "subject(s)."))
  }
  for (property in roles$subject_properties) {
    complete <- vapply(subjects, function(subject) {
      value <- data[[property]][!is.na(id) & id == subject]
      length(value) > 0L && all(!is.na(value))
    }, logical(1))
    constant <- vapply(subjects, function(subject) {
      value <- data[[property]][!is.na(id) & id == subject]
      length(unique(value[!is.na(value)])) == 1L
    }, logical(1))
    add(
      paste0("subject_property_", property),
      if (all(complete & constant)) "pass" else "error",
      if (all(complete & constant)) {
        paste(property, "is complete and constant within subject.")
      } else {
        paste(
          property, "is missing or varies within",
          sum(!(complete & constant)), "subject(s)."
        )
      }
    )
  }

  direct <- setdiff(.direct_identifier_names(names(data)), roles$id)
  add("direct_identifiers", if (length(direct)) "error" else "pass",
      if (length(direct)) paste("Possible direct identifiers:",
                                paste(direct, collapse = ", ")) else
        "No obvious direct-identifier column names remain.")
  endpoint_summary <- if (is.null(roles$dvid)) "DV" else
    sort(unique(as.character(data[[roles$dvid]][allowed])))
  summary <- list(
    rows = nrow(data), subjects = length(unique(id)),
    event_rows = sum(event), observation_rows = sum(allowed),
    observed_dv = sum(present), censored_observations = if (is.null(roles$cens))
      0L else sum(suppressWarnings(as.numeric(as.character(
        data[[roles$cens]][allowed]
      ))) != 0, na.rm = TRUE),
    endpoints = endpoint_summary
  )
  report <- .finish_validation(checks, summary)
  if (strict && !report$valid) {
    stop("PMX validation failed: ", paste(
      report$checks$message[report$checks$status == "error"], collapse = " "
    ), call. = FALSE)
  }
  report
}

#' @export
print.pmx_validation <- function(x, ...) {
  cat(if (isTRUE(x$valid)) "Valid PMX structure" else "Invalid PMX structure",
      "\n")
  print(x$checks, row.names = FALSE)
  invisible(x)
}
