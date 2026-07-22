#' Validate a pharmacometric event dataset
#'
#' Performs structural checks using explicit roles. Validation concerns event
#' coherence, schema usability, ordering, endpoint fields, and baseline
#' constancy; it does not assess scientific model fit or distributional
#' equivalence.
#'
#' @param data A data frame containing PMX event records.
#' @param roles A role mapping from [pmx_roles()].
#' @param strict If `TRUE`, stop when any error-level check fails.
#'
#' @return A `pmx_validation` list with `valid`, `checks`, and `summary`
#'   components.
#' @export
#'
#' @examples
#' data <- data.frame(
#'   ID = c(1L, 1L), TIME = c(0, 1), DV = c(0, 2.1),
#'   AMT = c(100, 0), EVID = c(1L, 0L), WT = c(70, 70)
#' )
#' roles <- pmx_roles("ID", "TIME", "DV", "AMT", "EVID",
#'                    covariates = "WT")
#' validate_pmx(data, roles)
validate_pmx <- function(data, roles, strict = FALSE) {
  checks <- list()
  add_check <- function(check, status, message) {
    checks[[length(checks) + 1L]] <<- data.frame(
      check = check,
      status = status,
      message = message,
      stringsAsFactors = FALSE
    )
  }

  if (!is.data.frame(data)) {
    add_check("data_frame", "error", "`data` is not a data frame or tibble.")
    report <- .finish_validation(checks, list())
    if (strict) stop(report$checks$message[1L], call. = FALSE)
    return(report)
  }
  if (!nrow(data)) {
    add_check("rows", "error", "The dataset has no rows.")
  } else {
    add_check("rows", "pass", paste(nrow(data), "rows found."))
  }
  if (anyDuplicated(names(data))) {
    add_check("column_names", "error", "Column names are not unique.")
  } else {
    add_check("column_names", "pass", "Column names are unique.")
  }

  role_error <- tryCatch({
    .assert_roles(data, roles)
    NULL
  }, error = function(error) conditionMessage(error))
  if (!is.null(role_error)) {
    add_check("roles", "error", role_error)
    report <- .finish_validation(checks, list(rows = nrow(data)))
    if (strict) stop(role_error, call. = FALSE)
    return(report)
  }
  add_check("roles", "pass", "All explicit role columns are available.")

  id <- data[[roles$id]]
  time <- data[[roles$time]]
  dv <- data[[roles$dv]]
  evid <- data[[roles$evid]]
  if (anyNA(id)) {
    add_check("id_missing", "error", "ID contains missing values.")
  } else {
    add_check("id_missing", "pass", "ID has no missing values.")
  }
  if (!is.numeric(time)) {
    add_check("time_type", "error", "TIME must be numeric.")
  } else if (any(!is.finite(time))) {
    add_check("time_finite", "error", "TIME contains non-finite values.")
  } else {
    add_check("time_finite", "pass", "TIME is numeric and finite.")
  }
  if (!is.numeric(dv)) {
    add_check("dv_type", "error", "DV must be numeric.")
  }
  if (anyNA(evid)) {
    add_check("evid_missing", "error", "EVID contains missing values.")
  } else {
    add_check("evid_missing", "pass", "EVID has no missing values.")
  }

  id_order <- .unique_in_order(id)
  badly_ordered <- vapply(id_order, function(subject) {
    if (is.na(subject)) return(FALSE)
    subject_time <- time[!is.na(id) & id == subject]
    is.numeric(subject_time) && any(diff(subject_time) < 0, na.rm = TRUE)
  }, logical(1))
  if (any(badly_ordered)) {
    add_check(
      "row_order", "error",
      paste("TIME decreases within", sum(badly_ordered), "subject(s).")
    )
  } else {
    add_check("row_order", "pass", "Rows are nondecreasing in TIME by subject.")
  }

  observation_candidate <- .is_zero(evid)
  observation_allowed <- .observation_rows(data, roles)
  observation_present <- observation_allowed & !is.na(dv)
  event_rows <- !is.na(evid) & !.is_zero(evid)
  if (!any(observation_allowed)) {
    add_check("observations", "error", "No genuine observation rows were found.")
  } else {
    add_check(
      "observations", "pass",
      paste(sum(observation_present), "non-missing observations found.")
    )
  }
  missing_observation <- observation_candidate & !observation_present
  if (any(missing_observation)) {
    add_check(
      "missing_observations", "warning",
      paste(sum(missing_observation),
            "EVID-zero rows carry missing DV or a nonzero MDV convention.")
    )
  } else {
    add_check("missing_observations", "pass",
              "No missing observation rows were found.")
  }
  if (!any(event_rows)) {
    add_check("events", "warning", "No nonzero-EVID event rows were found.")
  } else {
    add_check("events", "pass", paste(sum(event_rows), "event rows found."))
  }
  if (is.numeric(dv) && any(!is.finite(dv[observation_present]))) {
    add_check("dv_finite", "error", "Observed DV contains non-finite values.")
  } else {
    add_check("dv_finite", "pass", "All non-missing observed DVs are finite.")
  }

  if (!is.null(roles$dvid)) {
    missing_endpoint <- observation_allowed & is.na(data[[roles$dvid]])
    if (any(missing_endpoint)) {
      add_check(
        "endpoint", "error",
        paste(sum(missing_endpoint), "observation rows have missing DVID.")
      )
    } else {
      add_check("endpoint", "pass", "Every observation row has an endpoint.")
    }
  }

  if (!is.null(roles$amt)) {
    amt <- data[[roles$amt]]
    nonzero_observation_amt <- observation_allowed & !is.na(amt) & amt != 0
    if (any(nonzero_observation_amt)) {
      add_check(
        "observation_amount", "warning",
        paste(sum(nonzero_observation_amt),
              "observation rows have nonzero AMT; verify this convention.")
      )
    } else {
      add_check("observation_amount", "pass",
                "Observation rows do not carry nonzero AMT.")
    }
  }

  if (!is.null(roles$mdv)) {
    mdv <- data[[roles$mdv]]
    if (anyNA(mdv)) {
      add_check("mdv", "warning", "MDV contains missing values.")
    } else {
      add_check("mdv", "pass", "MDV is complete.")
    }
  }

  for (covariate in roles$covariates) {
    constant <- vapply(id_order, function(subject) {
      if (is.na(subject)) return(TRUE)
      value <- data[[covariate]][!is.na(id) & id == subject]
      length(unique(value[!is.na(value)])) <= 1L
    }, logical(1))
    if (all(constant)) {
      add_check(paste0("covariate_", covariate), "pass",
                paste(covariate, "is constant within subject."))
    } else {
      add_check(
        paste0("covariate_", covariate), "error",
        paste(covariate, "varies within", sum(!constant), "subject(s).")
      )
    }
  }

  endpoint <- .endpoint(data, roles)
  summary <- list(
    rows = nrow(data),
    subjects = length(unique(id)),
    event_rows = sum(event_rows),
    observation_rows = sum(observation_candidate),
    observed_dv = sum(observation_present),
    missing_observations = sum(missing_observation),
    endpoints = sort(unique(endpoint[observation_allowed]))
  )
  report <- .finish_validation(checks, summary)
  if (strict && !report$valid) {
    errors <- report$checks$message[report$checks$status == "error"]
    stop("PMX validation failed: ", paste(errors, collapse = " "),
         call. = FALSE)
  }
  report
}

.finish_validation <- function(checks, summary) {
  check_table <- if (length(checks)) do.call(rbind, checks) else data.frame()
  structure(
    list(
      valid = !nrow(check_table) || !any(check_table$status == "error"),
      checks = check_table,
      summary = summary
    ),
    class = "pmx_validation"
  )
}

#' @export
print.pmx_validation <- function(x, ...) {
  cat(if (isTRUE(x$valid)) "Valid PMX structure" else "Invalid PMX structure",
      "\n")
  if (nrow(x$checks)) {
    print(x$checks, row.names = FALSE)
  }
  invisible(x)
}
