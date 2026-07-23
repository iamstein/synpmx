.finish_validation <- function(checks, summary) {
  check_table <- if (length(checks)) do.call(rbind, checks) else data.frame()
  structure(list(
    valid = !nrow(check_table) || !any(check_table$status == "error"),
    checks = check_table, summary = summary
  ), class = "pmx_validation")
}

# Name a role together with the column the user mapped it to, so an error points
# at their data rather than at an abstract role. `time` mapped to "RFSTDTC"
# reads as: the `time` role (column 'RFSTDTC').
.role_col <- function(roles, role) {
  column <- roles[[role]]
  paste0("the `", role, "` role (column '", column, "')")
}

# A count of offending rows plus a concrete example: which row, and the value
# that failed. This is the difference between "TIME must be finite" and knowing
# that row 12 holds the string "2025-03-14T08:00".
.offenders <- function(mask, values, max_examples = 1L) {
  bad <- which(mask)
  n <- length(bad)
  if (!n) return("")
  examples <- utils::head(bad, max_examples)
  shown <- vapply(examples, function(i) {
    value <- values[[i]]
    rendered <- if (is.na(value)) "NA" else
      paste0("\"", format(value, trim = TRUE), "\"")
    paste0("row ", i, " = ", rendered)
  }, character(1))
  paste0(" (", n, " of ", length(values), " rows; e.g. ",
         paste(shown, collapse = ", "), ")")
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
      if (anyNA(id)) paste0(
        "ID is missing in some rows: ", .role_col(roles, "id"),
        .offenders(is.na(id), id), ". Every row must belong to a subject."
      ) else "ID has no missing values.")
  if (!is.numeric(time)) {
    add("time", "error", paste0(
      "TIME must be a numeric column, but ", .role_col(roles, "time"),
      " is ", class(time)[1L], ". A date/time string such as ",
      "\"2025-03-14T08:00\" must be converted to elapsed hours or days first."
    ))
  } else if (any(!is.finite(time))) {
    add("time", "error", paste0(
      "TIME must be finite, but ", .role_col(roles, "time"),
      " has non-finite values", .offenders(!is.finite(time), time), "."
    ))
  } else {
    add("time", "pass", "Actual time is numeric and finite.")
  }
  if (!is.numeric(dv)) {
    add("dv_type", "error", paste0(
      "DV must be a numeric column, but ", .role_col(roles, "dv"), " is ",
      class(dv)[1L], ". Text markers such as \"<LLOQ\", \"BLQ\", or \".\" ",
      "for missing must be converted to numbers (and, for below-limit rows, ",
      "handled through the `cens` role) before synthesis."
    ))
  } else add("dv_type", "pass", "DV is numeric.")
  add("evid_missing", if (anyNA(evid)) "error" else "pass",
      if (anyNA(evid)) paste0(
        "EVID is missing in some rows: ", .role_col(roles, "evid"),
        .offenders(is.na(evid), evid),
        ". Every row must be marked an observation (0) or an event (nonzero)."
      ) else "EVID has no missing values.")

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
    nonfinite <- present & !is.finite(dv)
    add("dv_finite", "error", paste0(
      "Observation rows must have a finite DV, but ", .role_col(roles, "dv"),
      " has non-finite values on observation rows",
      .offenders(nonfinite, dv),
      ". Mark rows with no measurement as events or missing (MDV/EVID) ",
      "rather than leaving DV as Inf or NaN."
    ))
  } else add("dv_finite", "pass", "Observed DVs are finite.")
  add("events", if (any(event)) "pass" else "warning",
      if (any(event)) paste(sum(event), "event rows found.") else
        "No nonzero-EVID event rows were found.")

  if (!is.null(roles$dvid)) {
    missing <- allowed & is.na(data[[roles$dvid]])
    add("endpoint", if (any(missing)) "error" else "pass",
        if (any(missing)) paste0(
          "Every observation must name its endpoint, but ",
          .role_col(roles, "dvid"), " is missing on observation rows",
          .offenders(missing, data[[roles$dvid]]), "."
        ) else "Every observation row has a DVID.")
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

    # Coherence between the flag and the value. The checks above accept any CENS
    # flag with any DV; they never ask whether the two agree. A threshold assay
    # cannot report a left-censored value (DV at the limit) above an ordinary
    # measurement of the same endpoint, nor a right-censored value below one, so
    # a crossing means CENS and DV were produced independently rather than from
    # one latent value. This is `REV-021`.
    limit <- if (is.null(roles$limit)) rep(NA_real_, nrow(data)) else
      suppressWarnings(as.numeric(data[[roles$limit]]))
    endpoint_id <- .endpoint(data, roles)
    incoherent <- FALSE
    for (name in unique(endpoint_id[allowed])) {
      here <- allowed & endpoint_id == name & is.finite(dv) & is.finite(cens)
      uncensored <- dv[here & cens == 0]
      if (!length(uncensored)) next
      # Point left-censoring reports the limit in DV; an uncensored value of the
      # same endpoint cannot sit below that limit without itself being censored.
      left <- dv[here & cens == 1 & !is.finite(limit)]
      if (length(left) && max(left) > min(uncensored)) incoherent <- TRUE
      # Right-censoring is the mirror image.
      right <- dv[here & cens == -1 & !is.finite(limit)]
      if (length(right) && min(right) < max(uncensored)) incoherent <- TRUE
    }
    add("cens_dv_coherence", if (incoherent) "error" else "pass",
        if (incoherent)
          paste("Censored and uncensored observations of an endpoint cross the",
                "censoring boundary: a CENS flag disagrees with its DV.") else
          "Censoring flags are consistent with reported DV values.")
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

  # An identifier-named column you have given a role to is one you have
  # consciously handled, so it is not flagged: a character `NAME` declared as
  # the `dvid` endpoint label is the common case. Only undeclared ones are
  # reported, and as a warning rather than an error, because `synpmx_avatar()`
  # drops undeclared columns by default -- the column will not reach the output
  # unless you also name it in `keep`.
  direct <- setdiff(.direct_identifier_names(names(data)),
                    .retained_role_columns(roles))
  add("direct_identifiers", if (length(direct)) "warning" else "pass",
      if (length(direct)) paste0(
        "Column name(s) look like direct identifiers: ",
        paste(direct, collapse = ", "),
        ". They are undeclared, so `synpmx_avatar()` drops them; declare one in ",
        "`keep` only if you intend to carry a real subject's value through."
      ) else "No undeclared direct-identifier column names remain.")
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
    messages <- report$checks$message[report$checks$status == "error"]
    stop("PMX validation failed with ", length(messages),
         if (length(messages) == 1L) " problem:\n" else " problems:\n",
         paste0("  ", seq_along(messages), ". ", messages, collapse = "\n"),
         "\nEach names the role and the column it maps to; fix the role ",
         "mapping in `pmx_roles()`, or the data.",
         call. = FALSE)
  }
  report
}

#' @export
print.pmx_validation <- function(x, ...) {
  if (isTRUE(x$valid)) {
    cat("Valid PMX structure.\n")
    return(invisible(x))
  }
  errors <- x$checks[x$checks$status == "error", , drop = FALSE]
  warnings <- x$checks[x$checks$status == "warning", , drop = FALSE]
  cat("Invalid PMX structure:", nrow(errors),
      if (nrow(errors) == 1L) "problem" else "problems", "to fix.\n\n")
  for (i in seq_len(nrow(errors))) {
    cat(i, ". ", errors$message[i], "\n", sep = "")
  }
  if (nrow(warnings)) {
    cat("\nWarnings (not fatal):\n")
    for (i in seq_len(nrow(warnings))) {
      cat("- ", warnings$message[i], "\n", sep = "")
    }
  }
  cat("\nEach problem names the role and the column it maps to. Fix the role",
      "mapping\nin `pmx_roles()`, or the data, and re-run.\n")
  invisible(x)
}
