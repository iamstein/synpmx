#' Declare pharmacometric column roles
#'
#' Column roles are always explicit. `pmxmock` never guesses critical event or
#' measurement columns from their names.
#'
#' @param id,time,dv,evid Single column names for subject ID, time, dependent
#'   variable, and event ID. These roles are required.
#' @param amt,cmt,dvid,mdv,rate Optional single column names for amount,
#'   compartment, endpoint, missing-DV indicator, and infusion rate.
#' @param covariates Character vector of baseline covariate column names, or
#'   `NULL`.
#'
#' @return A `pmx_roles` object consumed by [mock_pmx()], [validate_pmx()], and
#'   [compare_pmx()].
#' @export
#'
#' @examples
#' roles <- pmx_roles(
#'   id = "ID", time = "TIME", dv = "DV", amt = "AMT",
#'   evid = "EVID", cmt = "CMT", covariates = "WT"
#' )
pmx_roles <- function(id, time, dv, amt = NULL, evid, cmt = NULL,
                      dvid = NULL, mdv = NULL, rate = NULL,
                      covariates = NULL) {
  roles <- list(
    id = id,
    time = time,
    dv = dv,
    amt = amt,
    evid = evid,
    cmt = cmt,
    dvid = dvid,
    mdv = mdv,
    rate = rate,
    covariates = covariates
  )

  scalar_roles <- setdiff(names(roles), "covariates")
  for (role in scalar_roles) {
    value <- roles[[role]]
    if (!is.null(value) &&
        (!is.character(value) || length(value) != 1L || is.na(value) ||
         !nzchar(value))) {
      stop("`", role, "` must be one non-empty column name or NULL.",
           call. = FALSE)
    }
  }
  if (!is.null(covariates) &&
      (!is.character(covariates) || anyNA(covariates) ||
       any(!nzchar(covariates)))) {
    stop("`covariates` must be a character vector of column names or NULL.",
         call. = FALSE)
  }
  roles$covariates <- unique(covariates)

  used <- unlist(roles, use.names = FALSE)
  duplicated_roles <- unique(used[duplicated(used)])
  if (length(duplicated_roles)) {
    stop(
      "A column cannot have multiple roles: ",
      paste(duplicated_roles, collapse = ", "), ".",
      call. = FALSE
    )
  }

  structure(roles, class = "pmx_roles")
}

#' @export
print.pmx_roles <- function(x, ...) {
  cat("Pharmacometric column roles:\n")
  for (role in names(x)) {
    value <- x[[role]]
    cat("  ", role, ": ",
        if (length(value)) paste(value, collapse = ", ") else "<absent>",
        "\n", sep = "")
  }
  invisible(x)
}

.assert_roles <- function(data, roles) {
  if (!inherits(roles, "pmx_roles")) {
    stop("`roles` must be created by `pmx_roles()`.", call. = FALSE)
  }
  required <- c("id", "time", "dv", "evid")
  absent_required <- required[vapply(roles[required], is.null, logical(1))]
  if (length(absent_required)) {
    stop("Required roles are absent: ",
         paste(absent_required, collapse = ", "), ".", call. = FALSE)
  }
  columns <- unlist(roles, use.names = FALSE)
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns)) {
    stop("Role columns not found in `data`: ",
         paste(missing_columns, collapse = ", "), ".", call. = FALSE)
  }
  invisible(TRUE)
}
