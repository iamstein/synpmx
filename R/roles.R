#' Declare pharmacometric column roles
#'
#' Column roles are explicit: `synpmx` does not infer critical PMX
#' semantics from column names. Columns listed in `exclude` are removed before
#' fitting and do not appear in generated data.
#'
#' @param id,time,dv,evid Required single column names for subject ID, actual
#'   time, dependent variable, and event ID.
#' @param amt,cmt,dvid,mdv,rate Optional single column names for amount,
#'   compartment, endpoint, missing-DV indicator, and infusion rate.
#' @param nominal_time,tad,occasion Optional time metadata columns.
#' @param cens,limit Optional Monolix-style censoring indicator and other
#'   interval-boundary columns.
#' @param addl,ii Optional additional-dose and interdose-interval columns.
#' @param covariates Baseline covariate column names, or `NULL`.
#' @param subject_properties Subject-level assignment or grouping columns, such
#'   as `ACTARM`, `TRT`, or a nominal dose group. These are treated as
#'   categorical, must be constant and nonmissing within a subject, and are
#'   modeled jointly with that subject's regimen rather than as independent
#'   baseline covariates.
#' @param assigned_dose Optional nominal assigned-dose column. It may vary by
#'   occasion but must be constant within subject and occasion and agree with
#'   the positive event amount. Generated values are derived from the generated
#'   regimen rather than sampled independently.
#' @param exclude Columns explicitly excluded before private fitting, such as
#'   direct identifiers. An ID role is still required as the privacy unit.
#'
#' @return A `pmx_roles` object used by the fitting, generation, validation, and
#'   comparison functions.
#' @export
#'
#' @examples
#' roles <- pmx_roles(
#'   id = "ID", time = "TIME", dv = "DV", amt = "AMT",
#'   evid = "EVID", cmt = "CMT", tad = "TAD", covariates = "WT"
#' )
pmx_roles <- function(id, time, dv, amt = NULL, evid, cmt = NULL,
                      dvid = NULL, mdv = NULL, rate = NULL,
                      nominal_time = NULL, tad = NULL, occasion = NULL,
                      cens = NULL, limit = NULL, addl = NULL, ii = NULL,
                      covariates = NULL, subject_properties = NULL,
                      assigned_dose = NULL, exclude = NULL) {
  roles <- list(
    id = id, time = time, nominal_time = nominal_time, tad = tad,
    occasion = occasion, dv = dv, amt = amt, evid = evid, cmt = cmt,
    dvid = dvid, mdv = mdv, rate = rate, cens = cens, limit = limit,
    addl = addl, ii = ii, assigned_dose = assigned_dose,
    covariates = covariates, subject_properties = subject_properties,
    exclude = exclude
  )

  vector_roles <- c("covariates", "subject_properties", "exclude")
  scalar_roles <- setdiff(names(roles), vector_roles)
  for (role in scalar_roles) {
    value <- roles[[role]]
    if (!is.null(value) &&
        (!is.character(value) || length(value) != 1L || is.na(value) ||
         !nzchar(value))) {
      stop("`", role, "` must be one non-empty column name or NULL.",
           call. = FALSE)
    }
  }
  for (role in vector_roles) {
    value <- roles[[role]]
    if (!is.null(value) &&
        (!is.character(value) || anyNA(value) || any(!nzchar(value)))) {
      stop("`", role,
           "` must be a character vector of column names or NULL.",
           call. = FALSE)
    }
    roles[[role]] <- unique(value)
  }

  modeled <- unlist(roles[setdiff(names(roles), "exclude")],
                    use.names = FALSE)
  duplicated_roles <- unique(modeled[duplicated(modeled)])
  if (length(duplicated_roles)) {
    stop("A column cannot have multiple roles: ",
         paste(duplicated_roles, collapse = ", "), ".", call. = FALSE)
  }
  overlap <- intersect(modeled, roles$exclude)
  if (length(overlap)) {
    stop("Modeled role columns cannot also be excluded: ",
         paste(overlap, collapse = ", "), ".", call. = FALSE)
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
  columns <- unlist(roles[setdiff(names(roles), "exclude")],
                    use.names = FALSE)
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns)) {
    stop("Role columns not found in `data`: ",
         paste(missing_columns, collapse = ", "), ".", call. = FALSE)
  }
  invisible(TRUE)
}
