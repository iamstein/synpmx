#' Declare pharmacometric column roles
#'
#' Column roles are explicit: `synpmx` does not infer critical PMX
#' semantics from column names. The declaration is also the complete manifest of
#' what survives into synthetic data. [synpmx_avatar()] drops every column not
#' named by some role, so a column you forget is dropped rather than silently
#' copied out of a real subject. Name a column in `keep` to carry it through.
#'
#' @param id,time,dv,evid Required single column names for subject ID, actual
#'   time, dependent variable, and event ID.
#' @param amt,cmt,mdv,rate Optional single column names for amount,
#'   compartment, missing-DV indicator, and infusion rate.
#' @param dvid Endpoint-key column(s). Usually one column. A dataset that labels
#'   the same endpoint two ways — a numeric `YTYPE` beside a character `NAME` —
#'   may declare both, `dvid = c("YTYPE", "NAME")`. The first is the grouping
#'   key; validation checks the rest are a consistent 1:1 mapping with it and
#'   errors if they disagree, and [synpmx_avatar()] carries all of them through.
#'
#'   **`dvid` and `cmt` answer different questions.** `dvid` is which endpoint a
#'   measurement is; `cmt` is which compartment a dose enters, and it is read
#'   only on event rows. Nothing infers one from the other, so a source with
#'   more than one endpoint needs `dvid` — without it every measurement is
#'   treated as one endpoint, sharing a single value transform and a single
#'   censoring boundary. [synpmx_avatar()] refuses rather than let that happen
#'   silently.
#'
#'   A NONMEM `CMT` usually does both jobs, so **one column may be named as both
#'   roles**: `pmx_roles(..., cmt = "CMT", dvid = "CMT")`. This is the only
#'   permitted overlap; every other collision is an error.
#' @param nominal_time,tad,occasion Optional time metadata columns.
#' @param cens,limit Optional Monolix-style censoring indicator and other
#'   interval-boundary columns.
#' @param addl,ii Optional additional-dose and interdose-interval columns.
#' @param covariates Baseline covariate column names, or `NULL`.
#' @param subject_properties Treatment arm, dose group, cohort — any **assigned,
#'   subject-level stratum**, as opposed to a measured characteristic, which is a
#'   `covariate`. Must be constant within subject; subjects with no recorded
#'   value are grouped as their own stratum, with a warning rather than an
#'   error. Several columns
#'   may be named, and their combination defines the stratum.
#'
#'   [synpmx_avatar()] carries these verbatim from the subject that supplied the
#'   event skeleton and uses the stratum to group two things that are protocol
#'   properties rather than patient properties: the dose-to-covariate
#'   relationship, and the pool of attendance patterns an avatar may draw from.
#'   It is **not** a blending barrier — only route of administration is (see [synpmx_avatar()]), so donors
#'   are still borrowed across strata to reach the donor floor.
#'
#'   The differential-privacy engines model the same columns jointly with the
#'   regimen as a released category domain.
#' @param dose_covariate The covariate the dose is a fixed multiple of, for
#'   weight-based or body-surface-area dosing: name the **covariate**, not the
#'   amount. Declaring it tells [synpmx_avatar()] to recompute each avatar's
#'   `amt` from the avatar's own blended covariate at the anchor's own
#'   milligrams-per-unit, instead of copying the anchor's milligrams.
#'
#'   This matters twice over. An avatar's covariates are blended while its `amt`
#'   would otherwise be copied verbatim, so under proportional dosing the
#'   avatar's own mg/kg comes out wrong — every generated patient violates the
#'   protocol it claims to follow. And a copied amount is one real patient's
#'   real dose, which under proportional dosing discloses that patient's weight
#'   exactly.
#'
#'   Left `NULL`, [synpmx_avatar()] tries to *infer* the relationship, and that
#'   inference deliberately fails closed: it requires the dose-to-covariate
#'   ratio to collapse onto a handful of protocol levels, so a study that
#'   rounds doses to vial sizes, or escalates within a patient by ratios that
#'   are not quite equal, is refused and the amounts are left alone. Declaring
#'   the covariate skips inference entirely and holds each dose row's own
#'   ratio, so **intra-patient escalation is preserved exactly** — three doses
#'   at 1, 2 and 4 mg/kg stay at 1, 2 and 4 mg/kg. The run report says which
#'   path was taken and, when inference declined, why.
#' @param assigned_dose Differential-privacy engines only. A nominal
#'   assigned-dose column reconstructed from the generated regimen.
#'   [synpmx_avatar()] does not use this — carry the column with `keep`.
#' @param keep Columns to carry into synthetic data verbatim, copied from the
#'   same source subject that supplied the event skeleton, with no blending or
#'   synthesis. This is for assigned, subject-defining values you want kept
#'   faithful to a subject's dosing — a treatment arm, a dose group, a
#'   randomization sequence, or a redundant endpoint label such as a character
#'   `NAME` beside a numeric `dvid`. Because the value comes from the same
#'   anchor as the doses, it stays coherent with them. Contrast `covariates`,
#'   which are *blended* into new values across neighbours. A kept value is one
#'   real subject's real value, so use it only where the source data's own
#'   access controls and confidentiality obligations still apply.
#' @param exclude Differential-privacy engines only. Columns removed before
#'   private fitting, such as direct identifiers. [synpmx_avatar()] does not use
#'   this — it drops every undeclared column by default, so not naming a column
#'   is how you drop it.
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
#'
#' # Two columns for one endpoint, and a treatment arm carried through verbatim.
#' roles <- pmx_roles(
#'   id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
#'   dvid = c("YTYPE", "NAME"), covariates = "WT", keep = "ARM"
#' )
pmx_roles <- function(id, time, dv, amt = NULL, evid, cmt = NULL,
                      dvid = NULL, mdv = NULL, rate = NULL,
                      nominal_time = NULL, tad = NULL, occasion = NULL,
                      cens = NULL, limit = NULL, addl = NULL, ii = NULL,
                      covariates = NULL, subject_properties = NULL,
                      dose_covariate = NULL, assigned_dose = NULL,
                      keep = NULL, exclude = NULL) {
  roles <- list(
    id = id, time = time, nominal_time = nominal_time, tad = tad,
    occasion = occasion, dv = dv, amt = amt, evid = evid, cmt = cmt,
    dvid = dvid, mdv = mdv, rate = rate, cens = cens, limit = limit,
    addl = addl, ii = ii, assigned_dose = assigned_dose,
    dose_covariate = dose_covariate,
    covariates = covariates, subject_properties = subject_properties,
    keep = keep, exclude = exclude
  )

  vector_roles <- c("dvid", "covariates", "subject_properties", "keep",
                    "exclude")
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
  # One column may be both `cmt` and `dvid`. NONMEM's `CMT` routinely does both
  # jobs -- the dosing compartment on event rows, the endpoint key on
  # observation rows -- and the alternative is making the user copy the column
  # to itself under another name, which declares nothing extra and helps nobody.
  # Every other collision stays an error.
  shared <- intersect(roles$cmt, roles$dvid)
  if (length(shared)) modeled <- modeled[-match(shared, modeled)]
  # `dose_covariate` must name a column that is ALSO in `covariates`, and that
  # is the whole mechanism rather than an oversight: the avatar's amount is
  # rebuilt from its own *blended* value, so the column has to be blended, which
  # means it has to be declared as a covariate. Requiring the user to copy the
  # column to itself under a second name would declare nothing extra.
  dose_shared <- intersect(roles$dose_covariate, roles$covariates)
  if (length(dose_shared)) modeled <- modeled[-match(dose_shared, modeled)]
  if (!is.null(roles$dose_covariate) &&
      !roles$dose_covariate %in% roles$covariates) {
    stop("`dose_covariate` (\"", roles$dose_covariate, "\") must also be ",
         "named in `covariates`. The avatar's amount is recomputed from its ",
         "own blended value of that column, so the column has to be blended.",
         call. = FALSE)
  }
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
