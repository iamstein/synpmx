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
#'   roles**: `pmx_roles(..., cmt = "CMT", dvid = "CMT")`. `cmt` and `adm` may
#'   likewise be one column, for the same reason — a NONMEM dataset has no
#'   administration column, and the compartment a dose enters is what says how
#'   it was given and which drug it is. Those two are the only permitted
#'   overlaps; every other collision is an error.
#' @param nominal_time,occasion Optional time metadata columns.
#' @param tad Time-after-dose column. **It is an output, not an input.**
#'   [synpmx_avatar()] recomputes it from the generated times and the generated
#'   dose rows and overwrites whatever the source held, so declaring the role
#'   says which column to overwrite and to carry through — the source's values
#'   are never used to generate anything.
#'
#'   [validate_pmx()] does read them, and reports (as a non-fatal warning) where
#'   the declared column disagrees with time since the most recent dose row.
#'   That disagreement is worth knowing: a study may measure TAD from the end of
#'   an infusion, from a nominal dose time, or from an assigned dosing occasion
#'   rather than the most recent dose. `nlmixr2data::nimoData` disagrees on 45%
#'   of its observation rows. Where it does, the synthetic column follows the
#'   derivation rather than the source's convention.
#'
#'   Two limits. Samples taken before any dose are reported as 0, because
#'   [validate_pmx()] refuses a negative TAD — not because a baseline sample is
#'   genuinely zero hours after a dose it precedes. And where `addl` or `ii` is
#'   declared, the derivation cannot see the doses those imply, so the agreement
#'   check is skipped and says so.
#' @param cens,limit Optional Monolix-style censoring indicator and other
#'   interval-boundary columns.
#' @param addl,ii Optional additional-dose and interdose-interval columns, in
#'   the NONMEM sense: one record plus "and `ADDL` more like it, every `II`".
#'   Declared together or not at all. Every generator expands them on the way
#'   in, so each dose the patient received is a row to the dose skeleton, the
#'   derived time after dose and the population fit, and folds the output back
#'   on the way out, so the synthetic study keeps its source's encoding. See
#'   [pmx_expand_doses()].
#' @param adm,routes Optional administration-id column and what its values mean,
#'   as `adm = "ADM", routes = c("1" = "iv", "2" = "extravascular")`. Declared
#'   together or not at all: which id is which route is a convention of the
#'   dataset and cannot be read off the numbers, so reading it wrong would put
#'   doses in the wrong compartment without failing. Two routes are recognised.
#'   `iv` covers a bolus and an infusion, told apart per record by `rate` rather
#'   than declared twice; `extravascular` covers subcutaneous, oral and
#'   intramuscular, which are one first-order absorption with an absorption rate
#'   constant and a bioavailability. A study dosing both ways is fitted with one
#'   model that routes each dose record to its own compartment, so a patient may
#'   receive both.
#'
#'   `adm` may name the same column as `cmt`: a NONMEM dataset has no
#'   administration column, and the compartment a dose enters is the
#'   administration id. See `dose_endpoints` for the other thing an
#'   administration id says — which drug the dose is.
#' @param dose_endpoints Which endpoint each administration id doses, as
#'   `dose_endpoints = c("1" = "drug A PK", "2" = "drug B PK")`. Needs `adm`,
#'   and names endpoints by their `dvid` value.
#'
#'   Declare it for a study that gives **two different drugs** and measures a
#'   concentration of each. Undeclared, every dose record drives every
#'   concentration endpoint's fit and its generated profile, which is right for
#'   a parent and its metabolite -- one dose, two analytes -- and wrong for a
#'   combination, where each drug's model would be fitted against the other
#'   drug's doses as well as its own.
#'
#'   The two cases cannot be told apart from the data: a metabolite has no dose
#'   records of its own, so "this endpoint's doses are the ones marked `2`" and
#'   "this endpoint has no doses" are the same table. Hence a declaration, on
#'   the same argument `routes` makes.
#'
#'   Every administration id present on a dose record must be named, and every
#'   endpoint named must be one the fit puts a structural model on — both are
#'   checked against the study when the fit runs. A mapping that leaves a
#'   concentration with no dose record at all is refused, and that refusal is
#'   the signal that the endpoint is a metabolite and this is not the argument
#'   for it.
#'
#'   Declared, each drug gets its own dose schedule per arm — its own amounts,
#'   its own interval, its own ladder of reductions — each endpoint's fit and
#'   generated profile read only that drug's doses, each observation's time
#'   after dose is measured from that drug's last dose, and the generated study
#'   writes both drugs' dose records back with their own `adm` value and
#'   compartment.
#'
#'   In a NONMEM dataset the compartment does this job, so declare the column
#'   twice: `pmx_roles(cmt = "CMT", adm = "CMT", routes = ..., dose_endpoints =
#'   c("1" = "drug A PK", "3" = "drug B PK"))`.
#'
#'   Either way, [model_report()] says which dose records drove which fit
#'   whenever a study has more than one concentration.
#' @param covariates Baseline covariate column names, or `NULL`.
#' @param strata Treatment arm, dose group, cohort — any **assigned,
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
#' @param endpoint_types What kind of values an endpoint takes, as a named
#'   character vector keyed by the endpoint's `dvid` value — for example
#'   `c("PD - Binary" = "binary", "PD - Ordinal" = "ordinal")`. Use `"DV"` as
#'   the name when no `dvid` is declared. Each value is one of `"continuous"`,
#'   `"binary"`, `"ordinal"`, or `"integer"`; `"count"` is accepted for
#'   `"integer"`.
#'
#'   This exists because blending is a weighted mean, and a weighted mean of
#'   several patients' zeros and ones is a number between them. Without it a
#'   binary endpoint comes back continuous — measured on `xgxr::mad`, a 0/1
#'   endpoint came back as 600 distinct values from -0.13 to 1.08. Declared or
#'   inferred, [synpmx_avatar()] snaps a `"binary"` or `"ordinal"` endpoint onto
#'   the levels the source used and rounds an `"integer"` one to whole numbers.
#'
#'   Left `NULL`, the type is inferred per endpoint from whether every observed
#'   value is a whole number and how many distinct ones there are;
#'   [pmx_endpoint_types()] reports what that inference decided and why. Unlike
#'   `dose_covariate`, this inference is on by default, because the question it
#'   asks is answered outright by the data and the repair reproduces granularity
#'   the source already had. Declare the type where the data is misleading: an
#'   endpoint recorded as whole numbers that is genuinely continuous, or a scale
#'   with a level nobody in this cohort happened to hit.
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
#'
#' # A study dosed both ways, saying which administration id is which route.
#' roles <- pmx_roles(
#'   id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
#'   cmt = "CMT", adm = "ADM", routes = c("1" = "iv", "2" = "extravascular")
#' )
pmx_roles <- function(id, time, dv, amt = NULL, evid, cmt = NULL,
                      dvid = NULL, mdv = NULL, rate = NULL,
                      nominal_time = NULL, tad = NULL, occasion = NULL,
                      cens = NULL, limit = NULL, addl = NULL, ii = NULL,
                      adm = NULL, routes = NULL, dose_endpoints = NULL,
                      covariates = NULL, strata = NULL,
                      dose_covariate = NULL, assigned_dose = NULL,
                      endpoint_types = NULL,
                      keep = NULL, exclude = NULL) {
  roles <- list(
    id = id, time = time, nominal_time = nominal_time, tad = tad,
    occasion = occasion, dv = dv, amt = amt, evid = evid, cmt = cmt,
    dvid = dvid, mdv = mdv, rate = rate, cens = cens, limit = limit,
    addl = addl, ii = ii, adm = adm, routes = routes,
    dose_endpoints = dose_endpoints,
    assigned_dose = assigned_dose,
    dose_covariate = dose_covariate,
    covariates = covariates, strata = strata,
    keep = keep, exclude = exclude
  )
  roles$routes <- .validate_routes(adm, routes)
  roles$dose_endpoints <- .validate_dose_endpoints(adm, dose_endpoints)

  vector_roles <- c("dvid", "covariates", "strata", "keep",
                    "exclude")
  # `routes` and `dose_endpoints` are mappings rather than column names, so they
  # sit out of both loops below and out of the collision check.
  routes <- roles$routes
  dose_endpoints <- roles$dose_endpoints
  roles$routes <- NULL
  roles$dose_endpoints <- NULL
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
  # And one column may be both `cmt` and `adm`, for the same reason. A NONMEM
  # dataset has no administration column: the compartment a dose enters is the
  # administration id, and in a study giving two drugs it is what separates
  # them. Requiring a copy of the column under a second name would declare
  # nothing extra.
  shared_adm <- intersect(roles$cmt, roles$adm)
  if (length(shared_adm)) modeled <- modeled[-match(shared_adm, modeled)]
  # `dose_covariate` must name a column that is ALSO in `covariates`, and that
  # is the whole mechanism rather than an oversight: the avatar's amount is
  # rebuilt from its own *blended* value, so the column has to be blended, which
  # means it has to be declared as a covariate. Requiring the user to copy the
  # column to itself under a second name would declare nothing extra.
  dose_shared <- intersect(roles$dose_covariate, roles$covariates)
  if (length(dose_shared)) modeled <- modeled[-match(dose_shared, modeled)]
  if (!is.null(roles$dose_covariate) &&
      !roles$dose_covariate %in% roles$covariates) {
    stop(.condition_text(
      "`dose_covariate` (\"", roles$dose_covariate, "\") must also be named ",
      "in `covariates`.",
      why = paste("The avatar's amount is recomputed from its own blended",
                  "value of that column, so the column has to be blended.")),
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
  # Attached after the column checks above, because this one names endpoints
  # rather than columns: it must not join `modeled` and be tested for collisions
  # with real column names.
  roles$endpoint_types <- .validate_endpoint_types(endpoint_types)
  # Same reason as `endpoint_types`: a mapping from the values of a column, not
  # a column name.
  roles$routes <- routes
  roles$dose_endpoints <- dose_endpoints
  structure(roles, class = "pmx_roles")
}

# What each administration id means. The numbers in an `ADM`-style column are a
# convention of the dataset and nothing can be read off them -- Monolix's own
# library model ships with a header saying IV is 1 and code routing 2 to IV --
# so the mapping is declared or the column is not read at all. `iv` covers a
# bolus and an infusion, which are told apart per record by `rate` exactly as
# Monolix tells them apart; `extravascular` covers subcutaneous, oral and
# intramuscular, which are one first-order absorption with `ka` and `f`.
.pk_routes <- c("iv", "extravascular")

.validate_routes <- function(adm, routes) {
  if (is.null(adm) && is.null(routes)) return(NULL)
  if (is.null(adm)) {
    stop("`routes` needs `adm`: it says what the values of an administration ",
         "column mean, so there has to be a column for it to describe.",
         call. = FALSE)
  }
  if (is.null(routes)) {
    stop(.condition_text(
      "`adm` needs `routes`, as `routes = c(\"1\" = \"iv\", \"2\" = ",
      "\"extravascular\")`.",
      why = paste("Which administration id is which route is a convention of",
                  "the dataset, and reading it wrong puts the doses in the",
                  "wrong compartment without failing.")), call. = FALSE)
  }
  if (!is.character(routes) || !length(routes) || is.null(names(routes)) ||
      anyNA(routes) || any(!nzchar(names(routes)))) {
    stop("`routes` must be a named character vector, mapping each value of `",
         adm, "` to a route.", call. = FALSE)
  }
  unknown <- setdiff(routes, .pk_routes)
  if (length(unknown)) {
    stop(.condition_text(
      "`routes` names route(s) outside the set: ",
      paste(unique(unknown), collapse = ", "), ". Available: ",
      paste(.pk_routes, collapse = ", "), ".",
      why = paste("A bolus and an infusion are both `iv` and are told apart by",
                  "`rate`; subcutaneous, oral and intramuscular are all",
                  "`extravascular`.")), call. = FALSE)
  }
  if (anyDuplicated(names(routes))) {
    stop("`routes` names the same administration id twice: ",
         paste(unique(names(routes)[duplicated(names(routes))]),
               collapse = ", "), ".", call. = FALSE)
  }
  routes
}

# Which drug each administration id is. The same argument `routes` makes: an
# `ADM` value is a convention of the dataset, so the mapping is declared or the
# column is not read for this at all.
#
# Undeclared, every dose record drives every concentration endpoint, which is
# what a parent and its metabolite need and what two co-administered drugs must
# not have. Declaring it is the only way to tell those apart, because a
# metabolite has no dose records of its own and "this endpoint's doses are
# elsewhere" and "this endpoint has no doses" are the same table.
.validate_dose_endpoints <- function(adm, dose_endpoints) {
  if (is.null(dose_endpoints)) return(NULL)
  if (is.null(adm)) {
    stop(.condition_text(
      "`dose_endpoints` needs `adm`: it says which endpoint each value of an ",
      "administration column doses, so there has to be a column for it to ",
      "describe.",
      fix = paste("Name the column that separates the drugs as `adm` -- in a",
                  "NONMEM dataset that is usually `CMT`.")), call. = FALSE)
  }
  if (!is.character(dose_endpoints) || !length(dose_endpoints) ||
      is.null(names(dose_endpoints)) || anyNA(dose_endpoints) ||
      any(!nzchar(dose_endpoints)) || any(!nzchar(names(dose_endpoints)))) {
    stop("`dose_endpoints` must be a named character vector, mapping each ",
         "value of `", adm, "` to the endpoint it doses.", call. = FALSE)
  }
  if (anyDuplicated(names(dose_endpoints))) {
    stop("`dose_endpoints` names the same administration id twice: ",
         paste(unique(names(dose_endpoints)[duplicated(names(dose_endpoints))]),
               collapse = ", "), ".", call. = FALSE)
  }
  dose_endpoints
}

.validate_endpoint_types <- function(endpoint_types) {
  if (is.null(endpoint_types)) return(NULL)
  if (!is.character(endpoint_types) || anyNA(endpoint_types)) {
    stop("`endpoint_types` must be a named character vector, one of ",
         paste(.endpoint_type_choices, collapse = ", "),
         " per endpoint.", call. = FALSE)
  }
  keys <- names(endpoint_types)
  if (is.null(keys) || any(is.na(keys)) || any(!nzchar(keys))) {
    stop("`endpoint_types` must be named by endpoint: ",
         "c(\"PD - Binary\" = \"binary\"). Use \"DV\" when no `dvid` is ",
         "declared.", call. = FALSE)
  }
  if (anyDuplicated(keys)) {
    stop("`endpoint_types` names an endpoint more than once: ",
         paste(unique(keys[duplicated(keys)]), collapse = ", "), ".",
         call. = FALSE)
  }
  accepted <- c(.endpoint_type_choices, names(.endpoint_type_synonyms))
  unknown <- setdiff(endpoint_types, accepted)
  if (length(unknown)) {
    stop("`endpoint_types` must be one of ",
         paste(.endpoint_type_choices, collapse = ", "), "; got ",
         paste(unique(unknown), collapse = ", "), ".", call. = FALSE)
  }
  endpoint_types
}

#' @export
print.pmx_roles <- function(x, ...) {
  cat("Pharmacometric column roles:\n")
  for (role in names(x)) {
    value <- x[[role]]
    if (identical(role, "endpoint_types") && length(value)) {
      value <- paste0(names(value), " = ", value)
    }
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
  columns <- unlist(roles[setdiff(names(roles), .non_column_roles)],
                    use.names = FALSE)
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns)) {
    stop("Role columns not found in `data`: ",
         paste(missing_columns, collapse = ", "), ".", call. = FALSE)
  }
  invisible(TRUE)
}
