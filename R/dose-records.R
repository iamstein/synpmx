# Compressed dose records -----------------------------------------------------
#
# NONMEM writes a repeated regimen as one record plus `ADDL` and `II`: "and 364
# more like it, every day". Every reader in this package -- the dose skeleton,
# the derived time after dose, the dosing model, the route key, the population
# fit -- counts dose ROWS, and a study encoded this way shows it one dose where
# the patient received a year of them. `onc_sim` is the study that shows the
# cost: 237 dose rows standing for 71,180 doses, and a population fit that
# returned a clearance six orders of magnitude off because it was fitting a
# cohort that had received a three-hundredth of the drug. The fit had moved, so
# nothing flagged it.
#
# So a declared `addl`/`ii` is expanded on the way in and compressed again on
# the way out. Expansion is exact, because `ADDL` can only ever mean an equal
# interval. Compression is exact too, in the one direction that matters: it
# collapses a run of dose records that are identical in every column but time
# and equally spaced, which is precisely the set of records expansion produces,
# and leaves anything else alone. A synthetic study therefore comes back in the
# encoding its source used, and a source that wrote every dose out stays that
# way -- the compressor only runs where the roles declare the columns.

#' Expand or compress repeated dose records
#'
#' `pmx_expand_doses()` writes out the doses that `ADDL` and `II` imply, one row
#' each, so that every dose the patient received is a row. `pmx_compress_doses()`
#' does the reverse: a run of dose records that agree in every column except
#' time and sit at one interval becomes one record carrying `ADDL` and `II`.
#'
#' Every generator calls both itself -- expansion on the source it reads,
#' compression on the data it returns -- so the synthetic study comes back in
#' the encoding its source used. They are exported for the caller who wants to
#' see what the generator saw, or to hand `nlmixr2` the compact form directly.
#'
#' Neither does anything unless `pmx_roles()` declares `addl` and `ii`. A
#' source that writes every dose out is left exactly as it is.
#'
#' Expansion places dose `k` at `time + k * ii`, and moves a declared
#' `nominal_time` by the same amount. Compression requires the recorded and the
#' nominal gaps to agree, so a schedule whose recorded times drift off the plan
#' is not compressed -- `ADDL` can only place a dose on an exact interval, and
#' writing one would change the study.
#'
#' @param data A data frame of pharmacometric events.
#' @param roles A `pmx_roles` object declaring `addl` and `ii`.
#' @return `data`, with the repeated doses written out or folded back, and the
#'   `addl` and `ii` columns set accordingly. Rows are ordered by subject and
#'   time, with a dose before an observation at the same time.
#' @export
#' @examples
#' roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
#'                    evid = "EVID", addl = "ADDL", ii = "II")
#' compact <- data.frame(
#'   ID = 1, TIME = c(0, 12, 36), DV = c(NA, 3.1, 2.4),
#'   AMT = c(100, 0, 0), EVID = c(1, 0, 0), ADDL = c(2, 0, 0), II = c(24, 0, 0)
#' )
#' expanded <- pmx_expand_doses(compact, roles)
#' nrow(expanded)                        # three doses, two observations
#' identical(pmx_compress_doses(expanded, roles)$ADDL, compact$ADDL)
pmx_expand_doses <- function(data, roles) {
  out <- .expand_doses(data, roles)
  attr(out, "pmx_origin") <- NULL
  out
}

#' @rdname pmx_expand_doses
#' @export
pmx_compress_doses <- function(data, roles) {
  .check_dose_record_roles(roles)
  if (is.null(roles$addl)) return(data)
  data <- as.data.frame(data)
  n <- nrow(data)
  # The columns have to exist on the way out even where the generator never
  # wrote them, and they sit beside the amount, which is where a reader expects
  # them.
  for (role in c("addl", "ii")) {
    column <- roles[[role]]
    if (!column %in% names(data)) {
      data[[column]] <- if (identical(role, "addl")) rep(0L, n) else rep(0, n)
      # `addl` after the amount, `ii` after `addl`.
      anchor <- match(if (identical(role, "addl")) roles$amt %||% roles$evid
                      else roles$addl, names(data))
      if (is.finite(anchor) && anchor < ncol(data) - 1L) {
        data <- data[, append(setdiff(seq_len(ncol(data)), ncol(data)),
                              ncol(data), after = anchor), drop = FALSE]
      }
    }
  }
  if (!n) return(data)
  as_integer <- is.integer(data[[roles$addl]])
  dosed <- .event_rows(data, roles)
  time <- suppressWarnings(as.numeric(data[[roles$time]]))
  nominal <- if (is.null(roles$nominal_time)) NULL else
    suppressWarnings(as.numeric(data[[roles$nominal_time]]))
  ids <- as.character(data[[roles$id]])
  # Everything but the clock has to agree along a run. Amount, compartment,
  # rate, route, a covariate, a stratum: a change in any of them is a different
  # record, and the second record stays.
  compare <- setdiff(names(data), c(roles$time, roles$nominal_time,
                                    roles$addl, roles$ii))
  key <- do.call(paste, c(lapply(compare, function(column) {
    as.character(data[[column]])
  }), sep = "\r"))
  key[is.na(key)] <- ""

  addl_out <- rep(0, n)
  ii_out <- rep(0, n)
  drop <- rep(FALSE, n)
  for (rows in split(seq_len(n), factor(ids, levels = unique(ids)))) {
    d <- rows[dosed[rows] & is.finite(time[rows])]
    if (length(d) < 2L) next
    d <- d[order(time[d])]
    start <- 1L
    while (start < length(d)) {
      end <- start
      gap <- NA_real_
      while (end < length(d) && identical(key[d[end + 1L]], key[d[start]])) {
        step <- time[d[end + 1L]] - time[d[end]]
        if (!is.finite(step) || step <= 0) break
        if (is.na(gap)) gap <- step
        else if (abs(step - gap) > 1e-8 * gap) break
        if (!is.null(nominal)) {
          nominal_step <- nominal[d[end + 1L]] - nominal[d[end]]
          if (!is.finite(nominal_step) ||
              abs(nominal_step - gap) > 1e-8 * gap) break
        }
        end <- end + 1L
      }
      if (end > start) {
        addl_out[d[start]] <- end - start
        ii_out[d[start]] <- gap
        drop[d[(start + 1L):end]] <- TRUE
      }
      start <- end + 1L
    }
  }
  out <- data[!drop, , drop = FALSE]
  out[[roles$addl]] <- if (as_integer) as.integer(addl_out[!drop]) else
    addl_out[!drop]
  out[[roles$ii]] <- ii_out[!drop]
  rownames(out) <- NULL
  out
}

# The expansion, keeping which source row each output row came from as the
# `pmx_origin` attribute. `validate_pmx()` needs that map to compare a derived
# time after dose against the declared column row for row.
.expand_doses <- function(data, roles) {
  .check_dose_record_roles(roles)
  data <- as.data.frame(data)
  n <- nrow(data)
  if (is.null(roles$addl) || !n) {
    attr(data, "pmx_origin") <- seq_len(n)
    return(data)
  }
  addl <- suppressWarnings(as.numeric(data[[roles$addl]]))
  ii <- suppressWarnings(as.numeric(data[[roles$ii]]))
  addl[!is.finite(addl)] <- 0
  ii[!is.finite(ii)] <- 0
  dosed <- .event_rows(data, roles)
  repeated <- dosed & addl > 0
  bad <- which(repeated & (addl != round(addl) | ii <= 0))
  if (length(bad)) {
    stop(.condition_text(
      "`", roles$addl, "` and `", roles$ii, "` cannot be expanded on ",
      length(bad), " dose record(s): `", roles$addl, "` has to be a whole ",
      "number of additional doses and `", roles$ii, "` a positive interval ",
      "wherever `", roles$addl, "` is above zero.",
      items = utils::head(sprintf("row %d: %s = %s, %s = %s", bad,
                                  roles$addl, format(addl[bad]),
                                  roles$ii, format(ii[bad])), 10L)),
      call. = FALSE)
  }
  copies <- rep(1L, n)
  copies[repeated] <- as.integer(addl[repeated]) + 1L
  origin <- rep(seq_len(n), copies)
  k <- sequence(copies) - 1L
  out <- data[origin, , drop = FALSE]
  shift <- k * ii[origin]
  time <- suppressWarnings(as.numeric(out[[roles$time]])) + shift
  out[[roles$time]] <- time
  if (!is.null(roles$nominal_time)) {
    nominal <- suppressWarnings(as.numeric(out[[roles$nominal_time]]))
    moved <- is.finite(nominal)
    nominal[moved] <- nominal[moved] + shift[moved]
    out[[roles$nominal_time]] <- nominal
  }
  out[[roles$addl]] <- if (is.integer(data[[roles$addl]])) 0L else 0
  out[[roles$ii]] <- if (is.integer(data[[roles$ii]])) 0L else 0
  # Subjects in the order they arrived, time within subject, and a dose before
  # an observation at the same time. The last two keys hold source order among
  # ties, so a study already in order stays in it.
  ids <- as.character(out[[roles$id]])
  ordering <- order(match(ids, unique(ids)), time, !dosed[origin], origin, k)
  out <- out[ordering, , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "pmx_origin") <- origin[ordering]
  out
}

.check_dose_record_roles <- function(roles) {
  if (!inherits(roles, "pmx_roles")) {
    stop("`roles` must come from `pmx_roles()`.", call. = FALSE)
  }
  if (is.null(roles$addl) != is.null(roles$ii)) {
    stop("`addl` and `ii` are declared together or not at all: one says how ",
         "many more doses there are and the other how far apart, and neither ",
         "means anything alone.", call. = FALSE)
  }
  invisible(TRUE)
}

# The estimation table for `nlmixr2`, folded back to the compact form it was
# handed as. `nlmixr2` reads `ADDL`/`II` natively and is fastest on them, and
# on a study whose doses were compressed to begin with, folding them back is
# exact by construction. The columns of this table are fixed by
# `.model_estimation_data()`, so the roles are spelled out here rather than
# taken from the caller.
.compress_dose_runs <- function(frame) {
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID",
                     cmt = if ("CMT" %in% names(frame)) "CMT" else NULL,
                     rate = if ("RATE" %in% names(frame)) "RATE" else NULL,
                     addl = "ADDL", ii = "II")
  pmx_compress_doses(frame, roles)
}
