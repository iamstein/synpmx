#' Public PMX censoring fixture
#'
#' A tiny, fully simulated event table illustrating uncensored, left-censored,
#' right-censored, and interval-censored Monolix-style records. It contains no
#' patient or proprietary information.
#'
#' @return A data frame with `CENS = 0`, `1`, and `-1`; interval censoring uses
#'   `DV` as the reported upper boundary and `LIMIT` as the lower boundary.
#' @export
#'
#' @examples
#' fixture <- pmx_censoring_fixture()
#' roles <- pmx_roles(
#'   id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
#'   dvid = "DVID", mdv = "MDV", cens = "CENS", limit = "LIMIT"
#' )
#' validate_pmx(fixture, roles)$valid
pmx_censoring_fixture <- function() {
  data.frame(
    ID = rep(1L, 5L),
    TIME = c(0, 1, 2, 3, 4),
    DV = c(0, 5, 1, 10, 4),
    AMT = c(100, 0, 0, 0, 0),
    EVID = c(1L, 0L, 0L, 0L, 0L),
    DVID = factor(rep("marker", 5L), levels = "marker"),
    MDV = c(1L, 0L, 0L, 0L, 0L),
    CENS = c(0L, 0L, 1L, -1L, 1L),
    LIMIT = c(NA, NA, NA, NA, 2),
    stringsAsFactors = FALSE
  )
}

#' Fully simulated public repeated-dose fixture
#'
#' Creates a deterministic two-endpoint PMX event table for privacy-utility and
#' workflow tests. The data are generated from fixed formulas and contain no
#' patient or proprietary information, so no random seed is required.
#'
#' `TIME`, `NTIME`, and `TAD` are in hours; `AMT` is in arbitrary dose units;
#' `WT` is in kg; `AGE` is in years; `cp` uses arbitrary concentration units;
#' and `pd` uses arbitrary response units.
#'
#' @param n_subjects Positive number of fully simulated subjects. The default
#'   is 60 to provide a larger companion to six- and twelve-subject examples.
#'
#' @return A deterministic data frame with two doses per subject, dose-relative
#'   `cp`, global study-time `pd`, explicit nominal/TAD/occasion roles, and
#'   fixed public factor levels.
#' @export
#'
#' @examples
#' public_data <- pmx_simulated_fixture(12)
#' length(unique(public_data$ID))
#' table(public_data$DVID, public_data$EVID == 0)
pmx_simulated_fixture <- function(n_subjects = 60L) {
  n_subjects <- .positive_integer(n_subjects, "n_subjects")
  pieces <- lapply(seq_len(n_subjects), function(subject) {
    time <- c(0, 0, 0.25, 1, 2, 4, 6, 8, 12, 12, 12.25, 13, 14, 16, 18, 20)
    endpoint <- c(
      "cp", "pd", "cp", "cp", "cp", "pd", "cp", "pd",
      "cp", "pd", "cp", "cp", "cp", "pd", "cp", "pd"
    )
    evid <- c(1L, 0L, rep(0L, 6L), 1L, rep(0L, 7L))
    cp_base <- c(
      NA, NA, 1, 8, 4, NA, 1.5, NA,
      NA, NA, 1.2, 8.5, 4.2, NA, 1.6, NA
    )
    pd_base <- c(
      NA, 80, NA, NA, NA, 70, NA, 55,
      NA, 40, NA, NA, NA, 50, NA, 65
    )
    dv <- ifelse(endpoint == "cp", cp_base, pd_base)
    dv[evid != 0] <- 0
    subject_effect <- 1 + 0.12 * sin(
      2 * pi * subject / max(n_subjects, 6L)
    )
    dv[evid == 0] <- dv[evid == 0] * subject_effect
    dose_times <- c(0, 12)
    occasion <- pmax(1L, findInterval(time, dose_times))
    tad <- pmax(0, time - dose_times[pmin(occasion, 2L)])
    data.frame(
      ID = as.integer(subject), TIME = time, NTIME = time, TAD = tad,
      OCC = occasion, DV = dv,
      AMT = ifelse(evid != 0, 100 + 5 * sin(subject), 0), RATE = 0,
      EVID = evid,
      CMT = ifelse(evid != 0, 1L, ifelse(endpoint == "cp", 2L, 3L)),
      DVID = factor(endpoint, levels = c("cp", "pd")),
      MDV = ifelse(evid == 0, 0L, 1L), CENS = 0L, LIMIT = NA_real_,
      WT = rep(70 + 12 * sin(subject / 4), length(time)),
      AGE = rep(as.integer(25 + subject %% 45L), length(time)),
      SEX = factor(
        rep(ifelse(subject %% 2L, "female", "male"), length(time)),
        levels = c("female", "male")
      ),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}
