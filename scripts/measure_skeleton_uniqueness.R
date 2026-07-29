# How exposed is a source cohort's event skeleton, and how much of that does
# `coarsen_time = TRUE` actually remove?
#
# synpmx_avatar() copies each avatar's event skeleton verbatim from one anchor.
# A source subject whose observation time vector is shared with nobody therefore
# hands an identifying schedule to every avatar anchored on it -- the defect
# design/TEST_SIM.md records as SIM-014 and, until `coarsen_time`, gated only
# against the structural/DP engine.
#
# This script measures the before and after on every public dataset, so the
# registry in design/TEST_SIM.md can be kept honest and so the same table can be
# produced for a real study under scripts_private/. It reports three classes,
# because they behave differently and only one of them is coarsening's job:
#
#   obs_time_alone  the observation time vector is shared with nobody. What
#                   coarsening collapses. Expect ~100% under actual recorded
#                   times and ~0% under nominal ones.
#   n_obs_alone     the observation *count* is shared with nobody. Coarsening
#                   cannot change a count, so this is the residual it leaves --
#                   dropout, discontinuation, missed visits -- and it is the
#                   outlier screen's job, not the grid's.
#   signature_alone the pmx_roles() event signature -- dose structure, dose
#                   amount, endpoint set -- is shared with nobody. Weight-based
#                   dosing makes this unique regardless of schedule, and
#                   coarsening does not change it either.
#
# Reading it: a dataset where obs_time_alone falls to zero after coarsening and
# n_obs_alone is already low is fully served by the default. One where
# n_obs_alone stays high needs the screen as well. One where obs_time_alone
# stays high after coarsening has near-coincident samples that defeated the grid
# threshold, and should be inspected rather than trusted.
#
# Run with:  Rscript scripts/measure_skeleton_uniqueness.R

if (!requireNamespace("synpmx", quietly = TRUE)) {
  stop("Install synpmx before running this script: R CMD INSTALL .")
}
if (!"skeleton_uniqueness" %in% getNamespaceExports("synpmx")) {
  stop(
    "The installed synpmx predates skeleton_uniqueness(). Reinstall from ",
    "this repository with `R CMD INSTALL .` before running it."
  )
}
if (!requireNamespace("nlmixr2data", quietly = TRUE)) {
  stop("Install nlmixr2data to run the public-dataset measurements.")
}

library(synpmx)

# One row per dataset. `coarsen()` reaches for the same internal the generator
# uses, rather than reimplementing the snap, so this measures what actually
# happens rather than an approximation of it (AGENTS.md: share one metric
# implementation).
coarsen <- utils::getFromNamespace(".coarsen_source_time", "synpmx")

summarize <- function(data, roles, label, stage) {
  report <- skeleton_uniqueness(data, roles)
  n <- nrow(report)
  data.frame(
    dataset = label,
    stage = stage,
    subjects = n,
    obs_time_alone = attr(report, "n_alone"),
    n_obs_alone = attr(report, "n_alone_n_obs"),
    signature_alone = attr(report, "n_alone_signature"),
    pct_obs_time_alone = round(100 * attr(report, "n_alone") / n, 1),
    smallest_class = attr(report, "min_class"),
    stringsAsFactors = FALSE
  )
}

measure <- function(data, roles, label) {
  before <- summarize(data, roles, label, "source")
  coarsened <- coarsen(data, roles)
  after <- summarize(coarsened$source, roles, label, "coarsened")
  after$grid <- coarsened$grid
  before$grid <- NA_character_
  after$deviation_sd <- if (length(coarsened$deviations)) {
    round(stats::sd(coarsened$deviations), 4)
  } else 0
  before$deviation_sd <- NA_real_
  rbind(before, after)
}

load_dataset <- function(name) {
  environment <- new.env(parent = emptyenv())
  utils::data(list = name, package = "nlmixr2data", envir = environment)
  get(name, envir = environment)
}

# The same public datasets the simulation evaluator uses, so the two reports can
# be read side by side. Roles mirror design/TEST_SIM.md's dataset registry.
cases <- list(
  theo_md = function() list(
    data = load_dataset("theo_md"),
    roles = pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                      evid = "EVID", cmt = "CMT", covariates = "WT")
  ),
  warfarin = function() list(
    data = load_dataset("warfarin"),
    roles = pmx_roles(id = "id", time = "time", dv = "dv", amt = "amt",
                      evid = "evid", dvid = "dvid", covariates = "wt")
  ),
  nimoData = function() list(
    data = load_dataset("nimoData"),
    roles = pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                      evid = "EVID", cmt = "CMT", rate = "RATE",
                      covariates = "DOS")
  ),
  mavoglurant = function() list(
    data = load_dataset("mavoglurant"),
    roles = pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                      evid = "EVID", cmt = "CMT", rate = "RATE",
                      occasion = "OCC", covariates = c("WT", "AGE"))
  ),
  simulated_fixture = function() list(
    data = pmx_simulated_fixture(30),
    roles = pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                      evid = "EVID", cmt = "CMT", dvid = "DVID",
                      covariates = "WT")
  )
)

rows <- list()
for (label in names(cases)) {
  case <- tryCatch(cases[[label]](), error = function(e) NULL)
  if (is.null(case)) {
    message("skipped ", label, ": could not be loaded with the declared roles.")
    next
  }
  measured <- tryCatch(measure(case$data, case$roles, label),
                       error = function(e) {
                         message("skipped ", label, ": ", conditionMessage(e))
                         NULL
                       })
  if (!is.null(measured)) rows[[label]] <- measured
}

report <- do.call(rbind, rows)
rownames(report) <- NULL

cat("\nSkeleton uniqueness, before and after coarsen_time\n")
cat("==================================================\n\n")
print(report, row.names = FALSE)

cat("\nColumns count subjects that are ALONE in the named class.\n")
cat("obs_time_alone is what coarsening collapses; n_obs_alone is the residual\n")
cat("left for flag_identifiable_subjects(); signature_alone is dose structure\n")
cat("and amount, which coarsening does not change.\n")
cat("\nSource-derived. Not releasable unless the source is separately public.\n")

invisible(report)
