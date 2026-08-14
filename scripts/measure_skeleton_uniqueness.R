# How exposed is a source cohort's event skeleton, and how much of that does
# `coarsen_time = TRUE` actually remove?
#
# synpmx_avatar() copies each avatar's event skeleton verbatim from one anchor.
# A source subject whose observation time vector is shared with nobody therefore
# hands an identifying schedule to every avatar anchored on it. Until
# `coarsen_time`, that was gated only against the structural/DP engine.
#
# This script measures the before and after on every public dataset, so the
# recorded numbers can be kept honest and so the same table can be
# produced for a real study under scripts_private/. The columns match the ones
# the demo vignette explains in prose, under "How well did the obfuscation
# work?".
#
#   alone            subjects whose visit schedule NO other subject shares. An
#                    avatar built on one of them wears a schedule belonging to
#                    exactly one real person. The `coarsened` row is the number
#                    you would consider dropping. It splits by cause, and the
#                    two causes have different answers:
#
#     unique_moment    the subject was observed at a time nobody else was. The
#                      grid is meant to absorb this, so a nonzero count means
#                      the INFERRED grid found no shared schedule --
#                      declaring a `nominal_time` role is the fix.
#     unique_pattern   every one of the subject's times is shared, but the
#                      combination of visits they attended is theirs alone.
#                      Dropout and missed visits. No grid can fix this at any
#                      resolution; only dropping or remediating helps.
#
#   unique_n_samples subjects with an observation count nobody else has.
#   unique_dosing    subjects unique on the event signature -- dose structure and
#                    amount. Weight-based dosing or per-subject titration keeps
#                    this high no matter what the grid does.
#   smallest_group   the fewest subjects sharing any one schedule; the effective
#                    k on the schedule axis.
#
# Reading it: `unique_moment` falling to zero means the grid did its whole job
# and no tuning will improve it. A high `unique_pattern` is not a coarsening
# failure at all, and the remedy is `flag_identifiable_subjects()`.
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
# implementation). The grid kind and the deviations it removed are only
# available from the internal, so the `after` row still reaches for it; the
# counts themselves go through the public `coarsen_time` argument.
coarsen <- utils::getFromNamespace(".coarsen_source_time", "synpmx")

# `n_alone_n_obs` and `n_alone_signature` were never attribute names, so those
# two columns came back NULL and silently vanished from the table.
summarize <- function(data, roles, label, stage, coarsen_time = FALSE) {
  report <- skeleton_uniqueness(data, roles, coarsen_time = coarsen_time)
  data.frame(
    dataset = label,
    stage = stage,
    subjects = nrow(report),
    alone = attr(report, "n_unique_schedule"),
    unique_moment = attr(report, "n_unshared_time"),
    unique_pattern = attr(report, "n_unique_schedule") -
      attr(report, "n_unshared_time"),
    unique_n_samples = attr(report, "n_unique_obs_count"),
    unique_dosing = attr(report, "n_unique_dose_signature"),
    smallest_group = attr(report, "min_class"),
    stringsAsFactors = FALSE
  )
}

measure <- function(data, roles, label) {
  before <- summarize(data, roles, label, "source")
  after <- summarize(data, roles, label, "coarsened", coarsen_time = TRUE)
  coarsened <- coarsen(data, roles)
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
# be read side by side. Roles mirror tests/testthat/helper-simulation-evaluation.R.
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

cat("\n`alone` = subjects whose visit schedule no other subject shares; on the\n")
cat("coarsened row it is the number you would consider dropping. It splits into\n")
cat("unique_moment (seen at a time nobody else was -- the grid's job, fixed by\n")
cat("declaring nominal_time) and unique_pattern (which visits were attended --\n")
cat("dropout, which no grid can touch). unique_n_samples and unique_dosing are\n")
cat("likewise outside coarsening's reach.\n")
cat("\nSource-derived. Not releasable unless the source is separately public.\n")

invisible(report)
