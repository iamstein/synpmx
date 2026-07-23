# Try the AVATAR blending generator on a real modeling dataset.
#
# This is the DEFAULT method. Use it when the synthetic data STAYS INSIDE your
# trusted computing environment (you are the only consumer, governance and
# access controls apply). It resamples and blends whole source subject
# trajectories -- more faithful and far less ceremony than the differentially
# private path -- but makes NO formal privacy guarantee. If the data may cross a
# trust boundary, use try_dp_calibrated.R instead. See README.md and
# design/METHOD_DISCUSSION.md for the decision rule.
#
# Runs in the safe computing environment, in a gitignored folder so nothing here
# can be committed.
#
# Workflow:
#   1. fill in the CONFIG block (just the column roles)
#   2. run with DRY_RUN = TRUE to prove the plumbing on public data
#   3. set DRY_RUN = FALSE, point DATA_PATH at the real dataset, rerun
#   4. compare the synthetic table to the source

library(synpmx)
# If not installed here: devtools::load_all("/path/to/synpmx")

# ============================================================================
# CONFIG  --  everything you edit is in this block
# ============================================================================

DRY_RUN <- TRUE          # TRUE  = practice run on public stand-in data, no
                         #         real data read.
                         # FALSE = use the real dataset at DATA_PATH.
DATA_PATH <- "data/your_modeling_dataset.csv"
OUT_DIR   <- "output"
SEED      <- 1234        # reproducibility seed; the caller's RNG is untouched

# --- Column roles: map YOUR columns onto PMX roles --------------------------
# Only the roles the real dataset has. AVATAR needs no model, priors, design, or
# budget -- just the column meanings. List baseline covariates in `covariates`
# and they are blended from compatible donors like everything else.
ROLES <- pmx_roles(
  id       = "ID",
  time     = "TIME",
  dv       = "DV",
  amt      = "AMT",
  evid     = "EVID",
  dvid     = "DVID",         # NULL if a single endpoint; or c("YTYPE","NAME")
                             #   to declare two consistent endpoint-key columns
  cmt      = "CMT",
  mdv      = "MDV",
  rate     = NULL,           # set if you have infusion RATE rows
  occasion = NULL,           # set if TIME resets by occasion
  covariates = c("WT", "AGE", "SEX"),  # BLENDED into new values across donors
  keep     = NULL            # copied verbatim from the anchor: treatment arm,
                             #   dose group, anything kept faithful to its doses
  # Every column NOT named by a role above is dropped, and the run says which.
)

# Number of synthetic subjects. NULL keeps the source cohort size.
N_SUBJECTS <- NULL

# Optional per-dataset cleanup: fix units, drop screening rows, recode. Receives
# the raw data frame, returns a cleaned one.
PREP <- function(d) {
  d
}

# ============================================================================
# WORKFLOW  --  you should not need to edit below here
# ============================================================================

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

read_source <- function() {
  if (DRY_RUN) {
    message("DRY_RUN = TRUE: using a public simulated dataset, not your data.")
    return(pmx_simulated_fixture(40))
  }
  if (!file.exists(DATA_PATH)) {
    stop("DATA_PATH not found: ", normalizePath(DATA_PATH, mustWork = FALSE))
  }
  utils::read.csv(DATA_PATH, stringsAsFactors = FALSE)
}

raw   <- PREP(read_source())
roles <- if (DRY_RUN) {
  pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
            dvid = "DVID", cmt = "CMT", mdv = "MDV",
            covariates = c("WT", "AGE", "SEX"))
} else {
  ROLES
}

message("\n== Source structural validation ==")
report <- validate_pmx(raw, roles, strict = FALSE)
if (!report$valid) {
  print(report)   # names each problem, the role, and the column it maps to
  stop("Source failed structural validation. See the numbered problems above; ",
       "each says which role and column to fix in `pmx_roles()`.", call. = FALSE)
}
message("  valid; ", length(unique(raw[[roles$id]])), " subjects")

# --- Synthesize (AVATAR blending) -------------------------------------------
message("\n== AVATAR synthesis ==")
synthetic <- synpmx_avatar(raw, roles, n_subjects = N_SUBJECTS, seed = SEED)
stopifnot(validate_pmx(synthetic, roles)$valid)
message("  generated ", nrow(synthetic), " rows for ",
        length(unique(synthetic[[roles$id]])), " subjects; new identifiers: ",
        length(intersect(synthetic[[roles$id]], raw[[roles$id]])) == 0)

utils::write.csv(synthetic, file.path(OUT_DIR, "avatar_synthetic.csv"),
                 row.names = FALSE)

# --- Restricted diagnostic: source vs synthetic -----------------------------
# Source-derived. Keep it in the safe environment; do not export the figure.
message("\n== Restricted source-vs-synthetic diagnostic ==")
if (requireNamespace("ggplot2", quietly = TRUE)) {
  obs <- function(d, label) {
    keep <- d[[roles$evid]] == 0 & !is.na(d[[roles$dv]])
    data.frame(dataset = label,
               time = as.numeric(d[[roles$time]][keep]),
               dv   = as.numeric(d[[roles$dv]][keep]),
               endpoint = if (is.null(roles$dvid)) "DV" else
                 as.character(d[[roles$dvid]][keep]))
  }
  compare <- rbind(obs(raw, "source"), obs(synthetic, "synthetic"))
  p <- ggplot2::ggplot(compare, ggplot2::aes(time, dv, colour = dataset)) +
    ggplot2::geom_point(alpha = 0.4, size = 0.7) +
    ggplot2::facet_grid(dataset ~ endpoint, scales = "free_y") +
    ggplot2::labs(title = "RESTRICTED: source vs AVATAR synthetic",
                  subtitle = "do not export", x = "Time", y = "DV") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
  ggplot2::ggsave(file.path(OUT_DIR, "RESTRICTED_source_vs_synthetic.png"), p,
                  width = 7, height = 6, dpi = 110)
  message("  wrote ", file.path(OUT_DIR, "RESTRICTED_source_vs_synthetic.png"),
          " (restricted; keep in the safe environment)")
} else {
  message("  ggplot2 not available; skipping the diagnostic plot")
}

message("\nDone. The synthetic table carries no formal privacy guarantee; it is ",
        "built by blending real subject trajectories and is for ",
        "trusted-environment use only.")
