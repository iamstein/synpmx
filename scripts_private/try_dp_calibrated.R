# Try the DIFFERENTIALLY PRIVATE structural generator on a real modeling dataset.
#
# Use this when the synthetic data may REACH SOMEONE THE SOURCE DATA COULD NOT
# (shared externally, published, moved outside its access controls) and needs a
# formal (epsilon, delta) guarantee. If it stays under the source data's own
# controls and obligations -- including on your own workstation -- use
# try_avatar.qmd instead -- it is simpler and more faithful. See the
# "synpmx-privacy" vignette for the decision rule (which method when).
#
# Runs in the safe computing environment, in a gitignored folder so nothing here
# can be committed.
#
# It runs TWO versions and shows them side by side:
#   PRIOR      generated from public inputs alone. No data, no privacy budget.
#   CALIBRATED the same, with PK magnitude corrected by a private release.
#
# Workflow:
#   1. fill in the CONFIG block from the protocol and preclinical prediction
#   2. run with DRY_RUN = TRUE to prove the plumbing on public data
#   3. set DRY_RUN = FALSE, point DATA_PATH at the real dataset, rerun
#   4. read the pre-flight verdict, then compare the two versions

library(synpmx)
# If not installed here: devtools::load_all("/path/to/synpmx")

# ============================================================================
# CONFIG  --  everything you edit is in this block
# ============================================================================

DRY_RUN <- TRUE          # TRUE  = practice run on public stand-in data, no
                         #         real data read, no privacy budget spent.
                         # FALSE = use the real dataset at DATA_PATH.
DATA_PATH   <- "data/your_modeling_dataset.csv"
OUT_DIR     <- "output"
EPSILON     <- 1             # from governance, not chosen to look good
SEED        <- 1234          # ordinary generation seed; unrelated to DP noise

# --- 1. Trial design (from the protocol) ------------------------------------
# Dose escalation BOTH between cohorts AND within each patient: every patient
# gets three doses a week apart, and each cohort starts a step higher. Supply
# one increasing sequence per cohort; cohort_sizes gives the patients in each.
DESIGN <- pmx_trial_design(
  dose_escalation = list(c(1, 2, 4),       # cohort 1: 1 mg, then 2, then 4
                         c(2, 4, 8),        # cohort 2
                         c(4, 8, 16)),      # cohort 3
  cohort_sizes    = c(6, 6, 6),            # patients per cohort
  dose_times      = c(0, 7, 14),           # days (use consistent time units)
  sampling        = c(0, 1, 4, 24, 72, 167), # times after EACH dose
  duration        = 0,                     # 0 = oral/bolus; >0 = infusion hours
  source          = "FILL IN: protocol vX section Y"
)

# --- 2. PK structural model (from preclinical / FIH prediction) -------------
# Two-compartment model with oral absorption. `typical` are your PREDICTED
# parameter values, from allometric scaling -- NOT measured from this dataset.
#   cl = clearance, v = central volume, q = inter-compartmental clearance,
#   v2 = peripheral volume, ka = absorption rate constant.
MODEL <- pmx_structural_model(
  pk          = "2cmt_oral",
  typical     = c(cl = 10, v = 30, q = 15, v2 = 100, ka = 1),
  pd          = "none",                    # PK only for now; see PD block below
  source      = "FILL IN: allometric scaling memo vX",
  iiv         = c(cl = 0.3, v = 0.2),      # public between-subject CV, per param
  residual_cv = 0.15                       # public assay + model CV
)

# --- Prior: how wrong might the PREDICTION be? ------------------------------
# This is NOT a prior on cl, v, ka separately. The calibration releases ONE
# number: a multiplier on the whole exposure (clearance). The prior is the
# plausible range of that multiplier.
#
#   c(1/4, 4)  means "the true clearance could be anywhere from 1/4x to 4x my
#              prediction" -- i.e. the prediction might be off by up to 4-fold
#              in either direction. That is typical for allometric scaling.
#
# A tighter range (e.g. c(1/2, 2)) costs less privacy budget but is only honest
# if you genuinely believe the prediction that well. See the explanation printed
# by the pre-flight below.
PRIORS <- pmx_priors(
  pk = pmx_prior(c(1 / 4, 4),              # 4-fold either way
                 source = "FILL IN: allometric scaling accuracy")
)

# --- 3. PD (leave commented for the PK-only version) ------------------------
# For the next version, switch MODEL$pd to a simple time course and add a PD
# prior. PD is a plain shape with no exposure coupling: "constant", "linear"
# (needs slope), or "exponential" (needs plateau and rate).
#
# MODEL <- pmx_structural_model(
#   pk = "1cmt_oral", typical = c(cl = 10, v = 70, ka = 1,
#                                 baseline = 100, plateau = 60, rate = 0.05),
#   pd = "exponential", source = "...", iiv = c(cl = 0.3, v = 0.2),
#   residual_cv = 0.15
# )
# PRIORS <- pmx_priors(
#   pk = pmx_prior(c(1 / 4, 4), "..."),
#   pd = pmx_prior(c(1 / 4, 4), "baseline literature")   # baseline priors are tight
# )

# --- 4. Covariate columns (optional) ----------------------------------------
# Baseline covariates carried into the output so covariate-handling pipeline
# code has columns to run against. Set COVARIATES <- NULL to skip.
#
# Simplest: just list the column names. Values are resampled from the data --
# uniform over the (clipped) observed range for continuous columns, proportional
# resample for categorical -- the same approach as Novartis's synadam. The
# default 1st/99th-percentile clip keeps the two extreme patients from leaking;
# pass clip = NULL to expose the raw min/max like synadam does. This is NOT
# differentially private and is for trusted-environment use only.
COVARIATES <- pmx_covariates_auto(
  c("WT", "AGE", "SEX", "RACE", "EGFR")     # every covariate column you want
  # , clip = NULL                            # uncomment to match synadam exactly
)

# Alternative, fully DP but you cite a public range/levels per column, and each
# costs one privacy slice. Prefer this if the covariates may leave the safe
# environment. Uncomment to use instead:
# COVARIATES <- pmx_covariates(
#   WT  = pmx_covariate(range = c(40, 120), source = "inclusion criteria"),
#   SEX = pmx_covariate(levels = c("M", "F"), source = "protocol")
# )

# --- 5. Column roles: map YOUR columns onto PMX roles -----------------------
# Only the roles the real dataset has. Ignored when DRY_RUN = TRUE.
ROLES <- pmx_roles(
  id   = "ID",  time = "TIME", dv  = "DV", amt = "AMT", evid = "EVID",
  dvid = "DVID", cmt = "CMT",  mdv = "MDV",
  occasion = "OCC", exclude = NULL
)

# Optional per-dataset cleanup: recode DVID to "cp" (and later "pd"), fix units,
# drop screening rows. Receives the raw frame, returns a cleaned one.
PREP <- function(d) {
  # d$DVID <- ifelse(d$DVID == 1, "cp", NA)   # example
  d
}

# ============================================================================
# WORKFLOW  --  you should not need to edit below here
# ============================================================================

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Read the source data (real, or dry-run stand-in) ----------------------
read_source <- function() {
  if (DRY_RUN) {
    # A public stand-in dataset with a known true clearance, so we can confirm
    # the pipeline recovers it before ever touching real data. Same structure
    # as MODEL but with a clearance 2.2x higher than the prediction.
    truth <- pmx_structural_model(
      "2cmt_oral", c(cl = 22, v = 30, q = 15, v2 = 100, ka = 1),
      source = "dry-run stand-in truth"
    )
    stand_in <- synpmx_prior(truth, DESIGN, n_subjects = 40, seed = 1)
    # Attach fake covariate columns matching COVARIATES, so the calibrated fit
    # has something to summarize. Names containing SEX/RACE/SEXN are treated as
    # categorical; everything else as continuous.
    if (!is.null(COVARIATES)) {
      ids <- unique(stand_in$ID)
      set.seed(99)
      for (nm in names(COVARIATES)) {
        per_id <- if (grepl("SEX|RACE|ETHNIC|ARM", nm, ignore.case = TRUE)) {
          sample(c("A", "B"), length(ids), replace = TRUE)
        } else {
          round(stats::rnorm(length(ids), 70, 15))
        }
        stand_in[[nm]] <- per_id[match(stand_in$ID, ids)]
      }
    }
    return(stand_in)
  }
  if (!file.exists(DATA_PATH)) {
    stop("DATA_PATH not found: ", normalizePath(DATA_PATH, mustWork = FALSE))
  }
  utils::read.csv(DATA_PATH, stringsAsFactors = FALSE)
}

raw   <- PREP(read_source())
roles <- if (DRY_RUN) pmx_generated_roles() else ROLES

message("\n== Source structural validation (restricted) ==")
report <- validate_pmx(raw, roles, strict = FALSE)
if (!report$valid) {
  print(report)   # names each problem, the role, and the column it maps to
  stop("Source failed structural validation. See the numbered problems above; ",
       "each says which role and column to fix in `pmx_roles()`.", call. = FALSE)
}
message("  valid; ", length(unique(raw[[roles$id]])), " subjects")

n_subjects <- length(unique(raw[[roles$id]]))

# ---- PRIOR version: public inputs only, no data, no budget -----------------
# Matched to the source cohort size so the two versions are comparable. Prior
# mode reads no data, so bootstrap covariates (which resample the data) cannot
# be generated here; they are added only to the calibrated version below. DP
# covariates with a public range/levels work in both.
message("\n== PRIOR version (no data, no privacy budget) ==")
prior_synthetic <- tryCatch(
  synpmx_prior(MODEL, DESIGN, n_subjects = n_subjects, seed = SEED,
               covariates = COVARIATES),
  error = function(e) {
    message("  (bootstrap covariates need data; omitting them from PRIOR)")
    synpmx_prior(MODEL, DESIGN, n_subjects = n_subjects, seed = SEED)
  }
)
stopifnot(validate_pmx(prior_synthetic, pmx_generated_roles())$valid)
message("  generated ", nrow(prior_synthetic), " rows for ",
        length(unique(prior_synthetic$ID)), " subjects")

# ---- Pre-flight: is the calibrated release worth its budget? ---------------
message("\n== Pre-flight (no data read, no budget) ==")
print(pmx_preflight(PRIORS, epsilon = EPSILON, n_subjects = n_subjects,
                    covariates = COVARIATES))

# ---- CALIBRATED version: the only budget-spending step ---------------------
message("\n== CALIBRATED version ==")

# synpmx_calibrated() is one of the DP engines: complete and tested, but not
# under active development, carrying known open findings, and not
# independently privacy-audited (design/REVIEW_BACKLOG.md REV-023). It
# refuses to run against a real backend without this acknowledgment, which
# does not persist across runs of this script -- that is deliberate.
synpmx_enable_dp_engines()

calibrated_synthetic <- synpmx_calibrated(
  raw, roles, MODEL, DESIGN, PRIORS,
  epsilon       = EPSILON,
  covariates    = COVARIATES,
  seed          = SEED,
  backend       = if (DRY_RUN) "public" else "opendp",
  public_source = DRY_RUN          # NEVER TRUE for confidential data
)

# The release travels with the dataset, so the record is always to hand. Draw
# further datasets with synpmx_generate(calibrated_synthetic, seed = ...);
# calling synpmx_calibrated() again would spend the budget a second time.
release <- attr(calibrated_synthetic, "synpmx_release")
print(release)
message("\nProvenance (for the release record):")
print(release$provenance)
stopifnot(validate_pmx(calibrated_synthetic, pmx_generated_roles())$valid)

# ---- Save outputs (gitignored) ---------------------------------------------
utils::write.csv(prior_synthetic,
                 file.path(OUT_DIR, "synthetic_prior.csv"), row.names = FALSE)
utils::write.csv(calibrated_synthetic,
                 file.path(OUT_DIR, "synthetic_calibrated.csv"), row.names = FALSE)
saveRDS(fit, file.path(OUT_DIR, "calibrated_model.rds"))

# ---- RESTRICTED diagnostic: source vs both versions ------------------------
# Source-derived. Keep it in the safe environment; do not export these figures.
message("\n== Restricted comparison (do not export) ==")
if (requireNamespace("ggplot2", quietly = TRUE)) {
  obs <- function(d, r, label) {
    keep <- d[[r$evid]] == 0 & !is.na(d[[r$dv]])
    data.frame(dataset = label,
               time = as.numeric(d[[r$time]][keep]),
               dv   = as.numeric(d[[r$dv]][keep]))
  }
  gr <- pmx_generated_roles()
  compare <- rbind(
    obs(raw, roles, "source"),
    obs(prior_synthetic, gr, "prior"),
    obs(calibrated_synthetic, gr, "calibrated")
  )
  compare$dataset <- factor(compare$dataset,
                            levels = c("source", "prior", "calibrated"))
  p <- ggplot2::ggplot(compare, ggplot2::aes(time, dv, colour = dataset)) +
    ggplot2::geom_point(alpha = 0.4, size = 0.7) +
    ggplot2::facet_wrap(~dataset, ncol = 1) +
    ggplot2::labs(title = "RESTRICTED: source vs prior vs calibrated",
                  subtitle = "do not export", x = "Time", y = "DV") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
  ggplot2::ggsave(file.path(OUT_DIR, "RESTRICTED_comparison.png"), p,
                  width = 7, height = 8, dpi = 110)
  message("  wrote ", file.path(OUT_DIR, "RESTRICTED_comparison.png"))
} else {
  message("  ggplot2 not available; skipping the plot")
}

message("\nDone. Compare the prior and calibrated versions, and check the ",
        "pre-flight verdict and correction factor before trusting output.")
