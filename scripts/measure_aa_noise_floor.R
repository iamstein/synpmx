# How small a privacy loss can adversarial accuracy resolve?
#
# `compare_pmx_proximity()` reports a null interval: what AA does between two
# halves of the same real cohort, where nothing is wrong by construction. That
# interval is the noise floor. This script measures it on three public datasets
# and converts it into the quantity that matters, which is how many memorized
# patients it would take to move AA clear of the floor.
#
# The conversion. If a fraction p of the compared set has a synthetic near-copy,
# those subjects' indicators are 0 rather than 1/2, so AA ~ (1 - p)/2 against a
# clean 1/2, and the shift is p/2. Two thresholds are reported because two
# different tests are in use and they differ by about a factor of two:
#
#   detect_vs_null  the synpmx test, where observed AA must fall below
#                   `null_lower`, a drop of about half the interval width.
#                   Needs p > width.
#   detect_vs_diff  the Yale holdout test, where AA(C,S) - AA(T,S) must clear
#                   the noise in a difference of two AAs. Conservative: needs
#                   p > 2 * width.
#
# Numbers from this script appear in the "noise floor" table of
# vignettes/articles/synthetic-data-checking-review.Rmd. Rerun it and update
# that table if the profile space, the null, or the defaults change.
#
# Install synpmx first with: R CMD INSTALL .

if (!requireNamespace("synpmx", quietly = TRUE)) {
  stop("Install synpmx before running this script: R CMD INSTALL .")
}
if (!requireNamespace("nlmixr2data", quietly = TRUE)) {
  stop("Install nlmixr2data to run this script.")
}
if (!requireNamespace("xgxr", quietly = TRUE)) {
  stop("Install xgxr to run this script.")
}
library(synpmx)

noise_floor <- function(label, data, roles, seed = 1) {
  synthetic <- suppressWarnings(synpmx_avatar(data, roles, seed = seed))
  p <- compare_pmx_proximity(data, synthetic, roles)
  width <- p$null_upper - p$null_lower
  data.frame(
    dataset = label,
    source_n = length(unique(data[[roles$id]])),
    compared_n = p$n_compared,
    observed_aa = round(p$adversarial_accuracy, 3),
    null_width = round(width, 3),
    detect_vs_null = round(width, 2),
    detect_vs_diff = round(2 * width, 2),
    patients_vs_null = ceiling(width * p$n_compared),
    one_patient_delta = round(1 / (2 * p$n_compared), 4),
    verdict = p$verdict,
    stringsAsFactors = FALSE
  )
}

results <- rbind(
  noise_floor(
    "theo_md", nlmixr2data::theo_md,
    pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
              evid = "EVID", cmt = "CMT", covariates = "WT")
  ),
  noise_floor(
    "warfarin", nlmixr2data::warfarin,
    pmx_roles(id = "id", time = "time", dv = "dv", amt = "amt",
              evid = "evid", dvid = "dvid",
              covariates = c("wt", "age", "sex"))
  ),
  noise_floor(
    "case1_pkpd", xgxr::case1_pkpd,
    pmx_roles(id = "ID", time = "TIME", dv = "LIDV", amt = "AMT",
              evid = "EVID", cmt = "CMT", dvid = "NAME",
              nominal_time = "NOMTIME", strata = c("TRTACT", "DOSE"),
              covariates = "WEIGHTB", keep = "STUDY")
  )
)

print(results[, setdiff(names(results), "verdict")], row.names = FALSE)
cat("\nVerdicts:\n")
for (i in seq_len(nrow(results))) {
  cat(sprintf("  %-11s %s\n", results$dataset[i], results$verdict[i]))
}
