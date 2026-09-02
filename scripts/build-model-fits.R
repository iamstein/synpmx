# Build the stored `pmx_fitted_model` objects the model-generator vignettes
# knit against, so that `R CMD check` and the pkgdown site never compile a
# population model.
#
# Run this by hand on a machine whose R can link a compiled model, after any
# change to `synpmx_model_estimate()` that moves the numbers. The vignettes read
# what it writes and will otherwise document a fit that no longer exists.
#
#   Rscript scripts/build-model-fits.R
#
# On macOS, R's `FLIBS` points at `/opt/gfortran`. Where that directory is
# absent, every model compiles and none of them links, and the failure reports
# itself as a missing C compiler. Install the official gfortran build for your
# platform, or put an empty `FLIBS=` in `~/.R/Makevars`.

devtools::load_all(".", quiet = TRUE)

if (!requireNamespace("nlmixr2data", quietly = TRUE)) {
  stop("nlmixr2data is needed to build the stored fits.")
}

warfarin <- as.data.frame(nlmixr2data::warfarin)
# `warfarin`'s 16 recorded observation times are the protocol's, so declaring
# `nominal_time` from `time` asserts something true about this study rather than
# constructing a grid. The public-data survey makes the same declaration.
warfarin$ntime <- warfarin$time
warfarin_roles <- pmx_roles(
  id = "id", time = "time", nominal_time = "ntime", dv = "dv", amt = "amt",
  evid = "evid", dvid = "dvid", covariates = c("wt", "age", "sex")
)
warfarin_fit <- synpmx_model_estimate(warfarin, warfarin_roles, seed = 1)

dir.create("inst/extdata", showWarnings = FALSE, recursive = TRUE)
saveRDS(warfarin_fit, "inst/extdata/warfarin-model-fit.rds", version = 2)
message("wrote inst/extdata/warfarin-model-fit.rds")
print(model_report(warfarin_fit))

# The methods survey runs all six generators on `theo_md`, which has twelve
# subjects -- below the floor `synpmx_model_estimate()` sets for itself. The
# floor is lowered here deliberately, and the article says so where it reads
# this fit.
snap_to <- function(x, grid) {
  grid[max.col(-abs(outer(x, grid, "-")), ties.method = "first")]
}
theo_md <- as.data.frame(nlmixr2data::theo_md)
theo_doses <- seq(0, 144, by = 24)
theo_samples <- c(0, 0.25, 0.5, 1, 2, 3, 4, 5, 7, 9, 12, 24)
interval <- pmax(1L, findInterval(theo_md$TIME, theo_doses))
theo_md$NTIME <- ifelse(
  theo_md$EVID == 0,
  theo_doses[interval] +
    snap_to(theo_md$TIME - theo_doses[interval], theo_samples),
  theo_md$TIME
)
theo_roles <- pmx_roles(
  id = "ID", time = "TIME", nominal_time = "NTIME", dv = "DV", amt = "AMT",
  evid = "EVID", cmt = "CMT", covariates = "WT"
)
theo_fit <- synpmx_model_estimate(theo_md, theo_roles, seed = 1,
                                  min_subjects = 12L)

saveRDS(theo_fit, "inst/extdata/theo-md-model-fit.rds", version = 2)
message("wrote inst/extdata/theo-md-model-fit.rds")
print(model_report(theo_fit))

if (!requireNamespace("xgxr", quietly = TRUE)) {
  stop("xgxr is needed to build the demo fit.")
}

# The study the three demos share, declared exactly as they declare it: 180
# patients, six arms, a concentration endpoint and a continuous PD one. The
# `CENS` flag is cleared on the PD rows, where the source sets it on values
# that are not below any limit.
case1 <- as.data.frame(get(utils::data(list = "case1_pkpd", package = "xgxr")))
case1$CENS <- ifelse(case1$NAME == "PD - Continuous", 0, case1$CENS)
case1_roles <- pmx_roles(
  id = "ID", time = "TIME", nominal_time = "NOMTIME", dv = "LIDV",
  cens = "CENS", amt = "AMT", evid = "EVID", cmt = "CMT", dvid = "NAME",
  strata = c("TRTACT", "DOSE"), covariates = "WEIGHTB", keep = "STUDY"
)
case1_fit <- synpmx_model_estimate(case1, case1_roles, seed = 1)

saveRDS(case1_fit, "inst/extdata/case1-pkpd-model-fit.rds", version = 2)
message("wrote inst/extdata/case1-pkpd-model-fit.rds")
print(model_report(case1_fit))

# The public-data survey. One fit per study, including the three that need
# something declared before they will run: `theo_md` and `nimoData` are below
# the subject floor, and `nimoData` and `pheno_sd` dose every patient the same
# way, so nothing in either reads as dose-proportional and the concentration is
# named by hand. The article shows each refusal before the declaration.
snap_to_grid <- snap_to

mad <- as.data.frame(get(utils::data(list = "mad", package = "xgxr")))
mad_roles <- pmx_roles(
  id = "ID", time = "TIME", dv = "LIDV", amt = "AMT", evid = "EVID",
  cmt = "CMT", dvid = "NAME", mdv = "MDV", nominal_time = "NOMTIME",
  strata = c("TRTACT", "DOSE"), covariates = c("WEIGHTB", "SEX")
)
saveRDS(synpmx_model_estimate(mad, mad_roles, seed = 1),
        "inst/extdata/mad-model-fit.rds", version = 2)
message("wrote inst/extdata/mad-model-fit.rds")

wbcSim <- as.data.frame(nlmixr2data::wbcSim)
wbcSim$NTIME <- wbcSim$TIME
wbc_roles <- pmx_roles(
  id = "ID", time = "TIME", nominal_time = "NTIME", dv = "DV", amt = "AMT",
  evid = "EVID", cmt = "CMT", rate = "RATE"
)
saveRDS(synpmx_model_estimate(wbcSim, wbc_roles, seed = 1),
        "inst/extdata/wbcsim-model-fit.rds", version = 2)
message("wrote inst/extdata/wbcsim-model-fit.rds")

mavoglurant <- as.data.frame(nlmixr2data::mavoglurant)
mavo_design <- c(0, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 10, 12, 24, 36, 48)
mavoglurant$NTIME <- ifelse(mavoglurant$EVID == 0,
                            snap_to_grid(mavoglurant$TIME, mavo_design),
                            mavoglurant$TIME)
mavo_roles <- pmx_roles(
  id = "ID", time = "TIME", nominal_time = "NTIME", dv = "DV", amt = "AMT",
  evid = "EVID", cmt = "CMT", rate = "RATE", mdv = "MDV", occasion = "OCC",
  keep = "DOSE", covariates = c("AGE", "SEX", "WT", "HT")
)
saveRDS(synpmx_model_estimate(mavoglurant, mavo_roles, seed = 1),
        "inst/extdata/mavoglurant-model-fit.rds", version = 2)
message("wrote inst/extdata/mavoglurant-model-fit.rds")

nimoData <- as.data.frame(nlmixr2data::nimoData)
nimo_interval <- 168
last_occasion <- ave(nimoData$OCC, nimoData$ID, FUN = max)
nominal_tad <- round(nimoData$TAD / 24) * 24
pre_dose <- nimoData$EVID == 0 & nominal_tad >= nimo_interval &
  nimoData$OCC < last_occasion
nominal_tad[pre_dose] <- nimo_interval - 1
nimoData$NTIME <- (nimoData$OCC - 1) * nimo_interval +
  ifelse(nimoData$EVID == 0, nominal_tad, 0)
nimo_roles <- pmx_roles(
  id = "ID", time = "TIME", nominal_time = "NTIME", dv = "DV", amt = "AMT",
  evid = "EVID", rate = "RATE", mdv = "MDV", tad = "TAD", occasion = "OCC",
  covariates = c("BSA", "AGE", "HGT"), keep = "DOS"
)
saveRDS(synpmx_model_estimate(nimoData, nimo_roles, seed = 1,
                              min_subjects = 12L,
                              endpoint_roles = c(pk = "DV")),
        "inst/extdata/nimo-model-fit.rds", version = 2)
message("wrote inst/extdata/nimo-model-fit.rds")

pheno_sd <- as.data.frame(nlmixr2data::pheno_sd)
pheno_sd$NTIME <- pheno_sd$TIME
pheno_roles <- pmx_roles(
  id = "ID", time = "TIME", nominal_time = "NTIME", dv = "DV", amt = "AMT",
  evid = "EVID", mdv = "MDV", covariates = c("WT", "APGR")
)
saveRDS(synpmx_model_estimate(pheno_sd, pheno_roles, seed = 1,
                              endpoint_roles = c(pk = "DV")),
        "inst/extdata/pheno-model-fit.rds", version = 2)
message("wrote inst/extdata/pheno-model-fit.rds")
