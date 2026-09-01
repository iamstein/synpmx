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
