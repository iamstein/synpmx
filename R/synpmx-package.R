#' synpmx: structurally faithful synthetic pharmacometric data
#'
#' `synpmx` primarily uses AVATAR-style profile blending to generate
#' structurally coherent PMX synthetic datasets for use within the source data's
#' own access controls and confidentiality obligations — which follow the data
#' rather than the machine, so a local workstation under those controls is a
#' supported destination. The package also provides a separate calibrated
#' structural workflow with a formal differential-privacy backend for output
#' that will reach anyone the source data could not.
#'
#' The output is intended for model-workflow exploration. It preserves broad
#' magnitudes, scientific clocks, event conventions, and coarse variability;
#' it deliberately does not preserve parameter estimates, detailed
#' distributions, inferential validity, or scientific conclusions.
#'
#' @keywords internal
"_PACKAGE"

utils::globalVariables(c(
  "dataset_plot", "subject_plot", "time_plot", "dv_plot",
  "endpoint_plot"
))
