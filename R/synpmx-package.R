#' synpmx: structurally faithful synthetic pharmacometric data
#'
#' `synpmx` primarily uses AVATAR-style profile blending to generate
#' structurally coherent PMX synthetic datasets in a trusted environment. The
#' package also provides a separate calibrated structural workflow with a
#' formal differential-privacy backend when that boundary is required.
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
