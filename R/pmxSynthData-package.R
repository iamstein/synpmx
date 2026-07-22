#' pmxSynthData: source-calibrated private mock pharmacometric data
#'
#' `pmxSynthData` fits a low-dimensional, subject-level differentially private
#' population model once inside a restricted environment, then generates any
#' number of structurally coherent PMX mock datasets by post-processing. It
#' constructs new event tables rather than copying patient schedules or
#' retaining source-subject records.
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
