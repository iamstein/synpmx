#' pmxSynthData: Structurally Faithful Mock Pharmacometric Data
#'
#' `pmxSynthData` creates mock data for model-workflow exploration by combining a
#' sampled subject's coherent event template with AVATAR-like blends of
#' compatible subjects' covariates and longitudinal measurements. It does not
#' fit or call a PK, PD, or NLME model.
#'
#' The method is an AVATAR-inspired variant adapted to longitudinal PMX event
#' data; it is not an exact reproduction of published AVATAR software. The
#' output is not anonymous data and carries no formal privacy guarantee.
#'
#' @keywords internal
"_PACKAGE"

utils::globalVariables(c(
  "dataset_plot", "subject_plot", "time_plot", "dv_plot",
  "endpoint_plot"
))
