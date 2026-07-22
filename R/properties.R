#' Summarize fitted subject properties and associated regimens
#'
#' Subject properties are categorical treatment assignments or grouping
#' variables declared through [pmx_roles()], such as `ACTARM`, `TRT`, or a
#' nominal dose group. Their released probabilities and property-conditioned
#' regimen summaries are source-dependent and therefore part of the fitted
#' model's privacy accounting.
#'
#' @param private_model A fitted model from [fit_private_pmx()].
#'
#' @return A data frame with one row per declared property stratum. It has zero
#'   rows when no subject properties were declared.
#' @export
subject_property_summary <- function(private_model) {
  validate_private_model(private_model, strict = TRUE)
  properties <- private_model$population$subject_properties
  metric_names <- c(
    "probability", "released_count", "n_doses", "dose_interval",
    "dose_amount", "dose_rate", "infusion_duration", "observation_count"
  )
  if (is.null(properties) || !length(properties$strata)) {
    out <- stats::setNames(
      replicate(length(metric_names), numeric(), simplify = FALSE),
      metric_names
    )
    return(as.data.frame(out, stringsAsFactors = FALSE))
  }

  rows <- lapply(properties$strata, function(stratum) {
    values <- as.data.frame(stratum$values, stringsAsFactors = FALSE)
    event <- stratum$event
    cbind(
      values,
      data.frame(
        probability = stratum$probability,
        released_count = stratum$released_count,
        n_doses = event$n_doses,
        dose_interval = event$dose_interval,
        dose_amount = event$dose_amount,
        dose_rate = event$dose_rate,
        infusion_duration = event$infusion_duration,
        observation_count = event$observation_count,
        stringsAsFactors = FALSE
      )
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  schema <- private_model$public$design$schema
  for (name in properties$names) {
    out[[name]] <- .cast_public_column(
      out[[name]], .schema_column(schema, name)
    )
  }
  out
}
