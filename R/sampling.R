#' Summarize the fitted sampling design
#'
#' Returns post-processed, privacy-accounted sampling summaries from a fitted
#' model. For dose-relative and occasion endpoints, the summary separates the
#' probability that an occasion is sampled from the mean number of observations
#' conditional on sampling. It never consults the source data.
#'
#' @param private_model A model returned by [fit_private_pmx()].
#'
#' @return A data frame with one row per endpoint and possible dose occasion.
#' @export
sampling_summary <- function(private_model) {
  validate_private_model(private_model, strict = TRUE)
  regimen <- .resolved_regimen(private_model)
  rows <- list()
  for (name in names(private_model$public$endpoints)) {
    endpoint <- private_model$public$endpoints[[name]]
    if (!endpoint$alignment %in% c("dose_relative", "occasion")) next
    timing <- private_model$population$timing[[name]]
    schedule <- private_model$public$design$endpoint_occasion_grids[[name]]
    occasions <- if (regimen$n_doses > 0L) regimen$n_doses else max(
      length(timing$occasion_observation_count),
      if (is.null(schedule)) 0L else max(as.integer(names(schedule)))
    )
    if (occasions < 1L) next
    if (is.null(schedule)) {
      probability <- .inferred_occasion_presence(timing, occasions)
      count <- .inferred_occasion_counts(timing, occasions)
      basis <- "privacy_accounted_inference"
    } else {
      probability <- as.numeric(seq_len(occasions) %in%
                                  as.integer(names(schedule)))
      count <- vapply(seq_len(occasions), function(occasion) {
        length(schedule[[as.character(occasion)]])
      }, integer(1))
      basis <- "public_schedule"
    }
    rows[[length(rows) + 1L]] <- data.frame(
      endpoint = name,
      occasion = seq_len(occasions),
      sampling_probability = probability,
      observations_if_sampled = count,
      expected_observations = probability * count,
      basis = basis,
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) {
    return(data.frame(
      endpoint = character(), occasion = integer(),
      sampling_probability = numeric(), observations_if_sampled = numeric(),
      expected_observations = numeric(), basis = character(),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, rows)
}
