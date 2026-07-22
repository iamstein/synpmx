.structural_summary <- function(data, roles, label) {
  observation <- .is_zero(data[[roles$evid]])
  allowed <- .observation_rows(data, roles)
  present <- allowed & !is.na(data[[roles$dv]])
  endpoint <- .endpoint(data, roles)
  data.frame(
    dataset = label,
    rows = nrow(data),
    subjects = length(unique(data[[roles$id]])),
    event_rows = sum(!.is_zero(data[[roles$evid]])),
    observation_rows = sum(observation),
    observed_dv = sum(present),
    missing_observations = sum(observation & is.na(data[[roles$dv]])),
    endpoints = paste(sort(unique(endpoint[observation])), collapse = ", "),
    stringsAsFactors = FALSE
  )
}

.event_counts <- function(data, roles, label) {
  endpoint <- .endpoint(data, roles)
  event <- as.character(data[[roles$evid]])
  cmt <- if (is.null(roles$cmt)) "<absent>" else
    as.character(data[[roles$cmt]])
  event[is.na(event)] <- "<missing>"
  cmt[is.na(cmt)] <- "<missing>"
  table_data <- data.frame(
    dataset = label,
    endpoint = endpoint,
    evid = event,
    cmt = cmt,
    stringsAsFactors = FALSE
  )
  result <- stats::aggregate(
    rep(1L, nrow(table_data)), table_data, FUN = sum
  )
  names(result)[ncol(result)] <- "rows"
  result
}

.plot_data <- function(data, roles, label, observations_only = TRUE) {
  selected <- if (observations_only) {
    .observation_rows(data, roles, require_present = TRUE)
  } else {
    rep(TRUE, nrow(data))
  }
  data.frame(
    dataset_plot = label,
    subject_plot = as.character(data[[roles$id]][selected]),
    time_plot = as.numeric(data[[roles$time]][selected]),
    dv_plot = as.numeric(data[[roles$dv]][selected]),
    endpoint_plot = .endpoint(data, roles)[selected],
    row_type_plot = ifelse(
      .is_zero(data[[roles$evid]][selected]), "observation", "event"
    ),
    stringsAsFactors = FALSE
  )
}

.comparison_plots <- function(source, mock, roles) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(list())
  trajectories <- rbind(
    .plot_data(source, roles, "source"),
    .plot_data(mock, roles, "mock")
  )
  overlay <- ggplot2::ggplot(
    trajectories,
    ggplot2::aes(
      x = time_plot, y = dv_plot,
      group = interaction(dataset_plot, subject_plot),
      colour = dataset_plot
    )
  ) +
    ggplot2::geom_line(alpha = 0.35, linewidth = 0.45) +
    ggplot2::facet_wrap(~endpoint_plot, scales = "free_y") +
    ggplot2::labs(
      x = roles$time, y = roles$dv, colour = "Dataset",
      title = "Source and mock longitudinal trajectories",
      subtitle = "Structural exploration only; not a statistical equivalence check"
    ) +
    ggplot2::theme_minimal()

  faceted <- ggplot2::ggplot(
    trajectories,
    ggplot2::aes(
      x = time_plot, y = dv_plot, group = subject_plot
    )
  ) +
    ggplot2::geom_line(alpha = 0.4, linewidth = 0.45) +
    ggplot2::facet_grid(dataset_plot ~ endpoint_plot, scales = "free_y") +
    ggplot2::labs(
      x = roles$time, y = roles$dv,
      title = "Individual trajectories by dataset and endpoint"
    ) +
    ggplot2::theme_minimal()

  list(overlay = overlay, faceted = faceted)
}

#' Compare source and mock PMX structures
#'
#' Produces concise structural tables and, when `ggplot2` is installed,
#' exploratory trajectory plots. This function does not test or claim
#' distributional, scientific, or statistical equivalence.
#'
#' @param source Source PMX data.
#' @param mock Mock PMX data from [mock_pmx()].
#' @param roles Explicit roles from [pmx_roles()].
#'
#' @return A `pmx_comparison` list containing `summary`, `event_counts`,
#'   `column_classes`, validation reports, and a named `plots` list.
#' @export
#'
#' @examples
#' source <- data.frame(
#'   ID = rep(1:3, each = 3), TIME = rep(c(0, 1, 2), 3),
#'   DV = c(0, 2, 1, 0, 3, 1, 0, 4, 2),
#'   AMT = rep(c(100, 0, 0), 3), EVID = rep(c(1L, 0L, 0L), 3)
#' )
#' roles <- pmx_roles("ID", "TIME", "DV", "AMT", "EVID")
#' mock <- mock_pmx(source, roles, n_subjects = 2)
#' comparison <- compare_pmx(source, mock, roles)
#' comparison$summary
compare_pmx <- function(source, mock, roles) {
  .assert_roles(source, roles)
  .assert_roles(mock, roles)
  source_validation <- validate_pmx(source, roles)
  mock_validation <- validate_pmx(mock, roles)
  source_classes <- vapply(source, function(column) {
    paste(class(column), collapse = "/")
  }, character(1))
  mock_classes <- vapply(mock, function(column) {
    paste(class(column), collapse = "/")
  }, character(1))
  column_classes <- data.frame(
    column = names(source),
    source = unname(source_classes),
    mock = unname(mock_classes[names(source)]),
    matches = unname(source_classes == mock_classes[names(source)]),
    stringsAsFactors = FALSE
  )

  structure(
    list(
      summary = rbind(
        .structural_summary(source, roles, "source"),
        .structural_summary(mock, roles, "mock")
      ),
      event_counts = rbind(
        .event_counts(source, roles, "source"),
        .event_counts(mock, roles, "mock")
      ),
      column_classes = column_classes,
      validation = list(source = source_validation, mock = mock_validation),
      plots = .comparison_plots(source, mock, roles)
    ),
    class = "pmx_comparison"
  )
}

#' @export
print.pmx_comparison <- function(x, ...) {
  cat("PMX structural comparison\n")
  print(x$summary, row.names = FALSE)
  if (length(x$plots)) {
    cat("\nExploratory plots: ", paste(names(x$plots), collapse = ", "), "\n",
        sep = "")
  }
  invisible(x)
}
