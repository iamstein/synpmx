.structural_summary <- function(data, roles, label) {
  allowed <- .observation_rows(data, roles)
  endpoint <- .endpoint(data, roles)
  data.frame(
    dataset = label, rows = nrow(data),
    subjects = length(unique(data[[roles$id]])),
    event_rows = sum(.event_rows(data, roles)),
    observation_rows = sum(allowed),
    observed_dv = sum(allowed & !is.na(data[[roles$dv]])),
    endpoints = paste(sort(unique(endpoint[allowed])), collapse = ", "),
    stringsAsFactors = FALSE
  )
}

.event_counts <- function(data, roles, label) {
  table_data <- data.frame(
    dataset = label,
    endpoint = .endpoint(data, roles),
    evid = as.character(data[[roles$evid]]),
    cmt = if (is.null(roles$cmt)) "<absent>" else
      as.character(data[[roles$cmt]]),
    stringsAsFactors = FALSE
  )
  table_data[is.na(table_data)] <- "<missing>"
  result <- stats::aggregate(rep(1L, nrow(table_data)), table_data, sum)
  names(result)[ncol(result)] <- "rows"
  result
}

.plot_data <- function(data, roles, label) {
  selected <- .observation_rows(data, roles, require_present = TRUE)
  data.frame(
    dataset_plot = label,
    subject_plot = as.character(data[[roles$id]][selected]),
    time_plot = as.numeric(data[[roles$time]][selected]),
    dv_plot = as.numeric(data[[roles$dv]][selected]),
    endpoint_plot = .endpoint(data, roles)[selected],
    stringsAsFactors = FALSE
  )
}

.comparison_plots <- function(source, synthetic, roles) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(list())
  trajectories <- rbind(.plot_data(source, roles, "source"),
                        .plot_data(synthetic, roles, "synthetic"))
  overlay <- ggplot2::ggplot(
    trajectories,
    ggplot2::aes(x = time_plot, y = dv_plot,
                 group = interaction(dataset_plot, subject_plot),
                 colour = dataset_plot)
  ) +
    ggplot2::geom_line(alpha = 0.35, linewidth = 0.45) +
    ggplot2::facet_wrap(~endpoint_plot, scales = "free_y") +
    ggplot2::labs(
      x = roles$time, y = roles$dv, colour = "Dataset",
      title = "Restricted source-versus-synthetic trajectory diagnostic",
      subtitle = "Not releasable unless separately privatized and budgeted"
    ) + ggplot2::theme_minimal()
  faceted <- ggplot2::ggplot(
    trajectories,
    ggplot2::aes(x = time_plot, y = dv_plot, group = subject_plot)
  ) +
    ggplot2::geom_line(alpha = 0.4, linewidth = 0.45) +
    ggplot2::facet_grid(dataset_plot ~ endpoint_plot, scales = "free_y") +
    ggplot2::labs(x = roles$time, y = roles$dv,
                  title = "Restricted individual-trajectory diagnostic") +
    ggplot2::theme_minimal()
  list(overlay = overlay, faceted = faceted)
}

.mark_release <- function(x, status) {
  if (is.null(x)) return(NULL)
  attr(x, "release_status") <- status
  x
}

#' Compare source and generated PMX structures inside the restricted environment
#'
#' Any component that uses `source` is marked
#' `"restricted_not_releasable"`. A fitted private model does not make a new
#' source-derived comparison private; releasing such a diagnostic requires a
#' separate public justification or budgeted DP mechanism.
#'
#' @param source Source PMX data.
#' @param synthetic Generated synthetic PMX data.
#' @param roles Explicit roles from [pmx_roles()].
#' @param endpoints Optional endpoint declarations.
#'
#' @return A `pmx_comparison` containing component-level release metadata.
#' @export
compare_pmx <- function(source, synthetic, roles, endpoints = NULL) {
  .assert_roles(source, roles)
  .assert_roles(synthetic, roles)
  source_validation <- validate_pmx(source, roles, endpoints)
  synthetic_validation <- validate_pmx(synthetic, roles, endpoints)
  source_classes <- vapply(source, function(x) paste(class(x), collapse = "/"),
                           character(1))
  synthetic_classes <- vapply(
    synthetic, function(x) paste(class(x), collapse = "/"), character(1)
  )
  column_classes <- data.frame(
    column = names(source), source = unname(source_classes),
    synthetic = unname(synthetic_classes[names(source)]),
    matches = unname(source_classes == synthetic_classes[names(source)]),
    stringsAsFactors = FALSE
  )
  status <- data.frame(
    component = c("summary", "event_counts", "column_classes",
                  "validation.source", "validation.synthetic", "plots"),
    release_status = c(
      rep("restricted_not_releasable", 4L), "releasable_post_processing",
      "restricted_not_releasable"
    ),
    stringsAsFactors = FALSE
  )
  structure(list(
    summary = .mark_release(rbind(
      .structural_summary(source, roles, "source"),
      .structural_summary(synthetic, roles, "synthetic")
    ), "restricted_not_releasable"),
    event_counts = .mark_release(rbind(
      .event_counts(source, roles, "source"),
      .event_counts(synthetic, roles, "synthetic")
    ), "restricted_not_releasable"),
    column_classes = .mark_release(column_classes,
                                   "restricted_not_releasable"),
    validation = list(
      source = .mark_release(source_validation,
                             "restricted_not_releasable"),
      synthetic = .mark_release(synthetic_validation,
                                "releasable_post_processing")
    ),
    plots = .mark_release(.comparison_plots(source, synthetic, roles),
                          "restricted_not_releasable"),
    release_status = status
  ), class = "pmx_comparison")
}

#' @export
print.pmx_comparison <- function(x, ...) {
  cat("Restricted PMX source-versus-synthetic comparison\n")
  print(x$summary, row.names = FALSE)
  cat("Source-derived components are not releasable unless separately public or privately budgeted.\n")
  invisible(x)
}

# Distributional summaries -----------------------------------------------------
#
# compare_pmx() answers "is the structure the same?" (schema, event grammar, row
# counts). compare_pmx_distributions() answers "are the numbers in the same
# range?" -- the per-covariate and per-endpoint sanity check a user eyeballs
# right after generating data.

# n / mean / sd / min / quartiles / max for one numeric vector, in the long
# layout the rest of this file uses (a `dataset` column, rows rbind-ed).
.numeric_summary_row <- function(values, dataset, variable) {
  values <- as.numeric(values)
  values <- values[is.finite(values)]
  has <- length(values) > 0L
  q <- function(p) if (has) unname(stats::quantile(values, p)) else NA_real_
  data.frame(
    variable = variable, dataset = dataset,
    n = length(values),
    mean = if (has) mean(values) else NA_real_,
    sd = if (length(values) > 1L) stats::sd(values) else NA_real_,
    min = if (has) min(values) else NA_real_,
    q25 = q(0.25), median = q(0.5), q75 = q(0.75),
    max = if (has) max(values) else NA_real_,
    stringsAsFactors = FALSE
  )
}

# Per-level counts and proportions for a categorical covariate.
.categorical_summary_rows <- function(values, dataset, variable) {
  values <- as.character(values)
  values[is.na(values)] <- "<missing>"
  counts <- table(values)
  data.frame(
    variable = variable, dataset = dataset,
    level = names(counts),
    n = as.integer(counts),
    proportion = as.numeric(counts) / sum(counts),
    row.names = NULL, stringsAsFactors = FALSE
  )
}

# One baseline value per subject, keeping the column's own type (a factor stays
# a factor rather than collapsing to integer codes, as unlist() would do).
.subject_baseline_values <- function(data, roles, covariate) {
  subjects <- .unique_in_order(data[[roles$id]])
  idx <- vapply(subjects, function(s) {
    rows <- which(!is.na(data[[roles$id]]) & data[[roles$id]] == s)
    present <- rows[!is.na(data[[covariate]][rows])]
    if (length(present)) present[1L] else if (length(rows)) rows[1L] else NA_integer_
  }, integer(1))
  data[[covariate]][idx]
}

# n (observations) / n_subjects / distribution of DV per endpoint.
.endpoint_dv_summary <- function(data, roles, dataset) {
  selected <- .observation_rows(data, roles, require_present = TRUE)
  endpoint <- .endpoint(data, roles)[selected]
  dv <- as.numeric(data[[roles$dv]][selected])
  id <- as.character(data[[roles$id]][selected])
  labels <- sort(unique(endpoint))
  rows <- lapply(labels, function(lab) {
    keep <- endpoint == lab
    row <- .numeric_summary_row(dv[keep], dataset, lab)
    row$n_subjects <- length(unique(id[keep]))
    row[, c("variable", "dataset", "n", "n_subjects", "mean", "sd",
            "min", "q25", "median", "q75", "max")]
  })
  do.call(rbind, rows)
}

#' Compare per-covariate and per-endpoint distributions of source and synthetic
#'
#' A numeric sanity check to run right after generating data. For each baseline
#' covariate and each endpoint (`dvid`), it summarizes the distribution in the
#' source and in the synthetic dataset side by side. The dependent variable and
#' continuous covariates get n, mean, standard deviation, minimum, quartiles, and
#' maximum; categorical covariates get per-level counts and proportions.
#'
#' This is the distributional companion to [compare_pmx()]. That function answers
#' whether the *structure* matches — schema, event grammar, row and event counts;
#' this one answers whether the *numbers* land in the same range. It is a
#' diagnostic, not a validation of statistical fidelity: AVATAR and the
#' differentially private engines deliberately do not reproduce source
#' distributions exactly, so expect the summaries to be close in magnitude and
#' shape, not identical.
#'
#' Every table is source-derived, so each is marked
#' `"restricted_not_releasable"`: it reads real covariate and endpoint values and
#' stays under the source data's access controls like any other
#' source-versus-synthetic diagnostic.
#'
#' @param source Source PMX data.
#' @param synthetic Generated synthetic PMX data, or `NULL` to summarize `source`
#'   on its own.
#' @param roles Explicit roles from [pmx_roles()].
#'
#' @return A `pmx_distribution_summary`: a list of `endpoints`,
#'   `covariates_numeric`, and `covariates_categorical` data frames. Each is
#'   `NULL` when the dataset declares no columns of that kind.
#' @seealso [compare_pmx()] for the structural comparison.
#' @export
#' @examples
#' data <- pmx_simulated_fixture(20)
#' roles <- pmx_roles(
#'   id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
#'   cmt = "CMT", dvid = "DVID", covariates = c("WT", "SEX")
#' )
#' synthetic <- suppressWarnings(synpmx_avatar(data, roles, seed = 1))
#' compare_pmx_distributions(data, synthetic, roles)
compare_pmx_distributions <- function(source, synthetic = NULL, roles) {
  .assert_roles(source, roles)
  datasets <- list(source = source)
  if (!is.null(synthetic)) {
    .assert_roles(synthetic, roles)
    datasets$synthetic <- synthetic
  }

  numeric_rows <- list()
  categorical_rows <- list()
  for (covariate in roles$covariates) {
    for (label in names(datasets)) {
      values <- .subject_baseline_values(datasets[[label]], roles, covariate)
      if (is.numeric(values)) {
        numeric_rows[[length(numeric_rows) + 1L]] <-
          .numeric_summary_row(values, label, covariate)
      } else {
        categorical_rows[[length(categorical_rows) + 1L]] <-
          .categorical_summary_rows(values, label, covariate)
      }
    }
  }

  endpoint_rows <- lapply(names(datasets), function(label) {
    .endpoint_dv_summary(datasets[[label]], roles, label)
  })

  bind <- function(parts) if (length(parts)) do.call(rbind, parts) else NULL
  structure(
    list(
      endpoints = .mark_release(bind(endpoint_rows),
                                "restricted_not_releasable"),
      covariates_numeric = .mark_release(bind(numeric_rows),
                                         "restricted_not_releasable"),
      covariates_categorical = .mark_release(bind(categorical_rows),
                                             "restricted_not_releasable")
    ),
    class = "pmx_distribution_summary"
  )
}

# Round numeric columns for display without touching the stored exact values.
.round_for_print <- function(df, digits = 4L) {
  numeric_cols <- vapply(df, is.numeric, logical(1))
  df[numeric_cols] <- lapply(df[numeric_cols], signif, digits = digits)
  df
}

#' @export
print.pmx_distribution_summary <- function(x, ...) {
  cat("Restricted PMX source-versus-synthetic distribution summary\n")
  section <- function(title, df) {
    if (is.null(df)) return(invisible())
    cat("\n", title, ":\n", sep = "")
    print(.round_for_print(df), row.names = FALSE)
  }
  section("Endpoints (dependent variable on observation rows)", x$endpoints)
  section("Continuous covariates (baseline, per subject)", x$covariates_numeric)
  section("Categorical covariates (baseline, per subject)",
          x$covariates_categorical)
  cat("\nSource-derived; not releasable unless separately public or privately",
      "budgeted.\n")
  invisible(x)
}
