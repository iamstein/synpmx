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

  # Order every table primarily by variable, then by dataset, so a variable's
  # source and synthetic rows sit together for easy comparison. Variables keep
  # their order of first appearance; datasets keep source before synthetic.
  bind <- function(parts) {
    if (!length(parts)) return(NULL)
    df <- do.call(rbind, parts)
    df <- df[order(match(df$variable, unique(df$variable)),
                   match(df$dataset, names(datasets))), , drop = FALSE]
    rownames(df) <- NULL
    df
  }
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

# Post-generation outlier / identifiability check -----------------------------
#
# compare_pmx_distributions() compares whole distributions; this checks
# individuals. A subject is easy to single out -- and so to re-identify -- when
# its event structure or measurements are unlike anyone else's: an unusually
# long follow-up, an odd number of doses, a rare dose level, or an extreme
# value. synpmx_avatar() copies each avatar's event skeleton from one anchor, so
# a structurally unique source subject reappears structurally unique. This
# screens for exactly those subjects, one axis at a time, with a robust
# median/MAD outlier score.

.modified_z <- function(x) {
  z <- rep(NA_real_, length(x))
  finite <- is.finite(x)
  if (sum(finite) < 2L) return(z)
  centre <- stats::median(x[finite])
  spread <- stats::median(abs(x[finite] - centre))
  if (!is.finite(spread) || spread == 0) {
    # No robust spread: nearly everyone shares one value, so any departure from
    # it is the outlier signal.
    z[finite] <- ifelse(x[finite] == centre, 0,
                        sign(x[finite] - centre) * Inf)
  } else {
    z[finite] <- 0.6745 * (x[finite] - centre) / spread
  }
  z
}

#' Flag structurally unusual -- and so easily identifiable -- subjects
#'
#' A post-generation screen for subjects that stand out from the cohort and are
#' therefore easy to single out and re-identify: the per-subject counterpart to
#' [compare_pmx_distributions()], which compares whole distributions. Each
#' subject is scored, one axis at a time, on a robust median/MAD statistic across
#' four structural features:
#'
#' - **follow-up time** -- the last observation time (catches the lone
#'   long-followed subject);
#' - **number of doses** -- an unusual dosing-history length;
#' - **dose magnitude** -- a rare dose level (needs an `amt` role); and
#' - **DV value** -- an extreme peak measurement.
#'
#' A subject is flagged when it is an outlier on any axis. This matters because
#' [synpmx_avatar()] copies each avatar's event skeleton from a single anchor, so
#' a structurally unique source subject yields a structurally unique -- and
#' identifiable -- avatar even though its measurements are blended. Run it on the
#' synthetic data before the data leaves the source's access controls and drop or
#' regenerate the flagged subjects; it can also be run on the source itself to
#' see which real subjects are hardest to hide. It is a heuristic screen, not a
#' privacy guarantee, and is marked `"restricted_not_releasable"`.
#'
#' @param data A PMX dataset -- typically the synthetic output, or the source.
#' @param roles Explicit roles from [pmx_roles()].
#' @param threshold Absolute modified-z cutoff above which a subject is an
#'   outlier on an axis. Default 3.5, the Iglewicz--Hoaglin value.
#'
#' @return A `pmx_identifiability` data frame, most-unusual first, one row per
#'   subject: `subject_id`, the four axis values (`follow_up_time`, `n_doses`,
#'   `max_dose`, `max_dv`), `outlier_axes` (a comma-separated list of the axes on
#'   which it is unusual, empty if none), and `flagged`.
#' @seealso [compare_pmx_distributions()], [compare_pmx()].
#' @export
#' @examples
#' data <- pmx_simulated_fixture(30)
#' roles <- pmx_roles(
#'   id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
#'   cmt = "CMT", dvid = "DVID", covariates = "WT"
#' )
#' synthetic <- suppressWarnings(synpmx_avatar(data, roles, seed = 1))
#' flag_identifiable_subjects(synthetic, roles)
flag_identifiable_subjects <- function(data, roles, threshold = 3.5) {
  .assert_roles(data, roles)
  if (!is.numeric(threshold) || length(threshold) != 1L ||
      is.na(threshold) || threshold <= 0) {
    stop("`threshold` must be a single positive number.", call. = FALSE)
  }

  subjects <- .unique_in_order(data[[roles$id]])
  sub <- factor(as.character(data[[roles$id]]),
                levels = as.character(subjects))
  observed <- .observation_rows(data, roles, require_present = TRUE)
  dosed <- .dose_rows(data, roles)
  time <- suppressWarnings(as.numeric(data[[roles$time]]))
  dv <- suppressWarnings(as.numeric(data[[roles$dv]]))
  amt <- if (!is.null(roles$amt)) {
    suppressWarnings(as.numeric(data[[roles$amt]]))
  } else NULL

  safe_max <- function(v) {
    v <- v[is.finite(v)]
    if (length(v)) max(v) else NA_real_
  }
  by_subject <- function(values, keep) {
    grouped <- split(values[keep], sub[keep])
    vapply(grouped, safe_max, numeric(1))[as.character(subjects)]
  }

  follow_up_time <- by_subject(time, observed)
  max_dv <- by_subject(dv, observed)
  n_doses <- as.numeric(tapply(as.integer(dosed), sub, sum)[
    as.character(subjects)
  ])
  n_doses[is.na(n_doses)] <- 0
  max_dose <- if (!is.null(amt)) by_subject(amt, dosed) else
    rep(NA_real_, length(subjects))

  axes <- list(time = follow_up_time, doses = n_doses,
               dose = max_dose, dv = max_dv)
  axis_label <- c(time = "follow-up time", doses = "number of doses",
                  dose = "dose magnitude", dv = "DV value")
  outlier <- lapply(axes, function(v) {
    if (sum(is.finite(v)) < 2L) return(rep(FALSE, length(v)))
    z <- .modified_z(v)
    flagged <- (is.finite(z) & abs(z) > threshold) | is.infinite(z)
    flagged[is.na(flagged)] <- FALSE
    flagged
  })

  outlier_matrix <- do.call(cbind, outlier)
  outlier_axes <- apply(outlier_matrix, 1L, function(hit) {
    paste(axis_label[names(axes)[hit]], collapse = ", ")
  })
  order_key <- rowSums(outlier_matrix)

  out <- data.frame(
    subject_id = as.character(subjects),
    follow_up_time = follow_up_time,
    n_doses = n_doses,
    max_dose = max_dose,
    max_dv = max_dv,
    outlier_axes = outlier_axes,
    flagged = order_key > 0L,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  out <- out[order(-order_key, out$subject_id), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "n_flagged") <- sum(out$flagged)
  attr(out, "threshold") <- threshold
  .mark_release(
    structure(out, class = c("pmx_identifiability", "data.frame")),
    "restricted_not_releasable"
  )
}

# Skeleton uniqueness ---------------------------------------------------------
#
# flag_identifiable_subjects() finds subjects that are *extreme*. This finds
# subjects that are *alone*, which is a different property: a patient can sit
# squarely in the middle of every distribution and still be the only one with
# their exact visit schedule. synpmx_avatar() copies the anchor's event skeleton
# verbatim (`.jitter_skeleton_time()` cannot change that -- its clamp holds every
# time inside its own Voronoi cell), so a subject alone in its equivalence class
# hands its schedule to every avatar anchored on it.
#
# Two classes are scored because `coarsen_time` collapses one and not the other:
#
#   - `signature` -- the full `.event_signature()`: dose amounts, dose gaps,
#     endpoint set. Under actual recorded times this is near-universally unique
#     and coarsening is what collapses it.
#   - `n_obs` -- the observation count alone. Coarsening does not touch it, so
#     what survives is dropout and missed visits. That residual is the screen's
#     job, not the grid's.

.class_sizes <- function(key) {
  counts <- table(key)
  as.integer(counts[match(key, names(counts))])
}

#' Score how many subjects share each subject's event skeleton
#'
#' A source-side screen for subjects that are **alone**, the complement to
#' [flag_identifiable_subjects()], which finds subjects that are **extreme**. A
#' patient can be perfectly ordinary on every distribution and still hold the
#' only copy of their visit schedule, and [synpmx_avatar()] copies the anchor's
#' event skeleton verbatim, so such a subject hands an identifying schedule to
#' every avatar anchored on it.
#'
#' Three equivalence classes are scored per subject:
#'
#' - **`obs_time`** -- subjects sharing the exact observation time vector. This
#'   is the fingerprint, because [synpmx_avatar()] copies the anchor's event
#'   skeleton verbatim. Under nominal visit times the class is large, since the
#'   schedule is protocol-driven; under actual recorded times it is
#'   near-universally of size one. `coarsen_time = TRUE` collapses the second
#'   case into the first, and `alone` reports this class.
#' - **`n_obs`** -- subjects sharing the observation count. Coarsening cannot
#'   change a count, so this is what survives it: dropout, early
#'   discontinuation, and missed visits. That residual is the outlier screen's
#'   job rather than the grid's.
#' - **`signature`** -- subjects sharing the full [pmx_roles()] event signature:
#'   dose structure, dose amounts, and endpoint set. Note this does *not* include
#'   observation times -- it is the key donor compatibility uses. Weight-based
#'   dosing or per-subject titration makes it unique regardless of schedule, and
#'   coarsening does not change that either.
#'
#' Run it on the **source**, before generating, to decide whether coarsening is
#' needed and what is left over once it is applied. It is a heuristic screen, not
#' a privacy guarantee, and is marked `"restricted_not_releasable"`.
#'
#' @param data A PMX dataset -- normally the source.
#' @param roles Explicit roles from [pmx_roles()].
#'
#' @return A `pmx_skeleton_uniqueness` data frame, most-exposed first, one row
#'   per subject: `subject_id`, `n_obs`, `n_doses`, `signature_class`,
#'   `obs_time_class`, `n_obs_class` (each counting the subjects sharing that
#'   key, including itself), and `alone` (`TRUE` when `obs_time_class == 1`).
#'   Attributes `n_alone`, `n_alone_signature`, `n_alone_n_obs`, and `min_class`
#'   summarize the cohort.
#' @seealso [flag_identifiable_subjects()], [synpmx_avatar()].
#' @export
#' @examples
#' data <- pmx_simulated_fixture(30)
#' roles <- pmx_roles(
#'   id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
#'   cmt = "CMT", dvid = "DVID", covariates = "WT"
#' )
#' skeleton_uniqueness(data, roles)
skeleton_uniqueness <- function(data, roles) {
  .assert_roles(data, roles)
  subjects <- .unique_in_order(data[[roles$id]])
  sub <- factor(as.character(data[[roles$id]]),
                levels = as.character(subjects))
  observed <- .observation_rows(data, roles, require_present = TRUE)
  dosed <- .dose_rows(data, roles)

  subject_rows <- lapply(as.character(subjects), function(id) which(sub == id))
  observed_index <- observed
  signature <- vapply(subject_rows, function(rows) {
    .event_signature(data[rows, , drop = FALSE], roles)
  }, character(1))
  # `.event_signature()` covers dose structure, dose amounts, and the endpoint
  # set -- it does *not* include observation times. Donor compatibility does not
  # need them, but they are exactly what `synpmx_avatar()` copies verbatim from
  # the anchor, so the observation time vector is scored on its own. Under actual
  # recorded times this is the class that is universally of size one, and the one
  # `coarsen_time = TRUE` exists to collapse.
  time <- suppressWarnings(as.numeric(data[[roles$time]]))
  obs_time <- vapply(subject_rows, function(rows) {
    values <- sort(time[rows[observed_index[rows]]])
    paste(format(values, digits = 12, trim = TRUE), collapse = ",")
  }, character(1))

  count_by <- function(keep) {
    as.integer(tapply(as.integer(keep), sub, sum)[as.character(subjects)])
  }
  n_obs <- count_by(observed)
  n_obs[is.na(n_obs)] <- 0L
  n_doses <- count_by(dosed)
  n_doses[is.na(n_doses)] <- 0L

  signature_class <- .class_sizes(signature)
  obs_time_class <- .class_sizes(obs_time)
  n_obs_class <- .class_sizes(as.character(n_obs))

  out <- data.frame(
    subject_id = as.character(subjects),
    n_obs = n_obs,
    n_doses = n_doses,
    signature_class = signature_class,
    obs_time_class = obs_time_class,
    n_obs_class = n_obs_class,
    alone = obs_time_class == 1L,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  out <- out[order(out$obs_time_class, out$signature_class, out$n_obs_class,
                   out$subject_id), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "n_alone") <- sum(out$obs_time_class == 1L)
  attr(out, "n_alone_signature") <- sum(out$signature_class == 1L)
  attr(out, "n_alone_n_obs") <- sum(out$n_obs_class == 1L)
  attr(out, "min_class") <- if (nrow(out)) min(out$obs_time_class) else NA_integer_
  .mark_release(
    structure(out, class = c("pmx_skeleton_uniqueness", "data.frame")),
    "restricted_not_releasable"
  )
}

#' @export
print.pmx_skeleton_uniqueness <- function(x, ...) {
  n <- nrow(x)
  alone <- attr(x, "n_alone") %||% sum(x$signature_class == 1L)
  alone_sig <- attr(x, "n_alone_signature") %||% sum(x$signature_class == 1L)
  alone_obs <- attr(x, "n_alone_n_obs") %||% sum(x$n_obs_class == 1L)
  cat(sprintf(
    "Restricted PMX skeleton-uniqueness screen: %d of %d subject%s alone\n",
    alone, n, if (n == 1L) "" else "s"
  ))
  cat(sprintf(
    "Alone = the only subject with this observation time vector (%.0f%% here).\n\n",
    if (n) 100 * alone / n else 0
  ))
  cat(sprintf("  obs times alone: %3d  <- what `coarsen_time = TRUE` collapses\n",
              alone))
  cat(sprintf("  n_obs alone:     %3d  <- the residual it leaves, for the screen\n",
              alone_obs))
  cat(sprintf("  signature alone: %3d  <- dose structure/amount; coarsening does not change it\n",
              alone_sig))
  cat("\n")
  cat(if (n > 12L) "Twelve most exposed:\n" else "By class size:\n")
  print(.round_for_print(utils::head(as.data.frame(x), 12L)), row.names = FALSE)
  if (n > 12L) {
    cat(sprintf("... %d more row(s) in the returned table.\n", n - 12L))
  }
  cat("\nSource-derived; not releasable unless separately public or privately",
      "budgeted.\n")
  invisible(x)
}

#' @export
print.pmx_identifiability <- function(x, ...) {
  n <- nrow(x)
  flagged <- attr(x, "n_flagged") %||% sum(x$flagged)
  cat(sprintf(
    "Restricted PMX outlier / identifiability check: %d of %d subject%s flagged\n",
    flagged, n, if (n == 1L) "" else "s"
  ))
  cat("Flag = a robust outlier in follow-up time, dose count, dose magnitude,",
      "or DV value.\n\n")
  cat(if (n > 12L) "Twelve most unusual:\n" else "By outlier count:\n")
  print(.round_for_print(utils::head(as.data.frame(x), 12L)), row.names = FALSE)
  if (n > 12L) {
    cat(sprintf("... %d more row(s) in the returned table.\n", n - 12L))
  }
  cat("\nSource-derived; not releasable unless separately public or privately",
      "budgeted.\n")
  invisible(x)
}

# Apply the truncate/drop policy to one dataset (no regeneration). Shared by the
# public remediation function and its replacement loop.
.apply_remediation_policy <- function(data, roles, time, other, threshold) {
  report <- flag_identifiable_subjects(data, roles, threshold = threshold)
  time_flag <- grepl("follow-up time", report$outlier_axes, fixed = TRUE)
  time_only <- report$flagged & report$outlier_axes == "follow-up time"
  other_flagged <- report$flagged & !time_only

  # Truncate toward the longest follow-up that is not itself a time outlier.
  ordinary <- report$follow_up_time[!time_flag & is.finite(report$follow_up_time)]
  horizon <- if (length(ordinary)) max(ordinary) else NA_real_

  drop_ids <- character()
  if (other == "drop") drop_ids <- c(drop_ids, report$subject_id[other_flagged])

  truncate_ids <- character()
  if (any(time_only)) {
    if (time == "drop") {
      drop_ids <- c(drop_ids, report$subject_id[time_only])
    } else if (time == "truncate") {
      # Only a follow-up *longer* than the ordinary maximum can be truncated.
      # A subject flagged for an unusually *short* follow-up has nothing to
      # trim, so it is dropped instead.
      long <- time_only & is.finite(report$follow_up_time) &
        is.finite(horizon) & report$follow_up_time > horizon
      truncate_ids <- report$subject_id[long]
      drop_ids <- c(drop_ids, report$subject_id[time_only & !long])
    }
  }
  drop_ids <- unique(drop_ids)

  id <- as.character(data[[roles$id]])
  keep <- !(id %in% drop_ids)
  if (length(truncate_ids)) {
    times <- suppressWarnings(as.numeric(data[[roles$time]]))
    keep <- keep & !(id %in% truncate_ids & is.finite(times) & times > horizon)
  }
  out <- data[keep, , drop = FALSE]
  rownames(out) <- NULL
  list(data = out, dropped = drop_ids, truncated = truncate_ids, horizon = horizon)
}

#' Remove or shorten the subjects `flag_identifiable_subjects()` flags
#'
#' Applies a remediation policy to the outliers found by
#' [flag_identifiable_subjects()]. By default a subject flagged **only** for an
#' unusually long follow-up is *truncated* back to the cohort's longest ordinary
#' follow-up (its late rows dropped), and a subject flagged for **any other**
#' reason -- an extreme DV, a rare dose level, or an unusual dose count -- is
#' *dropped* entirely.
#'
#' The split is deliberate. Truncation is offered only for an unusually *long*
#' follow-up, the one structural outlier a value-level edit can genuinely fix:
#' shortening a long timeline leaves a shorter but ordinary subject. Everything
#' else is dropped -- an unusually *short* follow-up has nothing to trim, an
#' extreme-DV subject is elevated across its whole trajectory so removing points
#' would only mangle it, and a rare dose cannot be trimmed without breaking the
#' regimen. A subject flagged for both a long follow-up and another reason is
#' also dropped, since truncation would not resolve the other reason.
#'
#' When `source` is supplied, each dropped subject is **replaced**: fresh avatars
#' are generated from `source`, screened by the same policy, and appended (with
#' new ids) until the cohort is back to its original size. So the output keeps
#' the same number of subjects, minus any it could not refill within `max_tries`.
#' Truncation keeps its subject, so it never triggers a replacement.
#'
#' Detection is per subject, so one long-followed patient is truncated once and
#' one extreme patient dropped-and-replaced once -- there is no row-level outlier
#' spray. With replacement, this is a self-contained alternative to preventing
#' structural outliers at generation time (skeleton sampling, `REV-026`).
#'
#' @param data A PMX dataset, typically the synthetic output.
#' @param roles Explicit roles from [pmx_roles()].
#' @param source Optional source PMX data. When given, dropped subjects are
#'   replaced by fresh avatars generated from it, so the cohort size is
#'   preserved. When `NULL` (default), dropped subjects are simply removed.
#' @param time Action for a subject whose *only* outlier axis is follow-up time:
#'   `"truncate"` (default) to shorten it to the longest ordinary follow-up,
#'   `"drop"` to remove it, or `"keep"` to leave it.
#' @param other Action for a subject flagged for any non-time reason: `"drop"`
#'   (default) or `"keep"`.
#' @param threshold Passed to [flag_identifiable_subjects()].
#' @param seed Reproducibility seed for the replacement generation. The caller's
#'   random-number state is restored by [synpmx_avatar()].
#' @param max_tries Maximum regeneration batches when refilling dropped subjects.
#'
#' @return `data` with the policy applied, carrying attributes `dropped`,
#'   `truncated` (affected subject ids), `replaced` (count refilled), and
#'   `horizon` (the follow-up truncation used, or `NA`).
#' @seealso [flag_identifiable_subjects()].
#' @export
#' @examples
#' data <- pmx_simulated_fixture(30)
#' roles <- pmx_roles(
#'   id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
#'   cmt = "CMT", dvid = "DVID", covariates = "WT"
#' )
#' synthetic <- suppressWarnings(synpmx_avatar(data, roles, seed = 1))
#' cleaned <- remediate_identifiable_subjects(synthetic, roles, source = data)
remediate_identifiable_subjects <- function(data, roles, source = NULL,
                                            time = c("truncate", "drop", "keep"),
                                            other = c("drop", "keep"),
                                            threshold = 3.5, seed = NULL,
                                            max_tries = 20L) {
  time <- match.arg(time)
  other <- match.arg(other)
  res <- .apply_remediation_policy(data, roles, time, other, threshold)
  out <- res$data
  replaced <- 0L

  if (length(res$dropped) && !is.null(source)) {
    .assert_roles(source, roles)
    need <- length(res$dropped)
    base_seed <- if (is.null(seed)) {
      sample.int(.Machine$integer.max, 1L)
    } else as.integer(seed)
    # New ids must avoid every original id, including the dropped ones, so a
    # replacement never silently reuses a removed subject's label.
    avoid <- data[[roles$id]]
    tries <- 0L
    while (replaced < need && tries < max_tries) {
      tries <- tries + 1L
      batch <- suppressWarnings(suppressMessages(synpmx_avatar(
        source, roles, n_subjects = need - replaced, seed = base_seed + tries
      )))
      if (!setequal(names(batch), names(out))) {
        stop("Replacement schema does not match `data`; regenerate from the ",
             "same source and roles.", call. = FALSE)
      }
      batch <- batch[, names(out), drop = FALSE]
      clean <- .apply_remediation_policy(batch, roles, time, other,
                                         threshold)$data
      clean_ids <- .unique_in_order(clean[[roles$id]])
      if (!length(clean_ids)) next
      take <- clean_ids[seq_len(min(length(clean_ids), need - replaced))]
      chunk <- clean[as.character(clean[[roles$id]]) %in% as.character(take), ,
                     drop = FALSE]
      fresh <- .new_ids(avoid, length(take))
      avoid <- c(avoid, fresh)
      id_map <- stats::setNames(fresh, as.character(take))
      chunk[[roles$id]] <- id_map[as.character(chunk[[roles$id]])]
      out <- rbind(out, chunk)
      replaced <- replaced + length(take)
    }
    rownames(out) <- NULL
    if (replaced < need) {
      warning(sprintf(
        "Refilled only %d of %d dropped subject(s) in %d tries.",
        replaced, need, max_tries
      ), call. = FALSE)
    }
  }

  message(sprintf(
    "remediate_identifiable_subjects(): dropped %d, truncated %d%s%s.",
    length(res$dropped), length(res$truncated),
    if (length(res$truncated)) sprintf(" (to follow-up <= %.4g)", res$horizon)
      else "",
    if (replaced) sprintf(", replaced %d", replaced) else ""
  ))
  attr(out, "dropped") <- res$dropped
  attr(out, "truncated") <- res$truncated
  attr(out, "replaced") <- replaced
  attr(out, "horizon") <- res$horizon
  out
}
