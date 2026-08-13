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

# Three significant digits for display, without touching the stored exact
# values. Each cell is formatted on its own rather than the column rounded as a
# vector: these tables put an endpoint in the hundreds next to one around 0.05,
# and a column laid out in one pass pads every entry to the widest one's decimal
# places (134.9000). Integer columns -- the n counts -- are left alone, since
# signif() would turn 1234 subjects into 1230.
.format_for_print <- function(df, digits = 3L) {
  double_cols <- vapply(df, is.double, logical(1))
  df[double_cols] <- lapply(df[double_cols], function(x) {
    vapply(x, function(v) format(signif(v, digits), trim = TRUE), character(1))
  })
  df
}

#' @export
print.pmx_distribution_summary <- function(x, ...) {
  cat("Restricted PMX source-versus-synthetic distribution summary\n")
  section <- function(title, df) {
    if (is.null(df)) return(invisible())
    cat("\n", title, ":\n", sep = "")
    print(.format_for_print(df), row.names = FALSE)
  }
  section("Endpoints (dependent variable on observation rows)", x$endpoints)
  section("Continuous covariates (baseline, per subject)", x$covariates_numeric)
  section("Categorical covariates (baseline, per subject)",
          x$covariates_categorical)
  cat("\nSource-derived; not releasable unless separately public or privately",
      "budgeted.\n")
  invisible(x)
}

# Each of the four study templates under `scripts_private/` carried its own
# twenty-line `kable_distributions()` helper to get these out as tables rather
# than as a preformatted block. That is display code, it was copied four times,
# and it is the package's job.
#' @exportS3Method knitr::knit_print
knit_print.pmx_distribution_summary <- function(x, ...) {
  section <- function(df, caption) {
    if (is.null(df)) return(NULL)
    # The alignment has to be taken before formatting: .format_for_print() hands
    # back character columns, which kable would otherwise left-align.
    align <- ifelse(vapply(df, is.numeric, logical(1)), "r", "l")
    paste(knitr::kable(.format_for_print(df), row.names = FALSE,
                       align = align, caption = caption),
          collapse = "\n")
  }
  out <- c(
    section(x$endpoints,
            "RESTRICTED -- endpoints (dependent variable on observation rows)"),
    section(x$covariates_numeric,
            "RESTRICTED -- continuous covariates (baseline, per patient)"),
    section(x$covariates_categorical,
            "RESTRICTED -- categorical covariates (baseline, per patient)")
  )
  knitr::asis_output(paste(Filter(Negate(is.null), out), collapse = "\n\n"))
}

# The source-side rare-level census -------------------------------------------
#
# Numeric covariates are blended into a value nobody had. Categorical ones are
# NOT: they are `sample()`d from the donors' values, and `strata` are copied from
# the anchor by design, so a synthetic patient's category is always some real
# patient's actual category, copied. The disclosure that follows is not "this
# level is unique in the output" -- it is that a level too few REAL patients held
# appears in a table that may travel.
#
# There is a protection already, and it is geometric rather than designed: the
# sole holder of a level sits alone on its own one-hot axis, is nobody's nearest
# neighbour, and is rarely chosen as a donor. Measured, a level held by one
# patient tends not to reach the output and a level held by two does. Nothing
# enforces that and nothing reported it, which is what this closes.

#' Which rare source levels reached the synthetic output
#'
#' Censuses every categorical axis -- `strata` and each non-numeric covariate --
#' on both sides, and marks the levels **too few source patients held** to be
#' safely copied out. That floor is `min_pattern_share`, the same rule the
#' generator applies to visit sets and the scorecard applies to copied vectors,
#' so no new threshold is introduced.
#'
#' Read the `exposed` rows, and among them the ones that `reached` the output. A
#' level held by two real patients, appearing in a released table, says that
#' someone with that attribute was in this study; for a named trial with public
#' inclusion criteria that can be close to identifying on its own, and no cohort
#' size helps. The remedies are upstream of generation: drop the covariate from
#' `covariates`, or collapse its rare levels before generating.
#'
#' This reads real patient data on both sides and is marked
#' `"restricted_not_releasable"`.
#'
#' **What it cannot see.** Rarity *in the world*. If every living carrier of a
#' mutation is in this study, the source count is the whole population and looks
#' unremarkable. And it censuses each column on its own: with `d` covariates
#' there are `2^d` combinations that could single a patient out, and enumerating
#' them is not something this does.
#'
#' @param source Source PMX data.
#' @param synthetic Generated synthetic PMX data. When it carries a
#'   `"pmx_settings"` attribute, its `min_pattern_share` is used as the floor.
#' @param roles Explicit roles from [pmx_roles()].
#' @param floor Levels held by fewer than this many source patients are
#'   `exposed`. Left `NULL` it is taken from the run, and `2` otherwise -- the
#'   lowest value that means "more than one real patient".
#'
#' @return A `pmx_rare_levels` data frame, one row per categorical column and
#'   level, with `source_patients`, `synthetic_patients`, `exposed` and
#'   `reached`. Zero rows when the roles declare no categorical axis.
#' @seealso [synpmx_scorecard()], which reports this as row B5b,
#'   [compare_pmx_distributions()],
#'   `vignette("scorecard-synthetic-data-checks")`.
#' @export
#' @examples
#' data <- pmx_simulated_fixture(20)
#' roles <- pmx_roles(
#'   id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
#'   cmt = "CMT", dvid = "DVID", covariates = c("WT", "SEX")
#' )
#' synthetic <- suppressWarnings(synpmx_avatar(data, roles, seed = 1))
#' compare_pmx_rare_levels(data, synthetic, roles)
compare_pmx_rare_levels <- function(source, synthetic, roles, floor = NULL) {
  .assert_roles(source, roles)
  .assert_roles(synthetic, roles)
  if (is.null(floor)) {
    settings <- attr(synthetic, "pmx_settings")
    floor <- as.integer(settings$min_pattern_share %||% 2L)
  }
  if (!is.numeric(floor) || length(floor) != 1L || is.na(floor) || floor < 1) {
    stop("`floor` must be one integer of 1 or more.", call. = FALSE)
  }
  floor <- as.integer(floor)

  # Taken from the source, so a column the generator dropped is still censused.
  columns <- .scorecard_categorical(source, roles)
  rows <- lapply(columns, function(column) {
    held <- function(data) {
      if (is.null(data[[column]])) return(integer())
      table(.scorecard_holders(data, roles, column))
    }
    source_held <- held(source)
    synthetic_held <- held(synthetic)
    levels <- union(names(source_held), names(synthetic_held))
    # By position, since `[[` on a name never matches the empty string a blank
    # covariate cell becomes -- a level a real dataset does carry.
    count <- function(counts, level) {
      at <- match(level, names(counts))
      if (is.na(at)) 0L else as.integer(counts[[at]])
    }
    data.frame(
      column = column,
      level = levels,
      source_patients = vapply(levels, count, integer(1), counts = source_held),
      synthetic_patients = vapply(levels, count, integer(1),
                                  counts = synthetic_held),
      stringsAsFactors = FALSE
    )
  })

  out <- if (length(rows)) do.call(rbind, rows) else data.frame(
    column = character(), level = character(), source_patients = integer(),
    synthetic_patients = integer(), stringsAsFactors = FALSE
  )
  # A level with no source holder cannot expose a source patient, whatever it is
  # doing in the output -- and it should not be there at all, since categories
  # are copied rather than invented.
  out$exposed <- out$source_patients > 0L & out$source_patients < floor
  out$reached <- out$synthetic_patients > 0L
  out <- out[order(out$source_patients, out$column, out$level), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "floor") <- floor
  out <- structure(out, class = c("pmx_rare_levels", "data.frame"))
  .mark_release(out, "restricted_not_releasable")
}

.rare_levels_headline <- function(x) {
  exposed <- sum(x$exposed)
  leaked <- sum(x$exposed & x$reached)
  sprintf(paste("%d level(s) held by fewer than %d source patients;",
                "%d of them reached the output."),
          exposed, attr(x, "floor"), leaked)
}

#' @export
print.pmx_rare_levels <- function(x, ...) {
  plain <- as.data.frame(x)
  cat("Restricted PMX rare-level census (source against synthetic)\n\n")
  if (!nrow(plain)) {
    cat("No categorical axis is declared, so there is nothing to census.\n")
    return(invisible(x))
  }
  cat(.rare_levels_headline(x), "\n\n", sep = "")
  # The exposed rows are the point; the rest is the census they sit in, and on a
  # study with many levels printing all of it buries them. One level per line
  # rather than a data frame, because the level names are study labels -- "NON-
  # HISPANIC OR LATINO" -- and a wide frame wraps them into an unreadable block.
  shown <- plain[plain$exposed, , drop = FALSE]
  if (nrow(shown)) {
    for (i in seq_len(nrow(shown))) {
      cat(sprintf("  %s = %s\n    %d source patient(s), %d avatar(s)%s\n",
                  shown$column[[i]], shown$level[[i]],
                  shown$source_patients[[i]], shown$synthetic_patients[[i]],
                  if (shown$reached[[i]]) " -- REACHED THE OUTPUT" else ""))
    }
    cat("\nA level that reached the output is one real patient's attribute,\n",
        "copied. Drop the covariate or collapse its rare levels before\n",
        "generating.\n", sep = "")
  } else {
    cat("Every level in the output is one that at least ", attr(x, "floor"),
        " source patients held.\n", sep = "")
  }
  cat("\nSource-derived; not releasable.\n")
  invisible(x)
}

#' @exportS3Method knitr::knit_print
knit_print.pmx_rare_levels <- function(x, ...) {
  plain <- as.data.frame(x)
  if (!nrow(plain)) {
    return(knitr::asis_output("No categorical axis is declared, so there is ",
                              "nothing to census."))
  }
  knitr::knit_print(knitr::kable(
    plain, row.names = FALSE,
    caption = paste("RESTRICTED --", .rare_levels_headline(x))
  ))
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
  if (is.finite(spread) && spread > 0) {
    z[finite] <- 0.6745 * (x[finite] - centre) / spread
    return(z)
  }
  # The median absolute deviation is zero, which on trial data is the ordinary
  # case rather than a degenerate one: more than half a cohort completing the
  # protocol share one exact follow-up time, so the MAD collapses. Treating
  # every departure from the median as infinitely extreme -- which is what this
  # used to do -- then flags every patient who stopped early, however ordinary
  # their follow-up. Measured on a 21-patient study, 10 of 21 avatars were
  # flagged on follow-up time alone, including two agreeing to three decimal
  # places.
  #
  # The standard fallback (Iglewicz and Hoaglin) is the MEAN absolute deviation
  # from the median, which is zero only when every value is identical -- and
  # then there is nothing that departs from it to flag.
  average <- mean(abs(x[finite] - centre))
  if (!is.finite(average) || average == 0) {
    z[finite] <- 0
    return(z)
  }
  z[finite] <- (x[finite] - centre) / (1.253314 * average)
  z
}

# A stratum smaller than this is scored against the whole cohort instead: a
# median and a scale taken from four patients say more about the four than about
# the patient being screened.
.screen_stratum_floor <- 5L

# Which patients each subject should be compared against.
#
# "Does this patient stand out?" needs a comparison group, and the cohort is the
# wrong one as soon as a study assigns anything. On a six-arm dose-ranging study
# the top arm sits about 6.4 modified-z units from the cohort median dose purely
# by protocol, so a cohort-wide screen flags all thirty of its patients for
# receiving the dose they were assigned -- measured on `xgxr::case1_pkpd`, where
# it flagged 59 of 180 avatars, 31 of them on dose magnitude alone. Within the
# arm every dose is identical, `.modified_z()` returns 0, and nothing is
# flagged, which is the right answer; a patient who got 600 mg in a 300 mg arm
# still scores far past the threshold.
#
# Strata are the declared comparison group, so use them where they exist and are
# big enough to estimate a scale from. This is the same idea as
# `min_pattern_share` operating inside schedule groups rather than across the
# cohort.
.screen_strata <- function(data, roles, subjects) {
  if (!length(roles$strata)) return(NULL)
  key <- factor(as.character(data[[roles$id]]),
                levels = as.character(subjects))
  values <- lapply(roles$strata, function(column) {
    vapply(split(as.character(data[[column]]), key), function(v) {
      v <- v[!is.na(v)]
      if (length(v)) v[[1L]] else NA_character_
    }, character(1))[as.character(subjects)]
  })
  stratum <- do.call(paste, c(values, list(sep = "\r")))
  sizes <- table(stratum)
  # Anyone in a stratum too small to score within falls back to the cohort, and
  # they are pooled rather than compared with each other.
  stratum[stratum %in% names(sizes)[sizes < .screen_stratum_floor]] <- NA
  if (all(is.na(stratum))) return(NULL)
  stratum
}

# `.modified_z()` computed inside each stratum, with NA stratum scored against
# every subject that has one.
.modified_z_by <- function(v, stratum) {
  if (is.null(stratum)) return(.modified_z(v))
  z <- rep(NA_real_, length(v))
  pooled <- is.na(stratum)
  if (any(pooled)) z[pooled] <- .modified_z(v)[pooled]
  for (level in unique(stratum[!pooled])) {
    members <- which(!pooled & stratum == level)
    z[members] <- .modified_z(v[members])
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
#' Scores are computed **within each declared stratum** ([pmx_roles()]
#' `strata`), because "does this patient stand out?" needs a comparison group
#' and the whole cohort is the wrong one as soon as a study assigns anything. On
#' a six-arm dose-ranging study the top arm sits far from the cohort median dose
#' purely by protocol: scored cohort-wide, `xgxr::case1_pkpd` flags 59 of 180
#' avatars, 31 of them for receiving the dose their arm was assigned. Scored
#' within arm it flags 1, and a patient given twice their arm's dose is still
#' flagged. Strata holding fewer than five subjects are scored against the whole
#' cohort instead, since a scale estimated from four patients describes the four
#' rather than the one being screened. With no `strata` declared, every subject
#' is scored against the cohort, as before.
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
  stratum <- .screen_strata(data, roles, subjects)
  outlier <- lapply(axes, function(v) {
    if (sum(is.finite(v)) < 2L) return(rep(FALSE, length(v)))
    z <- .modified_z_by(v, stratum)
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
#     what survives is missing visits, whatever caused them. That residual is the screen's
#     job, not the grid's.

.class_sizes <- function(key) {
  counts <- table(key)
  as.integer(counts[match(key, names(counts))])
}

#' Score how many patients share each patient's event skeleton
#'
#' A source-side screen for patients that are **alone**, the complement to
#' [flag_identifiable_subjects()], which finds patients that are **extreme**. A
#' patient can be perfectly ordinary on every distribution and still hold the
#' only copy of their visit schedule, and [synpmx_avatar()] copies the anchor's
#' event skeleton verbatim, so such a patient hands an identifying schedule to
#' every avatar anchored on it.
#'
#' Three questions are asked of every patient, and the answer to each is a
#' count of how many patients share that property, the patient included. A
#' count of 1 means "nobody else":
#'
#' - **`n_share_schedule`** -- who else was observed at exactly this list of
#'   times? This is the fingerprint, because [synpmx_avatar()] copies the
#'   anchor's event skeleton verbatim. Under nominal visit times the count is
#'   large, since the schedule is protocol-driven; under actual recorded times
#'   it is near-universally 1. Coarsening collapses the second case into the
#'   first.
#' - **`n_share_obs_count`** -- who else has this many observations? Coarsening
#'   cannot change a count, so this is what survives it: missed visits, early
#'   discontinuation, and follow-up that has not reached the later visits.
#' - **`n_share_dosing`** -- who else has this dose structure and these dose
#'   amounts? This is the full [pmx_roles()] event signature and it does *not*
#'   include observation times; it is the key donor compatibility uses.
#'   Weight-based dosing or per-patient titration makes it unique regardless of
#'   schedule, and coarsening does not change that either.
#'
#' `n_share_rarest_time` splits the schedule count by cause, which matters
#' because the two causes have opposite remedies. A patient whose schedule is
#' unique **and** whose rarest single time was shared with nobody
#' (`n_share_rarest_time == 1`) was sampled at a one-off moment: a time grid is
#' meant to absorb that, and declaring `nominal_time` is the fix. A patient
#' whose schedule is unique while every individual time is shared
#' (`n_share_rarest_time >= 2`) has visits missing rather than moved, and no
#' grid at any resolution touches it. `why_unique` states which.
#'
#' # Before or after coarsening
#'
#' By default this scores the times **exactly as they appear in `data`**.
#' [synpmx_avatar()] snaps the source onto a shared visit grid first
#' (`coarsen_time = TRUE`, its default) and the numbers it records in
#' `pmx_settings` are therefore post-coarsening. Pass `coarsen_time = TRUE`
#' here to score the same grid the generator would build, and run it both ways
#' to see how much of the exposure coarsening actually removed. The printed
#' header always says which of the two you are looking at.
#'
#' Run it on the **source**, before generating. It is a heuristic screen, not a
#' privacy guarantee, and is marked `"restricted_not_releasable"`.
#'
#' @param data A PMX dataset -- normally the source.
#' @param roles Explicit roles from [pmx_roles()].
#' @param coarsen_time Score the coarsened visit grid [synpmx_avatar()] would
#'   build (`TRUE`) or the recorded times as given (`FALSE`, the default).
#'
#' @return A `pmx_skeleton_uniqueness` data frame, most-exposed first, one row
#'   per patient: `subject_id`, `n_obs`, `n_doses`, `n_share_dosing`,
#'   `n_share_schedule`, `n_share_rarest_time`, `n_share_obs_count` (each
#'   counting the patients sharing that property, including this one),
#'   `n_visits` and `nearest_set_diff` (how many visit slots -- one endpoint at
#'   one time -- separate this patient from the closest other one),
#'   `unique_schedule` (`TRUE` when `n_share_schedule == 1`), and `why_unique`.
#'   Attributes `n_unique_schedule`, `n_unique_dose_signature`,
#'   `n_unique_obs_count`, `n_unshared_time`, `min_class`, and `coarsened`
#'   summarize the cohort; `summary_table`, `sharing_table` and `by_endpoint`
#'   hold the tables `print()` shows.
#'
#' # Reading the count
#'
#' `n_share_schedule == 1` is exact-set equality, and on a real study that is a
#' harsh test: with forty visit slots and twenty patients, two patients who
#' differ by one missed sample score as "unique" exactly like two with nothing
#' in common. `nearest_set_diff` is what separates those cases, and the printed
#' output states it alongside the count. A cohort can read 15 of 21 unique while
#' every one of those 15 is a single missing sample away from somebody else.
#'
#' The `by_endpoint` table says *which* endpoint is responsible, since a
#' schedule is only as shared as its least shared part: a study measuring a
#' biomarker at every visit and PK at some of them is unique on the pooled
#' schedule the moment one PK sample is missing.
#'
#' None of this is something generation can lower — it is a property of the
#' source. What generation controls is whether an avatar ends up *carrying* one
#' of these schedules, which [pmx_masking_report()] reports as "avatars keeping
#' their anchor's own visit set".
#' @seealso [plot_pmx_schedule()] for the same information as a picture,
#'   [flag_identifiable_subjects()], [synpmx_avatar()].
#' @export
#' @examples
#' data <- pmx_simulated_fixture(30)
#' roles <- pmx_roles(
#'   id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
#'   cmt = "CMT", dvid = "DVID", covariates = "WT"
#' )
#' skeleton_uniqueness(data, roles)
#' skeleton_uniqueness(data, roles, coarsen_time = TRUE)
skeleton_uniqueness <- function(data, roles, coarsen_time = FALSE) {
  .assert_roles(data, roles)
  if (isTRUE(coarsen_time)) data <- .coarsen_source_time(data, roles)$source
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
  obs_times <- lapply(subject_rows, function(rows) {
    sort(time[rows[observed_index[rows]]])
  })
  obs_time <- vapply(obs_times, function(values) {
    paste(.time_key(values), collapse = ",")
  }, character(1))
  # Being alone on the whole vector has two very different causes, and they need
  # different remedies. Either the subject was observed at a moment nobody else
  # was -- a one-off visit no grid can hide, because the cell has one member --
  # or every individual time is shared and only the *pattern* of which visits
  # were attended is unique. The second is missing visits -- discontinuation, a
  # missed visit, or follow-up not yet that long -- and
  # `coarsen_time` cannot touch it however fine or coarse the grid.
  # `n_share_rarest_time` separates them: it is the smallest number of subjects
  # sharing any one of this subject's observation times.
  #
  # The keys are built per value by `.time_key()`, and that is the whole
  # correctness argument. `format(x, digits = 12)` picks ONE layout for the
  # whole vector it is handed, so hour 12 keys as "12" for a subject sampled
  # only on the hour and as "12.00000000000" for a subject who also has a
  # 1.2142857 sample. One real visit then splits across several keys, each with
  # a small count, and patients were reported as holding a one-off observation
  # time that eighteen others in fact shared -- which fired the
  # "unique observation times" alert, and the `nominal_time` advice with it, on
  # cohorts that were already on a shared grid.
  time_keys <- lapply(obs_times, function(values) .time_key(unique(values)))
  shared <- table(unlist(time_keys))
  min_time_share <- vapply(time_keys, function(keys) {
    if (!length(keys)) return(NA_integer_)
    min(as.integer(shared[keys]))
  }, integer(1))

  # How near the misses are. `n_share_schedule` is exact-set equality, which is
  # brutal on a real study: with forty visit slots and twenty patients, two
  # patients who differ by ONE missed sample are as "unique" as two with nothing
  # in common. A cohort that looks obviously similar on the schedule map can
  # still come out 15 of 21 unique, and the count alone cannot tell the reader
  # which of those two situations they are in. This can: it is the smallest
  # number of visit slots that separate this patient from any other.
  # Keyed on endpoint AND time, not on time alone. With two endpoints drawn at
  # overlapping visits, two patients can hold an identical set of distinct
  # *times* while differing in which endpoint was measured at them -- so a
  # pooled-time difference reads 0 for patients who are genuinely not alike.
  # A "visit slot" here is one endpoint measured at one time.
  cell_endpoint <- .endpoint(data, roles)
  visit_sets <- lapply(subject_rows, function(rows) {
    hit <- rows[observed_index[rows]]
    unique(paste0(cell_endpoint[hit], "@", .time_key(time[hit])))
  })
  nearest_diff <- vapply(seq_along(visit_sets), function(i) {
    others <- setdiff(seq_along(visit_sets), i)
    if (!length(others)) return(NA_integer_)
    min(vapply(others, function(j) {
      length(union(visit_sets[[i]], visit_sets[[j]])) -
        length(intersect(visit_sets[[i]], visit_sets[[j]]))
    }, integer(1)))
  }, integer(1))
  n_visits <- lengths(visit_sets)

  # Per endpoint, because a schedule is only as shared as its least shared part.
  # A study measuring PK at some visits and a biomarker at all of them is unique
  # on the pooled schedule the moment one PK sample is missing, and nothing in
  # the pooled counts says which endpoint did it.
  endpoint <- .endpoint(data, roles)
  by_endpoint <- lapply(sort(unique(endpoint[observed])), function(name) {
    keys <- vapply(subject_rows, function(rows) {
      hit <- rows[observed[rows] & endpoint[rows] == name]
      if (!length(hit)) return(NA_character_)
      paste(.time_key(sort(unique(time[hit]))), collapse = ",")
    }, character(1))
    present <- !is.na(keys)
    data.frame(
      endpoint = name,
      patients = sum(present),
      `distinct visit sets` = length(unique(keys[present])),
      `patients alone on theirs` = sum(present &
        keys %in% names(which(table(keys[present]) == 1L))),
      check.names = FALSE, stringsAsFactors = FALSE
    )
  })
  by_endpoint <- do.call(rbind, by_endpoint)

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

  # Every column is "how many patients share this with me, me included", so 1
  # always means "nobody else" and the reader needs one rule rather than four.
  # The old names (`obs_time_class`, `signature_class`, `min_time_share`) named
  # the equivalence class rather than the question, and nobody could read the
  # table without the help page open next to it.
  unique_schedule <- obs_time_class == 1L
  why_unique <- ifelse(
    !unique_schedule, "",
    ifelse(!is.na(min_time_share) & min_time_share == 1L,
           "one-off observation time", "set of visits attended")
  )
  out <- data.frame(
    subject_id = as.character(subjects),
    n_obs = n_obs,
    n_doses = n_doses,
    n_share_schedule = obs_time_class,
    n_share_rarest_time = min_time_share,
    n_share_obs_count = n_obs_class,
    n_share_dosing = signature_class,
    n_visits = n_visits,
    nearest_set_diff = nearest_diff,
    unique_schedule = unique_schedule,
    why_unique = why_unique,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  out <- out[order(out$n_share_schedule, out$n_share_dosing,
                   out$n_share_obs_count, out$subject_id), , drop = FALSE]
  rownames(out) <- NULL
  n_unique_schedule <- sum(out$n_share_schedule == 1L)
  n_unshared_time <- sum(out$n_share_rarest_time == 1L, na.rm = TRUE)
  attr(out, "n_unique_schedule") <- n_unique_schedule
  attr(out, "n_unique_dose_signature") <- sum(out$n_share_dosing == 1L)
  attr(out, "n_unique_obs_count") <- sum(out$n_share_obs_count == 1L)
  attr(out, "n_unshared_time") <- n_unshared_time
  attr(out, "min_class") <- if (nrow(out)) min(out$n_share_schedule) else
    NA_integer_
  attr(out, "coarsened") <- isTRUE(coarsen_time)
  attr(out, "by_endpoint") <- by_endpoint
  attr(out, "summary_table") <- .skeleton_summary_table(out)
  attr(out, "sharing_table") <- .skeleton_sharing_table(out)
  .mark_release(
    structure(out, class = c("pmx_skeleton_uniqueness", "data.frame")),
    "restricted_not_releasable"
  )
}

# The cohort answer, in the order a reader needs it: the headline count, then
# its two causes -- which have opposite remedies -- then the two properties
# that are unique for reasons a grid was never going to touch. Built as a data
# frame rather than printed with `cat()` so `knit_print()` can hand the same
# rows to `knitr::kable()` and a report gets a real table.
.skeleton_summary_table <- function(x) {
  n <- nrow(x)
  unique_schedule <- sum(x$n_share_schedule == 1L)
  unshared <- sum(x$n_share_rarest_time == 1L, na.rm = TRUE)
  pattern_only <- max(unique_schedule - unshared, 0L)
  rows <- list(
    c("Observation schedule nobody else has", unique_schedule,
      "the headline: an avatar anchored here wears one real patient's schedule"),
    c("... a one-off observation time", unshared,
      "sampled when nobody else was. A time grid can absorb this: declare `nominal_time`"),
    c("... the set of visits attended", pattern_only,
      "every time is shared. A missed visit, a discontinuation, or follow-up that has not reached the later visits. No grid touches this"),
    c("Observation count nobody else has", sum(x$n_share_obs_count == 1L),
      "survives any grid; the residual `flag_identifiable_subjects()` looks at"),
    c("Dosing nobody else has", sum(x$n_share_dosing == 1L),
      "dose amounts and gaps. Weight-based dosing makes this near-universal")
  )
  out <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(out) <- c("Patients whose ...", "n", "Meaning")
  out$n <- as.integer(out$n)
  out <- cbind(out[1L], `% of cohort` = if (n) round(100 * out$n / n) else 0,
               out[-1L])
  out[, c("Patients whose ...", "n", "% of cohort", "Meaning")]
}

# How crowded is each schedule? One row per group size, so a cohort where every
# patient sits in a group of 8 reads at a glance and does not need 60 patient
# rows to say so. This is the number that behaves like a k in k-anonymity on
# the schedule axis.
.skeleton_sharing_table <- function(x) {
  counts <- table(x$n_share_schedule)
  data.frame(
    `Patients sharing that schedule` = as.integer(names(counts)),
    `Patients` = as.integer(counts),
    `% of cohort` = if (nrow(x)) round(100 * as.integer(counts) / nrow(x)) else 0,
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

# Shared by print() and knit_print() so the two cannot drift. The verdict is
# stated in words first because the number on its own does not say whether the
# reader has a problem, and "44% unique" reads as alarming on a cohort where
# every one of those is an ordinary missed or not-yet-reached visit.
#
# Everything is recomputed from the columns rather than read off the cohort
# attributes. `[.data.frame` keeps the class and the attributes, so a reader who
# prints `screen[screen$unique_schedule, ]` -- the obvious thing to do -- would
# otherwise get the whole cohort's counts over a two-row table, and percentages
# above 100.
.skeleton_headline <- function(x) {
  n <- nrow(x)
  unique_schedule <- sum(x$n_share_schedule == 1L)
  unshared <- sum(x$n_share_rarest_time == 1L, na.rm = TRUE)
  pattern_only <- max(unique_schedule - unshared, 0L)
  scored <- if (isTRUE(attr(x, "coarsened"))) {
    paste("Scored AFTER coarsening, on the shared visit grid",
          "`synpmx_avatar()` builds. These are the numbers a run reports.")
  } else {
    paste("Scored on the recorded times AS GIVEN, before any coarsening.",
          "`synpmx_avatar()` coarsens first by default, so run this again",
          "with `coarsen_time = TRUE` to see what the grid removes.")
  }
  verdict <- if (!n) {
    "No patients."
  } else if (unique_schedule == 0L) {
    "Every patient shares their observation schedule with somebody. Nothing to do."
  } else if (unshared == 0L) {
    paste0(unique_schedule, " of ", n, " patients (",
           round(100 * unique_schedule / n),
           "%) have an observation schedule nobody else has, all of them ",
           "because of which visits they attended rather than when. No time ",
           "grid can change that; `min_pattern_share` is what stops those ",
           "sets being reused.")
  } else {
    paste0(unique_schedule, " of ", n, " patients (",
           round(100 * unique_schedule / n),
           "%) have an observation schedule nobody else has: ", unshared,
           " from a one-off observation time, ", pattern_only,
           " from which visits they attended. Declaring `nominal_time` ",
           "addresses the first group; nothing addresses the second.")
  }
  # Near-misses, stated in the same breath as the count. "15 of 21 are unique"
  # and "the typical one differs from its nearest neighbour by a single visit"
  # describe the same cohort and lead to opposite conclusions.
  alone <- x[x$n_share_schedule == 1L, , drop = FALSE]
  closeness <- if (!nrow(alone) || all(is.na(alone$nearest_set_diff))) {
    NULL
  } else {
    sprintf(paste("Those %d are not necessarily far apart. The typical one",
                  "differs from its nearest neighbour by %d of about %d visit",
                  "slots (range %d to %d), where a slot is one endpoint",
                  "measured at one time. A difference of one or two is a missed",
                  "sample, not a different schedule -- which is why the count",
                  "alone is a poor guide."),
            nrow(alone),
            stats::median(alone$nearest_set_diff, na.rm = TRUE),
            stats::median(alone$n_visits, na.rm = TRUE),
            min(alone$nearest_set_diff, na.rm = TRUE),
            max(alone$nearest_set_diff, na.rm = TRUE))
  }
  # The distinction that trips everybody: this counts SOURCE patients, and a
  # source is what it is. What the masking controls is whether any AVATAR ends
  # up carrying one of these, which is the run report's "avatars keeping their
  # anchor's own visit set". A cohort can be 15 of 21 unique here and still put
  # none of those schedules into the output.
  scope <- paste(
    "This is a property of the SOURCE, and nothing in generation can lower it.",
    "What generation controls is whether an avatar ends up wearing one of these",
    "schedules -- that is `pmx_masking_report()`'s \"avatars keeping their",
    "anchor's own visit set\", which should be near 0% however high the count",
    "above is."
  )
  list(scored = scored, verdict = verdict, closeness = closeness, scope = scope,
       summary = .skeleton_summary_table(x),
       sharing = .skeleton_sharing_table(x),
       by_endpoint = attr(x, "by_endpoint"))
}

.skeleton_footer <- paste(
  "One row per patient is in the returned data frame;",
  "`plot_pmx_schedule()` draws the same cohort. Source-derived;",
  "not releasable unless separately public or privately budgeted."
)

#' @export
print.pmx_skeleton_uniqueness <- function(x, ...) {
  parts <- .skeleton_headline(x)
  cat("Restricted PMX schedule-uniqueness screen\n")
  cat(.wrap_plain(parts$scored), "\n\n", sep = "")
  cat(.wrap_plain(parts$verdict), "\n\n", sep = "")
  if (!is.null(parts$closeness)) {
    cat(.wrap_plain(parts$closeness), "\n\n", sep = "")
  }
  cat(.wrap_plain(parts$scope), "\n\n", sep = "")
  print(parts$summary[, c("Patients whose ...", "n", "% of cohort")],
        row.names = FALSE)
  cat("\nHow crowded is each schedule (1 = nobody else has it):\n")
  print(parts$sharing, row.names = FALSE)
  if (!is.null(parts$by_endpoint) && nrow(parts$by_endpoint) > 1L) {
    cat("\nWhich endpoint is doing it:\n")
    print(parts$by_endpoint, row.names = FALSE)
  }
  cat("\n", .wrap_plain(.skeleton_footer), "\n", sep = "")
  invisible(x)
}

# A knitted report gets real tables instead of a preformatted block, which is
# the whole reason the summaries above are data frames. Registered lazily, so
# knitr stays a Suggests.
#' @exportS3Method knitr::knit_print
knit_print.pmx_skeleton_uniqueness <- function(x, ...) {
  parts <- .skeleton_headline(x)
  out <- c(
    paste0("**Schedule-uniqueness screen.** ", parts$scored),
    parts$verdict,
    parts$closeness,
    parts$scope,
    paste(knitr::kable(parts$summary, align = c("l", "r", "r", "l")),
          collapse = "\n"),
    # A plain heading rather than a `kable()` caption: pandoc binds a caption
    # to whichever table it decides is adjacent, and in a block of two it chose
    # the wrong one, labelling the summary above with the sharing table's title.
    "**How crowded is each schedule** (1 = nobody else has it):",
    paste(knitr::kable(parts$sharing, align = c("r", "r", "r")),
          collapse = "\n"),
    if (!is.null(parts$by_endpoint) && nrow(parts$by_endpoint) > 1L) {
      paste(c("**Which endpoint is doing it.** A schedule is only as shared as",
              "its least shared part.",
              paste(knitr::kable(parts$by_endpoint,
                                 align = c("l", "r", "r", "r")),
                    collapse = "\n")),
            collapse = "\n\n")
    },
    paste0("*", .skeleton_footer, "*")
  )
  knitr::asis_output(paste(out, collapse = "\n\n"))
}

# Schedule map ----------------------------------------------------------------
#
# [skeleton_uniqueness()] answers "how many patients are alone?" with a number,
# and a number does not say whether that is a protocol with two stragglers or a
# study where every patient was sampled ad hoc. The picture does, immediately:
# a coarsened cohort on a real protocol grid draws as vertical stripes with a
# ragged right edge (visits missing from the end), and an uncoarsenable one
# draws as scatter.
#
# Base graphics on purpose. This is a diagnostic a user runs mid-analysis on a
# restricted machine, so it should not depend on ggplot2 being installed, and a
# points-on-a-grid plot is exactly what base draws well.
# Deliberately red-free. Red means one thing in this figure -- "this is the bit
# that identifies somebody" -- and it meant three before: the second endpoint
# was `#e66101`, which at screen distance is the same colour as the `#d7191c`
# used for a unique schedule and for a singleton visit time. A reader could not
# tell whether a red dot was an endpoint or a warning.
.identifying_colour <- "#d7191c"

.schedule_palette <- function(n) {
  base <- c("#1f78b4", "#33a02c", "#6a3d9a", "#1b9e77", "#8c6d31", "#525252")
  if (n <= length(base)) base[seq_len(n)] else
    grDevices::hcl.colors(n, "Blue-Yellow")
}

#' Draw a cohort's dosing and observation schedule
#'
#' The picture behind [skeleton_uniqueness()]. One row per patient, one mark
#' per event: when they were dosed, and when each endpoint was observed. Read
#' it to decide whether a uniqueness count is a real problem or ordinary
#' an ordinary gap in follow-up.
#'
#' Two panels:
#'
#' - **the map** -- patients ordered by how long they were followed, so a
#'   ragged right-hand edge reads as a staircase. That edge is follow-up
#'   ending, whether because a patient discontinued or because the study has
#'   not reached their later visits yet. A patient whose observation
#'   schedule no other patient shares is marked in the margin, and their label
#'   is drawn in red.
#' - **the visit histogram** -- how many patients were observed at each time on
#'   the grid. A protocol grid gives tall bars at a handful of times. A bar of
#'   height one is a moment only one patient was sampled at, which is precisely
#'   what [synpmx_avatar()] would copy verbatim onto an avatar, and those bars
#'   are drawn in red.
#'
#' By default the times are **coarsened first**, so the picture shows the grid
#' [synpmx_avatar()] actually generates on. Pass `coarsen_time = FALSE` to see
#' the recorded times instead; drawing it both ways is the quickest way to see
#' what coarsening bought.
#'
#' Source-derived, like every diagnostic here: keep the figure inside the safe
#' environment.
#'
#' @param data A PMX dataset -- normally the source.
#' @param roles Explicit roles from [pmx_roles()].
#' @param coarsen_time Draw the coarsened visit grid (`TRUE`, the default) or
#'   the recorded times as given (`FALSE`).
#' @param max_patients Draw at most this many patients, evenly spread through
#'   the ordering so the shape of the cohort survives. Default 80.
#' @param main Plot title. Defaults to a description of what is drawn.
#'
#' @return The [skeleton_uniqueness()] table for the drawn data, invisibly.
#' @seealso [skeleton_uniqueness()], [synpmx_avatar()].
#' @export
#' @examples
#' data <- pmx_simulated_fixture(20)
#' roles <- pmx_roles(
#'   id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
#'   cmt = "CMT", dvid = "DVID", covariates = "WT"
#' )
#' plot_pmx_schedule(data, roles)
plot_pmx_schedule <- function(data, roles, coarsen_time = TRUE,
                              max_patients = 80L, main = NULL) {
  .assert_roles(data, roles)
  if (isTRUE(coarsen_time)) data <- .coarsen_source_time(data, roles)$source
  screen <- skeleton_uniqueness(data, roles)

  time <- suppressWarnings(as.numeric(data[[roles$time]]))
  id <- as.character(data[[roles$id]])
  observed <- .observation_rows(data, roles, require_present = TRUE) &
    is.finite(time)
  dosed <- .dose_rows(data, roles) & is.finite(time)
  endpoint <- .endpoint(data, roles)

  # Ordered by last observation, so the right-hand edge is the follow-up curve.
  # Ties broken by observation count, so two patients followed equally long sit
  # next to each other with the sparser one below.
  subjects <- .unique_in_order(id)
  last_seen <- vapply(subjects, function(s) {
    values <- time[observed & id == s]
    if (length(values)) max(values) else -Inf
  }, numeric(1))
  n_seen <- vapply(subjects, function(s) sum(observed & id == s), integer(1))
  subjects <- subjects[order(last_seen, n_seen, subjects)]
  if (length(subjects) > max_patients) {
    keep <- unique(round(seq(1, length(subjects), length.out = max_patients)))
    subjects <- subjects[keep]
    data_note <- sprintf("showing %d of %d patients", length(subjects),
                         length(unique(id)))
  } else {
    data_note <- sprintf("%d patients", length(subjects))
  }
  row <- match(id, subjects)
  drawn <- !is.na(row)

  unique_schedule <- screen$subject_id[screen$unique_schedule]
  endpoints <- sort(unique(endpoint[observed & drawn]))
  colours <- .schedule_palette(length(endpoints))
  names(colours) <- endpoints

  old <- graphics::par(no.readonly = TRUE)
  # `par()` does not undo `layout()` -- it is a separate mechanism -- so the
  # split has to be cleared explicitly or the caller's next plot lands in the
  # top panel of a two-panel screen. Registered first so it runs last, after
  # the saved parameters are back.
  on.exit(graphics::layout(1L), add = TRUE)
  on.exit(graphics::par(old), add = TRUE, after = FALSE)
  label_cex <- if (length(subjects) > 45L) 0.45 else 0.6
  graphics::layout(matrix(c(1L, 2L), ncol = 1L), heights = c(3, 1))
  graphics::par(mar = c(2.2, 7.5, 3.2, 1.2), cex.axis = 0.8)

  x_range <- range(time[drawn & (observed | dosed)], finite = TRUE)
  graphics::plot(NA, xlim = x_range, ylim = c(0.5, length(subjects) + 0.5),
                 xlab = "", ylab = "", yaxt = "n", bty = "n")
  # The rows worth looking at, marked by their background rather than by the
  # colour of their label: 15 red labels in a column of 21 is not a signal.
  flagged_rows <- which(subjects %in% unique_schedule)
  if (length(flagged_rows)) {
    graphics::rect(x_range[[1L]], flagged_rows - 0.5,
                   x_range[[2L]], flagged_rows + 0.5,
                   col = "#fdecea", border = NA)
  }
  graphics::abline(h = seq_along(subjects), col = "grey93")
  # Doses first and underneath: a vertical tick, so an observation drawn at the
  # same time still reads.
  graphics::points(time[dosed & drawn], row[dosed & drawn], pch = 124,
                   col = "grey35", cex = 0.9)
  # Endpoints are offset within the row rather than drawn on top of each other.
  # Where a study measures two endpoints at the same visits, overplotting hides
  # the one thing the picture is for -- seeing WHICH endpoint is missing where,
  # which is usually what makes a patient's schedule unique.
  offset <- if (length(endpoints) > 1L) {
    stats::setNames((seq_along(endpoints) - (length(endpoints) + 1) / 2) *
                      (0.6 / length(endpoints)), endpoints)
  } else stats::setNames(0, endpoints)
  for (name in endpoints) {
    hit <- observed & drawn & endpoint == name
    graphics::points(time[hit], row[hit] + offset[[name]], pch = 16,
                     cex = 0.55, col = colours[[name]])
  }
  graphics::axis(2, at = seq_along(subjects), labels = subjects, las = 1,
                 tick = FALSE, cex.axis = label_cex,
                 col.axis = "grey25", line = -0.6)
  for (i in flagged_rows) {
    graphics::axis(2, at = i, labels = subjects[i], las = 1, tick = FALSE,
                   cex.axis = label_cex, col.axis = .identifying_colour,
                   line = -0.6)
  }
  graphics::title(
    main = main %||% paste0(
      "Dosing and observation schedule, ",
      if (isTRUE(coarsen_time)) "coarsened onto the shared visit grid"
      else "recorded times as given"
    ),
    xlab = "", cex.main = 1
  )
  graphics::mtext(
    sprintf(paste("%s. Dot colour is the ENDPOINT (see legend). The %d shaded",
                  "rows are patients whose exact set of visits nobody else has"),
            data_note, sum(subjects %in% unique_schedule)),
    side = 3, line = 0.2, cex = 0.7, col = "grey30"
  )
  graphics::legend("bottomright", legend = c(endpoints, "dose"),
                   pch = c(rep(16, length(endpoints)), 124),
                   col = c(unname(colours), "grey35"), bty = "n",
                   cex = 0.7, horiz = TRUE)

  # How many *patients* -- not rows -- were observed at each time. A height of
  # one is a moment one patient alone was sampled at, and that is what the
  # verbatim skeleton copy would hand to every avatar anchored there.
  per_time <- table(unique(data.frame(
    t = time[observed], s = id[observed], stringsAsFactors = FALSE
  ))$t)
  times <- as.numeric(names(per_time))
  heights <- as.integer(per_time)
  graphics::par(mar = c(3.6, 7.5, 1.4, 1.2))
  graphics::plot(NA, xlim = x_range, ylim = c(0, max(heights, 1L)),
                 xlab = "", ylab = "", bty = "n")
  graphics::segments(times, 0, times, heights,
                     col = ifelse(heights == 1L, .identifying_colour, "#4d4d4d"),
                     lwd = 2)
  graphics::title(xlab = sprintf("Time (%s, source units)", roles$time),
                  ylab = "patients", cex.lab = 0.85, line = 2.2)
  graphics::mtext(sprintf(paste("How many patients were observed at each time.",
                                "%d bar(s) in red reach only 1 -- a moment one",
                                "patient alone was sampled at"),
                          sum(heights == 1L)),
                  side = 3, line = 0.1, cex = 0.7, col = "grey30")
  invisible(screen)
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
  print(.format_for_print(utils::head(as.data.frame(x), 12L)), row.names = FALSE)
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

# Masking report ---------------------------------------------------------------
#
# `attr(synthetic, "pmx_settings")` is a flat list of thirty-odd names, and
# every report that wanted to show it -- the demo vignette, four study scripts
# under scripts_private/ -- built its own `data.frame()` of labels by hand.
# They drifted, which is how a report ended up printing an unrounded
# `pattern_sampled_fraction` of 0.142857 under a label that did not describe
# it, and how "route arms / compatible event groups" ended up as one cell
# holding two unrelated numbers. One implementation, here, so the labels are
# maintained once and every study report says the same thing.
#
# Three columns rather than two on purpose: the value is meaningless without
# the sentence next to it, and that sentence is what a reader needs at the
# moment they meet the number, not in a paragraph underneath the table.

.grid_explanation <- function(settings) {
  derived <- settings$time_grid_derived_rows
  total <- settings$time_grid_rows
  share <- if (is.null(total) || is.na(total) || !total) NA_real_ else
    derived / total
  switch(
    as.character(settings$time_grid),
    nominal = paste("every visit was snapped to the declared `nominal_time`,",
                    "which is the protocol grid. This is the reliable case"),
    derived = paste("no usable `nominal_time`, so a grid was inferred from the",
                    "recorded times themselves. Declaring `nominal_time` is",
                    "better"),
    mixed = paste0(
      "some rows had a usable `nominal_time` and the rest did not",
      if (is.na(share)) "" else sprintf(
        "; %d of %d event rows (%.0f%%) were snapped to an *inferred* grid",
        derived, total, 100 * share),
      ". Fill in `nominal_time` on those rows to move this to `nominal`"),
    off = "`coarsen_time = FALSE`: recorded times were copied as they stand",
    none = "no time column could be coarsened; recorded times were copied",
    "unrecognised grid"
  )
}

.masking_rows <- function(settings, before = NULL) {
  # Every row that counts patients or avatars shows both the count and the
  # share, because neither is readable alone: "5%" of a 21-patient cohort is one
  # patient, and "15" means nothing without the cohort size next to it.
  n_source <- settings$source_subjects
  n_built <- settings$n_subjects
  both <- function(count, total) {
    if (is.null(count) || length(count) != 1L || is.na(count)) return("--")
    if (is.null(total) || is.na(total) || !total) return(format(count))
    sprintf("%s (%.0f%%)", format(count), 100 * count / total)
  }
  # Same, from the other direction: a fraction, shown with the count it implies.
  share <- function(x, total) {
    if (is.null(x) || is.na(x)) return("--")
    if (is.null(total) || is.na(total) || !total) {
      return(paste0(round(100 * x), "%"))
    }
    sprintf("%d of %d (%.0f%%)", round(x * total), total, 100 * x)
  }
  pct <- function(x) if (is.null(x) || is.na(x)) "--" else
    paste0(round(100 * x), "%")
  num <- function(x) if (is.null(x) || length(x) != 1L || is.na(x)) "--" else
    format(x)
  header <- function(text) c(paste0("**", text, "**"), "", "")

  rows <- list(
    header("Who was available to build on"),
    c("Patients in the source", num(settings$source_subjects), ""),
    c("&nbsp;&nbsp;excluded as structurally extreme",
      both(settings$anchors_screened_out, n_source),
      "`screen`: follow-up or dose count over twice the cohort's 90th percentile"),
    c("&nbsp;&nbsp;excluded, route arm too small",
      both(settings$anchors_route_excluded, n_source),
      "`on_donor_shortfall`: a route arm holding fewer than k + 1 patients"),
    c("&nbsp;&nbsp;left to anchor avatars on",
      both(settings$anchors_available, n_source),
      "an excluded patient still contributes as a donor"),
    c("Avatars built", num(settings$n_subjects),
      "cohort size is unaffected by the exclusions above"),

    header("Donor pools: who may be blended with whom"),
    c("Administration routes", num(settings$routes),
      "oral, infusion, and so on. Donors are NEVER blended across a route, so each is a separate pool"),
    c("Dose/schedule groups", num(settings$compatible_event_groups),
      "patients with an identical dose pattern and endpoint set. Donors are looked for here first; many small groups means the search falls back to the wider route pool"),

    header("How much of one real patient reaches one avatar"),
    c("Donor floor, k", num(settings$k),
      "real patients blended into each avatar"),
    c("Largest share one donor may hold", num(settings$max_donor_weight),
      "`max_donor_weight`"),
    c("&nbsp;&nbsp;that cap actually bound on",
      share(settings$cap_binding_fraction, n_built),
      "of avatars. Near 100% means the cap, not distance, is setting the weights"),
    c("Effective donors per avatar, mean",
      num(round(settings$mean_effective_donors, 2)),
      "1 / sum(w^2). This, not k, is how many patients an avatar is really made of"),

    header("Visit schedule: WHEN patients were observed"),
    c("Visit grid used", num(settings$time_grid), .grid_explanation(settings)),
    if (!is.null(before)) {
      c("Unique observation schedules, before coarsening",
        both(attr(before, "n_unique_schedule"), n_source),
        "patients whose list of observation times nobody else shares")
    },
    c("Unique observation schedules, after coarsening",
      both(settings$unique_schedule_n, n_source),
      "the count that matters: an avatar copies its anchor's times verbatim"),
    c("&nbsp;&nbsp;because of a one-off observation time",
      both(settings$unique_obs_time_n, n_source),
      "sampled when nobody else was. Declaring `nominal_time` is the fix"),
    c("&nbsp;&nbsp;because of which visits they attended",
      both(settings$unique_visit_set_n, n_source),
      "every time is shared. The visits themselves are missing -- a missed visit, a discontinuation, or follow-up that has not reached them -- and no grid can fix that"),

    header("Visit sets: WHICH of those visits each patient attended"),
    c("Distinct visit sets in the source", num(settings$patterns_total),
      "a visit set is which of the shared grid visits one patient actually had"),
    c(sprintf("&nbsp;&nbsp;held by fewer than %s patients, so not reused",
              num(settings$min_pattern_share)),
      both(settings$patterns_dropped, settings$patterns_total),
      "`min_pattern_share` is that threshold. These visit sets are lost, not approximated"),
    c("&nbsp;&nbsp;real patients holding those",
      both(settings$subjects_with_dropped_pattern, n_source),
      "those patients are NOT removed -- they still anchor avatars and still act as donors. Only their particular pattern of absences stops being copied"),
    c("Avatars given a visit set from the pool",
      share(settings$pattern_sampled_fraction, n_built),
      "drawn from the sets that cleared the threshold, or built from their shape -- never from their own anchor alone"),
    c("&nbsp;&nbsp;of those, misses placed fresh",
      share(settings$pattern_generated_fraction, n_built),
      "the kind of missingness was reused; exactly which visits were missed was invented"),
    c("&nbsp;&nbsp;of those, miss count moved",
      share(settings$pattern_shifted_fraction, n_built),
      "no arrangement at the wanted number of missing visits was free, so the count moved by a visit or two. Misses at the END of a record are the case that forces it, because for a given count there is exactly one such arrangement"),
    c("&nbsp;&nbsp;of those, a rare set swapped for a shared one",
      share(settings$pattern_substituted_fraction, n_built),
      "the anchor's own set was held by nobody else and no arrangement was free, so the group's most widely held set was used instead -- less faithful to that avatar, and it discloses nothing"),
    c("&nbsp;&nbsp;of those, moved to a different anchor",
      share(settings$pattern_reanchored_fraction, n_built),
      "the first anchor's own set was shared by nobody and nothing legal could be placed, so this avatar was anchored elsewhere. Every source patient stays a donor and stays available to anchor others"),
    c("Avatars keeping their anchor's own visit set",
      share(1 - (settings$pattern_sampled_fraction %||% NA_real_), n_built),
      "not a problem in itself: if several real patients share that set, copying it identifies nobody. Only the next row is a disclosure"),
    c("**Avatars carrying a visit set nobody else shares**",
      both(settings$identifying_visit_sets, n_built),
      "**this is the row that must be 0%.** That pattern of which visits have observations belongs to one real patient. It is non-zero only when the schedule group has no shared set to substitute; the run alerts when it happens"),

    c("&nbsp;&nbsp;of those, dosing re-truncated",
      share(settings$dose_truncated_fraction, n_built),
      "the anchor stopped dosing at a depth nobody else used, so the avatar stops at a different one -- shared, or used by nobody. Truncating a schedule to a real dose time is protocol-valid in a way that moving dose times is not"),
    c("Distinct dose schedules in the source",
      num(settings$dose_regimens_source), ""),
    c("&nbsp;&nbsp;represented in the synthetic cohort",
      both(settings$dose_regimens_represented, settings$dose_regimens_source),
      "a regimen only one patient received cannot be given to an avatar without pointing at them, so it is not represented at all. This is the cost of the guarantee below, and on a small cohort it is unavoidable rather than a setting to tune"),
    c("**Avatars carrying a dose schedule nobody else shares**",
      both(settings$identifying_dose_schedules, n_built),
      "**must also be 0%.** Dose events are copied from the anchor verbatim, so patients whose dose times nobody shares are not built upon. Non-zero only when EVERY patient is in that position, which individualised dosing can cause"),

    header("Dose"),
    c("Amounts recomputed from a covariate",
      if (is.na(settings$dose_basis)) "**no**" else
        paste0("**yes**, from `", settings$dose_basis, "`",
               if (isTRUE(settings$dose_basis_declared)) " (declared)"
               else " (inferred)"),
      settings$dose_basis_note %||%
        "inferred only where dose / covariate collapses onto a few protocol levels; fails closed. Declare `dose_covariate` in `pmx_roles()` to skip the inference"),
    if (is.na(settings$dose_basis)) {
      c("&nbsp;&nbsp;so `amt` is copied verbatim", "from the anchor",
        "each avatar's implied dose per kg is therefore its anchor's, not its own, and the amount still encodes one real patient's covariate. Declare `dose_covariate` if this study is weight- or BSA-based")
    },
    if (!is.na(settings$dose_basis) && !isTRUE(settings$dose_basis_declared)) {
      c("&nbsp;&nbsp;protocol levels found",
        paste(signif(settings$dose_levels, 4), collapse = ", "),
        paste0("dose per unit of `", settings$dose_basis,
               "`; every amount was snapped to the nearest of these"))
    }
  )
  do.call(rbind, Filter(Negate(is.null), rows))
}

#' Report what each masking mechanism did, and what it cost
#'
#' [synpmx_avatar()] records everything it removed on the
#' `"pmx_settings"` attribute of its result. This turns that flat list into the
#' table to read after a run: who was left to build on, how many real patients
#' reach one avatar, what the visit grid managed to collapse, which visit sets
#' were too rare to reuse, and whether dose amounts were recomputed.
#'
#' Every row carries a sentence saying what the number means, because none of
#' them mean anything on their own. The rows worth looking at hardest:
#'
#' - **Unique observation schedules, after coarsening** -- patients whose list
#'   of observation times nobody else shares. An avatar anchored on one wears a
#'   schedule belonging to one real person. Its two sub-rows have opposite
#'   remedies: a one-off observation time is what declaring `nominal_time`
#'   fixes, and a unique set of *attended* visits is missing visits, which no grid
#'   touches.
#' - **Shared by too few patients, so not reused** -- real patterns of missing
#'   visits and dose interruptions that will not appear in the synthetic data.
#'   Discarding them is what stops an avatar carrying a schedule traceable to
#'   one person. If this study's interruptions matter, lower
#'   `min_pattern_share` (2 is the lowest value that still guarantees no
#'   synthetic patient has a schedule unique to a real one).
#' - **Avatars carrying a visit set nobody else shares** -- the only row here
#'   that is a disclosure rather than a fidelity cost, and the one to drive to
#'   zero. Keeping the *anchor's own* set is fine whenever several real patients
#'   share it; it is a problem only when that set is unique to one of them.
#'   Where a shared set exists, one is substituted automatically, so this row is
#'   non-zero only when the whole schedule group has nothing shareable.
#' - **Amounts recomputed from a covariate** -- says outright whether
#'   weight-based or body-surface-area dosing was detected, and when it was
#'   not, why not. Detection is deliberately conservative: it fails closed and
#'   leaves amounts alone rather than rewriting a study that is not
#'   dose-proportional.
#'
#' At the default `min_pattern_share = 2`, "shared by too few patients" and
#' "real patients holding those" are necessarily equal -- a set is discarded
#' exactly when fewer than two patients share it, so every discarded set has
#' one holder. They diverge only at a floor of 3 or more.
#'
#' Marked `"restricted_not_releasable"` when `source` is supplied, since the
#' before-coarsening row then reads the source.
#'
#' @param synthetic A dataset from [synpmx_avatar()], carrying its
#'   `"pmx_settings"` attribute.
#' @param source Optionally the source dataset. Supplying it (with `roles`)
#'   adds the before-coarsening schedule count, so the table shows what
#'   coarsening removed rather than only what was left.
#' @param roles Explicit roles from [pmx_roles()]. Required with `source`.
#'
#' @return A `pmx_masking_report` data frame with columns `Quantity`, `Value`,
#'   and `What it means`. Section headers appear as rows whose `Quantity` is
#'   bold and whose other cells are empty.
#' @seealso [skeleton_uniqueness()], [compare_pmx_proximity()],
#'   [synpmx_avatar()].
#' @export
#' @examples
#' data <- pmx_simulated_fixture(30)
#' roles <- pmx_roles(
#'   id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
#'   cmt = "CMT", dvid = "DVID", covariates = "WT"
#' )
#' synthetic <- suppressWarnings(synpmx_avatar(data, roles, seed = 1))
#' pmx_masking_report(synthetic, data, roles)
pmx_masking_report <- function(synthetic, source = NULL, roles = NULL) {
  settings <- attr(synthetic, "pmx_settings")
  if (is.null(settings)) {
    stop("`synthetic` carries no `pmx_settings` attribute; it did not come ",
         "from `synpmx_avatar()`.", call. = FALSE)
  }
  if (!is.null(source) && is.null(roles)) {
    stop("`roles` is required when `source` is supplied.", call. = FALSE)
  }
  before <- if (is.null(source)) NULL else skeleton_uniqueness(source, roles)
  out <- as.data.frame(.masking_rows(settings, before),
                       stringsAsFactors = FALSE)
  names(out) <- c("Quantity", "Value", "What it means")
  rownames(out) <- NULL
  out <- structure(out, class = c("pmx_masking_report", "data.frame"))
  if (is.null(source)) out else
    .mark_release(out, "restricted_not_releasable")
}

# `**bold**` and `&nbsp;` are markdown, for the knitted table that is the point
# of this object. On a console they are noise, so the console path strips them
# and indents with spaces instead.
.plain_masking <- function(x) {
  x$Quantity <- gsub("&nbsp;&nbsp;", "  ", x$Quantity, fixed = TRUE)
  x$Quantity <- gsub("\\*\\*", "", x$Quantity)
  x$Value <- gsub("\\*\\*", "", x$Value)
  x
}

#' @export
print.pmx_masking_report <- function(x, ...) {
  plain <- .plain_masking(as.data.frame(x))
  cat("What the masking mechanisms did\n\n")
  for (i in seq_len(nrow(plain))) {
    if (!nzchar(plain$Value[[i]])) {
      cat(if (i > 1L) "\n" else "", plain$Quantity[[i]], "\n", sep = "")
      next
    }
    cat(sprintf("  %-48s %s\n", plain$Quantity[[i]], plain$Value[[i]]))
    if (nzchar(plain$`What it means`[[i]])) {
      cat(.wrap_plain(plain$`What it means`[[i]], initial = "      ",
                      prefix = "      "), "\n", sep = "")
    }
  }
  invisible(x)
}

#' @exportS3Method knitr::knit_print
knit_print.pmx_masking_report <- function(x, ...) {
  knitr::knit_print(knitr::kable(
    as.data.frame(x), row.names = FALSE, align = c("l", "r", "l"),
    caption = "What each masking mechanism did, and what it cost."
  ))
}

# Nearest-neighbour proximity ---------------------------------------------------
#
# Every masking mechanism has a measurement except the first one. Blending is
# what protects the *values* -- each avatar's covariates and trajectories are
# averaged across at least `k` real donors -- and `cap_binding_fraction` and
# `mean_effective_donors` describe how the blend was *built*, not whether the
# result landed too close to somebody real. `skeleton_uniqueness()` answers the
# structural question exactly, by counting who shares what; this answers the
# geometric one.
#
# Raw distance to the closest real record is the obvious statistic and a bad
# one. It has no natural scale, so it is only interpretable against a
# comparison, and against zero it mostly tracks cohort size and dimension: the
# nearest of N points gets closer as N grows, which would make a larger source
# look worse while actually making blending safer. Worse, if the synthetic
# marginals match the source -- the goal -- then those distances converge on the
# real-to-real distances anyway, so a small value is evidence the generator
# worked.
#
# The fix is to compare each point against its own neighbourhood instead. For
# every subject, ask whether its nearest neighbour lies in its own dataset or
# the other one. Under the ideal -- a synthetic subject indistinguishable from
# one more real draw -- that is a coin flip, and the statistic sits at one half.
# It falls toward zero when synthetic subjects sit closer to real ones than to
# each other, which is memorisation. It rises toward one when the two sets have
# separated, which is a utility failure rather than a privacy one. Because it
# compares distances measured in the same neighbourhood, the scale cancels.
#
# The null is built by running the identical statistic on two halves of the real
# cohort. Every small-sample artefact -- one set minimising over n - 1
# candidates while the other has n, the dependence between nearest-neighbour
# links, the small cohorts pharmacometrics actually has -- is then present in
# the null and the test alike, and cancels. That is worth more than any
# distributional assumption, and it is why the null is not simply asserted to be
# one half.

.nearest_within <- function(distance) {
  diag(distance) <- Inf
  apply(distance, 1L, min)
}

# One dataset's subjects against another's, both already in a common space.
.adversarial_accuracy <- function(a, b) {
  if (nrow(a) < 2L || nrow(b) < 2L) return(NA_real_)
  across <- as.matrix(stats::dist(rbind(a, b)))
  index_a <- seq_len(nrow(a))
  index_b <- nrow(a) + seq_len(nrow(b))
  own_a <- .nearest_within(across[index_a, index_a, drop = FALSE])
  own_b <- .nearest_within(across[index_b, index_b, drop = FALSE])
  other_a <- apply(across[index_a, index_b, drop = FALSE], 1L, min)
  other_b <- apply(across[index_b, index_a, drop = FALSE], 1L, min)
  mean(c(mean(other_a > own_a), mean(other_b > own_b)))
}

#' Are synthetic subjects sitting too close to real ones?
#'
#' The measurement for donor blending, the one masking mechanism
#' [synpmx_avatar()] applies to the *values* rather than the structure.
#' [skeleton_uniqueness()] answers the structural question by counting who
#' shares which schedule; this answers the geometric one, by asking whether each
#' subject's nearest neighbour lies in its own dataset or the other one.
#'
#' The reported statistic is a nearest-neighbour adversarial accuracy in
#' \eqn{[0, 1]}:
#'
#' - **near 0.5** — a synthetic subject is no more like a real subject than one
#'   real subject is like another. This is the target.
#' - **toward 0** — synthetic subjects sit closer to real subjects than to each
#'   other. That is memorisation, and it is the privacy failure.
#' - **toward 1** — the two sets have separated. Privacy is fine and utility is
#'   not.
#'
#' Raw distance to the closest real record is deliberately not the headline. It
#' has no natural scale, and measured against zero it mostly tracks cohort size —
#' the nearest of `N` points gets closer as `N` grows, so a larger source would
#' score worse while blending across more donors actually makes it safer. The
#' quantiles are still returned for context, alongside the real-to-real
#' quantiles they should be read against.
#'
#' The null interval comes from running the identical statistic on two halves of
#' the **source** cohort, so every small-sample artefact is present in the null
#' and the observed value alike and cancels. At the cohort sizes pharmacometrics
#' works with, that interval is wide: this will catch a blatant leak, not a
#' subtle one. Treat a value inside the interval as "nothing detected", never as
#' "nothing there".
#'
#' Marked `"restricted_not_releasable"`: it reads the source.
#'
#' @param source Source PMX data.
#' @param synthetic Generated synthetic PMX data.
#' @param roles Explicit roles from [pmx_roles()].
#' @param replicates Split-half replicates used to build the null. Default 50.
#' @param seed Seed for the subsampling and splits. The caller's RNG is left
#'   untouched.
#' @param pca_variance Variance retained when both datasets are projected into a
#'   common profile space. Default 0.90, matching [synpmx_avatar()].
#'
#' @return A one-row `pmx_proximity` data frame: `adversarial_accuracy`,
#'   `null_lower` / `null_upper` (the central 95% of the split-half null),
#'   `verdict`, `n_compared` (patients per side, the same on both arms and in
#'   the null), `n_null_replicates`, and the 5th-percentile nearest-neighbour
#'   distances `synthetic_to_source_q05` and `source_to_source_q05`.
#' @seealso [skeleton_uniqueness()], [compare_pmx_distributions()].
#' @export
#' @examples
#' data <- pmx_simulated_fixture(40)
#' roles <- pmx_roles(
#'   id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
#'   cmt = "CMT", dvid = "DVID", covariates = "WT"
#' )
#' synthetic <- suppressWarnings(synpmx_avatar(data, roles, seed = 1))
#' compare_pmx_proximity(data, synthetic, roles, replicates = 10)
compare_pmx_proximity <- function(source, synthetic, roles, replicates = 50L,
                                  seed = 1L, pca_variance = 0.90) {
  .assert_roles(source, roles)
  .assert_roles(synthetic, roles)
  if (!is.numeric(replicates) || length(replicates) != 1L ||
      is.na(replicates) || replicates < 1) {
    stop("`replicates` must be one positive integer.", call. = FALSE)
  }

  # Both datasets are projected together, so the coordinates are comparable by
  # construction rather than by a transform carried between two fits.
  marked_source <- source
  marked_synthetic <- synthetic
  marked_source[[roles$id]] <- paste0("S:", as.character(source[[roles$id]]))
  marked_synthetic[[roles$id]] <- paste0("G:",
                                         as.character(synthetic[[roles$id]]))
  # `synpmx_avatar()` drops undeclared columns, so the two tables need not share
  # a schema; the roles do, and only role-named columns feed the profile anyway.
  shared <- intersect(names(marked_source), names(marked_synthetic))
  combined <- rbind(marked_source[, shared, drop = FALSE],
                    marked_synthetic[, shared, drop = FALSE])
  profiles <- .build_profiles(combined, roles, pca_variance)
  coordinates <- profiles$coordinates
  is_synthetic <- startsWith(as.character(profiles$subjects), "G:")
  real <- coordinates[!is_synthetic, , drop = FALSE]
  fake <- coordinates[is_synthetic, , drop = FALSE]

  .with_local_seed(seed, {
    # Both arms of the comparison, and both arms of the null, are run at one
    # size. A statistic built from nearest neighbours depends on how many
    # candidates there were, so the sizes have to match or the null is measuring
    # something else.
    size <- min(nrow(real) %/% 2L, nrow(fake))
    if (size < 3L) {
      out <- data.frame(
        adversarial_accuracy = NA_real_, null_lower = NA_real_,
        null_upper = NA_real_, verdict = "too few subjects to compare",
        n_compared = size, n_null_replicates = as.integer(replicates),
        synthetic_to_source_q05 = NA_real_,
        source_to_source_q05 = NA_real_, stringsAsFactors = FALSE
      )
      return(.mark_release(
        structure(out, class = c("pmx_proximity", "data.frame")),
        "restricted_not_releasable"
      ))
    }

    take <- function(x, n) x[sample.int(nrow(x), n), , drop = FALSE]
    observed <- .adversarial_accuracy(take(fake, size), take(real, size))
    null <- vapply(seq_len(as.integer(replicates)), function(i) {
      split <- sample.int(nrow(real))
      .adversarial_accuracy(
        real[split[seq_len(size)], , drop = FALSE],
        real[split[size + seq_len(size)], , drop = FALSE]
      )
    }, numeric(1))
    null <- null[is.finite(null)]
    bounds <- if (length(null)) {
      unname(stats::quantile(null, c(0.025, 0.975)))
    } else c(NA_real_, NA_real_)

    # Context only, and reported at the 5th percentile because a leak lives in
    # the left tail rather than the middle.
    between <- as.matrix(stats::dist(rbind(fake, real)))[
      seq_len(nrow(fake)), nrow(fake) + seq_len(nrow(real)), drop = FALSE
    ]
    within <- .nearest_within(as.matrix(stats::dist(real)))
    q05 <- function(x) unname(stats::quantile(x, 0.05))

    verdict <- if (!is.finite(observed) || anyNA(bounds)) {
      "no null available"
    } else if (observed < bounds[[1L]]) {
      "below the null: synthetic subjects sit closer to real ones than real ones do to each other"
    } else if (observed > bounds[[2L]]) {
      "above the null: the two sets have separated, which is a utility concern"
    } else {
      "within the null: nothing detected at this cohort size"
    }

    out <- data.frame(
      adversarial_accuracy = observed,
      null_lower = bounds[[1L]], null_upper = bounds[[2L]],
      verdict = verdict, n_compared = size,
      n_null_replicates = as.integer(replicates),
      synthetic_to_source_q05 = q05(apply(between, 1L, min)),
      source_to_source_q05 = q05(within),
      stringsAsFactors = FALSE
    )
    .mark_release(
      structure(out, class = c("pmx_proximity", "data.frame")),
      "restricted_not_releasable"
    )
  })
}

# The statistic's name is jargon and its scale is not the usual one -- higher
# is not better, 0.5 is -- so the printout leads with the question in words and
# gives the number a sentence of its own. "(16 per side)" used to sit at the end
# of the headline with nothing saying what a side was.
.proximity_lines <- function(x) {
  if (is.na(x$adversarial_accuracy)) {
    return(c(
      "Question: is a synthetic patient closer to a real patient than real patients are to each other?",
      paste("Answer:", x$verdict)
    ))
  }
  reading <- if (x$adversarial_accuracy < x$null_lower) {
    "Too close. Synthetic patients sit nearer to real ones than real ones do to each other -- this is memorisation, and the privacy failure this check exists to find."
  } else if (x$adversarial_accuracy > x$null_upper) {
    "Too far apart. The two sets have separated, so a classifier could tell them apart. That is a utility problem, not a privacy one."
  } else {
    "Nothing detected. The value sits inside the null interval, which at this cohort size is wide -- read it as 'no blatant leak', never as 'no leak'."
  }
  c(
    "Question: is a synthetic patient closer to a real patient than real patients are to each other?",
    sprintf("Measured %.3f, on a scale where 0.5 means 'no closer' and is the target; 0 would mean every synthetic patient is glued to a real one.",
            x$adversarial_accuracy),
    sprintf("Expected %.3f to %.3f if nothing were wrong. That interval is not assumed -- it is the same statistic run %s times on two halves of the real cohort, %d patients per half, which is also how many synthetic patients were compared.",
            x$null_lower, x$null_upper, x$n_null_replicates %||% 50L,
            x$n_compared),
    paste("Verdict:", reading),
    sprintf("For context, distance to the nearest neighbour (5th percentile, so the closest pairs): synthetic-to-real %.3f versus real-to-real %.3f. These are only comparable to each other; the units are PCA profile space.",
            x$synthetic_to_source_q05, x$source_to_source_q05)
  )
}

#' @export
print.pmx_proximity <- function(x, ...) {
  cat("Restricted PMX nearest-neighbour proximity check\n\n")
  # `rbind()`-ing several reports to tabulate a set of datasets is the obvious
  # thing to do with these, and it keeps the class, so handle more than one row
  # rather than failing on a length-2 condition.
  if (nrow(x) != 1L) {
    print(.format_for_print(as.data.frame(x)), row.names = FALSE)
    cat("\n", .wrap_plain(paste(
      "0.5 is the target; toward 0 is memorisation. The null interval is wide",
      "at small cohorts, so a value inside it means nothing was detected.",
      "Source-derived; not releasable unless separately public or privately",
      "budgeted."
    )), "\n", sep = "")
    return(invisible(x))
  }
  for (line in .proximity_lines(x)) {
    cat(.wrap_plain(line, initial = "  ", prefix = "    "), "\n", sep = "")
  }
  cat("\n", .wrap_plain(paste(
    "Source-derived; not releasable unless separately public or privately",
    "budgeted."
  )), "\n", sep = "")
  invisible(x)
}

#' @exportS3Method knitr::knit_print
knit_print.pmx_proximity <- function(x, ...) {
  if (nrow(x) != 1L) return(knitr::knit_print(as.data.frame(x)))
  knitr::asis_output(paste(c(
    "**Nearest-neighbour proximity check.**",
    paste0("- ", .proximity_lines(x)),
    "",
    paste("*Source-derived; not releasable unless separately public or",
          "privately budgeted.*")
  ), collapse = "\n"))
}
