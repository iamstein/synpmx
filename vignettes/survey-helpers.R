# Helpers shared by the public-data surveys.
#
# `avatar-public-data-examples.Rmd` and `pca-public-data-examples.Rmd` draw the
# same figure from the same data in the same way, because the two surveys exist
# to be read against each other. Keeping one copy is what makes that true; two
# copies drift. A third generator's survey sources this file as well.
#
# Sourced from the vignette, which knits with its own directory as the working
# directory under `R CMD build` and may knit from the package root elsewhere.
# Not a vignette itself: it declares no engine, so the build does not try to
# knit it.

# Source and synthetic get the same two colours in every figure in both
# surveys.
comparison_colours <- c(Source = "#1B6CA8", Synthetic = "#D95F02")

observed_plot_data <- function(data, roles, dataset,
                               clock = "study_time") {
  observed <- as.character(data[[roles$evid]]) %in% c("0", "0.0")
  if (!is.null(roles$mdv)) {
    observed <- observed & as.character(data[[roles$mdv]]) %in% c("0", "0.0")
  }
  observed <- observed & !is.na(data[[roles$dv]])
  observation_rows <- which(observed)
  occasion <- rep(1L, length(observation_rows))
  tad <- rep(NA_real_, length(observation_rows))
  if (!is.null(roles$occasion)) {
    declared <- suppressWarnings(as.integer(
      data[[roles$occasion]][observation_rows]
    ))
    valid <- !is.na(declared) & declared >= 1L
    occasion[valid] <- declared[valid]
  }
  if (!is.null(roles$tad)) {
    declared <- suppressWarnings(as.numeric(data[[roles$tad]][observation_rows]))
    valid <- is.finite(declared)
    tad[valid] <- pmax(0, declared[valid])
  }
  subject_values <- data[[roles$id]]
  for (id in unique(subject_values[observation_rows])) {
    subject_rows <- which(!is.na(subject_values) & subject_values == id)
    events <- !(as.character(data[[roles$evid]][subject_rows]) %in%
                  c("0", "0.0"))
    if (!is.null(roles$amt)) {
      events <- events & as.numeric(data[[roles$amt]][subject_rows]) > 0
    }
    positions <- which(subject_values[observation_rows] == id)
    event_rows <- subject_rows[events]
    if (length(event_rows) && !is.null(roles$occasion)) {
      event_occasion <- suppressWarnings(as.integer(
        data[[roles$occasion]][event_rows]
      ))
      for (position in positions) {
        candidates <- event_rows[event_occasion == occasion[position]]
        if (length(candidates) && !is.finite(tad[position])) {
          origin <- min(as.numeric(data[[roles$time]][candidates]))
          tad[position] <-
            as.numeric(data[[roles$time]][observation_rows[position]]) - origin
        }
      }
    } else if (length(event_rows)) {
      dose_times <- sort(unique(as.numeric(data[[roles$time]][event_rows])))
      occasion[positions] <- pmax(1L, findInterval(
        as.numeric(data[[roles$time]][observation_rows[positions]]),
        dose_times
      ))
      occasion[positions] <- pmin(occasion[positions], length(dose_times))
      tad[positions] <-
        as.numeric(data[[roles$time]][observation_rows[positions]]) -
        dose_times[occasion[positions]]
    }
  }
  plotted_time <- if (identical(clock, "tad")) tad else
    as.numeric(data[[roles$time]][observation_rows])
  data.frame(
    dataset = factor(dataset, levels = c("Source", "Synthetic")),
    subject = as.character(data[[roles$id]][observation_rows]),
    time = plotted_time,
    dv = as.numeric(data[[roles$dv]][observation_rows]),
    occasion = occasion,
    endpoint = if (is.null(roles$dvid)) "DV" else
      as.character(data[[roles$dvid]][observation_rows]),
    stringsAsFactors = FALSE
  )
}

# Every dataset below gets the same figure: observation rows only, source beside
# synthetic on a shared y axis, one row per endpoint. `clock = "tad"` plots time
# after dose instead of study time, which is the readable view once a study
# doses more than once.
overlay_plot <- function(source, synthetic, roles, title,
                         clock = "study_time",
                         x_label = "Study time (hours)", y_label = "DV",
                         log_y = FALSE, alpha = 0.3) {
  plotted <- rbind(
    observed_plot_data(source, roles, "Source", clock),
    observed_plot_data(synthetic, roles, "Synthetic", clock)
  )
  # One line per patient on a study-time axis, and one line per occasion on a
  # dose-relative one, where the profiles are meant to lie on top of each other.
  plotted$series <- if (identical(clock, "tad")) {
    interaction(plotted$dataset, plotted$subject, plotted$occasion)
  } else {
    interaction(plotted$dataset, plotted$subject)
  }
  figure <- ggplot2::ggplot(
    plotted,
    ggplot2::aes(time, dv, group = series, colour = dataset)
  ) +
    ggplot2::geom_line(alpha = alpha) +
    ggplot2::geom_point(alpha = alpha, size = 0.7) +
    # Endpoint down the side, source and synthetic across. `facet_grid()` frees
    # a scale per ROW, so this orientation is the one that gives each endpoint
    # its own y while holding source and synthetic on a shared one -- which is
    # the only arrangement the eye can compare. The transpose does the opposite
    # on both counts: it lets the two panels drift onto different axes, and it
    # squeezes every endpoint onto one, so a PK concentration reading single
    # digits flattens to a line next to a PD score in the hundreds.
    #
    # A single-endpoint study needs no row strip; it would print the endpoint's
    # name down the side of the only row, beside the axis title already naming
    # the same thing.
    (if (length(unique(plotted$endpoint)) > 1L) {
      ggplot2::facet_grid(endpoint ~ dataset, scales = "free_y",
                          switch = "y")
    } else {
      ggplot2::facet_wrap(~dataset)
    }) +
    ggplot2::scale_colour_manual(values = comparison_colours) +
    ggplot2::labs(x = x_label, y = y_label, colour = "Dataset", title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none",
                   strip.placement = "outside")
  if (isTRUE(log_y)) figure + ggplot2::scale_y_log10() else figure
}

# Facet rows are endpoints now, so a five-endpoint study needs five times the
# height a one-endpoint study does. Chunks pass this to `fig.height` rather than
# every figure inheriting the one-endpoint default and arriving unreadable.
overlay_height <- function(data, roles, per_endpoint = 2.2, minimum = 3.4) {
  observed <- as.character(data[[roles$evid]]) %in% c("0", "0.0")
  endpoints <- if (is.null(roles$dvid)) 1L else {
    length(unique(data[[roles$dvid]][observed]))
  }
  max(minimum, per_endpoint * endpoints)
}

# A survey run is an attempt, because a generator may refuse a dataset it
# cannot model. These three keep a chunk option safe when the run object was
# never created, either because the generator refused or because the source
# package is not installed. `data` and `roles` are taken lazily, so an object
# a missing package never created is never forced.
ran <- function(name) {
  environment <- knitr::knit_global()
  exists(name, envir = environment) &&
    isTRUE(get(name, envir = environment)$ok)
}

run_overlay_height <- function(name, data, roles, default = 4) {
  if (ran(name)) overlay_height(data, roles) else default
}

run_distribution_height <- function(name, data, roles, default = 4) {
  if (ran(name)) compare_pmx_distributions_height(data, roles) else default
}
