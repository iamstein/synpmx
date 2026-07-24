# Version 2 end-to-end demonstrations for the public nlmixr2data datasets.
# Install synpmx first with: R CMD INSTALL .

if (!requireNamespace("synpmx", quietly = TRUE)) {
  stop("Install synpmx before running this script: R CMD INSTALL .")
}
required_api <- c(
  "synpmx_empirical", "sampling_summary", "subject_property_summary"
)
missing_api <- setdiff(required_api, getNamespaceExports("synpmx"))
if (length(missing_api)) {
  stop(
    "The installed synpmx is older than this demo. Reinstall the ",
    "package from this repository with `R CMD INSTALL .` before running it."
  )
}
if (!requireNamespace("nlmixr2data", quietly = TRUE)) {
  stop("Install nlmixr2data to run the integration demonstrations.")
}

# synpmx_empirical() is one of the DP engines: complete and tested, but not
# under active development, carrying known open findings, and not
# independently privacy-audited (design/REVIEW_BACKLOG.md REV-023). This
# script uses `backend = "public"` throughout, which is exempt from the
# session-level acknowledgment gate, but the call is here anyway so a reader
# who copies this pattern onto `backend = "opendp"` sees it demonstrated.
synpmx::synpmx_enable_dp_engines()

load_dataset <- function(name) {
  environment <- new.env(parent = emptyenv())
  utils::data(list = name, package = "nlmixr2data", envir = environment)
  get(name, envir = environment, inherits = FALSE)
}

demo_budget <- function() {
  synpmx::pmx_budget_allocation(
    subject_count = .10, event = .15, timing = .15,
    covariates = .10, endpoints = .50, censoring = 0
  )
}

record_kind <- function(data, roles) {
  event <- !is.na(data[[roles$evid]]) &
    !as.character(data[[roles$evid]]) %in% c("0", "0.0")
  observation <- !event & !is.na(data[[roles$dv]])
  if (!is.null(roles$mdv)) {
    observation <- observation &
      as.character(data[[roles$mdv]]) %in% c("0", "0.0")
  }
  ifelse(event, "dose/event",
         ifelse(observation, "observation", "non-observation"))
}

observed_plot_data <- function(data, roles, dataset,
                               clock = "study_time") {
  observed <- as.character(data[[roles$evid]]) %in% c("0", "0.0")
  if (!is.null(roles$mdv)) {
    observed <- observed &
      as.character(data[[roles$mdv]]) %in% c("0", "0.0")
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
      as.character(data[[roles$dvid]][observed]),
    stringsAsFactors = FALSE
  )
}

demo_design_summary <- function(data, roles, dataset, time_bounds,
                                clock = "study_time") {
  plotted <- observed_plot_data(data, roles, dataset, clock)
  if (!identical(clock, "tad")) {
    plotted <- plotted[
      plotted$time >= time_bounds[1L] & plotted$time <= time_bounds[2L],
      , drop = FALSE
    ]
  }
  cohort_ids <- as.character(unique(data[[roles$id]]))
  pieces <- lapply(sort(unique(plotted$endpoint)), function(endpoint) {
    rows <- plotted[plotted$endpoint == endpoint, , drop = FALSE]
    counts <- table(factor(rows$subject, levels = cohort_ids))
    data.frame(
      dataset = dataset,
      endpoint = endpoint,
      patients = length(cohort_ids),
      patients_with_endpoint = sum(counts > 0),
      observations = nrow(rows),
      mean_time_points_per_patient = mean(counts),
      median_time_points_per_patient = stats::median(counts),
      first_time = min(rows$time),
      last_time = max(rows$time),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, pieces)
}

check_demo_similarity <- function(source, synthetic, roles, time_bounds,
                                  label, clock = "study_time") {
  source_subjects <- length(unique(source[[roles$id]]))
  synthetic_subjects <- length(unique(synthetic[[roles$id]]))
  if (source_subjects != synthetic_subjects) {
    stop(label, ": source and synthetic patient counts differ (",
         source_subjects, " versus ", synthetic_subjects, ").")
  }
  source_summary <- demo_design_summary(
    source, roles, "Source", time_bounds, clock
  )
  synthetic_summary <- demo_design_summary(
    synthetic, roles, "Synthetic", time_bounds, clock
  )
  if (!setequal(source_summary$endpoint, synthetic_summary$endpoint)) {
    stop(label, ": source and synthetic endpoint sets differ.")
  }
  paired <- merge(
    source_summary, synthetic_summary, by = "endpoint",
    suffixes = c("_source", "_synthetic")
  )
  point_difference <- abs(
    paired$mean_time_points_per_patient_synthetic -
      paired$mean_time_points_per_patient_source
  )
  point_allowance <- pmax(
    1, 0.25 * paired$mean_time_points_per_patient_source
  )
  if (any(point_difference > point_allowance)) {
    stop(label, ": mean time points per patient differ materially for ",
         paste(paired$endpoint[point_difference > point_allowance],
               collapse = ", "), ".")
  }
  coverage_allowance <- 0.20 * diff(time_bounds)
  bad_coverage <-
    abs(paired$first_time_synthetic - paired$first_time_source) >
      coverage_allowance |
    abs(paired$last_time_synthetic - paired$last_time_source) >
      coverage_allowance
  if (any(bad_coverage)) {
    stop(label, ": source and synthetic time coverage differs materially for ",
         paste(paired$endpoint[bad_coverage], collapse = ", "), ".")
  }
  summary <- rbind(source_summary, synthetic_summary)
  summary$dataset <- factor(
    summary$dataset, levels = c("Source", "Synthetic")
  )
  summary[order(summary$dataset, summary$endpoint), , drop = FALSE]
}

overlay_plot <- function(source, synthetic, roles, name, clock = "study_time",
                         log_y = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  plot_data <- rbind(
    observed_plot_data(source, roles, "Source", clock),
    observed_plot_data(synthetic, roles, "Synthetic", clock)
  )
  grouping <- if (identical(clock, "tad")) {
    interaction(
      plot_data$dataset, plot_data$subject, plot_data$endpoint,
      plot_data$occasion
    )
  } else {
    interaction(plot_data$dataset, plot_data$subject, plot_data$endpoint)
  }
  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(time, dv, colour = dataset, group = grouping)
  ) +
    ggplot2::geom_line(alpha = 0.35) +
    ggplot2::geom_point(alpha = 0.65, size = 0.9) +
    ggplot2::facet_grid(dataset ~ endpoint, scales = "free_y") +
    ggplot2::scale_colour_manual(
      values = c(Source = "#1B6CA8", Synthetic = "#D95F02")
    ) +
    ggplot2::labs(
      x = if (identical(clock, "tad")) "Time after dose" else "Study time",
      y = if (log_y) "DV (log scale)" else "DV", colour = "Dataset",
      title = paste(name, "source and synthetic observations")
    ) +
    ggplot2::theme_minimal()
  if (log_y) {
    plot <- plot + if (requireNamespace("xgxr", quietly = TRUE)) {
      xgxr::xgx_scale_y_log10()
    } else {
      ggplot2::scale_y_log10()
    }
  }
  plot
}

run_public_demo <- function(name, roles, endpoints, bounds, design, limits,
                            seed, clock = "study_time",
                            comparison_bounds = NULL, log_y = FALSE) {
  source <- load_dataset(name)
  fit_bounds <- bounds(source)
  synthetic <- synpmx::synpmx_empirical(
    data = source, roles = roles, endpoints = endpoints,
    epsilon = 5, delta = 0, bounds = fit_bounds,
    public_design = design(source), contribution_limits = limits,
    budget_allocation = demo_budget(),
    seed = seed,
    backend = "public", public_source = TRUE
  )
  validation <- synpmx::validate_pmx(synthetic, roles, endpoints)
  comparison <- synpmx::compare_pmx(source, synthetic, roles, endpoints)
  overlay <- overlay_plot(source, synthetic, roles, name, clock, log_y = log_y)
  if (is.null(comparison_bounds)) comparison_bounds <- fit_bounds$time
  design_checks <- check_demo_similarity(
    source, synthetic, roles, comparison_bounds, name, clock
  )
  stopifnot(validation$valid)

  message("\n", name, ": ", nrow(synthetic), " generated rows; ",
          length(unique(synthetic[[roles$id]])), " generated subjects")
  message("Source data (first six rows):")
  print(utils::head(source))
  message("Synthetic data (first six rows):")
  print(utils::head(synthetic))
  source_kind <- record_kind(source, roles)
  synthetic_kind <- record_kind(synthetic, roles)
  message("Record counts (dose/event rows are not samples):")
  print(rbind(Source = table(source_kind), Synthetic = table(synthetic_kind)))
  message("Cohort and sampling-design checks by endpoint:")
  print(design_checks, row.names = FALSE)
  print(synpmx::privacy_report(synthetic))
  print(validation)
  print(comparison$release_status)
  if (interactive() && !is.null(overlay)) {
    print(overlay)
  }
  invisible(list(
    source = source,
    source_observations = source[source_kind == "observation", , drop = FALSE],
    synthetic = synthetic,
    synthetic_observations = synthetic[synthetic_kind == "observation", , drop = FALSE],
    overlay_plot = overlay, design_checks = design_checks,
    comparison = comparison, comparison_clock = clock
  ))
}

message("Production backend status:")
print(synpmx::dp_backend_status())
message(paste(
  "These five sources are public, so the demonstrations use the guarded",
  "public fixture backend and make no DP claim. For confidential fitting,",
  "install OpenDP and use the default backend."
))

theo_roles <- synpmx::pmx_roles(
  id = "ID", time = "TIME", dv = "DV", amt = "AMT",
  evid = "EVID", cmt = "CMT", covariates = "WT"
)
theo_endpoints <- list(cp = synpmx::pmx_endpoint(
  alignment = "dose_relative", transform = "log", shape = "occasion",
  cmt = 2
))
theophylline <- run_public_demo(
  "theo_md", theo_roles, theo_endpoints,
  bounds = function(source) synpmx::pmx_bounds(
    c(0, 170), list(cp = c(0, 30)), amt = c(0, 500),
    covariates = list(WT = c(40, 130))
  ),
  design = function(source) synpmx::pmx_public_design(
    synpmx::pmx_schema(source), dose_evid = 101, dose_cmt = 1
  ),
  limits = synpmx::pmx_contribution_limits(40, 8, 8, 30, 11),
  seed = 101
)
print(synpmx::sampling_summary(theophylline$synthetic))
message(paste(
  "Theophylline interpretation: the fit infers seven Q24H dose events and",
  "sampling concentrated around dose occasions 1, 2, and 7."
))

warfarin_roles <- synpmx::pmx_roles(
  id = "id", time = "time", dv = "dv", amt = "amt",
  evid = "evid", dvid = "dvid", covariates = c("wt", "age", "sex")
)
warfarin_endpoints <- list(
  cp = synpmx::pmx_endpoint(
    "cp", "dose_relative", "log", "occasion"
  ),
  pca = synpmx::pmx_endpoint(
    "pca", "study_time", "identity", "global"
  )
)
warfarin_result <- run_public_demo(
  "warfarin", warfarin_roles, warfarin_endpoints,
  bounds = function(source) synpmx::pmx_bounds(
    c(0, 144), list(cp = c(0, 25), pca = c(0, 120)),
    amt = c(0, 200),
    covariates = list(wt = c(40, 150), age = c(18, 100))
  ),
  design = function(source) synpmx::pmx_public_design(
    synpmx::pmx_schema(source), dose_evid = 1
  ),
  limits = synpmx::pmx_contribution_limits(
    30, 2, 2, c(cp = 20, pca = 12), 12
  ),
  seed = 202
)

wbc_roles <- synpmx::pmx_roles(
  id = "ID", time = "TIME", dv = "DV", amt = "AMT",
  evid = "EVID", cmt = "CMT", rate = "RATE",
  covariates = c("V2I", "V1I", "CLI")
)
wbc_endpoints <- list(wbc = synpmx::pmx_endpoint(
  alignment = "study_time", transform = "log", shape = "global", cmt = 3
))
wbc <- run_public_demo(
  "wbcSim", wbc_roles, wbc_endpoints,
  bounds = function(source) synpmx::pmx_bounds(
    c(0, 720), list(wbc = c(0, 30)), amt = c(-200, 200),
    rate = c(-200, 200),
    covariates = list(V2I = c(100, 1500), V1I = c(100, 1200),
                      CLI = c(100, 800))
  ),
  design = function(source) synpmx::pmx_public_design(
    synpmx::pmx_schema(source), dose_evid = 10101, dose_cmt = 1,
    endpoint_cmt = list(wbc = 3)
  ),
  limits = synpmx::pmx_contribution_limits(20, 2, 2, 12, 9),
  seed = 303
)

# NimoData uses DOS as a subject-level treatment-group property. The property
# conditions the inferred regimen, so generated 50, 100, 200, and 400 mg
# groups retain matching event amounts. WGT is deliberately excluded because
# it varies longitudinally in this source and is not a baseline covariate.
nimo_roles <- synpmx::pmx_roles(
  id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
  rate = "RATE", mdv = "MDV", tad = "TAD", occasion = "OCC",
  covariates = c("BSA", "AGE", "HGT"),
  subject_properties = "DOS", exclude = "WGT"
)
nimo_endpoints <- list(cp = synpmx::pmx_endpoint(
  alignment = "dose_relative", transform = "identity", shape = "occasion"
))
nimo <- run_public_demo(
  "nimoData", nimo_roles, nimo_endpoints,
  bounds = function(source) synpmx::pmx_bounds(
    c(0, 3000), list(cp = c(-1, 10)), amt = c(0, 500),
    rate = c(-1200, 1200),
    covariates = list(
      BSA = c(1, 2.5), AGE = c(18, 100), HGT = c(120, 210)
    )
  ),
  design = function(source) synpmx::pmx_public_design(
    synpmx::pmx_schema(source, exclude = "WGT"), dose_evid = 1,
    category_levels = list(DOS = c(50, 100, 200, 400))
  ),
  limits = synpmx::pmx_contribution_limits(60, 10, 10, 8, 12),
  seed = 404, clock = "tad", comparison_bounds = c(0, 3000)
)
message("NimoData released treatment groups and conditioned regimens:")
print(synpmx::subject_property_summary(nimo$synthetic), row.names = FALSE)

# Mavoglurant is a crossover dataset whose TIME clock resets within OCC. DOSE
# is not sampled as an independent covariate; it is regenerated from the
# positive AMT event for each generated subject/occasion.
mav_roles <- synpmx::pmx_roles(
  id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
  cmt = "CMT", rate = "RATE", mdv = "MDV", occasion = "OCC",
  assigned_dose = "DOSE", covariates = c("AGE", "SEX", "WT", "HT")
)
mav_endpoints <- list(cp = synpmx::pmx_endpoint(
  alignment = "dose_relative", transform = "log", shape = "occasion",
  cmt = 2
))
mavoglurant_result <- run_public_demo(
  "mavoglurant", mav_roles, mav_endpoints,
  bounds = function(source) synpmx::pmx_bounds(
    c(0, 120), list(cp = c(0, 2000)), amt = c(0, 60),
    rate = c(-350, 350),
    covariates = list(
      AGE = c(18, 90), WT = c(40, 150), HT = c(1.4, 2.1)
    )
  ),
  design = function(source) synpmx::pmx_public_design(
    synpmx::pmx_schema(source), dose_evid = 1, dose_cmt = 1,
    endpoint_cmt = list(cp = 2), category_levels = list(SEX = c(1, 2))
  ),
  limits = synpmx::pmx_contribution_limits(30, 2, 2, 15, 15),
  seed = 505, clock = "tad", comparison_bounds = c(0, 120), log_y = TRUE
)
message("Mavoglurant assigned-dose coherence by occasion:")
with(mavoglurant_result$synthetic, print(stats::aggregate(
  cbind(DOSE, AMT) ~ OCC,
  data = mavoglurant_result$synthetic[
    EVID != 0 & AMT > 0, , drop = FALSE
  ],
  FUN = function(value) paste(sort(unique(round(value, 4))), collapse = ", ")
), row.names = FALSE))
