# Version 2 end-to-end demonstrations for the public nlmixr2data datasets.
# Install pmxSynthData first with: R CMD INSTALL .

if (!requireNamespace("pmxSynthData", quietly = TRUE)) {
  stop("Install pmxSynthData before running this script: R CMD INSTALL .")
}
required_api <- c("fit_private_pmx", "sampling_summary")
missing_api <- setdiff(required_api, getNamespaceExports("pmxSynthData"))
if (length(missing_api)) {
  stop(
    "The installed pmxSynthData is older than this demo. Reinstall the ",
    "package from this repository with `R CMD INSTALL .` before running it."
  )
}
if (!requireNamespace("nlmixr2data", quietly = TRUE)) {
  stop("Install nlmixr2data to run the integration demonstrations.")
}

load_dataset <- function(name) {
  environment <- new.env(parent = emptyenv())
  utils::data(list = name, package = "nlmixr2data", envir = environment)
  get(name, envir = environment, inherits = FALSE)
}

demo_budget <- function() {
  pmxSynthData::pmx_budget_allocation(
    subject_count = .10, event = .15, timing = .15,
    covariates = .10, endpoints = .50, censoring = 0
  )
}

observed_plot_data <- function(data, roles, dataset) {
  observed <- as.character(data[[roles$evid]]) %in% c("0", "0.0") &
    !is.na(data[[roles$dv]])
  observation_rows <- which(observed)
  occasion <- rep(1L, length(observation_rows))
  subject_values <- data[[roles$id]]
  for (id in unique(subject_values[observation_rows])) {
    subject_rows <- which(!is.na(subject_values) & subject_values == id)
    events <- !(as.character(data[[roles$evid]][subject_rows]) %in%
                  c("0", "0.0"))
    if (!is.null(roles$amt)) {
      events <- events & as.numeric(data[[roles$amt]][subject_rows]) > 0
    }
    dose_times <- sort(unique(as.numeric(data[[roles$time]][
      subject_rows[events]
    ])))
    positions <- which(subject_values[observation_rows] == id)
    if (length(dose_times)) {
      occasion[positions] <- pmax(1L, findInterval(
        as.numeric(data[[roles$time]][observation_rows[positions]]),
        dose_times
      ))
    }
  }
  data.frame(
    dataset = factor(dataset, levels = c("Source", "Synthetic")),
    subject = as.character(data[[roles$id]][observation_rows]),
    time = as.numeric(data[[roles$time]][observation_rows]),
    dv = as.numeric(data[[roles$dv]][observation_rows]),
    occasion = occasion,
    endpoint = if (is.null(roles$dvid)) "DV" else
      as.character(data[[roles$dvid]][observed]),
    stringsAsFactors = FALSE
  )
}

demo_design_summary <- function(data, roles, dataset, time_bounds) {
  plotted <- observed_plot_data(data, roles, dataset)
  plotted <- plotted[
    plotted$time >= time_bounds[1L] & plotted$time <= time_bounds[2L],
    , drop = FALSE
  ]
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
                                  label) {
  source_subjects <- length(unique(source[[roles$id]]))
  synthetic_subjects <- length(unique(synthetic[[roles$id]]))
  if (source_subjects != synthetic_subjects) {
    stop(label, ": source and synthetic patient counts differ (",
         source_subjects, " versus ", synthetic_subjects, ").")
  }
  source_summary <- demo_design_summary(
    source, roles, "Source", time_bounds
  )
  synthetic_summary <- demo_design_summary(
    synthetic, roles, "Synthetic", time_bounds
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

overlay_plot <- function(source, mock, roles, name) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  plot_data <- rbind(
    observed_plot_data(source, roles, "Source"),
    observed_plot_data(mock, roles, "Synthetic")
  )
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(time, dv, colour = dataset,
                 group = interaction(dataset, subject))
  ) +
    ggplot2::geom_line(alpha = 0.35) +
    ggplot2::geom_point(alpha = 0.65, size = 0.9) +
    ggplot2::facet_wrap(
      ggplot2::vars(dataset, endpoint),
      ncol = max(1L, length(unique(plot_data$endpoint))),
      scales = "free_y"
    ) +
    ggplot2::scale_colour_manual(
      values = c(Source = "#1B6CA8", Synthetic = "#D95F02")
    ) +
    ggplot2::labs(
      x = "Study time", y = "DV", colour = "Dataset",
      title = paste(name, "source and synthetic observations")
    ) +
    ggplot2::theme_minimal()
}

run_public_demo <- function(name, roles, endpoints, bounds, design, limits,
                            seed) {
  source <- load_dataset(name)
  fit_bounds <- bounds(source)
  model <- pmxSynthData::fit_private_pmx(
    data = source, roles = roles, endpoints = endpoints,
    epsilon = 5, delta = 0, bounds = fit_bounds,
    public_design = design(source), contribution_limits = limits,
    budget_allocation = demo_budget(),
    backend = "public", public_source = TRUE
  )
  mock <- pmxSynthData::generate_pmx(model, seed = seed)
  validation <- pmxSynthData::validate_pmx(mock, roles, endpoints)
  comparison <- pmxSynthData::compare_pmx(source, mock, roles, endpoints)
  overlay <- overlay_plot(source, mock, roles, name)
  design_checks <- check_demo_similarity(
    source, mock, roles, fit_bounds$time, name
  )
  stopifnot(validation$valid)

  message("\n", name, ": ", nrow(mock), " generated rows; ",
          length(unique(mock[[roles$id]])), " generated subjects")
  message("Source data (first six rows):")
  print(utils::head(source))
  message("Mock data (first six rows):")
  print(utils::head(mock))
  source_kind <- ifelse(as.character(source[[roles$evid]]) %in% c("0", "0.0"),
                        "observation", "dose/event")
  mock_kind <- ifelse(as.character(mock[[roles$evid]]) %in% c("0", "0.0"),
                      "observation", "dose/event")
  message("Record counts (dose/event rows are not samples):")
  print(rbind(Source = table(source_kind), Mock = table(mock_kind)))
  message("Cohort and sampling-design checks by endpoint:")
  print(design_checks, row.names = FALSE)
  print(pmxSynthData::privacy_report(model))
  print(validation)
  print(comparison$release_status)
  if (interactive() && !is.null(overlay)) {
    print(overlay)
  }
  invisible(list(
    source = source,
    source_observations = source[source_kind == "observation", , drop = FALSE],
    model = model, mock = mock,
    mock_observations = mock[mock_kind == "observation", , drop = FALSE],
    overlay_plot = overlay, design_checks = design_checks,
    comparison = comparison
  ))
}

message("Production backend status:")
print(pmxSynthData::dp_backend_status())
message(paste(
  "These three sources are public, so the demonstrations use the guarded",
  "public fixture backend and make no DP claim. For confidential fitting,",
  "install OpenDP and use the default backend."
))

theo_roles <- pmxSynthData::pmx_roles(
  id = "ID", time = "TIME", dv = "DV", amt = "AMT",
  evid = "EVID", cmt = "CMT", covariates = "WT"
)
theo_endpoints <- list(cp = pmxSynthData::pmx_endpoint(
  alignment = "dose_relative", transform = "log", shape = "occasion",
  grid = c(0, .25, .5, 1, 2, 4, 5, 7, 9, 12, 23.75), cmt = 2
))
theophylline <- run_public_demo(
  "theo_md", theo_roles, theo_endpoints,
  bounds = function(source) pmxSynthData::pmx_bounds(
    c(0, 170), list(cp = c(0, 30)), amt = c(0, 500),
    covariates = list(WT = c(40, 130))
  ),
  design = function(source) pmxSynthData::pmx_public_design(
    pmxSynthData::pmx_schema(source), dose_evid = 101, dose_cmt = 1
  ),
  limits = pmxSynthData::pmx_contribution_limits(40, 8, 8, 30, 11),
  seed = 101
)
print(pmxSynthData::sampling_summary(theophylline$model))
message(paste(
  "Theophylline interpretation: the fit infers seven Q24H dose events and",
  "sampling concentrated around dose occasions 1, 2, and 7."
))

warfarin_roles <- pmxSynthData::pmx_roles(
  id = "id", time = "time", dv = "dv", amt = "amt",
  evid = "evid", dvid = "dvid", covariates = c("wt", "age", "sex")
)
warfarin_endpoints <- list(
  cp = pmxSynthData::pmx_endpoint(
    "cp", "dose_relative", "log", "occasion",
    grid = c(.5, 1, 1.5, 2, 3, 6, 9, 12, 24, 48, 72, 120)
  ),
  pca = pmxSynthData::pmx_endpoint(
    "pca", "study_time", "identity", "global",
    grid = c(0, 24, 36, 48, 72, 96, 120, 144)
  )
)
warfarin_result <- run_public_demo(
  "warfarin", warfarin_roles, warfarin_endpoints,
  bounds = function(source) pmxSynthData::pmx_bounds(
    c(0, 144), list(cp = c(0, 25), pca = c(0, 120)),
    amt = c(0, 200),
    covariates = list(wt = c(40, 150), age = c(18, 100))
  ),
  design = function(source) pmxSynthData::pmx_public_design(
    pmxSynthData::pmx_schema(source), dose_times = 0, n_doses = 1,
    dose_amount = 100, dose_evid = 1,
    endpoint_grids = lapply(warfarin_endpoints, `[[`, "grid")
  ),
  limits = pmxSynthData::pmx_contribution_limits(
    30, 2, 2, c(cp = 20, pca = 12), 12
  ),
  seed = 202
)

wbc_roles <- pmxSynthData::pmx_roles(
  id = "ID", time = "TIME", dv = "DV", amt = "AMT",
  evid = "EVID", cmt = "CMT", rate = "RATE",
  covariates = c("V2I", "V1I", "CLI")
)
wbc_endpoints <- list(wbc = pmxSynthData::pmx_endpoint(
  alignment = "study_time", transform = "log", shape = "global",
  grid = c(0, 72, 120, 168, 192, 240, 336, 504, 672), cmt = 3
))
wbc <- run_public_demo(
  "wbcSim", wbc_roles, wbc_endpoints,
  bounds = function(source) pmxSynthData::pmx_bounds(
    c(0, 720), list(wbc = c(0, 30)), amt = c(-200, 200),
    rate = c(-200, 200),
    covariates = list(V2I = c(100, 1500), V1I = c(100, 1200),
                      CLI = c(100, 800))
  ),
  design = function(source) pmxSynthData::pmx_public_design(
    pmxSynthData::pmx_schema(source), dose_times = 0, n_doses = 1,
    dose_amount = 120, dose_rate = 120, infusion_duration = 3,
    dose_evid = 10101, dose_cmt = 1,
    endpoint_grids = list(wbc = wbc_endpoints$wbc$grid),
    endpoint_cmt = list(wbc = 3)
  ),
  limits = pmxSynthData::pmx_contribution_limits(20, 2, 2, 12, 9),
  seed = 303
)
