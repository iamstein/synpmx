# Multi-seed simulation evaluation for public and fully simulated datasets.
#
# Example:
# Rscript scripts/evaluate_simulations.R \
#   --datasets=all --seeds=101:110 --backend=public \
#   --output=output/simulation-evaluation

sim_eval_script_args <- commandArgs()
sim_eval_file_arg <- grep("^--file=", sim_eval_script_args, value = TRUE)
sim_eval_root <- if (length(sim_eval_file_arg)) {
  normalizePath(file.path(dirname(sub("^--file=", "", sim_eval_file_arg[1L])), ".."))
} else {
  normalizePath(getwd())
}

if (file.exists(file.path(sim_eval_root, "DESCRIPTION")) &&
    requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(sim_eval_root, quiet = TRUE)
} else if (!requireNamespace("pmxSynthData", quietly = TRUE)) {
  stop(
    "Install pmxSynthData or run this script from the repository with pkgload installed.",
    call. = FALSE
  )
} else {
  library(pmxSynthData)
}

helper <- file.path(
  sim_eval_root, "tests", "testthat", "helper-simulation-evaluation.R"
)
if (!file.exists(helper)) stop("Evaluation helper not found: ", helper, call. = FALSE)
sys.source(helper, envir = .GlobalEnv)

parse_options <- function(arguments) {
  out <- list(
    datasets = "all", seeds = "101:110", backend = "public", epsilon = "5",
    output = file.path("output", "simulation-evaluation")
  )
  for (argument in arguments) {
    if (!grepl("^--[^=]+=", argument)) {
      stop("Arguments must use --name=value syntax: ", argument, call. = FALSE)
    }
    pieces <- strsplit(sub("^--", "", argument), "=", fixed = TRUE)[[1L]]
    name <- pieces[1L]
    if (!name %in% names(out)) stop("Unknown option --", name, call. = FALSE)
    out[[name]] <- paste(pieces[-1L], collapse = "=")
  }
  out
}

parse_seeds <- function(value) {
  if (grepl("^[0-9]+:[0-9]+$", value)) {
    limits <- as.integer(strsplit(value, ":", fixed = TRUE)[[1L]])
    seeds <- seq(limits[1L], limits[2L])
  } else {
    seeds <- suppressWarnings(as.integer(strsplit(value, ",", fixed = TRUE)[[1L]]))
  }
  if (!length(seeds) || anyNA(seeds) || any(seeds < 0L)) {
    stop("`--seeds` must be A:B or comma-separated nonnegative integers.",
         call. = FALSE)
  }
  unique(seeds)
}

options <- parse_options(commandArgs(trailingOnly = TRUE))
seeds <- parse_seeds(options$seeds)
epsilon <- suppressWarnings(as.numeric(options$epsilon))
if (length(epsilon) != 1L || !is.finite(epsilon) || epsilon <= 0) {
  stop("`--epsilon` must be one finite positive number.", call. = FALSE)
}
if (!options$backend %in% c("public", "opendp")) {
  stop("`--backend` must be public or opendp.", call. = FALSE)
}
if (options$backend == "opendp" && !isTRUE(dp_backend_status()$available)) {
  stop("The OpenDP backend is unavailable.", call. = FALSE)
}
opendp_status <- dp_backend_status()
message(
  "OpenDP production adapter available=", opendp_status$available,
  if (!is.na(opendp_status$version)) paste0("; version=", opendp_status$version)
  else ""
)

available <- names(sim_eval_registry(include_optional = TRUE))
datasets <- if (identical(options$datasets, "all")) available else
  strsplit(options$datasets, ",", fixed = TRUE)[[1L]]
unknown <- setdiff(datasets, available)
if (length(unknown)) {
  stop("Unavailable evaluation datasets: ", paste(unknown, collapse = ", "),
       call. = FALSE)
}

output <- if (grepl("^/", options$output)) options$output else
  file.path(sim_eval_root, options$output)
figures <- file.path(output, "figures")
dir.create(figures, recursive = TRUE, showWarnings = FALSE)

metric_rows <- list()
gate_rows <- list()
sampling_rows <- list()
regimen_rows <- list()
property_rows <- list()
first_synthetic <- list()
fitted_models <- list()

for (id in datasets) {
  message("Fitting ", id, " with backend=", options$backend,
          " and epsilon=", epsilon)
  case <- sim_eval_case(id)
  model <- sim_eval_fit(case, epsilon = epsilon, backend = options$backend)
  fitted_models[[id]] <- model
  regimen_rows[[id]] <- data.frame(
    dataset = id, epsilon = epsilon, backend = options$backend,
    n_doses = model$population$event$n_doses,
    dose_interval = model$population$event$dose_interval,
    dose_amount = model$population$event$dose_amount,
    dose_rate = model$population$event$dose_rate,
    infusion_probability = model$population$event$infusion_probability,
    infusion_duration = model$population$event$infusion_duration,
    fitted_subject_count = model$population$private_subject_count,
    stringsAsFactors = FALSE
  )
  fitted_properties <- subject_property_summary(model)
  if (nrow(fitted_properties)) {
    property_names <- model$population$subject_properties$names
    property_rows[[id]] <- data.frame(
      dataset = id,
      stratum = seq_len(nrow(fitted_properties)),
      property_values = apply(
        fitted_properties[, property_names, drop = FALSE], 1L,
        function(value) paste(
          paste(property_names, value, sep = "="), collapse = ";"
        )
      ),
      fitted_properties[
        , setdiff(names(fitted_properties), property_names), drop = FALSE
      ],
      epsilon = epsilon, backend = options$backend,
      stringsAsFactors = FALSE
    )
  }
  fitted_sampling <- sampling_summary(model)
  if (nrow(fitted_sampling)) {
    fitted_sampling$dataset <- id
    fitted_sampling$epsilon <- epsilon
    fitted_sampling$backend <- options$backend
    sampling_rows[[id]] <- fitted_sampling
  }
  cohort_size <- length(unique(case$source[[case$roles$id]]))
  for (seed in seeds) {
    synthetic <- generate_pmx(model, n_subjects = cohort_size, seed = seed)
    if (is.null(first_synthetic[[id]])) first_synthetic[[id]] <- synthetic
    gates <- sim_eval_gate_results(case, model, synthetic)
    gates$seed <- seed
    gates$epsilon <- epsilon
    gates$backend <- options$backend
    gate_rows[[paste(id, seed, sep = "_")]] <- gates
    metric_rows[[paste(id, seed, sep = "_")]] <- sim_eval_metric_rows(
      case, model, synthetic, seed, epsilon, options$backend
    )
  }
}

metrics <- do.call(rbind, metric_rows)
gates <- do.call(rbind, gate_rows)
regimens <- do.call(rbind, regimen_rows)
sampling <- if (length(sampling_rows)) do.call(rbind, sampling_rows) else
  data.frame()
properties <- if (length(property_rows)) do.call(rbind, property_rows) else
  data.frame()
failures <- gates[!gates$pass, , drop = FALSE]
rownames(metrics) <- rownames(gates) <- rownames(regimens) <- NULL
if (nrow(sampling)) rownames(sampling) <- NULL
if (nrow(properties)) rownames(properties) <- NULL
rownames(failures) <- NULL

utils::write.csv(metrics, file.path(output, "metrics-by-seed.csv"), row.names = FALSE)
utils::write.csv(gates, file.path(output, "gate-results.csv"), row.names = FALSE)
utils::write.csv(failures, file.path(output, "failures.csv"), row.names = FALSE)
utils::write.csv(regimens, file.path(output, "regimen-by-fit.csv"), row.names = FALSE)
utils::write.csv(sampling, file.path(output, "sampling-by-fit.csv"), row.names = FALSE)
utils::write.csv(
  properties, file.path(output, "subject-properties-by-fit.csv"),
  row.names = FALSE
)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  for (id in datasets) {
    case <- sim_eval_case(id)
    model <- fitted_models[[id]]
    plot_clock <- case$comparison_clock
    plot_suffix <- if (identical(plot_clock, "tad")) "tad" else "study-time"
    plot <- sim_eval_plot(case, first_synthetic[[id]], clock = plot_clock)
    ggplot2::ggsave(
      file.path(figures, paste0(id, "-", plot_suffix, ".png")), plot,
      width = max(8, 4.5 * length(case$endpoints)), height = 7, dpi = 140
    )
    if (id == "theo_md") {
      tad_plot <- sim_eval_plot(case, first_synthetic[[id]], clock = "tad")
      ggplot2::ggsave(
        file.path(figures, "theo_md-tad.png"), tad_plot,
        width = 9, height = 7, dpi = 140
      )
    }
    fitted_sampling <- sampling_summary(model)
    if (nrow(fitted_sampling)) {
      sampling_plot_data <- rbind(
        data.frame(
          endpoint = fitted_sampling$endpoint,
          occasion = fitted_sampling$occasion,
          metric = "Sampling probability",
          value = fitted_sampling$sampling_probability
        ),
        data.frame(
          endpoint = fitted_sampling$endpoint,
          occasion = fitted_sampling$occasion,
          metric = "Observations if sampled",
          value = fitted_sampling$observations_if_sampled
        )
      )
      sampling_plot <- ggplot2::ggplot(
        sampling_plot_data,
        ggplot2::aes(x = occasion, y = value, colour = endpoint,
                     group = endpoint)
      ) +
        ggplot2::geom_point() +
        ggplot2::facet_wrap(ggplot2::vars(metric), scales = "free_y") +
        ggplot2::labs(
          x = "Dose occasion", y = NULL, colour = "Endpoint",
          title = paste(id, "released sampling design")
        ) +
        ggplot2::theme_minimal()
      if (any(table(sampling_plot_data$metric,
                    sampling_plot_data$endpoint) > 1L)) {
        sampling_plot <- sampling_plot + ggplot2::geom_line()
      }
      ggplot2::ggsave(
        file.path(figures, paste0(id, "-sampling.png")), sampling_plot,
        width = 9, height = 4.5, dpi = 140
      )
    }
    curve_rows <- lapply(names(case$endpoints), function(endpoint_name) {
      endpoint <- model$public$endpoints[[endpoint_name]]
      trajectory <- model$population$trajectories[[endpoint_name]]
      working <- trajectory$mean_working
      if (length(working) >= 3L) {
        working[2:(length(working) - 1L)] <-
          (trajectory$mean_working[1:(length(working) - 2L)] +
             2 * trajectory$mean_working[2:(length(working) - 1L)] +
             trajectory$mean_working[3:length(working)]) / 4
      }
      dv <- if (identical(endpoint$transform_resolved, "log")) {
        exp(working) - endpoint$offset
      } else working
      data.frame(
        endpoint = endpoint_name, clock = trajectory$grid,
        released_mean = pmin(pmax(dv, endpoint$bound[1L]), endpoint$bound[2L])
      )
    })
    curve_data <- do.call(rbind, curve_rows)
    curve_plot <- ggplot2::ggplot(
      curve_data,
      ggplot2::aes(x = clock, y = released_mean, colour = endpoint)
    ) +
      ggplot2::geom_line(linewidth = .8) + ggplot2::geom_point() +
      ggplot2::facet_wrap(ggplot2::vars(endpoint), scales = "free") +
      ggplot2::labs(
        x = "Declared endpoint clock", y = "Released mean DV",
        colour = "Endpoint", title = paste(id, "released trajectory curves")
      ) +
      ggplot2::theme_minimal()
    ggplot2::ggsave(
      file.path(figures, paste0(id, "-released-curves.png")), curve_plot,
      width = max(7, 4 * length(case$endpoints)), height = 4.5, dpi = 140
    )
  }
}

git_value <- function(arguments) {
  tryCatch(
    paste(system2("git", c("-C", sim_eval_root, arguments), stdout = TRUE,
                  stderr = FALSE), collapse = "\n"),
    error = function(e) NA_character_
  )
}
dependencies <- c("ggplot2", "nlmixr2data", "opendp")
dependency_versions <- vapply(dependencies, function(package) {
  if (requireNamespace(package, quietly = TRUE)) {
    as.character(utils::packageVersion(package))
  } else "not installed"
}, character(1))
manifest_value <- function(value) {
  paste(capture.output(dput(unclass(value))), collapse = " ")
}
dataset_manifest <- unlist(lapply(datasets, function(id) {
  case <- sim_eval_case(id)
  c(
    paste0("dataset_", id, "_bounds: ", manifest_value(case$bounds)),
    paste0(
      "dataset_", id, "_contribution_limits: ",
      manifest_value(case$limits)
    ),
    paste0(
      "dataset_", id, "_budget_allocation: ",
      manifest_value(case$budget)
    )
  )
}), use.names = FALSE)
manifest <- c(
  paste("timestamp:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("package_version:", as.character(utils::packageVersion("pmxSynthData"))),
  paste("git_commit:", git_value(c("rev-parse", "HEAD"))),
  paste("dirty_worktree:", nzchar(git_value(c("status", "--porcelain")))),
  paste("R_version:", R.version.string),
  paste("platform:", R.version$platform),
  paste("backend:", options$backend),
  paste("opendp_available:", opendp_status$available),
  paste("opendp_version:", ifelse(
    is.na(opendp_status$version), "not installed", opendp_status$version
  )),
  paste("opendp_production:", opendp_status$production),
  paste("epsilon:", epsilon),
  paste("datasets:", paste(datasets, collapse = ",")),
  paste("seeds:", paste(seeds, collapse = ",")),
  paste0("dependency_", dependencies, ": ", dependency_versions),
  dataset_manifest
)
writeLines(manifest, file.path(output, "run-manifest.txt"))

html_escape <- function(x) {
  x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}
html_table <- function(data) {
  if (!nrow(data)) return("<p>None.</p>")
  header <- paste0("<th>", html_escape(names(data)), "</th>", collapse = "")
  rows <- apply(data, 1L, function(row) {
    paste0("<tr>", paste0("<td>", html_escape(row), "</td>", collapse = ""),
           "</tr>")
  })
  paste0("<table><thead><tr>", header, "</tr></thead><tbody>",
         paste(rows, collapse = "\n"), "</tbody></table>")
}
gate_summary <- stats::aggregate(
  pass ~ dataset, gates, function(x) paste0(sum(x), "/", length(x), " passed")
)
image_html <- character()
for (id in datasets) {
  case <- sim_eval_case(id)
  suffix <- if (identical(case$comparison_clock, "tad")) {
    "tad"
  } else "study-time"
  path <- paste0("figures/", id, "-", suffix, ".png")
  if (file.exists(file.path(output, path))) {
    image_html <- c(image_html, paste0(
      "<h2>", html_escape(id), "</h2><img src=\"", path,
      "\" alt=\"", html_escape(id), " comparison\">"
    ))
  }
  for (suffix in c("sampling", "released-curves")) {
    extra_path <- paste0("figures/", id, "-", suffix, ".png")
    if (file.exists(file.path(output, extra_path))) {
      image_html <- c(image_html, paste0(
        "<img src=\"", extra_path, "\" alt=\"", html_escape(id), " ",
        suffix, "\">"
      ))
    }
  }
}
summary_html <- paste0(
  "<!doctype html><html><head><meta charset=\"utf-8\"><title>",
  "pmxSynthData simulation evaluation</title><style>",
  "body{font-family:sans-serif;max-width:1200px;margin:auto;padding:2rem}",
  "table{border-collapse:collapse}th,td{border:1px solid #ccc;padding:.35rem}",
  "img{max-width:100%;height:auto}</style></head><body>",
  "<h1>pmxSynthData simulation evaluation</h1>",
  "<p>Backend: ", html_escape(options$backend), "; epsilon: ", epsilon,
  "; seeds: ", html_escape(paste(seeds, collapse = ", ")), ".</p>",
  "<h2>Gate summary</h2>", html_table(gate_summary),
  "<h2>Failures</h2>", html_table(failures),
  paste(image_html, collapse = "\n"), "</body></html>"
)
writeLines(summary_html, file.path(output, "summary.html"))

message("Evaluation outputs written to ", output)
message(nrow(gates) - nrow(failures), " gates passed; ",
        nrow(failures), " failed.")
if (nrow(failures)) quit(save = "no", status = 1L)
