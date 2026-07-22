sim_eval_budget <- function(censoring = 0) {
  pmx_budget_allocation(
    subject_count = 0.10, event = 0.15, timing = 0.15,
    covariates = 0.10, endpoints = 0.50 - censoring,
    censoring = censoring
  )
}

sim_eval_load_nlmixr2 <- function(name) {
  environment <- new.env(parent = emptyenv())
  utils::data(list = name, package = "nlmixr2data", envir = environment)
  get(name, envir = environment, inherits = FALSE)
}

sim_eval_censoring_source <- function() {
  one <- pmx_censoring_fixture()
  out <- do.call(rbind, lapply(seq_len(8L), function(id) {
    subject <- one
    subject$ID <- as.integer(id)
    selected <- subject$CENS == 0 & subject$EVID == 0
    subject$DV[selected] <- subject$DV[selected] + id / 20
    subject
  }))
  rownames(out) <- NULL
  out
}

sim_eval_case <- function(id) {
  if (!id %in% c(
    "censoring", "theo_md", "warfarin", "wbcSim", "nimoData",
    "mavoglurant"
  )) {
    stop("Unknown simulation-evaluation dataset `", id, "`.", call. = FALSE)
  }
  if (id != "censoring" && !requireNamespace("nlmixr2data", quietly = TRUE)) {
    stop("Package `nlmixr2data` is required for `", id, "`.", call. = FALSE)
  }

  if (id == "censoring") {
    source <- sim_eval_censoring_source()
    roles <- pmx_roles(
      id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
      dvid = "DVID", mdv = "MDV", cens = "CENS", limit = "LIMIT"
    )
    endpoints <- list(marker = pmx_endpoint(
      dvid = "marker", alignment = "study_time", transform = "identity",
      shape = "global",
      censoring = list(left = 1, right = 10, interval = c(2, 4))
    ))
    bounds <- pmx_bounds(
      c(0, 5), list(marker = c(0, 12)), amt = c(0, 200),
      limit = list(marker = c(0, 12))
    )
    design <- pmx_public_design(pmx_schema(source), dose_evid = 1)
    limits <- pmx_contribution_limits(8, 2, 2, 5, 4)
    budget <- pmx_budget_allocation(.10, .10, .15, 0, .45, .20)
  } else if (id == "theo_md") {
    source <- sim_eval_load_nlmixr2("theo_md")
    roles <- pmx_roles(
      id = "ID", time = "TIME", dv = "DV", amt = "AMT",
      evid = "EVID", cmt = "CMT", covariates = "WT"
    )
    endpoints <- list(cp = pmx_endpoint(
      alignment = "dose_relative", transform = "log", shape = "occasion",
      cmt = 2
    ))
    bounds <- pmx_bounds(
      c(0, 170), list(cp = c(0, 30)), amt = c(0, 500),
      covariates = list(WT = c(40, 130))
    )
    design <- pmx_public_design(
      pmx_schema(source), dose_evid = 101, dose_cmt = 1
    )
    limits <- pmx_contribution_limits(40, 8, 8, 30, 11)
    budget <- sim_eval_budget()
  } else if (id == "warfarin") {
    source <- sim_eval_load_nlmixr2("warfarin")
    roles <- pmx_roles(
      id = "id", time = "time", dv = "dv", amt = "amt",
      evid = "evid", dvid = "dvid", covariates = c("wt", "age", "sex")
    )
    endpoints <- list(
      cp = pmx_endpoint("cp", "dose_relative", "log", "occasion"),
      pca = pmx_endpoint("pca", "study_time", "identity", "global")
    )
    bounds <- pmx_bounds(
      c(0, 144), list(cp = c(0, 25), pca = c(0, 120)),
      amt = c(0, 200),
      covariates = list(wt = c(40, 150), age = c(18, 100))
    )
    design <- pmx_public_design(pmx_schema(source), dose_evid = 1)
    limits <- pmx_contribution_limits(
      30, 2, 2, c(cp = 20, pca = 12), 12
    )
    budget <- sim_eval_budget()
  } else if (id == "wbcSim") {
    source <- sim_eval_load_nlmixr2("wbcSim")
    roles <- pmx_roles(
      id = "ID", time = "TIME", dv = "DV", amt = "AMT",
      evid = "EVID", cmt = "CMT", rate = "RATE",
      covariates = c("V2I", "V1I", "CLI")
    )
    endpoints <- list(wbc = pmx_endpoint(
      alignment = "study_time", transform = "log", shape = "global", cmt = 3
    ))
    bounds <- pmx_bounds(
      c(0, 720), list(wbc = c(0, 30)), amt = c(-200, 200),
      rate = c(-200, 200),
      covariates = list(
        V2I = c(100, 1500), V1I = c(100, 1200), CLI = c(100, 800)
      )
    )
    design <- pmx_public_design(
      pmx_schema(source), dose_evid = 10101, dose_cmt = 1,
      endpoint_cmt = list(wbc = 3)
    )
    limits <- pmx_contribution_limits(20, 2, 2, 12, 9)
    budget <- sim_eval_budget()
  } else if (id == "nimoData") {
    source <- sim_eval_load_nlmixr2("nimoData")
    roles <- pmx_roles(
      id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
      rate = "RATE", mdv = "MDV", tad = "TAD", occasion = "OCC",
      covariates = c("BSA", "AGE", "HGT"),
      subject_properties = "DOS", exclude = "WGT"
    )
    endpoints <- list(cp = pmx_endpoint(
      alignment = "dose_relative", transform = "identity", shape = "occasion"
    ))
    bounds <- pmx_bounds(
      c(0, 3000), list(cp = c(-1, 10)), amt = c(0, 500),
      rate = c(-1200, 1200),
      covariates = list(
        BSA = c(1, 2.5), AGE = c(18, 100), HGT = c(120, 210)
      )
    )
    design <- pmx_public_design(
      pmx_schema(source, exclude = "WGT"), dose_evid = 1,
      category_levels = list(DOS = c(50, 100, 200, 400))
    )
    limits <- pmx_contribution_limits(60, 10, 10, 8, 12)
    budget <- sim_eval_budget()
  } else {
    source <- sim_eval_load_nlmixr2("mavoglurant")
    roles <- pmx_roles(
      id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
      cmt = "CMT", rate = "RATE", mdv = "MDV", occasion = "OCC",
      assigned_dose = "DOSE", covariates = c("AGE", "SEX", "WT", "HT")
    )
    endpoints <- list(cp = pmx_endpoint(
      alignment = "dose_relative", transform = "log", shape = "occasion",
      cmt = 2
    ))
    bounds <- pmx_bounds(
      c(0, 120), list(cp = c(0, 2000)), amt = c(0, 60),
      rate = c(-350, 350),
      covariates = list(
        AGE = c(18, 90), WT = c(40, 150), HT = c(1.4, 2.1)
      )
    )
    design <- pmx_public_design(
      pmx_schema(source), dose_evid = 1, dose_cmt = 1,
      endpoint_cmt = list(cp = 2), category_levels = list(SEX = c(1, 2))
    )
    limits <- pmx_contribution_limits(30, 2, 2, 15, 15)
    budget <- sim_eval_budget()
  }

  structure(list(
    id = id, source = source, roles = roles, endpoints = endpoints,
    bounds = bounds, design = design, limits = limits, budget = budget,
    comparison_clock = if (id %in% c("nimoData", "mavoglurant")) {
      "tad"
    } else "study_time"
  ), class = "pmx_sim_eval_case")
}

sim_eval_registry <- function(include_optional = TRUE) {
  ids <- "censoring"
  if (include_optional && requireNamespace("nlmixr2data", quietly = TRUE)) {
    ids <- c(
      ids, "theo_md", "warfarin", "wbcSim", "nimoData", "mavoglurant"
    )
  }
  stats::setNames(lapply(ids, sim_eval_case), ids)
}

sim_eval_fit <- function(case, epsilon = 5, backend = "public") {
  suppressWarnings(fit_private_pmx(
    case$source, case$roles, case$endpoints,
    epsilon = epsilon, delta = 0, bounds = case$bounds,
    public_design = case$design,
    contribution_limits = case$limits,
    budget_allocation = case$budget,
    backend = backend, public_source = identical(backend, "public")
  ))
}

sim_eval_observation_rows <- function(data, roles) {
  observed <- as.character(data[[roles$evid]]) %in% c("0", "0.0")
  if (!is.null(roles$mdv)) {
    observed <- observed & as.character(data[[roles$mdv]]) %in% c("0", "0.0")
  }
  observed & !is.na(data[[roles$dv]])
}

sim_eval_event_rows <- function(data, roles) {
  !is.na(data[[roles$evid]]) &
    !as.character(data[[roles$evid]]) %in% c("0", "0.0")
}

sim_eval_endpoint_name <- function(data, rows, roles, endpoints) {
  if (is.null(roles$dvid)) return(rep(names(endpoints)[1L], length(rows)))
  declared <- vapply(endpoints, function(x) as.character(x$dvid), character(1))
  names(declared) <- names(endpoints)
  names(declared)[match(as.character(data[[roles$dvid]][rows]), declared)]
}

sim_eval_plot_data <- function(data, roles, endpoints, dataset,
                               clock = "study_time", time_bounds = NULL) {
  rows <- which(sim_eval_observation_rows(data, roles))
  ids <- data[[roles$id]]
  occasion <- rep(1L, length(rows))
  tad <- rep(NA_real_, length(rows))
  if (!is.null(roles$occasion)) {
    declared <- suppressWarnings(as.integer(data[[roles$occasion]][rows]))
    valid <- !is.na(declared) & declared >= 1L
    occasion[valid] <- declared[valid]
  }
  if (!is.null(roles$tad)) {
    declared <- suppressWarnings(as.numeric(data[[roles$tad]][rows]))
    tad[is.finite(declared)] <- pmax(0, declared[is.finite(declared)])
  }
  for (id in unique(ids[rows])) {
    subject_rows <- which(!is.na(ids) & ids == id)
    event <- sim_eval_event_rows(data[subject_rows, , drop = FALSE], roles)
    if (!is.null(roles$amt)) {
      event <- event & as.numeric(data[[roles$amt]][subject_rows]) > 0
    }
    positions <- which(ids[rows] == id)
    event_rows <- subject_rows[event]
    if (length(event_rows) && !is.null(roles$occasion)) {
      event_occasion <- suppressWarnings(as.integer(
        data[[roles$occasion]][event_rows]
      ))
      for (position in positions) {
        candidates <- event_rows[event_occasion == occasion[position]]
        if (length(candidates) && !is.finite(tad[position])) {
          origin <- min(as.numeric(data[[roles$time]][candidates]))
          tad[position] <- as.numeric(data[[roles$time]][rows[position]]) - origin
        }
      }
    } else if (length(event_rows)) {
      dose_times <- sort(unique(as.numeric(data[[roles$time]][event_rows])))
      occasion[positions] <- pmax(1L, findInterval(
        as.numeric(data[[roles$time]][rows[positions]]), dose_times
      ))
      occasion[positions] <- pmin(occasion[positions], length(dose_times))
      tad[positions] <- as.numeric(data[[roles$time]][rows[positions]]) -
        dose_times[occasion[positions]]
    }
  }
  time <- if (identical(clock, "tad")) tad else
    as.numeric(data[[roles$time]][rows])
  out <- data.frame(
    dataset = factor(dataset, levels = c("Source", "Synthetic")),
    subject = as.character(ids[rows]),
    endpoint = sim_eval_endpoint_name(data, rows, roles, endpoints),
    occasion = occasion, time = time,
    dv = as.numeric(data[[roles$dv]][rows]),
    stringsAsFactors = FALSE
  )
  if (!is.null(time_bounds) && !identical(clock, "tad")) {
    out <- out[
      out$time >= time_bounds[1L] & out$time <= time_bounds[2L],
      , drop = FALSE
    ]
  }
  out
}

sim_eval_design_summary <- function(data, roles, endpoints, dataset,
                                    time_bounds = NULL,
                                    clock = "study_time") {
  plotted <- sim_eval_plot_data(
    data, roles, endpoints, dataset, clock = clock,
    time_bounds = time_bounds
  )
  cohort <- as.character(unique(data[[roles$id]]))
  pieces <- lapply(names(endpoints), function(endpoint) {
    selected <- plotted[plotted$endpoint == endpoint, , drop = FALSE]
    counts <- table(factor(selected$subject, levels = cohort))
    finite_time <- selected$time[is.finite(selected$time)]
    data.frame(
      dataset = dataset, endpoint = endpoint,
      patients = length(cohort), patients_with_endpoint = sum(counts > 0),
      observations = nrow(selected), mean_points = mean(counts),
      median_points = stats::median(counts),
      first_time = if (length(finite_time)) min(finite_time) else NA_real_,
      last_time = if (length(finite_time)) max(finite_time) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, pieces)
}

sim_eval_sampling_summary <- function(data, roles, endpoints) {
  plotted <- sim_eval_plot_data(data, roles, endpoints, "Synthetic")
  if (!nrow(plotted)) return(data.frame())
  cohort <- as.character(unique(data[[roles$id]]))
  pieces <- list()
  for (endpoint in names(endpoints)) {
    selected <- plotted[plotted$endpoint == endpoint, , drop = FALSE]
    occasions <- sort(unique(selected$occasion))
    for (occasion in occasions) {
      rows <- selected[selected$occasion == occasion, , drop = FALSE]
      count <- table(factor(rows$subject, levels = cohort))
      active <- count > 0
      pieces[[length(pieces) + 1L]] <- data.frame(
        endpoint = endpoint, occasion = occasion,
        activation = mean(active),
        conditional_count = if (any(active)) mean(count[active]) else 0,
        expected_count = mean(count), stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, pieces)
}

sim_eval_directional_peaks <- function(data, roles, endpoints) {
  plotted <- sim_eval_plot_data(data, roles, endpoints, "Synthetic", "tad")
  local <- names(endpoints)[vapply(endpoints, function(endpoint) {
    endpoint$alignment %in% c("dose_relative", "occasion")
  }, logical(1))]
  plotted <- plotted[plotted$endpoint %in% local & is.finite(plotted$time), ]
  if (!nrow(plotted)) return(data.frame())
  groups <- split(
    plotted,
    interaction(plotted$subject, plotted$endpoint, plotted$occasion, drop = TRUE)
  )
  do.call(rbind, lapply(groups, function(profile) {
    profile <- profile[order(profile$time), , drop = FALSE]
    direction <- sign(diff(profile$dv))
    direction <- direction[direction != 0]
    data.frame(
      subject = profile$subject[1L], endpoint = profile$endpoint[1L],
      occasion = profile$occasion[1L], observations = nrow(profile),
      directional_peaks = sum(diff(direction) < 0),
      stringsAsFactors = FALSE
    )
  }))
}

sim_eval_timing_vector_copied <- function(source, synthetic, roles) {
  source_vectors <- split(as.numeric(source[[roles$time]]), source[[roles$id]])
  synthetic_vectors <- split(
    as.numeric(synthetic[[roles$time]]), synthetic[[roles$id]]
  )
  any(vapply(synthetic_vectors, function(candidate) {
    any(vapply(source_vectors, identical, logical(1), y = candidate))
  }, logical(1)))
}

sim_eval_public_overrides_absent <- function(design) {
  scalars <- c(
    "dose_times", "dose_interval", "n_doses", "dose_amount", "dose_rate",
    "infusion_duration"
  )
  all(vapply(design[scalars], is.null, logical(1))) &&
    !length(design$endpoint_grids) && !length(design$endpoint_occasion_grids)
}

sim_eval_gate_results <- function(case, model, synthetic) {
  source <- case$source
  roles <- case$roles
  endpoints <- case$endpoints
  retained_source <- source[
    , setdiff(names(source), roles$exclude), drop = FALSE
  ]
  result <- list()
  add <- function(gate, pass, value = "", threshold = "", detail = "") {
    result[[length(result) + 1L]] <<- data.frame(
      dataset = case$id, gate = gate, pass = isTRUE(pass),
      value = as.character(value), threshold = as.character(threshold),
      detail = detail, stringsAsFactors = FALSE
    )
  }

  validation <- validate_pmx(synthetic, roles, endpoints)
  add("valid_pmx", validation$valid, validation$valid, "TRUE")
  source_ids <- unique(as.character(source[[roles$id]]))
  synthetic_ids <- unique(as.character(synthetic[[roles$id]]))
  add("new_ids", !any(synthetic_ids %in% source_ids),
      sum(synthetic_ids %in% source_ids), "0 reused IDs")
  add("cohort_size", length(synthetic_ids) == length(source_ids),
      length(synthetic_ids), length(source_ids))
  add("column_names", identical(names(synthetic), names(retained_source)),
      paste(names(synthetic), collapse = ","),
      paste(names(retained_source), collapse = ","))
  source_class <- vapply(
    retained_source, function(x) paste(class(x), collapse = "/"),
    character(1)
  )
  synthetic_class <- vapply(
    synthetic, function(x) paste(class(x), collapse = "/"), character(1)
  )
  add("column_classes", identical(synthetic_class, source_class),
      paste(synthetic_class, collapse = ","), paste(source_class, collapse = ","))
  add("no_source_timing_vector",
      !sim_eval_timing_vector_copied(source, synthetic, roles), "none", "none")
  add("no_source_derived_design",
      sim_eval_public_overrides_absent(model$public$design), "absent", "absent")

  source_summary <- sim_eval_design_summary(
    source, roles, endpoints, "Source", case$bounds$time,
    case$comparison_clock
  )
  synthetic_summary <- sim_eval_design_summary(
    synthetic, roles, endpoints, "Synthetic", case$bounds$time,
    case$comparison_clock
  )
  paired <- merge(
    source_summary, synthetic_summary, by = "endpoint",
    suffixes = c("_source", "_synthetic")
  )
  add("endpoint_set", nrow(paired) == length(endpoints) &&
        all(paired$patients_with_endpoint_synthetic > 0),
      paste(paired$endpoint, collapse = ","), paste(names(endpoints), collapse = ","))
  point_difference <- abs(paired$mean_points_synthetic - paired$mean_points_source)
  point_allowance <- pmax(1, .25 * paired$mean_points_source)
  add("mean_points_per_subject", all(point_difference <= point_allowance),
      paste(round(point_difference, 3), collapse = ","),
      paste(round(point_allowance, 3), collapse = ","))
  coverage_difference <- pmax(
    abs(paired$first_time_synthetic - paired$first_time_source),
    abs(paired$last_time_synthetic - paired$last_time_source)
  )
  coverage_allowance <- .20 * diff(case$bounds$time)
  add("time_coverage", all(coverage_difference <= coverage_allowance),
      paste(round(coverage_difference, 3), collapse = ","), coverage_allowance)

  if (case$id == "theo_md") {
    add("theo_dose_count", model$population$event$n_doses == 7L,
        model$population$event$n_doses, 7)
    add("theo_dose_interval",
        abs(model$population$event$dose_interval - 24) <= .1,
        round(model$population$event$dose_interval, 3), "24 +/- 0.1")
    occasion_counts <- lapply(split(synthetic, synthetic[[roles$id]]), function(x) {
      doses <- sort(as.numeric(x[[roles$time]][sim_eval_event_rows(x, roles)]))
      observations <- as.numeric(x[[roles$time]][sim_eval_observation_rows(x, roles)])
      table(factor(pmax(1L, findInterval(observations, doses)), levels = 1:7))
    })
    inactive <- all(vapply(occasion_counts, function(x) all(x[3:6] == 0),
                           logical(1)))
    add("theo_inactive_occasions", inactive, inactive, "occasions 3-6 empty")
    peaks <- sim_eval_directional_peaks(synthetic, roles, endpoints)
    intensive <- peaks$observations >= 3
    add("theo_single_peak",
        any(intensive) && all(peaks$directional_peaks[intensive] <= 1L),
        if (any(intensive)) max(peaks$directional_peaks[intensive]) else NA,
        "<= 1")
    inside <- vapply(split(synthetic, synthetic[[roles$id]]), function(x) {
      event <- sim_eval_event_rows(x, roles)
      if (!is.null(roles$amt)) event <- event & as.numeric(x[[roles$amt]]) > 0
      doses <- sort(as.numeric(x[[roles$time]][event]))
      observations <- as.numeric(x[[roles$time]][sim_eval_observation_rows(x, roles)])
      occasion <- pmax(1L, findInterval(observations, doses))
      next_dose <- ifelse(occasion < length(doses), doses[occasion + 1L], Inf)
      all(observations < next_dose)
    }, logical(1))
    add("theo_within_occasion", all(inside), all(inside), "TRUE")
  } else if (case$id == "warfarin") {
    plotted <- sim_eval_plot_data(synthetic, roles, endpoints, "Synthetic")
    cp <- plotted[plotted$endpoint == "cp", , drop = FALSE]
    cp_last <- vapply(split(cp$time, cp$subject), max, numeric(1))
    add("warfarin_cp_after_24", any(cp$time > 24), max(cp$time), "> 24")
    add("warfarin_cp_median_last", stats::median(cp_last) >= 72,
        round(stats::median(cp_last), 3), ">= 72")
    cp_pair <- paired[paired$endpoint == "cp", , drop = FALSE]
    add("warfarin_cp_points",
        abs(cp_pair$mean_points_synthetic - cp_pair$mean_points_source) <= 1,
        round(cp_pair$mean_points_synthetic, 3),
        paste0(round(cp_pair$mean_points_source, 3), " +/- 1"))
    add("warfarin_endpoint_patients",
        all(paired$patients_with_endpoint_synthetic == paired$patients_source),
        paste(paired$patients_with_endpoint_synthetic, collapse = ","),
        paste(paired$patients_source, collapse = ","))
  } else if (case$id == "wbcSim") {
    event <- sim_eval_event_rows(synthetic, roles)
    amount <- as.numeric(synthetic[[roles$amt]][event])
    rate <- as.numeric(synthetic[[roles$rate]][event])
    add("wbc_infusion_pairs", any(amount > 0) && any(amount < 0),
        paste(range(amount), collapse = ","), "positive and negative")
    add("wbc_amount_rate", all(abs(amount - rate) <= 1e-8),
        max(abs(amount - rate)), "0")
    add("wbc_no_4580", !any(as.numeric(synthetic[[roles$time]]) == 4580),
        any(as.numeric(synthetic[[roles$time]]) == 4580), "FALSE")
    plotted <- sim_eval_plot_data(synthetic, roles, endpoints, "Synthetic")
    recovery <- vapply(split(plotted, plotted$subject), function(x) {
      x <- x[order(x$time), , drop = FALSE]
      if (nrow(x) < 3L) return(FALSE)
      nadir <- which.min(x$dv)
      nadir > 1L && nadir < nrow(x) && x$dv[nadir] < x$dv[1L] &&
        x$dv[nrow(x)] > x$dv[nadir]
    }, logical(1))
    add("wbc_decline_recovery", mean(recovery) >= .80,
        round(mean(recovery), 3), ">= 0.80")
  } else if (case$id == "nimoData") {
    property <- subject_property_summary(model)
    add(
      "nimo_dose_strata",
      identical(as.integer(property$DOS), c(50L, 100L, 200L, 400L)),
      paste(property$DOS, collapse = ","), "50,100,200,400"
    )
    add(
      "nimo_stratum_regimens",
      all(abs(property$dose_amount - property$DOS) <= 1e-8) &&
        all(property$n_doses == 10L) &&
        all(property$infusion_duration > 0),
      paste(round(property$dose_amount, 3), collapse = ","),
      "dose amount equals DOS; 10 doses; positive duration"
    )
    coherent <- vapply(split(synthetic, synthetic[[roles$id]]), function(x) {
      property_value <- unique(x$DOS)
      amount <- as.numeric(x[[roles$amt]])[
        sim_eval_event_rows(x, roles) & as.numeric(x[[roles$amt]]) > 0
      ]
      length(property_value) == 1L && length(amount) == 10L &&
        all(abs(amount - property_value) <= 1e-8)
    }, logical(1))
    add("nimo_property_regimen_coherence", all(coherent), mean(coherent), "1")
    amount <- as.numeric(synthetic[[roles$amt]])
    add(
      "nimo_infusion_pairs", any(amount > 0) && any(amount < 0),
      paste(round(range(amount), 3), collapse = ","),
      "positive starts and negative stops"
    )
    observed <- sim_eval_observation_rows(synthetic, roles)
    nonterminal <- observed & synthetic[[roles$occasion]] < 10L
    inside <- vapply(which(nonterminal), function(row) {
      subject <- synthetic[[roles$id]][row]
      next_occasion <- synthetic[[roles$occasion]][row] + 1L
      candidate <- synthetic[[roles$id]] == subject &
        synthetic[[roles$occasion]] == next_occasion &
        sim_eval_event_rows(synthetic, roles) &
        as.numeric(synthetic[[roles$amt]]) > 0
      any(candidate) &&
        as.numeric(synthetic[[roles$time]][row]) <
          min(as.numeric(synthetic[[roles$time]][candidate]))
    }, logical(1))
    add("nimo_nonterminal_boundary", all(inside), mean(inside), "1")
    terminal_tad <- as.numeric(synthetic[[roles$tad]])[
      observed & synthetic[[roles$occasion]] == 10L
    ]
    positive_times <- split(
      as.numeric(synthetic[[roles$time]])[
        sim_eval_event_rows(synthetic, roles) &
          as.numeric(synthetic[[roles$amt]]) > 0
      ],
      synthetic[[roles$id]][
        sim_eval_event_rows(synthetic, roles) &
          as.numeric(synthetic[[roles$amt]]) > 0
      ]
    )
    typical_interval <- stats::median(unlist(lapply(positive_times, diff)))
    add(
      "nimo_terminal_washout",
      length(terminal_tad) && max(terminal_tad) > typical_interval,
      if (length(terminal_tad)) round(max(terminal_tad), 3) else NA,
      paste(">", round(typical_interval, 3))
    )
  } else if (case$id == "mavoglurant") {
    positive_event <- sim_eval_event_rows(synthetic, roles) &
      as.numeric(synthetic[[roles$amt]]) > 0
    assigned <- as.numeric(synthetic[[roles$assigned_dose]])
    amount <- as.numeric(synthetic[[roles$amt]])
    add(
      "mav_assigned_dose_event",
      all(abs(assigned[positive_event] - amount[positive_event]) <= 1e-8),
      max(abs(assigned[positive_event] - amount[positive_event])), "0"
    )
    group <- interaction(
      synthetic[[roles$id]], synthetic[[roles$occasion]], drop = TRUE
    )
    constant <- vapply(split(assigned, group), function(x) {
      length(unique(x)) == 1L
    }, logical(1))
    add("mav_assigned_dose_constant", all(constant), mean(constant), "1")
    doses <- table(synthetic[[roles$id]][positive_event])
    add("mav_two_occasions", all(doses == 2L), paste(range(doses), collapse = ","),
        "2 positive dose events per subject")
    sex_spec <- model$population$covariates$SEX
    sex_values <- as.character(unique(synthetic$SEX))
    add(
      "mav_numeric_category",
      identical(sex_spec$type, "categorical") &&
        all(sex_values %in% c("1", "2")),
      paste(sort(sex_values), collapse = ","), "categorical levels 1,2"
    )
  } else if (case$id == "censoring") {
    observed <- sim_eval_observation_rows(synthetic, roles)
    states <- unique(as.integer(synthetic[[roles$cens]][observed]))
    source_states <- unique(as.integer(source[[roles$cens]]))
    add("censoring_source_states", all(c(-1L, 0L, 1L) %in% source_states),
        paste(sort(source_states), collapse = ","), "-1,0,1")
    add("censoring_generated_states",
        all(states %in% c(-1L, 0L, 1L)) && any(states != 0L),
        paste(sort(states), collapse = ","), "valid and at least one censored")
    interval <- observed & synthetic[[roles$cens]] == 1L &
      is.finite(synthetic[[roles$limit]])
    add("censoring_interval_limits", any(interval) &&
          all(synthetic[[roles$limit]][interval] <= synthetic[[roles$dv]][interval]),
        sum(interval), "> 0 coherent intervals")
  }

  do.call(rbind, result)
}

sim_eval_metric_rows <- function(case, model, synthetic, seed, epsilon,
                                 backend) {
  source <- sim_eval_design_summary(
    case$source, case$roles, case$endpoints, "Source", case$bounds$time,
    case$comparison_clock
  )
  generated <- sim_eval_design_summary(
    synthetic, case$roles, case$endpoints, "Synthetic", case$bounds$time,
    case$comparison_clock
  )
  paired <- merge(source, generated, by = "endpoint",
                  suffixes = c("_source", "_synthetic"))
  gates <- sim_eval_gate_results(case, model, synthetic)
  data.frame(
    dataset = case$id, seed = as.integer(seed), epsilon = epsilon,
    backend = backend, endpoint = paired$endpoint,
    source_subjects = paired$patients_source,
    synthetic_subjects = paired$patients_synthetic,
    source_observations = paired$observations_source,
    synthetic_observations = paired$observations_synthetic,
    source_mean_points = paired$mean_points_source,
    synthetic_mean_points = paired$mean_points_synthetic,
    source_first_time = paired$first_time_source,
    synthetic_first_time = paired$first_time_synthetic,
    source_last_time = paired$last_time_source,
    synthetic_last_time = paired$last_time_synthetic,
    gates_passed = sum(gates$pass), gates_failed = sum(!gates$pass),
    stringsAsFactors = FALSE
  )
}

sim_eval_plot <- function(case, synthetic, clock = case$comparison_clock) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  data <- rbind(
    sim_eval_plot_data(
      case$source, case$roles, case$endpoints, "Source", clock,
      case$bounds$time
    ),
    sim_eval_plot_data(
      synthetic, case$roles, case$endpoints, "Synthetic", clock,
      case$bounds$time
    )
  )
  grouping <- if (identical(clock, "tad")) {
    interaction(data$dataset, data$subject, data$endpoint, data$occasion)
  } else {
    interaction(data$dataset, data$subject, data$endpoint)
  }
  ggplot2::ggplot(
    data, ggplot2::aes(x = time, y = dv, colour = dataset, group = grouping)
  ) +
    ggplot2::geom_line(alpha = .30) +
    ggplot2::geom_point(alpha = .55, size = .7) +
    ggplot2::facet_grid(dataset ~ endpoint, scales = "free_y") +
    ggplot2::scale_colour_manual(
      values = c(Source = "#1B6CA8", Synthetic = "#D95F02")
    ) +
    ggplot2::labs(
      x = if (identical(clock, "tad")) "Time after dose" else "Study time",
      y = "DV", colour = "Dataset",
      title = paste(case$id, "simulation evaluation")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
}
