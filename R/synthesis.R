.validate_generator_options <- function(n_subjects, source_n, event_method,
                                        dv_method, k, pca_variance,
                                        subject_noise_sd, residual_noise_sd,
                                        residual_phi, time_jitter) {
  if (is.null(n_subjects)) n_subjects <- source_n
  if (length(n_subjects) != 1L || is.na(n_subjects) ||
      n_subjects < 1 || n_subjects != as.integer(n_subjects)) {
    stop("`n_subjects` must be NULL or one positive integer.", call. = FALSE)
  }
  if (!identical(event_method, "template")) {
    stop("The prototype supports only `event_method = \"template\"`.",
         call. = FALSE)
  }
  if (!identical(dv_method, "avatar_blend")) {
    stop("The prototype supports only `dv_method = \"avatar_blend\"`.",
         call. = FALSE)
  }
  if (length(k) != 1L || is.na(k) || k < 1 || k != as.integer(k)) {
    stop("`k` must be one positive integer.", call. = FALSE)
  }
  if (length(pca_variance) != 1L || is.na(pca_variance) ||
      pca_variance <= 0 || pca_variance > 1) {
    stop("`pca_variance` must be in (0, 1].", call. = FALSE)
  }
  for (argument in c("subject_noise_sd", "residual_noise_sd", "time_jitter")) {
    value <- get(argument)
    if (length(value) != 1L || is.na(value) || !is.finite(value) || value < 0) {
      stop("`", argument, "` must be one finite nonnegative number.",
           call. = FALSE)
    }
  }
  if (length(residual_phi) != 1L || is.na(residual_phi) ||
      !is.finite(residual_phi) || abs(residual_phi) >= 1) {
    stop("`residual_phi` must be finite and strictly between -1 and 1.",
         call. = FALSE)
  }
  as.integer(n_subjects)
}

.jitter_skeleton_time <- function(skeleton, roles, time_jitter) {
  if (time_jitter == 0) return(skeleton)
  time <- skeleton[[roles$time]]
  unique_time <- sort(unique(time))
  if (!length(unique_time)) return(skeleton)
  jittered <- unique_time + stats::rnorm(length(unique_time), sd = time_jitter)
  if (length(unique_time) > 1L) {
    midpoint <- (unique_time[-1L] + unique_time[-length(unique_time)]) / 2
    first_lower <- if (min(unique_time) >= 0) 0 else -Inf
    lower <- c(first_lower, midpoint)
    upper <- c(midpoint, Inf)
    jittered <- pmin(pmax(jittered, lower), upper)
  } else {
    jittered <- max(0, jittered)
  }
  skeleton[[roles$time]] <- jittered[match(time, unique_time)]
  skeleton
}

.subject_value <- function(data, rows, column) {
  .first_present(data[[column]][rows])
}

.synthesize_covariates <- function(skeleton, data, roles, donor_indices,
                                   weights, profiles, subject_noise_sd) {
  for (covariate in roles$covariates) {
    template <- data[[covariate]]
    values <- lapply(donor_indices, function(donor) {
      .subject_value(data, profiles$subject_rows[[donor]], covariate)
    })

    if (is.numeric(template) && !is.factor(template)) {
      numeric_values <- vapply(values, function(value) {
        if (!length(value) || is.na(value)) NA_real_ else as.numeric(value)
      }, numeric(1))
      okay <- is.finite(numeric_values)
      if (!any(okay)) {
        value <- .first_present(skeleton[[covariate]])
      } else {
        available <- numeric_values[okay]
        available_weights <- weights[okay] / sum(weights[okay])
        positive <- all(available > 0)
        skewed <- positive && length(available) > 1L &&
          max(available) / max(stats::median(available), .Machine$double.eps) > 3
        working <- if (skewed) log(available) else available
        center <- sum(working * available_weights)
        spread <- if (length(working) > 1L) stats::sd(working) else 0
        if (!is.finite(spread) || spread <= 0) {
          spread <- max(abs(center) * 0.05, 0.01)
        }
        value <- center + stats::rnorm(1L, sd = subject_noise_sd * spread)
        if (skewed) value <- exp(value)
        if (positive) value <- max(value, sqrt(.Machine$double.eps))
      }
    } else {
      character_values <- vapply(values, function(value) {
        if (!length(value) || is.na(value)) NA_character_ else as.character(value)
      }, character(1))
      okay <- !is.na(character_values)
      if (!any(okay)) {
        value <- .first_present(skeleton[[covariate]])
      } else {
        available_weights <- weights[okay] / sum(weights[okay])
        value <- sample(character_values[okay], 1L, prob = available_weights)
        if (is.logical(template)) value <- identical(value, "TRUE")
      }
    }
    skeleton[[covariate]][] <- value
  }
  skeleton
}

# Censoring -------------------------------------------------------------------
#
# AVATAR makes no formal privacy guarantee, so it may read the censoring
# boundary straight off the source data. A differentially private engine may
# not: there the limit is a declared public input, because a boundary inferred
# from confidential records is itself a release. That asymmetry is why this
# lives here rather than being shared with `fit_private`'s path.

# One endpoint's censoring boundaries, recovered from how the source reports
# censored rows under the Monolix convention: DV sits at the boundary, and
# LIMIT carries the other end of an interval when there is one.
.source_censoring <- function(data, roles, endpoint_name) {
  if (is.null(roles$cens)) return(NULL)
  observed <- .observation_rows(data, roles, require_present = TRUE) &
    .endpoint(data, roles) == endpoint_name
  if (!any(observed)) return(NULL)
  cens <- suppressWarnings(as.numeric(as.character(data[[roles$cens]][observed])))
  dv <- suppressWarnings(as.numeric(data[[roles$dv]][observed]))
  limit <- if (is.null(roles$limit)) rep(NA_real_, sum(observed)) else
    suppressWarnings(as.numeric(data[[roles$limit]][observed]))

  out <- list()
  left <- is.finite(cens) & cens == 1 & is.finite(dv) & !is.finite(limit)
  interval <- is.finite(cens) & cens == 1 & is.finite(dv) & is.finite(limit)
  right <- is.finite(cens) & cens == -1 & is.finite(dv)
  # A study with more than one assay limit collapses to the most conservative
  # boundary rather than inventing a per-row rule the source cannot support.
  if (any(left)) out$left <- max(dv[left])
  if (any(right)) out$right <- min(dv[right])
  if (any(interval)) {
    out$interval <- c(min(limit[interval]), max(dv[interval]))
  }
  if (!length(out)) NULL else out
}

# Replace reported values on censored rows with a latent draw inside the
# censoring region, so that blending sees plausible values rather than a stack
# of identical boundary substitutions. Without this every censored donor drags
# the blend toward the limit, and the synthetic data inherits a floor the real
# study does not have.
.impute_censored <- function(data, roles) {
  if (is.null(roles$cens)) return(data)
  observed <- .observation_rows(data, roles, require_present = TRUE)
  if (!any(observed)) return(data)
  endpoint <- .endpoint(data, roles)
  cens <- suppressWarnings(as.numeric(as.character(data[[roles$cens]])))
  dv <- suppressWarnings(as.numeric(data[[roles$dv]]))
  limit <- if (is.null(roles$limit)) rep(NA_real_, nrow(data)) else
    suppressWarnings(as.numeric(data[[roles$limit]]))

  for (name in unique(endpoint[observed])) {
    rows <- which(observed & endpoint == name & is.finite(cens) & cens != 0)
    if (!length(rows)) next
    for (i in rows) {
      if (cens[i] == 1 && is.finite(limit[i])) {
        # Interval-censored: the value lies between LIMIT and DV.
        dv[i] <- stats::runif(1L, min(limit[i], dv[i]), max(limit[i], dv[i]))
      } else if (cens[i] == 1) {
        # Left-censored: uniform below the limit. A uniform draw rather than a
        # fixed LLOQ/2 avoids replacing one artificial spike with another.
        dv[i] <- stats::runif(1L, 0, dv[i])
      } else if (cens[i] == -1) {
        # Right-censored: no upper bound is reported, so extend by the same
        # relative width a left-censored draw would span.
        dv[i] <- dv[i] * stats::runif(1L, 1, 2)
      }
    }
  }
  data[[roles$dv]] <- dv
  data
}

.donor_trajectory <- function(data, roles, rows, endpoint_name, transform) {
  subject_data <- data[rows, , drop = FALSE]
  observed <- .observation_rows(subject_data, roles, require_present = TRUE)
  endpoint <- .endpoint(subject_data, roles)
  selected <- observed & endpoint == endpoint_name
  time <- .aligned_time(subject_data, roles)[selected]
  value <- .transform_dv(subject_data[[roles$dv]][selected], transform)
  okay <- is.finite(time) & is.finite(value)
  list(time = time[okay], value = value[okay])
}

.interpolate_trajectory <- function(trajectory, target_time) {
  time <- trajectory$time
  value <- trajectory$value
  if (!length(time)) return(rep(NA_real_, length(target_time)))
  if (length(unique(time)) == 1L) return(rep(mean(value), length(target_time)))

  absolute <- stats::approx(time, value, xout = target_time,
                            ties = mean, rule = 1)$y
  missing <- !is.finite(absolute)
  if (any(missing) && length(unique(target_time)) > 1L) {
    target_range <- range(target_time)
    donor_range <- range(time)
    target_fraction <- (target_time[missing] - target_range[1L]) /
      diff(target_range)
    mapped_time <- donor_range[1L] + target_fraction * diff(donor_range)
    absolute[missing] <- stats::approx(
      time, value, xout = mapped_time, ties = mean, rule = 2
    )$y
  }
  absolute
}

.endpoint_noise_scale <- function(data, roles, endpoint_name, transform) {
  observed <- .observation_rows(data, roles, require_present = TRUE)
  endpoint <- .endpoint(data, roles)
  values <- .transform_dv(
    data[[roles$dv]][observed & endpoint == endpoint_name], transform
  )
  if (identical(transform$method, "log_offset")) return(1)
  scale <- stats::sd(values, na.rm = TRUE)
  if (!is.finite(scale) || scale <= 0) {
    scale <- max(abs(stats::median(values, na.rm = TRUE)) * 0.1, 0.01)
  }
  scale
}

.synthesize_trajectories <- function(skeleton, data, roles, donor_indices,
                                     weights, profiles, subject_noise_sd,
                                     residual_noise_sd, residual_phi,
                                     warnings,
                                     censoring_by_endpoint = list()) {
  allowed <- .observation_rows(skeleton, roles)
  present <- allowed & !is.na(skeleton[[roles$dv]])
  endpoint <- .endpoint(skeleton, roles)
  source_observed <- .observation_rows(data, roles, require_present = TRUE)
  source_endpoint <- .endpoint(data, roles)

  for (endpoint_name in unique(endpoint[present])) {
    target_rows <- which(present & endpoint == endpoint_name)
    target_time <- .aligned_time(skeleton, roles)[target_rows]
    transform <- profiles$transforms[[endpoint_name]]
    if (is.null(transform)) {
      transform <- .choose_transform(
        data[[roles$dv]][source_observed & source_endpoint == endpoint_name]
      )
    }
    donor_matrix <- vapply(donor_indices, function(donor) {
      trajectory <- .donor_trajectory(
        data, roles, profiles$subject_rows[[donor]], endpoint_name, transform
      )
      .interpolate_trajectory(trajectory, target_time)
    }, numeric(length(target_time)))
    if (is.null(dim(donor_matrix))) {
      donor_matrix <- matrix(donor_matrix, ncol = length(donor_indices))
    }
    if (subject_noise_sd == 0 && residual_noise_sd == 0 &&
        ncol(donor_matrix) > 1L) {
      donor_spread <- apply(donor_matrix, 1L, function(values) {
        values <- values[is.finite(values)]
        if (length(values) < 2L) Inf else diff(range(values))
      })
      if (length(donor_spread) &&
          all(donor_spread <= sqrt(.Machine$double.eps))) {
        warnings$add(paste0(
          "Endpoint `", endpoint_name,
          "` had indistinguishable donor trajectories with noise disabled; ",
          "an exact source-shaped trajectory was mathematically unavoidable."
        ))
      }
    }
    blended <- apply(donor_matrix, 1L, .weighted_available, weights = weights)
    missing_blend <- !is.finite(blended)
    if (any(missing_blend)) {
      fallback <- .transform_dv(
        data[[roles$dv]][source_observed & source_endpoint == endpoint_name],
        transform
      )
      fallback <- stats::median(fallback[is.finite(fallback)], na.rm = TRUE)
      if (!is.finite(fallback)) {
        stop("No usable DV values exist for endpoint `", endpoint_name, "`.",
             call. = FALSE)
      }
      blended[missing_blend] <- fallback
      warnings$add(paste0(
        "Endpoint `", endpoint_name,
        "` required a dataset-median interpolation fallback."
      ))
    }

    scale <- .endpoint_noise_scale(data, roles, endpoint_name, transform)
    shift <- stats::rnorm(1L, sd = subject_noise_sd * scale)
    residual <- .ar1_noise(
      length(blended), residual_phi, residual_noise_sd * scale
    )
    # The back-transformed blend is the latent value: what the subject would
    # have measured with no assay limit. DV, CENS, and LIMIT are then
    # reconstructed from it together, so the three always agree.
    generated <- .inverse_dv(blended + shift + residual, transform)
    skeleton[[roles$dv]][target_rows] <- generated
    censoring <- censoring_by_endpoint[[endpoint_name]]
    if (!is.null(censoring)) {
      skeleton <- .censor_latent(skeleton, target_rows, generated, roles,
                                 public = censoring)
    }
  }
  skeleton
}

.source_uses_standard_mdv <- function(data, roles) {
  if (is.null(roles$mdv)) return(FALSE)
  expected_observed <- .is_zero(data[[roles$evid]]) & !is.na(data[[roles$dv]])
  actual_observed <- .is_zero(data[[roles$mdv]])
  all(expected_observed == actual_observed)
}

.derive_standard_mdv <- function(skeleton, roles) {
  observed <- .is_zero(skeleton[[roles$evid]]) &
    !is.na(skeleton[[roles$dv]])
  template <- skeleton[[roles$mdv]]
  values <- ifelse(observed, 0, 1)
  if (is.factor(template)) {
    skeleton[[roles$mdv]] <- factor(
      as.character(values), levels = levels(template),
      ordered = is.ordered(template)
    )
  } else if (is.character(template)) {
    skeleton[[roles$mdv]] <- as.character(values)
  } else if (is.integer(template)) {
    skeleton[[roles$mdv]] <- as.integer(values)
  } else {
    skeleton[[roles$mdv]] <- as.numeric(values)
  }
  skeleton
}

.select_donors <- function(anchor, profiles, k, warnings) {
  compatible <- which(profiles$signatures == profiles$signatures[anchor])
  candidates <- setdiff(compatible, anchor)
  if (!length(candidates)) {
    warnings$add(
      "A compatible event-pattern group contained only its anchor; the anchor was used as the sole measurement donor and randomized noise supplied the only trajectory perturbation."
    )
    return(list(indices = anchor, distances = 0, weights = 1))
  }
  distances <- .neighbor_distances(profiles$coordinates, anchor, candidates)
  order_index <- order(distances, candidates)
  keep <- min(as.integer(k), length(candidates))
  if (keep < k) {
    warnings$add(paste0(
      "`k` was reduced to ", keep,
      " in at least one compatible event-pattern group."
    ))
  }
  chosen <- candidates[order_index[seq_len(keep)]]
  chosen_distances <- distances[order_index[seq_len(keep)]]
  if (length(chosen) < 2L) {
    warnings$add(
      "A compatible event-pattern group supplied fewer than two non-anchor donors."
    )
  }
  if (all(chosen_distances <= sqrt(.Machine$double.eps))) {
    warnings$add(
      "Duplicated subject profiles produced zero neighbor distances; epsilon-stabilized randomized weights were used."
    )
  }
  list(
    indices = chosen,
    distances = chosen_distances,
    weights = .randomized_weights(chosen_distances)
  )
}

#' Synthesize a structurally faithful PMX dataset (AVATAR-style)
#'
#' Samples complete subject event templates and fills them with AVATAR-like,
#' endpoint-specific blends of compatible subjects' baseline covariates and
#' longitudinal measurements. Event-control fields such as EVID, AMT, RATE,
#' CMT, and DVID are never averaged or independently generated.
#'
#' This is an AVATAR-inspired adaptation, not an exact implementation of
#' published AVATAR software. It creates synthetic data for model-workflow
#' exploration. It does not provide formal anonymization or preserve scientific
#' parameter or covariate-response relationships.
#'
#' @details
#' For selected compatible donors, randomized raw weights are
#' `Exp(1) / max(distance, epsilon) * 2^(-randomized_rank)`. They are normalized
#' and, when multiple donors are available, a dominant weight is capped at
#' 0.80 with its excess redistributed. The
#' same subject weights are used for covariates and all endpoints; weights are
#' renormalized locally when a donor lacks a requested endpoint/time value.
#'
#' Positive-like endpoints use an offset log scale and are constrained to be
#' nonnegative after back-transformation. Other endpoints use the identity
#' scale. Transform choices and interpolation alignment are recorded in the
#' returned `pmx_settings` attribute.
#'
#' @param data A source PMX data frame or tibble.
#' @param roles Explicit roles from [pmx_roles()]. Columns listed in
#'   `roles$exclude` are omitted from the generated output.
#' @param n_subjects Number of synthetic subjects. `NULL` retains the source count.
#' @param seed Reproducibility seed. The caller's random-number state is
#'   restored on exit.
#' @param event_method Event generation method. The prototype supports
#'   `"template"`.
#' @param dv_method Measurement method. The prototype supports
#'   `"avatar_blend"`.
#' @param k Maximum number of compatible non-anchor donors.
#' @param pca_variance Fraction of usable profile variance retained for
#'   neighborhood distances.
#' @param subject_noise_sd Nonnegative subject perturbation multiplier.
#' @param residual_noise_sd Nonnegative within-trajectory noise multiplier.
#' @param residual_phi AR(1) correlation in observation order, strictly between
#'   -1 and 1.
#' @param time_jitter Standard deviation for coherent tied-time jitter. Zero,
#'   the default, leaves the event template's times unchanged.
#'
#' @return An ordinary data frame or tibble with retained source columns, order,
#'   and practical classes. A lightweight `pmx_settings` attribute records the
#'   generator choices and endpoint transformations.
#' @export
#'
#' @examples
#' source <- data.frame(
#'   ID = rep(1:3, each = 4),
#'   TIME = rep(c(0, 0, 1, 2), 3),
#'   DV = c(0, 0.2, 2, 1, 0, 0.3, 3, 1.5, 0, 0.4, 4, 2),
#'   AMT = rep(c(100, 0, 0, 0), 3),
#'   EVID = rep(c(1L, 0L, 0L, 0L), 3),
#'   CMT = rep(c(1L, 2L, 2L, 2L), 3),
#'   WT = rep(c(60, 70, 80), each = 4)
#' )
#' roles <- pmx_roles("ID", "TIME", "DV", "AMT", "EVID", "CMT", NULL,
#'                    NULL, NULL, "WT")
#' synthetic <- synpmx_avatar(source, roles, n_subjects = 2, seed = 123)
#' validate_pmx(synthetic, roles)$valid
synpmx_avatar <- function(data, roles, n_subjects = NULL, seed = 123,
                     event_method = "template",
                     dv_method = "avatar_blend", k = 5,
                     pca_variance = 0.90, subject_noise_sd = 0.15,
                     residual_noise_sd = 0.05, residual_phi = 0.6,
                     time_jitter = 0) {
  if (!is.data.frame(data)) stop("`data` must be a data frame or tibble.",
                                 call. = FALSE)
  .assert_roles(data, roles)
  # Allowlist, not blocklist: keep only what a role names, so a column the user
  # never mentioned -- a secondary identifier, a site, a randomization date --
  # cannot ride out of a real subject into synthetic data by being forgotten.
  # `exclude` is still honoured for anyone who names a column there explicitly.
  retained_names <- setdiff(intersect(names(data),
                                      .retained_role_columns(roles)),
                            roles$exclude)
  dropped <- setdiff(names(data), retained_names)
  if (length(dropped)) {
    message("synpmx_avatar(): dropped ", length(dropped),
            " undeclared column(s): ", paste(dropped, collapse = ", "),
            ".\n  Declare a column in `keep` to carry it through verbatim.")
  }
  source <- data[, retained_names, drop = FALSE]
  source_roles <- roles
  source_roles$exclude <- NULL
  class(source_roles) <- "pmx_roles"
  validate_pmx(source, source_roles, strict = TRUE)
  subjects <- .unique_in_order(source[[source_roles$id]])
  n_subjects <- .validate_generator_options(
    n_subjects, length(subjects), event_method, dv_method, k, pca_variance,
    subject_noise_sd, residual_noise_sd, residual_phi, time_jitter
  )

  .with_local_seed(seed, {
    warnings <- .warning_collector()
    # Boundaries come from the reported data, before imputation overwrites the
    # substituted values they are read from.
    censoring_by_endpoint <- list()
    if (!is.null(source_roles$cens)) {
      endpoint_names <- unique(.endpoint(
        source, source_roles
      )[.observation_rows(source, source_roles, require_present = TRUE)])
      for (name in endpoint_names) {
        censoring_by_endpoint[[name]] <- .source_censoring(
          source, source_roles, name
        )
      }
      censored_rows <- sum(suppressWarnings(as.numeric(as.character(
        source[[source_roles$cens]]
      ))) != 0, na.rm = TRUE)
      if (censored_rows && !length(Filter(Negate(is.null),
                                          censoring_by_endpoint))) {
        warnings$add(paste0(
          "A `cens` role is declared and ", censored_rows, " row(s) are ",
          "flagged, but no censoring boundary could be read from them; ",
          "CENS was carried through without being reconstructed."
        ))
      }
      source <- .impute_censored(source, source_roles)
    }
    profiles <- .build_profiles(source, source_roles, pca_variance)
    new_ids <- .new_ids(source[[source_roles$id]], n_subjects)
    anchors <- sample.int(length(subjects), n_subjects, replace = TRUE)
    standard_mdv <- .source_uses_standard_mdv(source, source_roles)
    generated <- vector("list", n_subjects)

    for (synthetic_index in seq_len(n_subjects)) {
      anchor <- anchors[synthetic_index]
      skeleton <- source[profiles$subject_rows[[anchor]], , drop = FALSE]
      original_order <- seq_len(nrow(skeleton))
      skeleton <- .jitter_skeleton_time(skeleton, source_roles, time_jitter)
      donors <- .select_donors(anchor, profiles, k, warnings)
      skeleton <- .synthesize_covariates(
        skeleton, source, source_roles, donors$indices, donors$weights, profiles,
        subject_noise_sd
      )
      skeleton <- .synthesize_trajectories(
        skeleton, source, source_roles, donors$indices, donors$weights, profiles,
        subject_noise_sd, residual_noise_sd, residual_phi, warnings,
        censoring_by_endpoint
      )
      if (standard_mdv) skeleton <- .derive_standard_mdv(skeleton, source_roles)
      new_id <- new_ids[synthetic_index]
      if (is.factor(skeleton[[source_roles$id]])) {
        new_label <- as.character(new_id)
        levels(skeleton[[source_roles$id]]) <- unique(c(
          levels(skeleton[[source_roles$id]]), new_label
        ))
        skeleton[[source_roles$id]][] <- new_label
      } else {
        skeleton[[source_roles$id]][] <- new_id
      }
      row_order <- order(skeleton[[source_roles$time]], original_order)
      generated[[synthetic_index]] <- skeleton[row_order, , drop = FALSE]
    }

    result <- do.call(rbind, generated)
    rownames(result) <- NULL
    result <- .restore_schema(result, source, source_roles)
    settings <- list(
      seed = as.integer(seed),
      n_subjects = n_subjects,
      event_method = event_method,
      dv_method = dv_method,
      k = as.integer(k),
      pca_variance = pca_variance,
      subject_noise_sd = subject_noise_sd,
      residual_noise_sd = residual_noise_sd,
      residual_phi = residual_phi,
      time_jitter = time_jitter,
      roles = unclass(source_roles),
      endpoint_transforms = profiles$transforms,
      alignment = paste(
        "time relative to first positive dose within compatible schedules;",
        "normalized observation-window fallback"
      ),
      compatible_event_groups = length(unique(profiles$signatures)),
      max_donor_weight = 0.80,
      warnings = warnings$messages
    )
    attr(result, "pmx_settings") <- settings
    validate_pmx(result, source_roles, strict = TRUE)
    if (length(warnings$messages)) {
      warning(
        "Synthetic generation used documented small-group/profile fallbacks:\n- ",
        paste(warnings$messages, collapse = "\n- "),
        call. = FALSE
      )
    }
    result
  })
}
