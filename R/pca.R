# PCA generation ------------------------------------------------------------
#
# Documented in `pca-algorithm.Rmd`. Change the behaviour of a step here and
# change its section there in the same commit.
#
# The contrast with `synpmx_avatar()` is what the principal components are used
# for. AVATAR uses them as a distance metric: the coordinates decide who donates
# to whom, and the values written out are weighted averages of real patients'
# numbers. Here they are a basis. Subject trajectories become scores, the scores
# become the response of a fitted linear model, and new subjects come from
# drawing that model's residual. No donor index exists in this file, and no
# number a patient measured reaches the output.

# One subject-by-feature matrix: every retained grid cell of every endpoint,
# beside every baseline covariate. Covariates and trajectories share one matrix
# so that a single draw carries the relationship between them.
#
# Documented in `pca-algorithm.Rmd`, Step 2.
.pca_features <- function(source, roles, min_column_patients,
                          transform_source = NULL) {
  subjects <- .unique_in_order(source[[roles$id]])
  n <- length(subjects)
  subject_rows <- lapply(subjects, function(subject) {
    !is.na(source[[roles$id]]) & source[[roles$id]] == subject
  })
  observed <- .observation_rows(source, roles, require_present = TRUE)
  endpoint <- .endpoint(source, roles)
  aligned <- .aligned_time(source, roles)
  endpoints <- sort(unique(endpoint[observed]))

  # The log transform's offset is half the smallest positive value, and it exists
  # to keep log() finite near zero. It must be read from what the assay actually
  # reported, not from the imputed values: `.impute_censored()` draws uniformly
  # on the raw scale, so a value censored at an LLOQ of 0.05 can land at 5e-05
  # and drag the offset three orders of magnitude down with it. The log range
  # then runs to -9.5 instead of -2.6, and 46% of the observations occupy seven
  # log units of what is purely uniform noise. Generation adds a Gaussian
  # residual on that scale, and exponentiating a wide Gaussian has a heavy upper
  # tail, so a low-dose arm that is entirely below quantification comes back with
  # detectable concentrations that rise with time.
  basis_for_transform <- transform_source %||% source
  transforms <- stats::setNames(lapply(endpoints, function(ep) {
    .choose_transform(suppressWarnings(as.numeric(
      basis_for_transform[[roles$dv]][observed & endpoint == ep]
    )))
  }), endpoints)

  # The whole nominal grid, with no cap on its width. AVATAR caps at fifteen
  # points because more columns buy a distance metric nothing; a basis has no
  # such ceiling, and `prcomp()` is bounded by the subject count either way.
  grids <- stats::setNames(lapply(endpoints, function(ep) {
    sort(unique(aligned[observed & endpoint == ep & is.finite(aligned)]))
  }), endpoints)

  features <- list()
  kinds <- character()
  members <- list()

  for (covariate in roles$covariates) {
    template <- source[[covariate]]
    values <- lapply(subject_rows, function(rows) .first_present(template[rows]))
    if (is.numeric(template) && !is.factor(template)) {
      name <- paste0("cov_", covariate)
      features[[name]] <- as.numeric(unlist(values))
      kinds[name] <- "covariate_continuous"
      members[[name]] <- list(covariate = covariate)
    } else {
      characters <- vapply(values, function(value) {
        if (!length(value) || is.na(value)) NA_character_ else as.character(value)
      }, character(1))
      levels_present <- if (is.factor(template)) levels(template) else
        sort(unique(characters[!is.na(characters)]))
      for (level in levels_present) {
        name <- paste0("cov_", covariate, "__", make.names(level))
        indicator <- as.numeric(characters == level)
        indicator[is.na(characters)] <- NA_real_
        features[[name]] <- indicator
        kinds[name] <- "covariate_level"
        members[[name]] <- list(covariate = covariate, level = level)
      }
    }
  }

  for (ep in endpoints) {
    grid <- grids[[ep]]
    if (!length(grid)) next
    trajectories <- t(vapply(subject_rows, function(rows) {
      selected <- rows & observed & endpoint == ep
      values <- .transform_dv(
        suppressWarnings(as.numeric(source[[roles$dv]][selected])),
        transforms[[ep]]
      )
      .trajectory_on_grid(aligned[selected], values, grid)
    }, numeric(length(grid))))
    held <- colSums(is.finite(trajectories))
    # A grid cell only a handful of patients ever reached is filled from the
    # median for everybody else, so it carries those patients and describes
    # nobody. Dropping it is cheaper than modelling it.
    for (column in which(held >= min_column_patients)) {
      name <- paste0("dv_", make.names(ep), "__", column)
      features[[name]] <- trajectories[, column]
      kinds[name] <- "endpoint_cell"
      members[[name]] <- list(endpoint = ep, time = grid[column],
                              patients = held[[column]])
    }
  }

  if (!length(features)) {
    stop("No usable features: every endpoint grid cell was held by fewer than ",
         min_column_patients, " patients.", call. = FALSE)
  }
  matrix_out <- as.matrix(as.data.frame(features, check.names = FALSE))
  rownames(matrix_out) <- as.character(subjects)
  list(matrix = matrix_out, kinds = kinds[colnames(matrix_out)],
       members = members[colnames(matrix_out)], transforms = transforms,
       grids = grids, endpoints = endpoints, subjects = subjects, n = n)
}

# Total amount administered, per subject. The one design fact that enters the
# model as a predictor rather than as a response.
.pca_subject_dose <- function(source, roles, subjects) {
  if (is.null(roles$amt)) return(stats::setNames(rep(0, length(subjects)),
                                                 as.character(subjects)))
  amount <- suppressWarnings(as.numeric(source[[roles$amt]]))
  amount[!is.finite(amount)] <- 0
  dose_rows <- .event_rows(source, roles)
  totals <- vapply(subjects, function(subject) {
    rows <- !is.na(source[[roles$id]]) & source[[roles$id]] == subject
    sum(amount[rows & dose_rows])
  }, numeric(1))
  stats::setNames(totals, as.character(subjects))
}

# Median fill, standardize, decompose, then model the retained scores against
# the arm.
#
# Two ways to carry dose. `"factor"` gives each arm its own mean score vector
# and its own residual covariance; `"log"` regresses the scores on `log1p` of
# the total dose, which is one coefficient instead of one mean per arm and
# extrapolates to doses the study did not run. `"factor"` is the default because
# neither of the things `"log"` assumes holds on a typical dose-ranging study:
# a lower limit of quantification (LLOQ) puts the low arms on a floor, so the
# dose-response is not a straight line through the origin, and the arms on that
# floor have far less between-subject spread than the top arm, so a single
# pooled covariance smears them upward and pulls the top arm in. Measured on a
# six-arm study, `"log"` compressed a source PK range of 0.05--0.44 into
# 0.10--0.18. Use `"log"` where the arms are genuinely log-linear and there are
# too many dose levels to spend a mean on each.
#
# Documented in `pca-algorithm.Rmd`, Steps 3 and 4.
.pca_fit <- function(features, dose, group, pca_variance, n_components,
                     dose_term) {
  raw <- features$matrix
  n <- nrow(raw)
  centers <- apply(raw, 2L, stats::median, na.rm = TRUE)
  filled <- raw
  for (column in seq_len(ncol(filled))) {
    missing <- !is.finite(filled[, column])
    filled[missing, column] <- centers[[column]]
  }
  centers <- colMeans(filled)
  scales <- apply(filled, 2L, stats::sd)
  keep <- is.finite(scales) & scales > sqrt(.Machine$double.eps)
  if (!any(keep)) {
    stop("Every feature is constant across subjects; there is nothing to fit.",
         call. = FALSE)
  }
  filled <- filled[, keep, drop = FALSE]
  centers <- centers[keep]
  scales <- scales[keep]
  standardized <- sweep(sweep(filled, 2L, centers, "-"), 2L, scales, "/")

  pca <- stats::prcomp(standardized, center = FALSE, scale. = FALSE,
                       rank. = min(n - 1L, ncol(standardized)))
  variance <- pca$sdev^2
  explained <- cumsum(variance) / sum(variance)
  # Capped at a fifth of the cohort as well as by explained variance. With many
  # components against few subjects the basis interpolates, and a drawn subject
  # can land on a real one.
  ceiling_k <- max(1L, min(length(pca$sdev), floor(n / 5)))
  k <- if (!is.null(n_components)) {
    min(as.integer(n_components), length(pca$sdev))
  } else {
    min(which(explained >= pca_variance)[1L], ceiling_k)
  }
  if (!is.finite(k) || k < 1L) k <- 1L
  scores <- pca$x[, seq_len(k), drop = FALSE]
  rotation <- pca$rotation[, seq_len(k), drop = FALSE]

  chol_of <- function(covariance) {
    ridge <- max(max(diag(covariance)), 1e-12) * 1e-06
    chol(covariance + diag(ridge, nrow(covariance)))
  }

  log_dose <- log1p(as.numeric(dose))
  has_dose <- stats::sd(log_dose) > sqrt(.Machine$double.eps)
  group <- as.character(group)
  groups <- unique(group)

  if (identical(dose_term, "log") && has_dose) {
    design <- cbind(intercept = 1, log_dose = log_dose)
    coefficients <- qr.solve(design, scores)
    residuals <- scores - design %*% coefficients
    covariance <- crossprod(residuals) / max(1L, n - ncol(design))
    model <- list(kind = "log", coefficients = coefficients,
                  chol = chol_of(covariance), covariance = covariance)
  } else {
    pooled_residual <- scores
    means <- list()
    for (value in groups) {
      rows <- group == value
      means[[value]] <- colMeans(scores[rows, , drop = FALSE])
      pooled_residual[rows, ] <- sweep(scores[rows, , drop = FALSE], 2L,
                                       means[[value]], "-")
    }
    pooled <- crossprod(pooled_residual) / max(1L, n - length(groups))
    # Each arm keeps its own spread, shrunk toward the pooled one by its size.
    # An arm on the LLOQ floor is genuinely tighter than the top arm and must
    # not be given the cohort's variability; an arm of twenty cannot support a
    # ten-by-ten covariance on its own either.
    chols <- list()
    covariances <- list()
    for (value in groups) {
      rows <- group == value
      size <- sum(rows)
      own <- crossprod(pooled_residual[rows, , drop = FALSE]) / max(1L, size - 1L)
      weight <- size / (size + k)
      covariance <- weight * own + (1 - weight) * pooled
      covariances[[value]] <- covariance
      chols[[value]] <- chol_of(covariance)
    }
    model <- list(kind = "factor", means = means, chols = chols,
                  covariance = pooled, covariances = covariances)
  }

  list(
    columns = colnames(filled), kinds = features$kinds[colnames(filled)],
    members = features$members[colnames(filled)],
    centers = centers, scales = scales, rotation = rotation,
    model = model, dose_term = model$kind, groups = groups,
    min_group = if (identical(model$kind, "log")) n else
      as.integer(min(table(group))),
    variance = variance, explained = explained, k = k, has_dose = has_dose,
    n_source = n, transforms = features$transforms, grids = features$grids,
    endpoints = features$endpoints
  )
}

# Draw scores for new subjects and invert the standardization, giving one
# feature vector each.
#
# Documented in `pca-algorithm.Rmd`, Step 5.
.pca_draw <- function(fit, doses, groups) {
  n <- length(groups)
  model <- fit$model
  noise <- matrix(stats::rnorm(n * fit$k), nrow = n)
  scores <- matrix(0, nrow = n, ncol = fit$k)
  if (identical(model$kind, "log")) {
    design <- cbind(intercept = 1, log_dose = log1p(as.numeric(doses)))
    scores <- design %*% model$coefficients + noise %*% model$chol
  } else {
    for (i in seq_len(n)) {
      value <- groups[[i]]
      scores[i, ] <- model$means[[value]] +
        noise[i, , drop = FALSE] %*% model$chols[[value]]
    }
  }
  standardized <- scores %*% t(fit$rotation)
  drawn <- sweep(sweep(standardized, 2L, fit$scales, "*"), 2L, fit$centers, "+")
  colnames(drawn) <- fit$columns
  drawn
}

# The dosing model: a planned schedule per arm, plus three rates that say how
# patients departed from it.
#
# A dose-ranging study gives everyone the same schedule and the interesting
# quantity is that schedule. An oncology study gives everyone the same *planned*
# schedule and then almost nobody follows it: doses are reduced for toxicity,
# cycles are skipped, and patients come off treatment at different times. The
# relative dose intensity that results is often the thing the dataset exists to
# analyse, so a generator that hands every patient the modal schedule has
# removed it.
#
# Both cases are the same model here. Each arm carries:
#
#   planned       one row per cycle: the nominal time, and the modal amount at
#                 that time. Within-patient escalation lives here, because the
#                 planned amount is read per cycle rather than once.
#   levels        the dose ladder, as ratios of the planned amount, descending
#                 from 1. Read from ratios several patients share.
#   discontinuation   P(this cycle is the last | still on treatment)
#   interruption      P(this cycle is skipped | still on treatment)
#   reduction         P(drop one level | still on, a lower level exists)
#
# A study where nobody reduces, skips or stops early has all three rates at zero
# and one level, and the model then reproduces the planned schedule exactly for
# every patient. That is the point: no detector decides which kind of study this
# is, because the rates already say it.
#
# What leaves the source is three probabilities and a ladder, per arm. No
# patient's own sequence of reductions is copied, and the planned grid stops at
# the last cycle `min_arm_patients` patients reached, so a schedule length only
# one patient has cannot be generated -- the same reasoning as `SIM-047` on the
# observation side.
#
# Documented in `pca-algorithm.Rmd`, Step 6.
.pca_dose_model <- function(times, amounts, member_rows, aligned, amount,
                            dose_rows, floor) {
  n <- length(member_rows)
  per_patient <- lapply(seq_len(n), function(i) {
    selected <- member_rows[[i]] & dose_rows
    order_by <- order(aligned[selected])
    data.frame(time = aligned[selected][order_by],
               amt = amount[selected][order_by])
  })

  # The planned grid: nominal dose times enough of the arm reached. Cycles only
  # one or two patients got to are dropped rather than modelled, so generated
  # follow-up cannot run past what the arm shares.
  all_times <- sort(unique(unlist(lapply(per_patient, function(p) p$time))))
  reached <- vapply(all_times, function(t) {
    sum(vapply(per_patient, function(p) any(abs(p$time - t) < 1e-8), logical(1)))
  }, integer(1))
  planned_times <- all_times[reached >= floor]
  if (!length(planned_times)) planned_times <- all_times[which.max(reached)]

  # Each patient's amounts over the planned cycles, and their ratio to their
  # OWN starting dose. Clinical practice reduces to a fraction of what this
  # patient started on, and reading the ratio that way also survives a study
  # dosed by body weight, where no two patients share an amount at all.
  n_cycles <- length(planned_times)
  amounts_at <- lapply(per_patient, function(p) {
    vapply(planned_times, function(t) {
      hit <- which(abs(p$time - t) < 1e-8)
      if (length(hit)) p$amt[hit[1L]] else NA_real_
    }, numeric(1))
  })
  first_amt <- vapply(per_patient, function(p) {
    if (!nrow(p)) NA_real_ else p$amt[which.min(p$time)]
  }, numeric(1))
  ratios <- lapply(seq_len(n), function(i) {
    if (!is.finite(first_amt[i]) || first_amt[i] <= 0) {
      return(rep(NA_real_, n_cycles))
    }
    amounts_at[[i]] / first_amt[i]
  })

  # The planned amount at a cycle is the modal amount among patients who are
  # still on their starting dose there. Taking the modal amount over everyone
  # would drift downward as reductions accumulate, and a patient at half dose
  # would then read as two thirds of a plan that had already fallen. Reading it
  # per cycle is also what lets a protocol-prescribed escalation stay part of
  # the plan rather than register as a departure from it.
  planned_amt <- vapply(seq_len(n_cycles), function(k) {
    unreduced <- vapply(seq_len(n), function(i) {
      r <- ratios[[i]][k]
      is.finite(r) && abs(r - 1) < 0.05
    }, logical(1))
    values <- vapply(which(unreduced), function(i) amounts_at[[i]][k], numeric(1))
    if (!length(values)) {
      values <- unlist(lapply(amounts_at, function(a) a[k]))
      values <- values[is.finite(values)]
    }
    if (!length(values)) return(0)
    counts <- table(sprintf("%.10g", values))
    as.numeric(names(counts)[which.max(counts)])
  }, numeric(1))
  planned <- data.frame(cycle = seq_len(n_cycles), time = planned_times,
                        amt = planned_amt)

  span <- integer(n)
  skipped <- integer(n)
  for (i in seq_len(n)) {
    p <- per_patient[[i]]
    if (!nrow(p)) { span[i] <- 0L; next }
    within <- planned$time <= max(p$time) + 1e-8
    span[i] <- sum(within)
    skipped[i] <- sum(!is.finite(amounts_at[[i]][within]))
    ratios[[i]][!within] <- NA_real_
  }

  # The ladder, and it is built from WITHIN-patient decreases rather than from
  # the spread of ratios. A study dosed by body weight gives every patient a
  # slightly different amount, so pooled ratios scatter continuously around 1
  # and any threshold on them invents a ladder out of ordinary between-patient
  # variation. A dose reduction is something else: the same patient receiving
  # less than they did at the previous cycle. Where no patient's amount ever
  # falls, there are no reductions and the ladder is a single level.
  drop_tolerance <- 0.05
  dropped_to <- unlist(lapply(ratios, function(r) {
    finite <- which(is.finite(r))
    if (length(finite) < 2L) return(NULL)
    values <- r[finite]
    fell <- which(diff(values) < -drop_tolerance * utils::head(values, -1L))
    if (!length(fell)) return(NULL)
    values[fell + 1L]
  }))
  levels <- 1
  if (length(dropped_to)) {
    counts <- table(round(dropped_to, 2))
    holders <- vapply(names(counts), function(value) {
      sum(vapply(ratios, function(r) {
        finite <- which(is.finite(r))
        if (length(finite) < 2L) return(FALSE)
        values <- r[finite]
        fell <- which(diff(values) < -drop_tolerance * utils::head(values, -1L))
        length(fell) > 0 &&
          any(abs(round(values[fell + 1L], 2) - as.numeric(value)) < 1e-8)
      }, logical(1)))
    }, integer(1))
    shared <- as.numeric(names(counts))[holders >= floor]
    shared <- shared[shared > 0 & shared < 1 - drop_tolerance]
    levels <- sort(unique(c(1, shared)), decreasing = TRUE)
  }

  level_of <- function(r) {
    out <- rep(NA_integer_, length(r))
    finite <- which(is.finite(r))
    for (j in finite) {
      out[j] <- which.min(abs(levels - r[j]))
    }
    out
  }

  # Discrete-time hazards, pooled over the arm.
  at_risk_stop <- 0L; stops <- 0L
  at_risk_skip <- 0L; skips <- 0L
  at_risk_drop <- 0L; drops <- 0L
  for (i in seq_len(n)) {
    if (!span[i]) next
    at_risk_stop <- at_risk_stop + span[i]
    # A patient who reached the last planned cycle was never observed to stop.
    if (span[i] < n_cycles) stops <- stops + 1L
    at_risk_skip <- at_risk_skip + span[i]
    skips <- skips + skipped[i]
    if (length(levels) > 1L) {
      values <- ratios[[i]][is.finite(ratios[[i]])]
      if (length(values) > 1L) {
        idx <- level_of(values)
        at_risk_drop <- at_risk_drop +
          sum(utils::head(idx, -1L) < length(levels))
        drops <- drops +
          sum(diff(values) < -drop_tolerance * utils::head(values, -1L))
      }
    }
  }
  rate <- function(events, at_risk) {
    if (!at_risk) return(0)
    min(1, max(0, events / at_risk))
  }

  list(
    planned = planned,
    levels = levels,
    discontinuation = rate(stops, at_risk_stop),
    interruption = rate(skips, at_risk_skip),
    reduction = rate(drops, at_risk_drop),
    patients = n,
    distinct = length(unique(vapply(per_patient, function(p) {
      paste(sprintf("%.6g", p$time), sprintf("%.6g", p$amt), collapse = "|")
    }, character(1)))),
    source_doses = mean(vapply(per_patient, nrow, integer(1)))
  )
}

# Simulate one patient's schedule from an arm's dosing model. Reduction is
# decided before the cycle is dosed, discontinuation after it, and an
# interruption skips the cycle without ending treatment.
#
# Documented in `pca-algorithm.Rmd`, Step 6.
.pca_draw_schedule <- function(dosing) {
  planned <- dosing$planned
  levels <- dosing$levels
  level <- 1L
  keep <- logical(nrow(planned))
  amounts <- numeric(nrow(planned))
  for (i in seq_len(nrow(planned))) {
    if (level < length(levels) && stats::runif(1) < dosing$reduction) {
      level <- level + 1L
    }
    if (stats::runif(1) >= dosing$interruption) {
      keep[i] <- TRUE
      amounts[i] <- planned$amt[i] * levels[level]
    }
    if (stats::runif(1) < dosing$discontinuation) break
  }
  data.frame(time = planned$time[keep], amt = amounts[keep])
}

# The dosing model and the visit model, one of each per arm.
#
# Both are summaries of the arm rather than facts about a patient. The dosing
# model is built above: a planned schedule and three rates. The visit model is,
# per endpoint and per retained nominal time, the fraction of the arm that has
# an observation there, so attendance is drawn per visit rather than a real
# patient's set of attended visits being reused.
#
# Documented in `pca-algorithm.Rmd`, Step 6.
.pca_arm_models <- function(source, roles, fit, subject_group, floor) {
  subjects <- .unique_in_order(source[[roles$id]])
  aligned <- .aligned_time(source, roles)
  observed <- .observation_rows(source, roles, require_present = TRUE)
  endpoint <- .endpoint(source, roles)
  # Every dosing event, not only the ones carrying drug. A placebo arm records
  # its administrations with `AMT = 0`, and `.dose_rows()` drops those whenever
  # any positive amount exists in the study, which would leave the placebo arm
  # with no dosing events at all.
  dose_rows <- .event_rows(source, roles)
  amount <- if (is.null(roles$amt)) rep(0, nrow(source)) else
    suppressWarnings(as.numeric(source[[roles$amt]]))
  amount[!is.finite(amount)] <- 0

  cells <- which(fit$kinds == "endpoint_cell")
  arms <- unique(subject_group)
  dosing <- list()
  visits <- list()
  sizes <- integer()

  for (arm in arms) {
    members <- subjects[subject_group == arm]
    member_rows <- lapply(members, function(subject) {
      !is.na(source[[roles$id]]) & source[[roles$id]] == subject
    })
    sizes[arm] <- length(members)

    dosing[[arm]] <- .pca_dose_model(NULL, NULL, member_rows, aligned, amount,
                                     dose_rows, floor)

    probability <- vapply(cells, function(column) {
      member <- fit$members[[column]]
      reached <- vapply(member_rows, function(rows) {
        selected <- rows & observed & endpoint == member$endpoint
        any(is.finite(aligned[selected]) &
              abs(aligned[selected] - member$time) < sqrt(.Machine$double.eps))
      }, logical(1))
      mean(reached)
    }, numeric(1))
    visits[[arm]] <- list(cells = cells, probability = probability)
  }
  list(dosing = dosing, visits = visits, sizes = sizes, arms = arms,
       cells = cells)
}

# Everything else the generator needs to write a table in the source's shape,
# read once, here, so that generation touches no patient row.
#
# The column prototypes are zero-length vectors: they carry class and factor
# levels and no values. `id_offset` is one number, the largest source ID, so a
# synthetic ID cannot collide with a real one.
#
# Documented in `pca-algorithm.Rmd`, Step 8.
.pca_schema <- function(source, roles, fit, subject_group) {
  observed <- .observation_rows(source, roles, require_present = TRUE)
  endpoint <- .endpoint(source, roles)
  dose_rows <- .event_rows(source, roles)

  mode_of <- function(rows, column) {
    values <- source[[column]][rows]
    values <- values[!is.na(values)]
    if (!length(values)) return(NA)
    counts <- table(as.character(values))
    values[match(names(counts)[which.max(counts)], as.character(values))]
  }
  cmt_dose <- if (is.null(roles$cmt)) NULL else mode_of(dose_rows, roles$cmt)
  cmt_obs <- stats::setNames(lapply(fit$endpoints, function(ep) {
    if (is.null(roles$cmt)) NULL else mode_of(observed & endpoint == ep,
                                              roles$cmt)
  }), fit$endpoints)

  carried <- intersect(c(roles$strata, roles$keep), names(source))
  subjects <- .unique_in_order(source[[roles$id]])
  arm_values <- list()
  for (arm in unique(subject_group)) {
    subject <- subjects[subject_group == arm][1L]
    rows <- which(!is.na(source[[roles$id]]) & source[[roles$id]] == subject)
    arm_values[[arm]] <- lapply(stats::setNames(carried, carried),
                                function(column) {
      value <- .first_present(source[[column]][rows])
      if (is.factor(value)) as.character(value) else value
    })
  }

  # The assay limit, per endpoint: one or two numbers read from the source, and
  # the only way the generator can put a value back on the boundary. Without it
  # every value below the limit is emitted as itself and an arm that is entirely
  # below quantification comes back as a spread of small numbers rather than the
  # flat line the study recorded.
  censoring <- stats::setNames(lapply(fit$endpoints, function(ep) {
    .source_censoring(source, roles, ep)
  }), fit$endpoints)

  identifiers <- source[[roles$id]]
  list(
    censoring = censoring,
    columns = names(source),
    prototypes = lapply(stats::setNames(names(source), names(source)),
                        function(column) source[[column]][0L]),
    id_class = class(identifiers)[[1L]],
    id_offset = if (is.numeric(identifiers) && any(!is.na(identifiers))) {
      max(identifiers, na.rm = TRUE)
    } else 0,
    id_levels = if (is.factor(identifiers)) levels(identifiers) else NULL,
    cmt_dose = cmt_dose, cmt_obs = cmt_obs,
    carried = carried, arm_values = arm_values,
    endpoint_specs = .endpoint_value_types(source, roles)
  )
}

# The nominal grid is not optional here, and it is not inferred.
#
# `synpmx_avatar()` derives a grid when none is declared, and that derivation is
# a masking mechanism: it decides where the visits are so that no patient's own
# time vector identifies them. In this algorithm the grid is the model's axis.
# Every feature is a cell on it, dose times and observation times are placed on
# it together, and a grid inferred from elapsed time puts a sample on the wrong
# side of a dose as soon as the doses were recorded as actuals. Deriving it
# would be a statement about the protocol, which the caller is the only one in a
# position to make.
#
# Documented in `pca-algorithm.Rmd`, Step 2.
.pca_require_nominal_time <- function(source, roles) {
  if (is.null(roles$nominal_time)) {
    stop("`synpmx_pca()` requires `nominal_time` in `pmx_roles()`. The nominal ",
         "grid is the axis every feature sits on, and inferring it from ",
         "recorded times is a statement about the protocol that only you can ",
         "make. Add the protocol's planned times as a column and declare it.",
         call. = FALSE)
  }
  nominal <- suppressWarnings(as.numeric(source[[roles$nominal_time]]))
  relevant <- .observation_rows(source, roles, require_present = TRUE) |
    .event_rows(source, roles)
  missing <- relevant & !is.finite(nominal)
  if (any(missing)) {
    stop("`nominal_time` is missing on ", sum(missing), " of ", sum(relevant),
         " dose and observation rows. Every row the model reads needs a ",
         "nominal time; fill them in or drop those rows before calling.",
         call. = FALSE)
  }
  source[[roles$time]] <- nominal
  source
}

# An arm of one or two has no between-subject spread to model: its mean score
# vector is that patient, and its covariance is noise around them. Refusing is
# the only honest answer, and it is loud rather than a silent pooling the caller
# never asked for.
.pca_require_arms <- function(group, minimum) {
  minimum <- as.integer(minimum)
  if (!is.finite(minimum) || minimum < 1L) {
    stop("`min_arm_patients` must be one positive integer.", call. = FALSE)
  }
  sizes <- table(group)
  short <- sizes[sizes < minimum]
  if (length(short)) {
    stop("`synpmx_pca()` needs at least ", minimum,
         " patients in every arm. Short: ",
         paste(sprintf("%s (%d)", names(short), as.integer(short)),
               collapse = ", "),
         ". Pool the arm, drop the column from `strata`, or exclude those ",
         "patients before calling.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Summarize a trial into the quantities a synthetic copy is built from
#'
#' The only stage that reads patient data. Reduces each subject's trajectories
#' and baseline covariates to principal-component scores, models those scores
#' against the arm, and fits a dosing model and a visit model per arm. The
#' returned object holds nothing but summaries: no patient row survives it.
#'
#' Run this, look at what it produced with [pca_report()],
#' [pca_dosing()], [pca_visits()] and [pca_components()], then pass
#' it to [synpmx_pca_generate()]. Generation reads the model and nothing else,
#' so what those four functions show is the whole of what the synthetic data was
#' built from.
#'
#' `nominal_time` is required. The grid it names is the axis every feature sits
#' on, and dose rows and observation rows are placed on it together, so it is
#' the caller's statement about the protocol rather than something inferred from
#' recorded times.
#'
#' No formal privacy claim is made.
#'
#' @param data Source PMX event data.
#' @param roles Explicit column roles from [pmx_roles()], including
#'   `nominal_time`.
#' @param seed Seed for the one random step in summarizing: censored values are
#'   replaced by a draw inside the censoring region before the basis is fitted,
#'   so that a column where most subjects sit at the assay limit describes the
#'   patients rather than the assay. The boundary is reapplied at generation.
#' @param dose_term How dose enters the score model. `"factor"` gives each arm
#'   its own mean score vector and its own residual covariance. `"log"`
#'   regresses the scores on `log1p()` of the total dose, which spends one
#'   coefficient rather than one mean per arm and extrapolates to doses the
#'   study did not run, at the cost of assuming the dose-response is log-linear
#'   and that every arm has the same between-subject spread. A lower limit of
#'   quantification breaks both assumptions, which is why `"factor"` is the
#'   default.
#' @param pca_variance Cumulative variance the retained components must reach.
#' @param n_components Number of components, overriding `pca_variance`.
#' @param min_column_patients Minimum distinct patients holding an observation
#'   in a grid cell for that cell to be modelled. Defaults to the larger of 3
#'   and a tenth of the cohort.
#' @param min_arm_patients Minimum patients in every arm. An arm below it has
#'   no spread of its own to model, so its mean score vector and its covariance
#'   would describe the one or two patients in it. The function refuses rather
#'   than summarizing them; pool the arm, drop the column from `strata`, or
#'   exclude those patients before calling.
#'
#' @return A `pmx_trial_summary`.
#' @seealso [synpmx_pca_generate()], [synpmx_pca()], [pca_report()].
#' @export
#' @examples
#' data <- pmx_simulated_fixture(60)
#' roles <- pmx_roles(
#'   id = "ID", time = "TIME", nominal_time = "NTIME", dv = "DV", amt = "AMT",
#'   evid = "EVID", cmt = "CMT", dvid = "DVID", mdv = "MDV"
#' )
#' trial_summary <- synpmx_pca_summarize(data, roles)
#' trial_summary
#' pca_report(trial_summary)
synpmx_pca_summarize <- function(data, roles, seed = NULL,
                                 dose_term = c("factor", "log"),
                                 pca_variance = 0.9, n_components = NULL,
                                 min_column_patients = NULL,
                                 min_arm_patients = 3L) {
  dose_term <- match.arg(dose_term)
  if (!inherits(roles, "pmx_roles")) {
    stop("`roles` must come from `pmx_roles()`.", call. = FALSE)
  }
  data <- as.data.frame(data)
  source <- data[, intersect(.retained_role_columns(roles), names(data)),
                 drop = FALSE]
  n_source <- length(.unique_in_order(source[[roles$id]]))
  if (n_source < 10L) {
    stop("`synpmx_pca_summarize()` needs at least 10 subjects to fit a basis; ",
         "this source has ", n_source, ".", call. = FALSE)
  }
  min_column_patients <- as.integer(
    min_column_patients %||% max(3L, ceiling(0.10 * n_source))
  )

  source <- .pca_require_nominal_time(source, roles)

  # The censoring boundary is read from the source before anything replaces the
  # reported values, because after imputation there is no boundary left to find.
  censoring_source <- source

  # Fit on latent values rather than on a stack of identical boundary
  # substitutions. A column where most subjects sit exactly at the limit has
  # almost no variance, so the basis would learn the assay rather than the
  # patients, and the generated data would inherit a floor slightly above the
  # real one. `.impute_censored()` draws inside the censoring region; the
  # boundary is put back at emit. Same reasoning as `synpmx_avatar()`, Step 11.
  #
  # This is the one random step in summarizing, which is why this function takes
  # a seed of its own.
  source <- if (is.null(seed)) {
    .impute_censored(source, roles)
  } else {
    .with_local_seed(seed, .impute_censored(source, roles))
  }

  features <- .pca_features(source, roles, min_column_patients,
                            transform_source = censoring_source)
  dose <- .pca_subject_dose(source, roles, features$subjects)
  strata_key <- as.character(.subject_strata(source, roles))
  subject_group <- vapply(features$subjects, function(subject) {
    rows <- which(!is.na(source[[roles$id]]) & source[[roles$id]] == subject)
    strata_key[rows[1L]]
  }, character(1))
  .pca_require_arms(subject_group, min_arm_patients)

  fit <- .pca_fit(features, dose, subject_group, pca_variance, n_components,
                  dose_term)
  fit$min_column_patients <- min_column_patients
  arm_models <- .pca_arm_models(source, roles, fit, subject_group,
                                min_arm_patients)

  structure(list(
    basis = fit,
    scores = fit$model,
    arms = list(arms = arm_models$arms, sizes = arm_models$sizes),
    dosing = arm_models$dosing,
    visits = arm_models$visits,
    schema = .pca_schema(censoring_source, roles, fit, subject_group),
    roles = roles,
    settings = list(dose_term = dose_term, pca_variance = pca_variance,
                    n_components = n_components,
                    min_column_patients = min_column_patients,
                    min_arm_patients = min_arm_patients),
    n_source = n_source
  ), class = "pmx_trial_summary")
}

#' Generate a synthetic PMX dataset from a trial summary
#'
#' Draws new subjects from a [synpmx_pca_summarize()] trial summary. This stage reads no
#' patient data: its arguments are the model and a subject count, so everything
#' the synthetic dataset is built from is visible in the model itself.
#'
#' Each generated subject is assigned an arm, keeping each arm's share of the
#' cohort. Its scores are that arm's mean plus a fresh residual, its dose
#' schedule is the one the arm holds in common, and it attends each visit with
#' the frequency the arm attended it. No individual's schedule and no
#' individual's visit set exists in the model to be copied.
#'
#' @param trial_summary A `pmx_trial_summary` from
#'   [synpmx_pca_summarize()], or a dataset generated from one.
#' @param n_subjects Number of synthetic subjects. Defaults to the number the
#'   model was fitted on.
#' @param seed Generation seed.
#'
#' @return A data frame in the source's shape, carrying the trial summary as an
#'   attribute.
#' @seealso [synpmx_pca_summarize()], [synpmx_pca()].
#' @export
#' @examples
#' data <- pmx_simulated_fixture(60)
#' roles <- pmx_roles(
#'   id = "ID", time = "TIME", nominal_time = "NTIME", dv = "DV", amt = "AMT",
#'   evid = "EVID", cmt = "CMT", dvid = "DVID", mdv = "MDV"
#' )
#' trial_summary <- synpmx_pca_summarize(data, roles)
#' synthetic <- synpmx_pca_generate(trial_summary, seed = 1)
#' nrow(synthetic) > 0
synpmx_pca_generate <- function(trial_summary, n_subjects = NULL,
                                seed = NULL) {
  trial_summary <- .pca_trial_summary(trial_summary)
  n_subjects <- as.integer(n_subjects %||% trial_summary$n_source)
  if (!is.finite(n_subjects) || n_subjects < 1L) {
    stop("`n_subjects` must be one positive integer.", call. = FALSE)
  }
  out <- if (is.null(seed)) {
    .pca_generate(trial_summary, n_subjects)
  } else {
    .with_local_seed(seed, .pca_generate(trial_summary, n_subjects))
  }
  attr(out, "pmx_trial_summary") <- trial_summary
  attr(out, "pmx_source") <- "pca"
  out
}

#' Summarize a PMX dataset and generate a synthetic one from the summary
#'
#' A single call for [synpmx_pca_summarize()] followed by
#' [synpmx_pca_generate()]. Use the two separately to look at what the summary
#' contains before generating from it; it is on the result either way, as the
#' `pmx_trial_summary` attribute.
#'
#' Where [synpmx_avatar()] blends values from real neighbouring patients, this
#' writes out no number a patient measured. What it carries out of the source is
#' a mean, a scale, a set of principal-component loadings, one mean score vector
#' per arm, a residual covariance, and a dosing and visit model per arm.
#' [pca_report()] inventories all of it.
#'
#' `nominal_time` is required. No formal privacy claim is made.
#'
#' @param data Source PMX event data.
#' @param roles Explicit column roles from [pmx_roles()], including
#'   `nominal_time`.
#' @param n_subjects Number of synthetic subjects. Defaults to the source count.
#' @param seed Generation seed.
#' @param ... Passed to [synpmx_pca_summarize()].
#'
#' @return A data frame in the source's shape, carrying the trial summary as an
#'   attribute.
#' @seealso [synpmx_pca_summarize()], [synpmx_pca_generate()],
#'   [synpmx_avatar()], [pca_report()].
#' @export
#' @examples
#' data <- pmx_simulated_fixture(60)
#' roles <- pmx_roles(
#'   id = "ID", time = "TIME", nominal_time = "NTIME", dv = "DV", amt = "AMT",
#'   evid = "EVID", cmt = "CMT", dvid = "DVID", mdv = "MDV"
#' )
#' synthetic <- synpmx_pca(data, roles, seed = 1)
#' nrow(synthetic) > 0
synpmx_pca <- function(data, roles, n_subjects = NULL, seed = NULL, ...) {
  synpmx_pca_generate(
    synpmx_pca_summarize(data, roles, seed = seed, ...),
    n_subjects = n_subjects, seed = seed
  )
}

# Arms keep their source share of the cohort, rounded to the requested total.
.pca_assign_arms <- function(trial_summary, n_subjects) {
  sizes <- trial_summary$arms$sizes
  weights <- sizes / sum(sizes)
  counts <- as.integer(round(weights * n_subjects))
  short <- n_subjects - sum(counts)
  if (short != 0L) {
    order_index <- order(weights, decreasing = short > 0)
    for (i in seq_len(abs(short))) {
      j <- order_index[(i - 1L) %% length(counts) + 1L]
      counts[j] <- max(0L, counts[j] + sign(short))
    }
  }
  rep(trial_summary$arms$arms, times = counts)[seq_len(n_subjects)]
}

.pca_new_ids <- function(schema, n) {
  width <- max(3L, nchar(as.character(n)))
  labels <- sprintf(paste0("syn_%0", width, "d"), seq_len(n))
  switch(
    schema$id_class,
    factor = factor(labels, levels = labels),
    integer = as.integer(schema$id_offset + seq_len(n)),
    numeric = as.numeric(schema$id_offset + seq_len(n)),
    labels
  )
}

# Generation. Every argument is a summary; no patient row is in scope.
#
# Documented in `pca-algorithm.Rmd`, Steps 5 to 8.
.pca_generate <- function(trial_summary, n_subjects) {
  fit <- trial_summary$basis
  roles <- trial_summary$roles
  schema <- trial_summary$schema

  assignment <- .pca_assign_arms(trial_summary, n_subjects)
  # Each subject's schedule is simulated before their scores are drawn, because
  # the total dose it comes to is what the score model is conditioned on. On a
  # study with reductions two patients in one arm no longer receive the same
  # amount, and the exposure that follows should reflect that.
  schedules <- lapply(assignment, function(arm) {
    .pca_draw_schedule(trial_summary$dosing[[arm]])
  })
  doses <- vapply(schedules, function(schedule) {
    if (nrow(schedule)) sum(schedule$amt) else 0
  }, numeric(1))
  drawn <- .pca_draw(fit, doses, assignment)

  pieces <- vector("list", n_subjects)
  for (i in seq_len(n_subjects)) {
    arm <- assignment[[i]]
    schedule <- schedules[[i]]
    visits <- trial_summary$visits[[arm]]

    rows <- list()
    if (nrow(schedule)) {
      rows[[length(rows) + 1L]] <- data.frame(
        TIME = schedule$time, DV = NA_real_, AMT = schedule$amt, EVID = 1L,
        CMT = if (is.null(schema$cmt_dose)) NA else schema$cmt_dose,
        DVID = NA_character_, stringsAsFactors = FALSE
      )
    }
    attended <- stats::runif(length(visits$cells)) < visits$probability
    for (index in which(attended)) {
      column <- visits$cells[[index]]
      member <- fit$members[[column]]
      value <- .inverse_dv(drawn[i, column], fit$transforms[[member$endpoint]])
      value <- .snap_endpoint_values(value,
                                     schema$endpoint_specs[[member$endpoint]])
      rows[[length(rows) + 1L]] <- data.frame(
        TIME = member$time, DV = value, AMT = 0, EVID = 0L,
        CMT = if (is.null(schema$cmt_obs[[member$endpoint]])) NA else
          schema$cmt_obs[[member$endpoint]],
        DVID = member$endpoint, stringsAsFactors = FALSE
      )
    }
    if (!length(rows)) next
    piece <- do.call(rbind, rows)
    piece <- piece[order(piece$TIME, piece$EVID == 0L), , drop = FALSE]
    piece$.subject <- i
    piece$.arm <- arm
    pieces[[i]] <- piece
  }

  frame <- do.call(rbind, pieces[!vapply(pieces, is.null, logical(1))])
  rownames(frame) <- NULL

  out <- data.frame(row.names = seq_len(nrow(frame)))
  ids <- .pca_new_ids(schema, n_subjects)
  out[[roles$id]] <- ids[frame$.subject]
  out[[roles$time]] <- frame$TIME
  out[[roles$nominal_time]] <- frame$TIME
  out[[roles$dv]] <- frame$DV
  if (!is.null(roles$amt)) out[[roles$amt]] <- frame$AMT
  out[[roles$evid]] <- frame$EVID
  if (!is.null(roles$cmt)) out[[roles$cmt]] <- frame$CMT
  if (!is.null(roles$dvid)) {
    for (column in roles$dvid) out[[column]] <- frame$DVID
  }
  if (!is.null(roles$mdv)) out[[roles$mdv]] <- as.integer(frame$EVID != 0L)
  if (!is.null(roles$cens)) out[[roles$cens]] <- 0L
  if (!is.null(roles$limit)) out[[roles$limit]] <- NA_real_
  if (!is.null(roles$rate)) out[[roles$rate]] <- 0
  if (!is.null(roles$occasion)) {
    out[[roles$occasion]] <- unlist(lapply(
      split(frame, frame$.subject),
      function(part) pmax(1L, cumsum(part$EVID != 0L))
    ), use.names = FALSE)
  }
  if (!is.null(roles$tad)) out[[roles$tad]] <- .derived_tad(out, roles)
  if (!is.null(roles$assigned_dose)) {
    out[[roles$assigned_dose]] <- doses[frame$.subject]
  }
  for (column in schema$carried) {
    lookup <- stats::setNames(
      lapply(trial_summary$arms$arms,
             function(arm) schema$arm_values[[arm]][[column]]),
      trial_summary$arms$arms
    )
    out[[column]] <- .restore_column(
      unlist(lookup[frame$.arm], use.names = FALSE), schema$prototypes[[column]]
    )
  }
  for (covariate in roles$covariates) {
    prototype <- schema$prototypes[[covariate]]
    continuous <- paste0("cov_", covariate)
    if (continuous %in% fit$columns) {
      out[[covariate]] <- drawn[frame$.subject, continuous]
    } else {
      level_columns <- fit$columns[
        vapply(fit$members, function(m) identical(m$covariate, covariate) &&
                 !is.null(m$level), logical(1))
      ]
      if (!length(level_columns)) next
      levels_named <- vapply(fit$members[level_columns],
                             function(m) m$level, character(1))
      picked <- levels_named[max.col(drawn[, level_columns, drop = FALSE],
                                     ties.method = "first")]
      out[[covariate]] <- .restore_column(picked[frame$.subject], prototype)
    }
  }
  # The drawn value is the latent one: what the subject would have measured with
  # no assay limit. DV, CENS and LIMIT are reconstructed from it together, after
  # any discrete endpoint has been snapped, so the three always agree.
  #
  # Documented in `pca-algorithm.Rmd`, Step 7.
  if (!is.null(roles$cens)) {
    for (endpoint_name in fit$endpoints) {
      censoring <- schema$censoring[[endpoint_name]]
      if (is.null(censoring)) next
      rows <- which(frame$EVID == 0L & frame$DVID == endpoint_name)
      if (!length(rows)) next
      out <- .censor_latent(out, rows, out[[roles$dv]][rows], roles,
                            public = censoring)
    }
  }

  out <- out[, intersect(schema$columns, names(out)), drop = FALSE]
  out <- out[order(out[[roles$id]], out[[roles$time]],
                   out[[roles$evid]] == 0L), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# `.subject_strata()` joins the strata columns with a control character, which
# is right for a key and wrong in a table someone reads. The arm keeps its key
# internally and is labelled with the columns joined readably.
.pca_arm_label <- function(arm) {
  gsub("\r", " / ", arm, fixed = TRUE)
}

# Accept a generated dataset or the trial summary itself.
.pca_trial_summary <- function(x) {
  out <- attr(x, "pmx_trial_summary") %||% x
  if (!inherits(out, "pmx_trial_summary")) {
    stop("Expected a dataset from `synpmx_pca()` or a `pmx_trial_summary` ",
         "from `synpmx_pca_summarize()`.", call. = FALSE)
  }
  out
}

.pca_basis <- function(x) .pca_trial_summary(x)$basis

#' The planned dose schedule each arm was generated from
#'
#' One row per arm and cycle, giving the nominal time and the amount the arm was
#' planned to receive there. The planned amount is read per cycle, so a
#' protocol-prescribed escalation is part of the plan rather than a departure
#' from it.
#'
#' This is what every generated patient in the arm starts from.
#' [pca_dose_rates()] gives the three hazards that then move them off it:
#' reductions, interruptions and discontinuation. On a study with no dose
#' modifications those rates are zero and this schedule is what every patient
#' receives.
#'
#' The grid stops at the last cycle `min_arm_patients` of the arm's patients
#' reached, so a treatment duration only one patient had cannot be generated.
#'
#' @param x A dataset from [synpmx_pca()], or its trial summary.
#'
#' @return A data frame with `arm`, `cycle`, `time` and `planned_amt`.
#' @seealso [synpmx_pca()], [pca_dose_rates()], [pca_visits()].
#' @export
#' @examples
#' data <- pmx_simulated_fixture(60)
#' roles <- pmx_roles(
#'   id = "ID", time = "TIME", nominal_time = "NTIME", dv = "DV", amt = "AMT",
#'   evid = "EVID", cmt = "CMT", dvid = "DVID", mdv = "MDV"
#' )
#' head(pca_dosing(synpmx_pca(data, roles, seed = 1)))
pca_dosing <- function(x) {
  trial_summary <- .pca_trial_summary(x)
  out <- do.call(rbind, lapply(trial_summary$arms$arms, function(arm) {
    entry <- trial_summary$dosing[[arm]]
    planned <- entry$planned
    if (!nrow(planned)) return(NULL)
    data.frame(arm = .pca_arm_label(arm), cycle = planned$cycle,
               time = planned$time, planned_amt = planned$amt,
               stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}

#' The dose-modification rates each arm was generated from
#'
#' One row per arm, giving the three discrete-time hazards that turn a planned
#' schedule into the schedule a patient actually received, and the dose ladder
#' reductions move down.
#'
#' A study where nobody reduces, skips or stops early has all three rates at
#' zero and a single level, and every generated patient then receives the
#' planned schedule exactly. A study with dose modifications --- oncology being
#' the usual case --- has non-zero rates, and generated patients differ from one
#' another in the same way and to the same degree the source patients did.
#'
#' `source_doses` and `distinct` describe what the arm actually contained, so a
#' generated dataset can be checked against them: `source_doses` is the mean
#' number of dosing events per patient, and `distinct` is how many different
#' schedules the arm held.
#'
#' @param x A dataset from [synpmx_pca()], or its trial summary.
#'
#' @return A data frame with `arm`, `planned_cycles`, `levels`,
#'   `discontinuation`, `interruption`, `reduction`, `patients`,
#'   `source_doses` and `distinct`.
#' @seealso [synpmx_pca()], [pca_dosing()], [pca_report()].
#' @export
#' @examples
#' data <- pmx_simulated_fixture(60)
#' roles <- pmx_roles(
#'   id = "ID", time = "TIME", nominal_time = "NTIME", dv = "DV", amt = "AMT",
#'   evid = "EVID", cmt = "CMT", dvid = "DVID", mdv = "MDV"
#' )
#' pca_dose_rates(synpmx_pca(data, roles, seed = 1))
pca_dose_rates <- function(x) {
  trial_summary <- .pca_trial_summary(x)
  out <- do.call(rbind, lapply(trial_summary$arms$arms, function(arm) {
    entry <- trial_summary$dosing[[arm]]
    data.frame(
      arm = .pca_arm_label(arm),
      planned_cycles = nrow(entry$planned),
      levels = paste(format(entry$levels, trim = TRUE), collapse = ", "),
      discontinuation = round(entry$discontinuation, 4),
      interruption = round(entry$interruption, 4),
      reduction = round(entry$reduction, 4),
      patients = entry$patients,
      source_doses = round(entry$source_doses, 1),
      distinct = entry$distinct,
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

#' The visit model each arm was generated from
#'
#' One row per arm, endpoint and modelled nominal time, giving the probability
#' that a generated subject in that arm has an observation there. It is the
#' fraction of the arm's patients who did, so attendance is drawn per visit
#' rather than a real patient's visit set being reused.
#'
#' @param x A dataset from [synpmx_pca()], or its trial summary.
#'
#' @return A data frame with `arm`, `endpoint`, `time`, `probability` and
#'   `patients`, the last being how many patients across the study hold that
#'   cell at all.
#' @seealso [synpmx_pca()], [pca_dosing()], [pca_report()].
#' @export
#' @examples
#' data <- pmx_simulated_fixture(60)
#' roles <- pmx_roles(
#'   id = "ID", time = "TIME", nominal_time = "NTIME", dv = "DV", amt = "AMT",
#'   evid = "EVID", cmt = "CMT", dvid = "DVID", mdv = "MDV"
#' )
#' head(pca_visits(synpmx_pca(data, roles, seed = 1)))
pca_visits <- function(x) {
  trial_summary <- .pca_trial_summary(x)
  fit <- trial_summary$basis
  out <- do.call(rbind, lapply(trial_summary$arms$arms, function(arm) {
    entry <- trial_summary$visits[[arm]]
    members <- fit$members[entry$cells]
    data.frame(
      arm = .pca_arm_label(arm),
      endpoint = vapply(members, function(m) m$endpoint, character(1)),
      time = vapply(members, function(m) as.numeric(m$time), numeric(1)),
      probability = entry$probability,
      patients = vapply(members, function(m) as.numeric(m$patients),
                        numeric(1)),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

#' @export
print.pmx_trial_summary <- function(x, ...) {
  fit <- x$basis
  cat("A trial summary, from synpmx_pca_summarize()\n\n")
  cat("  fitted on   ", x$n_source, "patients,",
      length(x$arms$arms), "arm(s):",
      paste(sprintf("%s (%d)", .pca_arm_label(x$arms$arms),
                    as.integer(x$arms$sizes[x$arms$arms])),
            collapse = ", "), "\n")
  cat("  endpoints   ",
      paste(sprintf("%s (%d visits modelled)", fit$endpoints,
                    vapply(fit$endpoints, function(ep) {
                      sum(vapply(fit$members[fit$kinds == "endpoint_cell"],
                                 function(m) identical(m$endpoint, ep),
                                 logical(1)))
                    }, integer(1))), collapse = ", "), "\n")
  covariates <- x$roles$covariates
  cat("  covariates  ",
      if (length(covariates)) paste(covariates, collapse = ", ") else "none",
      "\n")
  share <- cumsum(fit$variance / sum(fit$variance))[fit$k]
  cat("  components  ", fit$k, sprintf("(%.0f%% of variance)", 100 * share),
      "\n")
  cat("  dose term   ", x$settings$dose_term, "\n")
  cycles <- vapply(x$dosing, function(d) nrow(d$planned), integer(1))
  varies <- vapply(x$dosing, function(d) {
    d$discontinuation > 0 || d$interruption > 0 || d$reduction > 0
  }, logical(1))
  cat("  dosing      ",
      sprintf("%d planned cycle(s) per arm", stats::median(cycles)),
      if (any(varies)) {
        sprintf("| %d of %d arm(s) reduce, skip or stop early",
                sum(varies), length(varies))
      } else "| no reductions, interruptions or early stops", "\n\n")
  cat("synpmx_pca_generate() reads this object and nothing else.",
      "To look inside it:\n")
  cat("  pca_report()      what it read out of the source data\n")
  cat("  pca_dosing()      the planned dose schedule, per arm\n")
  cat("  pca_dose_rates()  reduction, interruption and discontinuation\n")
  cat("  pca_visits()      the probability of a visit, per arm\n")
  cat("  pca_components()  the loadings, over time\n")
  invisible(x)
}

.pca_score_numbers <- function(fit) {
  if (identical(fit$dose_term, "log")) {
    nrow(fit$model$coefficients) * fit$k
  } else {
    length(fit$model$means) * fit$k
  }
}

#' What the PCA fit read out of the source data
#'
#' One row per released quantity: what it is, how many numbers it holds, and
#' the smallest number of patients standing behind any one of them. That last
#' column is where disclosure risk sits. A grid cell or a covariate mean is
#' backed by the whole cohort, while a rare covariate level can be backed by a
#' single patient.
#'
#' @param x A dataset from [synpmx_pca()], or the fit itself.
#'
#' @return A `pca_report` data frame.
#' @seealso [synpmx_pca()], [pca_components()].
#' @export
#' @examples
#' data <- pmx_simulated_fixture(60)
#' pca_report(synpmx_pca(data, pmx_generated_roles(), seed = 1))
pca_report <- function(x) {
  trial_summary <- .pca_trial_summary(x)
  fit <- trial_summary$basis
  dosing_numbers <- sum(vapply(trial_summary$dosing, function(d) {
    2L * nrow(d$planned) + length(d$levels) + 3L
  }, integer(1)))
  visit_numbers <- sum(vapply(trial_summary$visits,
                              function(v) length(v$probability), integer(1)))
  arm_numbers <- length(trial_summary$arms$arms) *
    length(trial_summary$schema$carried)
  censoring_numbers <- sum(vapply(trial_summary$schema$censoring,
                                  function(c) length(unlist(c)), integer(1)))
  cells <- which(fit$kinds == "endpoint_cell")
  cell_patients <- if (length(cells)) {
    min(vapply(fit$members[cells], function(m) as.numeric(m$patients),
               numeric(1)))
  } else NA_real_
  p <- length(fit$columns)
  rows <- data.frame(
    quantity = c("visit grid", "feature centers", "feature scales",
                 "loadings", "score means", "score covariance",
                 "endpoint transforms", "assay limits", "dosing model",
                 "visit model", "arm constants"),
    what = c(
      "Nominal times modelled, per endpoint",
      "Mean of each grid cell and covariate",
      "Standard deviation of the same",
      "Component loadings on each feature",
      if (identical(fit$dose_term, "log")) "Score regression on log dose" else
        "Mean score vector, per arm",
      "Residual covariance between components",
      "Log or identity, per endpoint",
      "Censoring boundary, per endpoint",
      "Planned cycles, the dose ladder, and three rates, per arm",
      "Probability of a visit, per arm, endpoint and time",
      "Strata and kept columns, one value per arm"
    ),
    numbers = c(
      length(cells), p, p, p * fit$k, .pca_score_numbers(fit),
      fit$k * fit$k, length(fit$transforms), censoring_numbers,
      dosing_numbers, visit_numbers, arm_numbers
    ),
    min_patients = c(
      cell_patients, cell_patients, cell_patients,
      fit$n_source, fit$min_group, fit$min_group, fit$n_source,
      fit$n_source, fit$min_group, fit$min_group, fit$min_group
    ),
    stringsAsFactors = FALSE
  )
  structure(rows, class = c("pca_report", "data.frame"),
            components = fit$k, subjects = fit$n_source)
}

#' @export
print.pca_report <- function(x, ...) {
  cat("What the PCA fit read out of the source data\n\n")
  cat("  subjects:", attr(x, "subjects"),
      " components retained:", attr(x, "components"), "\n\n")
  print(as.data.frame(x), row.names = FALSE)
  invisible(x)
}

#' Every feature the components are built on
#'
#' One row per column of the matrix `synpmx_pca_summarize()` decomposed: one
#' per baseline covariate, and one per endpoint per retained nominal time. This
#' is the grid the whole method sits on, so it is where to look first when a
#' generated dataset is missing a visit or a covariate.
#'
#' `center` and `scale` are the column's mean and standard deviation on the
#' modelling scale, which is the log scale for a positive endpoint --- see
#' `transform` in the same row. `patients` is how many subjects hold an
#' observation in that cell; cells held by fewer than `min_column_patients` were
#' dropped rather than modelled, so they do not appear here at all.
#'
#' @param x A dataset from [synpmx_pca()], or its trial summary.
#'
#' @return A data frame with `feature`, `kind`, `endpoint`, `time`, `covariate`,
#'   `level`, `patients`, `center`, `scale` and `transform`. Marked
#'   `"restricted_not_releasable"`: the counts and moments are read from real
#'   data.
#' @seealso [synpmx_pca_summarize()], [pca_components()], [pca_report()].
#' @export
#' @examples
#' data <- pmx_simulated_fixture(60)
#' roles <- pmx_roles(
#'   id = "ID", time = "TIME", nominal_time = "NTIME", dv = "DV", amt = "AMT",
#'   evid = "EVID", cmt = "CMT", dvid = "DVID", mdv = "MDV"
#' )
#' head(pca_features(synpmx_pca_summarize(data, roles)))
pca_features <- function(x) {
  fit <- .pca_basis(x)
  field <- function(name, empty) {
    vapply(fit$members, function(member) {
      value <- member[[name]]
      if (is.null(value)) empty else value[[1L]]
    }, empty)
  }
  endpoints <- field("endpoint", NA_character_)
  out <- data.frame(
    feature = fit$columns,
    kind = unname(fit$kinds),
    endpoint = endpoints,
    time = field("time", NA_real_),
    covariate = field("covariate", NA_character_),
    level = field("level", NA_character_),
    patients = field("patients", NA_real_),
    center = unname(fit$centers),
    scale = unname(fit$scales),
    transform = vapply(endpoints, function(endpoint) {
      if (is.na(endpoint)) NA_character_ else
        fit$transforms[[endpoint]]$method
    }, character(1)),
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  .mark_release(out, "restricted_not_releasable")
}

#' The score model each arm is generated from
#'
#' One row per arm and retained component, giving the mean score that arm is
#' centred on and the residual standard deviation it is scattered by. A
#' generated subject's scores are their arm's `mean` plus a draw whose spread is
#' `sd`, so these two columns are the whole of the between-subject variability
#' the synthetic data will have.
#'
#' The `sd` values differ by arm on purpose. An arm sitting on an assay limit is
#' genuinely tighter than one well above it, and giving every arm the pooled
#' spread would smear the low arms upward.
#'
#' Both `dose_term` settings report the same shape. Under `"factor"` the mean is
#' the arm's own; under `"log"` it is the regression evaluated at that arm's
#' planned total dose, and `sd` is the shared residual.
#'
#' @param x A dataset from [synpmx_pca()], or its trial summary.
#'
#' @return A data frame with `arm`, `component`, `mean` and `sd`. Marked
#'   `"restricted_not_releasable"`.
#' @seealso [synpmx_pca_summarize()], [pca_components()], [pca_report()].
#' @export
#' @examples
#' data <- pmx_simulated_fixture(60)
#' roles <- pmx_roles(
#'   id = "ID", time = "TIME", nominal_time = "NTIME", dv = "DV", amt = "AMT",
#'   evid = "EVID", cmt = "CMT", dvid = "DVID", mdv = "MDV"
#' )
#' pca_scores(synpmx_pca_summarize(data, roles))
pca_scores <- function(x) {
  trial_summary <- .pca_trial_summary(x)
  fit <- trial_summary$basis
  scores <- trial_summary$scores
  arms <- trial_summary$arms$arms

  out <- do.call(rbind, lapply(arms, function(arm) {
    if (identical(scores$kind, "log")) {
      planned <- trial_summary$dosing[[arm]]$planned
      dose <- if (nrow(planned)) sum(planned$amt) else 0
      design <- cbind(intercept = 1, log_dose = log1p(dose))
      mean <- as.numeric(design %*% scores$coefficients)
      spread <- sqrt(diag(scores$covariance))
    } else {
      mean <- as.numeric(scores$means[[arm]])
      spread <- sqrt(diag(scores$covariances[[arm]]))
    }
    data.frame(
      arm = .pca_arm_label(arm),
      component = paste0("PC", seq_len(fit$k)),
      mean = mean, sd = as.numeric(spread),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  .mark_release(out, "restricted_not_releasable")
}

#' Component loadings over time, and the variance each component explains
#'
#' The loadings are what makes a principal component readable. Plotted against
#' time rather than tabulated, a component that is flat and positive is overall
#' magnitude and one that crosses zero separates early from late.
#'
#' @param x A dataset from [synpmx_pca()], or the fit itself.
#'
#' @return A data frame with one row per component and retained grid cell,
#'   carrying `variance_explained` as an attribute.
#' @seealso [synpmx_pca()], [pca_report()].
#' @export
#' @examples
#' data <- pmx_simulated_fixture(60)
#' head(pca_components(synpmx_pca(data, pmx_generated_roles(), seed = 1)))
pca_components <- function(x) {
  fit <- .pca_basis(x)
  features <- as.data.frame(pca_features(x))
  out <- do.call(rbind, lapply(seq_len(fit$k), function(j) {
    cbind(features[, c("feature", "kind", "endpoint", "time", "covariate",
                       "patients")],
          component = paste0("PC", j),
          loading = unname(fit$rotation[, j]),
          stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  share <- fit$variance / sum(fit$variance)
  attr(out, "variance_explained") <- data.frame(
    component = paste0("PC", seq_len(fit$k)),
    variance_explained = share[seq_len(fit$k)],
    cumulative = cumsum(share)[seq_len(fit$k)],
    stringsAsFactors = FALSE
  )
  out
}
