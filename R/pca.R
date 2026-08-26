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
.pca_features <- function(source, roles, min_column_patients) {
  subjects <- .unique_in_order(source[[roles$id]])
  n <- length(subjects)
  subject_rows <- lapply(subjects, function(subject) {
    !is.na(source[[roles$id]]) & source[[roles$id]] == subject
  })
  observed <- .observation_rows(source, roles, require_present = TRUE)
  endpoint <- .endpoint(source, roles)
  aligned <- .aligned_time(source, roles)
  endpoints <- sort(unique(endpoint[observed]))

  transforms <- stats::setNames(lapply(endpoints, function(ep) {
    .choose_transform(suppressWarnings(as.numeric(
      source[[roles$dv]][observed & endpoint == ep]
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

# The dosing model and the visit model, one of each per arm.
#
# Both are summaries of the arm rather than facts about a patient. The dosing
# model is the schedule the arm holds in common -- the modal set of (time,
# amount) pairs among its patients -- reported with the share of the arm that
# holds it, so a schedule chosen by a bare plurality is visible as one. The
# visit model is, per endpoint and per retained nominal time, the fraction of
# the arm that has an observation there.
#
# Nothing here is an individual's schedule or an individual's visit set, which
# is why neither can leave. The cost is the variety: an arm generates one
# schedule, however many its patients had.
#
# Documented in `pca-algorithm.Rmd`, Step 6.
.pca_arm_models <- function(source, roles, fit, subject_group) {
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

    schedules <- lapply(member_rows, function(rows) {
      selected <- rows & dose_rows
      order_by <- order(aligned[selected])
      data.frame(time = aligned[selected][order_by],
                 amt = amount[selected][order_by])
    })
    keys <- vapply(schedules, function(schedule) {
      paste(sprintf("%.6g", schedule$time), sprintf("%.6g", schedule$amt),
            collapse = "|")
    }, character(1))
    counts <- sort(table(keys), decreasing = TRUE)
    modal <- names(counts)[1L]
    schedule <- schedules[[which(keys == modal)[1L]]]
    dosing[[arm]] <- list(
      schedule = schedule,
      patients = as.integer(counts[[1L]]),
      distinct = length(counts),
      share = as.numeric(counts[[1L]]) / length(members)
    )

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

#' Summarize a PMX dataset into a principal-component model
#'
#' The only stage that reads patient data. Reduces each subject's trajectories
#' and baseline covariates to principal-component scores, models those scores
#' against the arm, and fits a dosing model and a visit model per arm. The
#' returned object holds nothing but summaries: no patient row survives it.
#'
#' Run this, look at what it produced with [pmx_pca_report()],
#' [pmx_pca_dosing()], [pmx_pca_visits()] and [pmx_pca_components()], then pass
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
#' @return A `pmx_pca_model`.
#' @seealso [synpmx_pca_generate()], [synpmx_pca()], [pmx_pca_report()].
#' @export
#' @examples
#' data <- pmx_simulated_fixture(60)
#' roles <- pmx_roles(
#'   id = "ID", time = "TIME", nominal_time = "NTIME", dv = "DV", amt = "AMT",
#'   evid = "EVID", cmt = "CMT", dvid = "DVID", mdv = "MDV"
#' )
#' model <- synpmx_pca_summarize(data, roles)
#' model
#' pmx_pca_report(model)
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

  features <- .pca_features(source, roles, min_column_patients)
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
  arm_models <- .pca_arm_models(source, roles, fit, subject_group)

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
  ), class = "pmx_pca_model")
}

#' Generate a synthetic PMX dataset from a principal-component model
#'
#' Draws new subjects from a [synpmx_pca_summarize()] model. This stage reads no
#' patient data: its arguments are the model and a subject count, so everything
#' the synthetic dataset is built from is visible in the model itself.
#'
#' Each generated subject is assigned an arm, keeping each arm's share of the
#' cohort. Its scores are that arm's mean plus a fresh residual, its dose
#' schedule is the one the arm holds in common, and it attends each visit with
#' the frequency the arm attended it. No individual's schedule and no
#' individual's visit set exists in the model to be copied.
#'
#' @param model A `pmx_pca_model` from [synpmx_pca_summarize()], or a dataset
#'   generated from one.
#' @param n_subjects Number of synthetic subjects. Defaults to the number the
#'   model was fitted on.
#' @param seed Generation seed.
#'
#' @return A data frame in the source's shape, carrying the model as an
#'   attribute.
#' @seealso [synpmx_pca_summarize()], [synpmx_pca()].
#' @export
#' @examples
#' data <- pmx_simulated_fixture(60)
#' roles <- pmx_roles(
#'   id = "ID", time = "TIME", nominal_time = "NTIME", dv = "DV", amt = "AMT",
#'   evid = "EVID", cmt = "CMT", dvid = "DVID", mdv = "MDV"
#' )
#' model <- synpmx_pca_summarize(data, roles)
#' synthetic <- synpmx_pca_generate(model, seed = 1)
#' nrow(synthetic) > 0
synpmx_pca_generate <- function(model, n_subjects = NULL, seed = NULL) {
  model <- .pca_model(model)
  n_subjects <- as.integer(n_subjects %||% model$n_source)
  if (!is.finite(n_subjects) || n_subjects < 1L) {
    stop("`n_subjects` must be one positive integer.", call. = FALSE)
  }
  out <- if (is.null(seed)) {
    .pca_generate(model, n_subjects)
  } else {
    .with_local_seed(seed, .pca_generate(model, n_subjects))
  }
  attr(out, "pmx_pca_model") <- model
  attr(out, "pmx_pca_fit") <- model$basis
  attr(out, "pmx_source") <- "pca"
  out
}

#' Summarize a PMX dataset and generate a synthetic one from the summary
#'
#' A single call for [synpmx_pca_summarize()] followed by
#' [synpmx_pca_generate()]. Use the two separately to look at what the summary
#' contains before generating from it; the model is on the result either way, as
#' the `pmx_pca_model` attribute.
#'
#' Where [synpmx_avatar()] blends values from real neighbouring patients, this
#' writes out no number a patient measured. What it carries out of the source is
#' a mean, a scale, a set of principal-component loadings, one mean score vector
#' per arm, a residual covariance, and a dosing and visit model per arm.
#' [pmx_pca_report()] inventories all of it.
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
#' @return A data frame in the source's shape, carrying the model as an
#'   attribute.
#' @seealso [synpmx_pca_summarize()], [synpmx_pca_generate()],
#'   [synpmx_avatar()], [pmx_pca_report()].
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
.pca_assign_arms <- function(model, n_subjects) {
  sizes <- model$arms$sizes
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
  rep(model$arms$arms, times = counts)[seq_len(n_subjects)]
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
.pca_generate <- function(model, n_subjects) {
  fit <- model$basis
  roles <- model$roles
  schema <- model$schema

  assignment <- .pca_assign_arms(model, n_subjects)
  doses <- vapply(assignment, function(arm) {
    schedule <- model$dosing[[arm]]$schedule
    if (nrow(schedule)) sum(schedule$amt) else 0
  }, numeric(1))
  drawn <- .pca_draw(fit, doses, assignment)

  pieces <- vector("list", n_subjects)
  for (i in seq_len(n_subjects)) {
    arm <- assignment[[i]]
    schedule <- model$dosing[[arm]]$schedule
    visits <- model$visits[[arm]]

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
    out[[roles$assigned_dose]] <- unname(vapply(frame$.arm, function(arm) {
      schedule <- model$dosing[[arm]]$schedule
      if (nrow(schedule)) sum(schedule$amt) else 0
    }, numeric(1)))
  }
  for (column in schema$carried) {
    lookup <- stats::setNames(
      lapply(model$arms$arms,
             function(arm) schema$arm_values[[arm]][[column]]),
      model$arms$arms
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

# Accept a generated dataset, the model, or the basis itself.
.pca_model <- function(x) {
  model <- attr(x, "pmx_pca_model") %||% x
  if (!inherits(model, "pmx_pca_model")) {
    stop("Expected a dataset from `synpmx_pca()` or a `pmx_pca_model`.",
         call. = FALSE)
  }
  model
}

.pca_basis <- function(x) {
  attr(x, "pmx_pca_fit") %||%
    (if (inherits(x, "pmx_pca_model")) x$basis else x)
}

#' The dosing model each arm was generated from
#'
#' One row per arm and dose, giving the time and the amount. This is the
#' schedule the arm holds in common: the modal set of times and amounts among
#' its patients, never an individual's. `share` is the fraction of the arm
#' holding it, so a schedule picked by a bare plurality is visible as one, and
#' `distinct` is how many schedules that arm actually contained.
#'
#' A study recording its dose times as actuals rather than as planned times
#' holds close to one schedule per patient. The generated data then holds one
#' per arm, which is the mechanism working and a real loss of variety.
#'
#' @param x A dataset from [synpmx_pca()], or its model.
#'
#' @return A data frame with `arm`, `dose`, `time`, `amt`, `share`,
#'   `patients` and `distinct`.
#' @seealso [synpmx_pca()], [pmx_pca_visits()], [pmx_pca_report()].
#' @export
#' @examples
#' data <- pmx_simulated_fixture(60)
#' roles <- pmx_roles(
#'   id = "ID", time = "TIME", nominal_time = "NTIME", dv = "DV", amt = "AMT",
#'   evid = "EVID", cmt = "CMT", dvid = "DVID", mdv = "MDV"
#' )
#' pmx_pca_dosing(synpmx_pca(data, roles, seed = 1))
pmx_pca_dosing <- function(x) {
  model <- .pca_model(x)
  out <- do.call(rbind, lapply(model$arms$arms, function(arm) {
    entry <- model$dosing[[arm]]
    schedule <- entry$schedule
    if (!nrow(schedule)) return(NULL)
    data.frame(arm = .pca_arm_label(arm), dose = seq_len(nrow(schedule)),
               time = schedule$time, amt = schedule$amt,
               share = entry$share, patients = entry$patients,
               distinct = entry$distinct, stringsAsFactors = FALSE)
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
#' @param x A dataset from [synpmx_pca()], or its model.
#'
#' @return A data frame with `arm`, `endpoint`, `time`, `probability` and
#'   `patients`, the last being how many patients across the study hold that
#'   cell at all.
#' @seealso [synpmx_pca()], [pmx_pca_dosing()], [pmx_pca_report()].
#' @export
#' @examples
#' data <- pmx_simulated_fixture(60)
#' roles <- pmx_roles(
#'   id = "ID", time = "TIME", nominal_time = "NTIME", dv = "DV", amt = "AMT",
#'   evid = "EVID", cmt = "CMT", dvid = "DVID", mdv = "MDV"
#' )
#' head(pmx_pca_visits(synpmx_pca(data, roles, seed = 1)))
pmx_pca_visits <- function(x) {
  model <- .pca_model(x)
  fit <- model$basis
  out <- do.call(rbind, lapply(model$arms$arms, function(arm) {
    entry <- model$visits[[arm]]
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
print.pmx_pca_model <- function(x, ...) {
  fit <- x$basis
  cat("A synpmx PCA model\n\n")
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
  doses <- vapply(x$dosing, function(d) nrow(d$schedule), integer(1))
  shares <- vapply(x$dosing, function(d) d$share, numeric(1))
  cat("  dosing      ", sprintf("%d dose(s) per arm", stats::median(doses)),
      sprintf("| shared by %.0f%%-%.0f%% of each arm",
              100 * min(shares), 100 * max(shares)), "\n\n")
  cat("Generation reads this object and nothing else. To look inside it:\n")
  cat("  pmx_pca_report(model)      what it read out of the source data\n")
  cat("  pmx_pca_dosing(model)      the dose schedule each arm shares\n")
  cat("  pmx_pca_visits(model)      the probability of a visit, per arm\n")
  cat("  pmx_pca_components(model)  the loadings, over time\n")
  invisible(x)
}

.pca_model_numbers <- function(fit) {
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
#' @return A `pmx_pca_report` data frame.
#' @seealso [synpmx_pca()], [pmx_pca_components()].
#' @export
#' @examples
#' data <- pmx_simulated_fixture(60)
#' pmx_pca_report(synpmx_pca(data, pmx_generated_roles(), seed = 1))
pmx_pca_report <- function(x) {
  fit <- .pca_basis(x)
  model <- tryCatch(.pca_model(x), error = function(e) NULL)
  dosing_numbers <- if (is.null(model)) NA_integer_ else
    2L * sum(vapply(model$dosing, function(d) nrow(d$schedule), integer(1)))
  visit_numbers <- if (is.null(model)) NA_integer_ else
    sum(vapply(model$visits, function(v) length(v$probability), integer(1)))
  arm_numbers <- if (is.null(model)) NA_integer_ else
    length(model$arms$arms) * length(model$schema$carried)
  censoring_numbers <- if (is.null(model)) NA_integer_ else
    sum(vapply(model$schema$censoring, function(c) length(unlist(c)), integer(1)))
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
      "Dose times and amounts each arm shares",
      "Probability of a visit, per arm, endpoint and time",
      "Strata and kept columns, one value per arm"
    ),
    numbers = c(
      length(cells), p, p, p * fit$k, .pca_model_numbers(fit),
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
  structure(rows, class = c("pmx_pca_report", "data.frame"),
            components = fit$k, subjects = fit$n_source)
}

#' @export
print.pmx_pca_report <- function(x, ...) {
  cat("What the PCA fit read out of the source data\n\n")
  cat("  subjects:", attr(x, "subjects"),
      " components retained:", attr(x, "components"), "\n\n")
  print(as.data.frame(x), row.names = FALSE)
  invisible(x)
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
#' @seealso [synpmx_pca()], [pmx_pca_report()].
#' @export
#' @examples
#' data <- pmx_simulated_fixture(60)
#' head(pmx_pca_components(synpmx_pca(data, pmx_generated_roles(), seed = 1)))
pmx_pca_components <- function(x) {
  fit <- .pca_basis(x)
  cells <- which(fit$kinds == "endpoint_cell")
  grid <- do.call(rbind, lapply(cells, function(column) {
    member <- fit$members[[column]]
    data.frame(feature = fit$columns[[column]], endpoint = member$endpoint,
               time = member$time, patients = member$patients,
               stringsAsFactors = FALSE)
  }))
  out <- do.call(rbind, lapply(seq_len(fit$k), function(j) {
    cbind(grid, component = paste0("PC", j),
          loading = fit$rotation[cells, j],
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
