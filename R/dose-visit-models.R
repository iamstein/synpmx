# The dosing model, the visit model, and the schedule draw ------------------
#
# Shared by every generator that has to reproduce a study's dosing and
# attendance rather than model it. `synpmx_pca()` was the first caller and
# these functions moved out of `R/pca.R` unchanged, so that a second generator
# calls the same code rather than a copy of it. Nothing here reads a fitted
# basis, a set of components, or anything else specific to one method: the
# inputs are the source rows, the arm each subject belongs to, and the grid
# cells attendance is measured on.
#
# Documented in `pca-algorithm.Rmd`, Step 6. A second generator's algorithm
# document describes the same functions from its own side.

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
.dose_model <- function(member_rows, aligned, amount, dose_rows, floor) {
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
.draw_schedule <- function(dosing) {
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
# `cells` is the grid attendance is measured on: one row per endpoint and
# nominal time, carrying the `index` the caller knows that cell by. Passing the
# grid rather than a fitted object is what keeps this callable by a generator
# that has no components.
#
# Both are summaries of the arm rather than facts about a patient. The dosing
# model is built above: a planned schedule and three rates. The visit model is,
# per endpoint and per retained nominal time, the fraction of the arm that has
# an observation there, so attendance is drawn per visit rather than a real
# patient's set of attended visits being reused.
.arm_models <- function(source, roles, cells, subject_group, floor) {
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

  arms <- unique(subject_group)
  index <- .named(cells$index, cells$name)
  dosing <- list()
  visits <- list()
  sizes <- integer()

  for (arm in arms) {
    members <- subjects[subject_group == arm]
    member_rows <- lapply(members, function(subject) {
      !is.na(source[[roles$id]]) & source[[roles$id]] == subject
    })
    sizes[arm] <- length(members)

    dosing[[arm]] <- .dose_model(member_rows, aligned, amount, dose_rows,
                                 floor)

    probability <- vapply(seq_len(nrow(cells)), function(row) {
      reached <- vapply(member_rows, function(rows) {
        selected <- rows & observed & endpoint == cells$endpoint[row]
        any(is.finite(aligned[selected]) &
              abs(aligned[selected] - cells$time[row]) <
                sqrt(.Machine$double.eps))
      }, logical(1))
      mean(reached)
    }, numeric(1))
    visits[[arm]] <- list(cells = index, probability = .named(probability, cells$name))
  }
  list(dosing = dosing, visits = visits, sizes = sizes, arms = arms,
       cells = index)
}

# `cells` may carry a `name` column, and a caller that has names for its grid
# cells gets them back on everything keyed by one. A caller that does not is
# not given invented ones.
.named <- function(x, labels) {
  if (is.null(labels)) x else stats::setNames(x, labels)
}

# An arm of one or two has no between-subject spread to model: whatever the arm
# summary is -- a mean score vector, a dose ladder, a per-visit attendance rate
# -- it is that patient, and its spread is noise around them. Refusing is the
# only honest answer, and it is loud rather than a silent pooling the caller
# never asked for. Shared by every generator that summarizes an arm.
.require_arms <- function(group, minimum, what) {
  minimum <- as.integer(minimum)
  if (!is.finite(minimum) || minimum < 1L) {
    stop("`min_arm_patients` must be one positive integer.", call. = FALSE)
  }
  sizes <- table(group)
  short <- sizes[sizes < minimum]
  if (length(short)) {
    stop("`", what, "` needs at least ", minimum,
         " patients in every arm. Short: ",
         paste(sprintf("%s (%d)", names(short), as.integer(short)),
               collapse = ", "),
         ". Pool the arm, drop the column from `strata`, or exclude those ",
         "patients before calling.", call. = FALSE)
  }
  invisible(TRUE)
}


# What the generated table has to look like to be the same study: the columns
# and their types, the compartment numbers, the assay limit per endpoint, the
# values carried verbatim per arm, and how subject identifiers were written.
# Read from the source once, by whichever generator is summarizing it.
.source_schema <- function(source, roles, endpoints, subject_group) {
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
  cmt_obs <- stats::setNames(lapply(endpoints, function(ep) {
    if (is.null(roles$cmt)) NULL else mode_of(observed & endpoint == ep,
                                              roles$cmt)
  }), endpoints)

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
  censoring <- stats::setNames(lapply(endpoints, function(ep) {
    .source_censoring(source, roles, ep)
  }), endpoints)

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

.assign_arms <- function(arms, sizes, n_subjects) {
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
  rep(arms, times = counts)[seq_len(n_subjects)]
}

.new_subject_ids <- function(schema, n) {
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
