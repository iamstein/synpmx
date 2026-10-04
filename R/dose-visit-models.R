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

# The planned value at a cycle: the most common one, where at least `floor`
# patients share it. Where none is shared -- a study dosed per kilogram gives
# every patient their own amount, and a fixed infusion time then gives each their
# own rate -- the mode of values held once is one patient's dose (REV-060), so the
# arm's mean is used instead, rounded to two significant figures.
.shared_mode <- function(values, floor) {
  counts <- table(sprintf("%.10g", values))
  if (max(counts) >= floor) return(as.numeric(names(counts)[which.max(counts)]))
  signif(mean(values), 2L)
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
.dose_model <- function(member_rows, aligned, amount, dose_rows, floor,
                        route = NULL, rate = NULL, adm = NULL) {
  n <- length(member_rows)
  if (is.null(route)) route <- rep(NA_character_, length(amount))
  if (is.null(rate)) rate <- rep(0, length(amount))
  # The administration id, carried so that the generated dose record can be
  # written back with the one the study used. Reversing the route mapping is
  # not enough for a study whose ids separate two drugs given the same way:
  # both are `extravascular`, and the reverse lookup would stamp every dose
  # with whichever id came first.
  if (is.null(adm)) adm <- rep(NA_character_, length(amount))
  per_patient <- lapply(seq_len(n), function(i) {
    selected <- member_rows[[i]] & dose_rows
    order_by <- order(aligned[selected])
    data.frame(time = aligned[selected][order_by],
               amt = amount[selected][order_by],
               route = route[selected][order_by],
               rate = rate[selected][order_by],
               adm = adm[selected][order_by],
               stringsAsFactors = FALSE)
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
    .shared_mode(values, floor)
  }, numeric(1))
  # How each planned dose is given, alongside how much. The route and the rate
  # are properties of the administration rather than of the patient, so each is
  # the modal value among the patients dosed at that cycle -- the same reading
  # the amount gets. A schedule that switches route partway through, or runs an
  # infusion at one cycle and a bolus at the next, is carried as that.
  modal <- function(values) {
    values <- values[!is.na(values)]
    if (!length(values)) return(NA)
    counts <- table(as.character(values))
    key <- names(counts)[which.max(counts)]
    values[match(key, as.character(values))]
  }
  at_cycle <- function(column) {
    lapply(seq_len(n_cycles), function(k) {
      unlist(lapply(per_patient, function(p) {
        hit <- which(abs(p$time - planned_times[k]) < 1e-8)
        if (length(hit)) p[[column]][hit[1L]] else NULL
      }))
    })
  }
  planned_route <- vapply(at_cycle("route"),
                          function(v) as.character(modal(v)), character(1))
  planned_rate <- vapply(at_cycle("rate"), function(v) {
    v <- suppressWarnings(as.numeric(v))
    v <- v[is.finite(v)]
    if (!length(v)) 0 else .shared_mode(v, floor)
  }, numeric(1))
  planned_adm <- vapply(at_cycle("adm"),
                        function(v) as.character(modal(v)), character(1))
  planned <- data.frame(cycle = seq_len(n_cycles), time = planned_times,
                        amt = planned_amt, route = planned_route,
                        rate = planned_rate, adm = planned_adm,
                        stringsAsFactors = FALSE)

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
  #
  # A level is the fraction dropped to, to one significant figure (REV-073):
  # 0.74, 0.75 and 0.78 are one level, 0.8, rather than three that each one or
  # two patients hold. Fewer, coarser levels are fewer released numbers, and
  # more patients stand behind each.
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
    counts <- table(signif(dropped_to, 1L))
    holders <- vapply(names(counts), function(value) {
      sum(vapply(ratios, function(r) {
        finite <- which(is.finite(r))
        if (length(finite) < 2L) return(FALSE)
        values <- r[finite]
        fell <- which(diff(values) < -drop_tolerance * utils::head(values, -1L))
        length(fell) > 0 &&
          any(abs(signif(values[fell + 1L], 1L) - as.numeric(value)) < 1e-8)
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
  dropped <- integer(n)
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
        dropped[i] <- sum(diff(values) <
                            -drop_tolerance * utils::head(values, -1L))
        drops <- drops + dropped[i]
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
    # How many patients, events and patient-cycles stand behind each rate.
    # Read by `.pool_rates()` and removed by it, so that the dose model
    # that is stored and released carries the rates alone.
    support = data.frame(
      rate = c("discontinuation", "interruption", "reduction"),
      patients = c(sum(span > 0 & span < n_cycles), sum(skipped > 0),
                   sum(dropped > 0)),
      events = c(stops, skips, drops),
      at_risk = c(at_risk_stop, at_risk_skip, at_risk_drop),
      stringsAsFactors = FALSE),
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
  # A schedule summarised before the route and the rate were carried has
  # neither column; it is one route given as a bolus, which is what these
  # defaults say.
  if (is.null(planned$route)) planned$route <- rep(NA_character_, nrow(planned))
  if (is.null(planned$rate)) planned$rate <- rep(0, nrow(planned))
  if (is.null(planned$adm)) planned$adm <- rep(NA_character_, nrow(planned))
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
  data.frame(time = planned$time[keep], amt = amounts[keep],
             route = planned$route[keep], adm = planned$adm[keep],
             # The rate follows the amount. Where a dose was reduced, an
             # infusion given over the same duration runs at a lower rate, which
             # is what a reduction means at the bedside; keeping the rate and
             # shortening the infusion would be a different decision and not one
             # anything here saw made.
             rate = ifelse(planned$amt[keep] > 0,
                           planned$rate[keep] * amounts[keep] /
                             pmax(planned$amt[keep], .Machine$double.eps),
                           0),
             stringsAsFactors = FALSE)
}

# The dosing model and the visit model, one of each per arm, with everything
# read from patients pooled over the arms (REV-071, REV-072).
#
# `cells` is the grid attendance is measured on: one row per endpoint and
# nominal time, carrying the `index` the caller knows that cell by. Passing the
# grid rather than a fitted object is what keeps this callable by a generator
# that has no components.
#
# Both are summaries of the study rather than facts about a patient. The dosing
# model is built above: a planned schedule per arm, which is the protocol, and
# three rates and a ladder, pooled over the arms. The visit model is the
# visits each arm has and one attendance rate for the study, so attendance is
# drawn per visit rather than a real patient's set of attended visits being
# reused.
.arm_models <- function(source, roles, cells, subject_group, floor,
                        dose_groups = NULL) {
  subjects <- .unique_in_order(source[[roles$id]])
  aligned <- .aligned_time(source, roles)
  observed <- .observation_rows(source, roles, require_present = TRUE)
  endpoint <- .endpoint(source, roles)
  # Every dosing event, not only the ones carrying drug. A placebo arm records
  # its administrations with `AMT = 0`, and `.dose_rows()` drops those whenever
  # any positive amount exists in the study, which would leave the placebo arm
  # with no dosing events at all.
  # A negative amount is not a dose. NONMEM writes the end of an infusion as a
  # mirror record -- `AMT` and `RATE` both negated at the stop time -- and
  # `wbcSim` is written that way, so reading every event as an administration
  # planned a schedule of alternating doses and anti-doses, and the generated
  # data carried them (`SIM-072`).
  amount <- if (is.null(roles$amt)) rep(0, nrow(source)) else
    suppressWarnings(as.numeric(source[[roles$amt]]))
  amount[!is.finite(amount)] <- 0
  dose_rows <- .event_rows(source, roles) & amount >= 0
  route <- .dose_routes(source, roles)
  admin <- if (is.null(roles$adm)) rep(NA_character_, nrow(source)) else
    as.character(source[[roles$adm]])
  rate <- if (is.null(roles$rate)) rep(0, nrow(source)) else
    suppressWarnings(as.numeric(source[[roles$rate]]))
  rate[!is.finite(rate) | rate < 0] <- 0

  arms <- unique(subject_group)
  index <- .named(cells$index, cells$name)
  dosing <- list()
  visits <- list()
  attenders <- list()
  sizes <- integer()

  for (arm in arms) {
    members <- subjects[subject_group == arm]
    member_rows <- lapply(members, function(subject) {
      !is.na(source[[roles$id]]) & source[[roles$id]] == subject
    })
    sizes[arm] <- length(members)

    # One dose model, or one per drug where `dose_groups` says which dose
    # record is whose. A combination study's arm has two schedules -- two
    # amounts, two intervals, two ladders of reduction -- and pooling them
    # plans a regimen neither drug was given: at a nominal time both drugs were
    # dosed at, only the first record would survive into the plan.
    dosing[[arm]] <- if (is.null(dose_groups)) {
      .dose_model(member_rows, aligned, amount, dose_rows, floor, route, rate,
                  admin)
    } else {
      stats::setNames(lapply(.unique_in_order(dose_groups[!is.na(dose_groups)]),
                             function(group) {
        .dose_model(member_rows, aligned, amount,
                    dose_rows & !is.na(dose_groups) & dose_groups == group,
                    floor, route, rate, admin)
      }), .unique_in_order(dose_groups[!is.na(dose_groups)]))
    }

    # Which visits each patient was observed at: one row per patient.
    attenders[[arm]] <- matrix(vapply(seq_len(nrow(cells)), function(row) {
      vapply(member_rows, function(rows) {
        selected <- rows & observed & endpoint == cells$endpoint[row]
        any(is.finite(aligned[selected]) &
              abs(aligned[selected] - cells$time[row]) <
                sqrt(.Machine$double.eps))
      }, logical(1))
    }, logical(length(members))), nrow = length(members))
  }

  # One attendance rate for the study (REV-071, REV-072), as the PD shapes are
  # one time course. An arm has a visit where at least `floor` of its patients
  # were observed there, which is the protocol's schedule for that arm and keeps
  # a placebo arm from being sampled for drug; the rate is the share of those
  # scheduled visits, over every arm, at which a patient was observed. One
  # number rests on every patient, where a share per visit gave an attacker who
  # knows a candidate's attendance one more thing to match; what is lost is
  # when in the study visits were missed. The threshold rule applies to the
  # patients on each side: fewer than `floor` patients who missed any visit
  # leaves the rate at 1.
  counts <- vapply(arms, function(arm) colSums(attenders[[arm]]),
                   numeric(nrow(cells)))
  if (!is.matrix(counts)) {
    counts <- matrix(counts, nrow = nrow(cells), dimnames = list(NULL, arms))
  }
  has <- counts >= floor
  seen <- sum(vapply(arms, function(arm) {
    sum(attenders[[arm]][, has[, arm], drop = FALSE])
  }, numeric(1)))
  scheduled <- sum(vapply(arms, function(arm) {
    nrow(attenders[[arm]]) * sum(has[, arm])
  }, numeric(1)))
  missed_any <- sum(vapply(arms, function(arm) {
    sum(rowSums(!attenders[[arm]][, has[, arm], drop = FALSE]) > 0)
  }, numeric(1)))
  seen_any <- sum(vapply(arms, function(arm) {
    sum(rowSums(attenders[[arm]][, has[, arm], drop = FALSE]) > 0)
  }, numeric(1)))
  rate <- if (scheduled > 0) seen / scheduled else 0
  thin <- min(missed_any, seen_any) > 0 && min(missed_any, seen_any) < floor
  if (thin) rate <- as.numeric(rate >= 0.5)
  for (arm in arms) {
    visits[[arm]] <- list(cells = index, probability = .named(
      ifelse(has[, arm], rate, 0), cells$name))
  }
  attendance <- data.frame(
    arm = "pooled", cell = NA_integer_, attenders = seen_any,
    misses = missed_any, rounded = thin,
    # Arm visits whose one or two patients were left out of them.
    masked = sum(counts > 0 & !has), stringsAsFactors = FALSE)
  pooled <- .pool_rates(dosing, floor)
  list(dosing = pooled$dosing, visits = visits, sizes = sizes, arms = arms,
       cells = index,
       audit = list(rates = pooled$audit, attendance = attendance))
}

# pmxmodel-algorithm.Rmd, Step 4, and pca-algorithm.Rmd, Step 6 (REV-065,
# REV-071). The three dose-change rates are pooled over every arm dosing the
# same drug, as the PD shapes are: each is the events over the patient-cycles
# at risk across those arms, and zero where fewer than `floor` patients had the
# change, because a rate one or two patients decide is their dosing history.
# The ladder is pooled with them, as every level some arm's patients share,
# so a reduction can step down in any arm. Pooling puts every arm's patients
# behind each rate, and a dose-dependent reduction or discontinuation becomes
# the study's. The counts behind each rate are removed from the dose model
# afterwards.
.pool_rates <- function(dosing, floor) {
  entries <- list()
  for (arm in names(dosing)) {
    models <- .dose_group_models(dosing[[arm]])
    groups <- names(models) %||% rep("", length(models))
    for (k in seq_along(models)) {
      entries[[length(entries) + 1L]] <- list(arm = arm, group = groups[[k]],
                                              k = k, model = models[[k]])
    }
  }
  audit <- list()
  for (group in unique(vapply(entries, function(e) e$group, character(1)))) {
    members <- which(vapply(entries, function(e) identical(e$group, group),
                            logical(1)))
    ladder <- sort(unique(unlist(lapply(entries[members], function(e) {
      e$model$levels
    }))), decreasing = TRUE)
    for (rate in c("discontinuation", "interruption", "reduction")) {
      support <- do.call(rbind, lapply(entries[members], function(e) {
        e$model$support[e$model$support$rate == rate, , drop = FALSE]
      }))
      if (is.null(support) || !nrow(support)) next
      patients <- sum(support$patients)
      enough <- patients >= floor && sum(support$at_risk) > 0
      value <- if (enough) min(1, sum(support$events) / sum(support$at_risk)) else 0
      for (e in members) entries[[e]]$model[[rate]] <- value
      audit[[length(audit) + 1L]] <- data.frame(
        arm = "pooled", group = group, rate = rate, patients = patients,
        action = if (enough) "pooled over arms" else if (patients > 0)
          "set to zero" else "none",
        patients_after = if (enough) patients else 0L,
        stringsAsFactors = FALSE)
    }
    for (e in members) entries[[e]]$model$levels <- ladder
  }
  for (entry in entries) {
    entry$model$support <- NULL
    if (!is.null(dosing[[entry$arm]]$planned)) {
      dosing[[entry$arm]] <- entry$model
    } else {
      dosing[[entry$arm]][[entry$k]] <- entry$model
    }
  }
  list(dosing = dosing,
       audit = if (length(audit)) do.call(rbind, audit) else NULL)
}

# `cells` may carry a `name` column, and a caller that has names for its grid
# cells gets them back on everything keyed by one. A caller that does not is
# not given invented ones.
.named <- function(x, labels) {
  if (is.null(labels)) x else stats::setNames(x, labels)
}

# An arm of one or two has no between-subject spread to model: whatever the arm
# summary is -- a mean score vector, a dose ladder, a per-visit attendance rate
# -- it is that patient, and its spread is noise around them. So the arm is left
# out of the summary rather than modelled from one or two people, and the caller
# is told which patients went and why: the synthetic cohort is missing an arm the
# source has, which is a fact about the output rather than a detail of the run.
# Returns one logical per element of `group`. Shared by every generator that
# summarizes an arm.
.drop_short_arms <- function(group, minimum, what) {
  minimum <- as.integer(minimum)
  if (!is.finite(minimum) || minimum < 1L) {
    stop("`min_arm_patients` must be one positive integer.", call. = FALSE)
  }
  sizes <- table(group)
  short <- sizes[sizes < minimum]
  if (!length(short)) return(rep(TRUE, length(group)))
  warning(.condition_text(
    "`", what, "` dropped ", sum(as.integer(short)), " patient(s) in ",
    length(short), " arm(s) below `min_arm_patients` = ", minimum, ":",
    items = sprintf("%s (%d)", .arm_label(names(short)), as.integer(short)),
    why = paste("An arm that size is deemed too small to include, so it is",
                "absent from the synthetic data."),
    fix = paste("Pool the arm, drop the column from `strata`, or lower",
                "`min_arm_patients` to keep it.")), call. = FALSE)
  !(group %in% names(short))
}


# What the generated table has to look like to be the same study: the columns
# and their types, the compartment numbers, the assay limit per endpoint, the
# values carried verbatim per arm, and how subject identifiers were written.
# Read from the source once, by whichever generator is summarizing it.
.source_schema <- function(source, roles, endpoints, subject_group,
                           min_patients = 3L) {
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
  # Where the study doses more than one way, the dose compartment is a property
  # of the route: an intravenous dose and a subcutaneous one are recorded in
  # different compartments and have to be written back that way, or the
  # generated table says every dose went to the same place.
  routes <- .dose_routes(source, roles)
  cmt_dose_route <- if (is.null(roles$cmt) || is.null(roles$adm)) NULL else {
    present <- .unique_in_order(routes[dose_rows & !is.na(routes)])
    stats::setNames(lapply(present, function(route) {
      mode_of(dose_rows & !is.na(routes) & routes == route, roles$cmt)
    }), present)
  }
  adm_class <- if (is.null(roles$adm)) NULL else
    class(source[[roles$adm]])[[1L]]
  # And per administration id, which is what separates two drugs dosed the same
  # way: both are `extravascular`, so the route lookup above cannot tell their
  # compartments apart and a generated dose record would put both drugs in one.
  cmt_dose_adm <- if (is.null(roles$cmt) || is.null(roles$adm)) NULL else {
    ids <- as.character(source[[roles$adm]])
    present <- .unique_in_order(ids[dose_rows & !is.na(ids)])
    stats::setNames(lapply(present, function(id) {
      mode_of(dose_rows & !is.na(ids) & ids == id, roles$cmt)
    }), present)
  }
  cmt_obs <- stats::setNames(lapply(endpoints, function(ep) {
    if (is.null(roles$cmt)) NULL else mode_of(observed & endpoint == ep,
                                              roles$cmt)
  }), endpoints)

  # Each arm's value of every carried column, copied from the arm's first
  # patient as `keep` declares. A `keep` value fewer than `min_patients`
  # patients in the arm hold is one or two patients' value, and is written as
  # missing instead, whatever its type (REV-066, REV-069): a numeric column that
  # varies within an arm, an age or a date, would otherwise put the first
  # patient's own number on every synthetic patient in the arm. A `strata`
  # column is exempt: it defines the arm, whose own floor is
  # `min_arm_patients`, and a caller who lowers that floor below this one has
  # said an arm that small is still an arm.
  carried <- intersect(c(roles$strata, roles$keep), names(source))
  subjects <- .unique_in_order(source[[roles$id]])
  arm_values <- list()
  removed <- character()
  for (arm in unique(subject_group)) {
    members <- subjects[subject_group == arm]
    member_rows <- lapply(members, function(subject) {
      which(!is.na(source[[roles$id]]) & source[[roles$id]] == subject)
    })
    arm_values[[arm]] <- lapply(stats::setNames(carried, carried),
                                function(column) {
      value <- .first_present(source[[column]][member_rows[[1L]]])
      if (is.factor(value)) value <- as.character(value)
      if (!column %in% roles$strata && !is.na(value)) {
        holders <- sum(vapply(member_rows, function(rows) {
          own <- .first_present(source[[column]][rows])
          !is.na(own) && identical(as.character(own), as.character(value))
        }, logical(1)))
        if (holders < min_patients) {
          removed <<- c(removed, sprintf("`%s` in arm %s", column,
                                         .arm_label(arm)))
          # Missing in the value's own class, so a date stays a date.
          value <- value[NA_integer_]
        }
      }
      value
    })
  }
  if (length(removed)) {
    warning(.condition_text(
      "A carried value held by fewer than ", min_patients,
      " patients in its arm is written as missing:",
      items = removed,
      why = paste("`keep` copies each arm's value from its first patient, and",
                  "a value one or two patients hold is theirs."),
      fix = "Declare in `keep` only values that are constant within an arm."),
      call. = FALSE)
  }

  # The assay limit, per endpoint: one or two numbers read from the source, and
  # the only way the generator can put a value back on the boundary. Without it
  # every value below the limit is emitted as itself and an arm that is entirely
  # below quantification comes back as a spread of small numbers rather than the
  # flat line the study recorded.
  censoring <- stats::setNames(lapply(endpoints, function(ep) {
    .source_censoring(source, roles, ep)
  }), endpoints)

  # No source identifier is kept (REV-058). A zero-length factor still carries
  # its levels, so the ID column's prototype is stored without them: a study
  # whose IDs are a factor otherwise wrote every patient's label into the
  # schema. A numeric ID keeps only an offset, the largest ID rounded up to a
  # power of ten, so that a synthetic ID cannot collide with a real one and
  # the largest real ID is not stored.
  #
  # Nor is any other level fewer than `min_patients` patients hold
  # (pmxmodel-algorithm.Rmd, Step 5; REV-067).
  # Every level of the source factor otherwise travelled in its prototype and
  # into the generated column: a `strata` column kept the label of an arm
  # dropped for having too few patients, and a `keep` column a site one patient
  # came from. A carried column keeps the values its arms carry, which already
  # meet their floors; any other factor, a covariate included (REV-070), keeps
  # the levels at least `min_patients` patients hold. The model generator's
  # covariate model then narrows a covariate to the levels it draws from.
  identifiers <- source[[roles$id]]
  prototypes <- lapply(stats::setNames(names(source), names(source)),
                       function(column) {
    prototype <- source[[column]][0L]
    if (!is.factor(prototype) || identical(column, roles$id)) {
      return(prototype)
    }
    kept <- if (column %in% carried) {
      unlist(lapply(arm_values, function(values) {
        as.character(values[[column]])
      }), use.names = FALSE)
    } else {
      .held_levels(source[[column]], identifiers, min_patients)
    }
    factor(character(), levels = levels(prototype)[levels(prototype) %in% kept],
           ordered = is.ordered(prototype))
  })
  if (is.factor(prototypes[[roles$id]])) {
    prototypes[[roles$id]] <- factor(character())
  }
  list(
    censoring = censoring,
    columns = names(source),
    prototypes = prototypes,
    id_class = class(identifiers)[[1L]],
    id_offset = .id_offset(identifiers),
    cmt_dose = cmt_dose, cmt_dose_route = cmt_dose_route,
    cmt_dose_adm = cmt_dose_adm,
    adm_class = adm_class, cmt_obs = cmt_obs,
    carried = carried, arm_values = arm_values,
    endpoint_specs = .released_endpoint_specs(source, roles, min_patients)
  )
}

# The values of `x` at least `minimum` distinct patients hold, on any row.
.held_levels <- function(x, ids, minimum) {
  present <- !is.na(x) & !is.na(ids)
  if (!any(present)) return(character())
  holders <- tapply(as.character(ids[present]), as.character(x[present]),
                    function(h) length(unique(h)))
  names(holders)[holders >= minimum]
}

# Each endpoint's value type as the schema stores it. Generation reads the
# type, the levels and the sign, to snap a generated value onto the scale.
# pmxmodel-algorithm.Rmd, Step 5.
#
# Without the reason text (REV-060): an inferred reason quotes the observed
# range -- "53 whole-number levels, from 9 to 100" -- and a minimum and a
# maximum are each one patient's value. And without a level fewer than
# `min_patients` patients recorded (REV-068): a binary or ordinal endpoint's
# levels are every value recorded, so a grade one patient reached once was
# stored, and with it the highest grade. The discrete model already folds such
# a level away, so generation never draws it. An inferred binary or ordinal
# endpoint left with fewer than two levels is stored as whole numbers, because
# its type alone would say the rarer level had been seen; a declared type is
# the caller's and stays.
.released_endpoint_specs <- function(source, roles, min_patients) {
  observed <- .observation_rows(source, roles, require_present = TRUE)
  endpoint <- .endpoint(source, roles)
  ids <- as.character(source[[roles$id]])
  dv <- suppressWarnings(as.numeric(source[[roles$dv]]))
  specs <- .endpoint_value_types(source, roles)
  stats::setNames(lapply(names(specs), function(name) {
    spec <- specs[[name]]
    if (length(spec$levels)) {
      at <- which(observed & endpoint == name & is.finite(dv) & !is.na(ids))
      held <- vapply(spec$levels, function(level) {
        length(unique(ids[at][abs(dv[at] - level) <= 1e-8]))
      }, integer(1))
      spec$levels <- spec$levels[held >= min_patients]
      if (length(spec$levels)) spec$nonnegative <- min(spec$levels) >= 0
      if (!isTRUE(spec$declared) && length(spec$levels) < 2L) {
        spec$type <- "integer"
        spec$levels <- NULL
      }
    }
    spec$reason <- if (isTRUE(spec$declared)) {
      "declared in `pmx_roles(endpoint_types = )`"
    } else {
      paste(spec$type, "values, inferred from the data")
    }
    spec
  }), names(specs))
}

# The offset synthetic numeric IDs count up from: the largest source ID rounded
# up to a power of ten, so no synthetic ID equals a real one and the largest
# real ID itself is not stored (REV-058). A non-numeric ID needs none.
.id_offset <- function(identifiers) {
  if (!is.numeric(identifiers) || !any(is.finite(identifiers))) return(0)
  largest <- max(identifiers[is.finite(identifiers)])
  if (largest < 1) return(0)
  10^ceiling(log10(largest + 1))
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
