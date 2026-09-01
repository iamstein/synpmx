# What the study is: which endpoint is the drug, and what design produced it ---
#
# Two questions answered from the source before anything is fitted. Neither one
# decides the model: the first decides what the model is *of*, and the second
# prunes the candidate set the fit then chooses within. Both are inference with
# a declaration available as an override, which is the fork `dose_covariate`
# already answers the same way -- infer where the data settles the question, and
# let the caller say so where it does not.

# One row per observation, carrying what every signal below needs: which dose
# interval the row falls in, how long after that dose it was drawn, and what
# amount the subject started on. Computed once because four signals and the
# richness counts all read it.
#
# Two time axes, and the signals read the nominal one. A shape, a peak position
# and a count of distinct sampling times are all statements about the protocol,
# and on recorded times they are statements about the clock instead: no two
# patients share a time, so a "median profile" is one point per time and a
# study sampled six times looks like a study sampled two hundred. Estimation
# reads `time` -- a population fit is about the dose that was actually given --
# and everything here reads `nominal_time`.
.model_dose_relative <- function(source, roles, time, dosed) {
  id <- as.character(source[[roles$id]])
  interval <- rep(NA_integer_, nrow(source))
  tad <- rep(NA_real_, nrow(source))
  first_time <- rep(NA_real_, nrow(source))
  for (rows in split(seq_len(nrow(source)), id)) {
    dose_at <- rows[dosed[rows] & is.finite(time[rows])]
    if (!length(dose_at)) next
    dose_time <- sort(time[dose_at])
    index <- findInterval(time[rows], dose_time)
    interval[rows] <- ifelse(index >= 1L, index, NA_integer_)
    tad[rows] <- ifelse(index >= 1L, time[rows] - dose_time[pmax(index, 1L)],
                        NA_real_)
    first_time[rows] <- dose_time[1L]
  }
  list(interval = interval, tad = tad, first_time = first_time)
}

.model_observations <- function(source, roles) {
  time <- suppressWarnings(as.numeric(source[[roles$time]]))
  nominal <- if (is.null(roles$nominal_time)) time else
    suppressWarnings(as.numeric(source[[roles$nominal_time]]))
  dv <- suppressWarnings(as.numeric(source[[roles$dv]]))
  endpoint <- .endpoint(source, roles)
  dosed <- .dose_rows(source, roles)
  observed <- .observation_rows(source, roles, require_present = TRUE) &
    is.finite(nominal) & is.finite(dv)
  amount <- if (is.null(roles$amt)) rep(NA_real_, nrow(source)) else
    suppressWarnings(as.numeric(source[[roles$amt]]))

  planned <- .model_dose_relative(source, roles, nominal, dosed)
  actual <- .model_dose_relative(source, roles, time, dosed)
  first_amt <- rep(NA_real_, nrow(source))
  for (rows in split(seq_len(nrow(source)), as.character(source[[roles$id]]))) {
    dose_at <- rows[dosed[rows] & is.finite(time[rows])]
    if (!length(dose_at)) next
    first_amt[rows] <- amount[dose_at[which.min(time[dose_at])]]
  }

  data.frame(
    subject = as.character(source[[roles$id]])[observed],
    endpoint = endpoint[observed],
    time = time[observed], ntime = nominal[observed], dv = dv[observed],
    interval = planned$interval[observed], tad = planned$tad[observed],
    actual_tad = actual$tad[observed],
    first_dose_time = planned$first_time[observed],
    first_dose_amt = first_amt[observed],
    stringsAsFactors = FALSE
  )
}

# The dose interval the cohort sampled most densely. Every shape and richness
# statement is made inside one interval, because a profile read across a dose is
# two profiles superposed and its peak is wherever the doses happened to fall.
.model_richest_interval <- function(observations) {
  present <- observations[is.finite(observations$interval), , drop = FALSE]
  if (!nrow(present)) return(NA_integer_)
  counts <- table(present$interval)
  as.integer(names(counts)[which.max(counts)])
}

# The cohort's median profile within one dose interval, on the nominal grid.
# Three things read it and all three need the same curve: the shape signal, the
# route, and the richness counts.
.model_median_profile <- function(observations, endpoint, interval) {
  rows <- observations[observations$endpoint == endpoint &
                         !is.na(observations$interval) &
                         observations$interval == interval, , drop = FALSE]
  if (!nrow(rows)) return(NULL)
  median_by_time <- tapply(rows$dv, sprintf("%.10g", rows$tad), stats::median)
  order_by <- order(as.numeric(names(median_by_time)))
  times <- as.numeric(names(median_by_time))[order_by]
  values <- as.numeric(median_by_time)[order_by]
  list(rows = rows, time = times, value = values,
       peak_time = times[which.max(values)])
}

# Step 1: which continuous endpoint is the drug concentration ----------------
#
# Four signals. Two of them are required and two break ties, which is a
# statement about how much each one is worth: being absent before the first dose
# and scaling with the dose are properties only a drug concentration has, while
# a compartment number is a convention the dataset's author chose and a rise-and
# -fall shape is one many biomarkers also have.

# Signal 1. The `cmt` role puts the endpoint where the doses go, or one above a
# dosing compartment nobody observes -- the depot-and-central convention, where
# drug is administered into 1 and measured in 2.
.model_signal_compartment <- function(source, roles, observations, endpoint) {
  if (is.null(roles$cmt)) return(NA)
  modal <- function(rows) {
    values <- as.character(source[[roles$cmt]][rows])
    values <- values[!is.na(values)]
    if (!length(values)) return(NA_character_)
    names(which.max(table(values)))
  }
  dose_cmt <- modal(.dose_rows(source, roles))
  observed <- .observation_rows(source, roles, require_present = TRUE)
  endpoint_cmt <- modal(observed & .endpoint(source, roles) == endpoint)
  if (is.na(dose_cmt) || is.na(endpoint_cmt)) return(NA)
  if (identical(dose_cmt, endpoint_cmt)) return(TRUE)
  observed_cmts <- unique(as.character(source[[roles$cmt]][observed]))
  numeric_dose <- suppressWarnings(as.numeric(dose_cmt))
  numeric_endpoint <- suppressWarnings(as.numeric(endpoint_cmt))
  if (!is.finite(numeric_dose) || !is.finite(numeric_endpoint)) return(FALSE)
  !dose_cmt %in% observed_cmts && numeric_endpoint == numeric_dose + 1
}

# Signal 2, required. Drug is absent before it is given. An observation drawn
# before the first dose is either missing from the record or sits at the assay
# limit, and a baseline biomarker is exactly the endpoint that has one.
.model_signal_post_dose <- function(source, roles, observations, endpoint,
                                    censoring) {
  rows <- observations[observations$endpoint == endpoint, , drop = FALSE]
  if (!nrow(rows)) return(FALSE)
  limit <- censoring[[endpoint]]$left %||% -Inf
  by_subject <- split(rows, rows$subject)
  clean <- vapply(by_subject, function(part) {
    before <- part$time < part$first_dose_time[1L] - 1e-8 |
      !is.finite(part$interval)
    !any(before) || all(part$dv[before] <= limit + 1e-8)
  }, logical(1))
  mean(clean) > 0.5
}

# Signal 3. Rises to a maximum and comes down again, once. Read off the median
# profile rather than any subject's, because one patient's noise is not a shape.
.model_signal_shape <- function(observations, endpoint, interval) {
  profile <- .model_median_profile(observations, endpoint, interval)
  if (is.null(profile) || length(profile$value) < 3L) return(NA)
  peak <- which.max(profile$value)
  after <- profile$value[seq(peak, length(profile$value))]
  # No rise after the peak beyond the noise the profile already shows.
  tolerance <- 0.05 * max(abs(profile$value))
  all(diff(after) <= tolerance)
}

# The dose levels the study actually has, or none.
#
# A level is a dose several patients were assigned, which is what makes the
# median of an arm a statement about that arm. A study dosed by body weight
# gives every patient a slightly different amount, and reading each of those as
# its own level would compare a median of one against a median of one over a
# dose ratio near 1 -- which almost any endpoint passes, including the biomarker
# this signal exists to reject. Where the amounts are a continuum there are no
# levels to compare and the signal is not computable, exactly as on a study with
# one arm.
.model_dose_levels <- function(observations) {
  amounts <- vapply(split(observations$first_dose_amt, observations$subject),
                    function(x) x[1L], numeric(1))
  amounts <- amounts[is.finite(amounts) & amounts > 0]
  if (!length(amounts)) return(numeric())
  levels <- sort(unique(amounts))
  if (length(levels) < 2L || length(levels) > 0.5 * length(amounts)) {
    return(numeric())
  }
  levels
}

# Signal 4, required. Concentration scales with the dose given; a biomarker does
# not, or does so through an exposure-response relationship far weaker than
# proportional. Undefined on a single-level study and on one dosed by body
# weight, where the other three signals decide instead.
.model_signal_proportional <- function(observations, endpoint) {
  rows <- observations[observations$endpoint == endpoint &
                         is.finite(observations$first_dose_amt) &
                         observations$first_dose_amt > 0, , drop = FALSE]
  if (!nrow(rows)) return(NA)
  levels <- .model_dose_levels(observations)
  if (length(levels) < 2L) return(NA)
  peak <- stats::aggregate(
    cbind(dv, first_dose_amt) ~ subject, data = rows,
    FUN = function(x) max(x, na.rm = TRUE)
  )
  peak <- peak[peak$first_dose_amt %in% levels, , drop = FALSE]
  if (!nrow(peak)) return(NA)
  low <- stats::median(peak$dv[peak$first_dose_amt == levels[1L]])
  high <- stats::median(peak$dv[peak$first_dose_amt == levels[length(levels)]])
  if (!is.finite(low) || low <= 0 || !is.finite(high) || high <= 0) return(NA)
  dose_ratio <- levels[length(levels)] / levels[1L]
  abs(log((high / low) / dose_ratio)) <= log(2)
}

# The four signals, as a table, for every continuous endpoint. Returned whether
# or not the classification succeeds, because a caller told "no endpoint passed"
# needs to see which signal refused each one.
.model_endpoint_signals <- function(source, roles, observations) {
  specs <- .endpoint_value_types(source, roles)
  continuous <- .model_fittable_endpoints(specs)
  censoring <- stats::setNames(lapply(continuous, function(ep) {
    .source_censoring(source, roles, ep)
  }), continuous)
  interval <- .model_richest_interval(observations)

  data.frame(
    endpoint = continuous,
    compartment = vapply(continuous, function(ep) {
      .model_signal_compartment(source, roles, observations, ep)
    }, logical(1)),
    post_dose = vapply(continuous, function(ep) {
      .model_signal_post_dose(source, roles, observations, ep, censoring)
    }, logical(1)),
    shape = vapply(continuous, function(ep) {
      .model_signal_shape(observations, ep, interval)
    }, logical(1)),
    proportional = vapply(continuous, function(ep) {
      .model_signal_proportional(observations, ep)
    }, logical(1)),
    row.names = NULL, stringsAsFactors = FALSE
  )
}

# The classification. `endpoint_roles` is the override, and it is a declaration
# rather than a hint: a named endpoint is accepted whatever the signals say,
# because the caller knows what the assay measured and the signals are guessing.
.model_classify_endpoints <- function(source, roles, observations,
                                      endpoint_roles = NULL) {
  signals <- .model_endpoint_signals(source, roles, observations)
  specs <- .endpoint_value_types(source, roles)
  continuous <- signals$endpoint
  discrete <- setdiff(names(specs), continuous)

  if (!is.null(endpoint_roles)) {
    declared <- unname(endpoint_roles[["pk"]] %||% endpoint_roles[[1L]])
    if (!declared %in% names(specs)) {
      stop("`endpoint_roles` names `", declared, "`, which is not an endpoint ",
           "in this data. Endpoints here are: ",
           paste(names(specs), collapse = ", "), ".", call. = FALSE)
    }
    if (!declared %in% continuous) {
      stop("`endpoint_roles` names `", declared, "`, which is a ",
           specs[[declared]]$type, " endpoint. A drug concentration has to be ",
           "a continuous time course to be fitted; binary and ordinal ",
           "endpoints are generated from their per-visit marginals instead.",
           call. = FALSE)
    }
    return(list(pk = declared, pd = setdiff(continuous, declared),
                discrete = discrete, signals = signals, decided_by = "declared"))
  }

  required <- signals$post_dose &
    (is.na(signals$proportional) | signals$proportional)
  required[is.na(required)] <- FALSE
  passing <- signals$endpoint[required]

  if (!length(passing)) {
    stop("No endpoint looks like a drug concentration: none is both absent ",
         "before the first dose and dose-proportional. Name the concentration ",
         "with `endpoint_roles = c(pk = \"...\")`, or use `synpmx_avatar()` or ",
         "`synpmx_pca()`, which fit no structural model. Signals read: ",
         .model_signal_summary(signals), call. = FALSE)
  }
  if (length(passing) > 1L) {
    breaks <- signals[required, , drop = FALSE]
    score <- rowSums(cbind(!is.na(breaks$compartment) & breaks$compartment,
                           !is.na(breaks$shape) & breaks$shape))
    best <- which(score == max(score))
    if (length(best) > 1L) {
      stop(length(best), " endpoints look equally like a drug concentration (",
           paste(breaks$endpoint[best], collapse = ", "), "): the required ",
           "signals pass for all of them and the compartment and shape ",
           "signals do not separate them. Name the concentration with ",
           "`endpoint_roles = c(pk = \"...\")`.", call. = FALSE)
    }
    passing <- breaks$endpoint[best]
  }

  list(pk = passing, pd = setdiff(continuous, passing), discrete = discrete,
       signals = signals, decided_by = "inferred")
}

# What a time course can be fitted to. `integer` joins `continuous` because a
# prothrombin activity or a cell count recorded as whole numbers is a continuous
# quantity that was rounded, and rounding it back is what `.snap_endpoint_values()`
# already does at emit. `binary` and `ordinal` are not time courses and come
# from per-visit marginals instead.
.model_fittable_endpoints <- function(specs) {
  names(specs)[vapply(specs, function(s) {
    s$type %in% c("continuous", "integer")
  }, logical(1))]
}

.model_signal_summary <- function(signals) {
  paste(vapply(seq_len(nrow(signals)), function(i) {
    sprintf("%s (post-dose %s, proportional %s)", signals$endpoint[i],
            .model_yes_no(signals$post_dose[i]),
            .model_yes_no(signals$proportional[i]))
  }, character(1)), collapse = "; ")
}

.model_yes_no <- function(x) {
  if (is.na(x)) "not computable" else if (x) "yes" else "no"
}

# Step 2: what design produced it --------------------------------------------
#
# Detection prunes; the fit chooses. Nothing here decides a model on its own,
# and where the source does not separate two cases both survive into the
# candidate set and AIC settles it.

# A `rate` role carrying a nonzero value, or a dose record with a duration, is
# an infusion and there is nothing to infer.
#
# Otherwise there is one property that separates the two remaining routes, and
# it is not how many patients peak early. A drug given by mouth cannot be in the
# blood at the moment it is swallowed, so a concentration observed at time zero
# after a dose is intravenous; a profile that rises before it falls is oral. A
# profile that declines from a first sample drawn some time after the dose is
# both of those things at once -- an intravenous bolus, or an oral dose whose
# absorption finished before anybody looked -- and no amount of reading the
# source separates them. There both routes are offered and AIC settles it, which
# is what "detection prunes, the fit chooses" means.
#
# Counting subjects instead is the mistake this replaces. On `warfarin`, 22 of
# 32 patients are first sampled at 24 hours and every one of their profiles
# declines from its first point, so a per-subject vote reads 69% intravenous
# from a study that is oral -- it counts the sampling schedule, not the drug.
.model_detect_route <- function(source, roles, observations, pk_endpoint,
                                interval) {
  dosed <- .dose_rows(source, roles)
  if (!is.null(roles$rate)) {
    rate <- suppressWarnings(as.numeric(source[[roles$rate]][dosed]))
    if (any(is.finite(rate) & rate != 0)) {
      return(list(route = "infusion", rising = NA_real_,
                  reason = "a nonzero `rate` on the dose records"))
    }
  }
  profile <- .model_median_profile(observations, pk_endpoint, interval)
  if (is.null(profile) || length(profile$value) < 2L) {
    return(list(route = "both", rising = NA_real_,
                reason = "too few distinct sampling times to place a peak"))
  }

  # The share of subjects whose own profile rises before it falls. Reported
  # rather than decisive: it is the number that says how much of the cohort the
  # median profile speaks for.
  by_subject <- split(profile$rows, profile$rows$subject)
  rising <- mean(vapply(by_subject, function(part) {
    part <- part[order(part$tad), , drop = FALSE]
    nrow(part) >= 2L && which.max(part$dv) > 1L
  }, logical(1)))

  if (profile$peak_time > min(profile$time) + 1e-8) {
    return(list(route = "oral", rising = rising,
                reason = sprintf(
                  "the median profile rises to a peak at %.4g before declining, and %.0f%% of subjects do too",
                  profile$peak_time, 100 * rising)))
  }
  if (min(profile$time) <= 1e-8) {
    return(list(route = "iv", rising = rising,
                reason = "the concentration is at its highest at the moment of the dose, which oral absorption cannot produce"))
  }
  list(route = "both", rising = rising,
       reason = sprintf(
         "the median profile declines from its first sample at %.4g after the dose, which an intravenous bolus and a fully absorbed oral dose both produce",
         min(profile$time)))
}

# Whether the sampling would support a two-compartment model. Reported rather
# than acted on: the candidate set is one-compartment only, and this is what
# tells a reader that asking for `pk = "2cmt_oral"` is worth their time.
#
# The reading itself is the standard one. A distribution phase cannot be
# identified from troughs, so it needs four distinct times in one dose interval
# with two after the peak. A sample *before* the peak is required only where
# there is an ascending limb to sample -- an intravenous bolus peaks at the
# dose, and asking for a sample before it would refuse every such study.
.model_sampling_richness <- function(observations, pk_endpoint, interval) {

  profile <- .model_median_profile(observations, pk_endpoint, interval)
  if (is.null(profile)) {
    return(list(per_subject = 0, before_peak = 0, after_peak = 0, rich = FALSE))
  }
  peak_time <- profile$peak_time
  by_subject <- split(profile$rows, profile$rows$subject)
  distinct <- vapply(by_subject, function(part) length(unique(part$tad)),
                     integer(1))
  before <- vapply(by_subject, function(part) {
    length(unique(part$tad[part$tad < peak_time - 1e-8]))
  }, integer(1))
  after <- vapply(by_subject, function(part) {
    length(unique(part$tad[part$tad > peak_time + 1e-8]))
  }, integer(1))

  list(
    per_subject = stats::median(distinct),
    before_peak = stats::median(before),
    after_peak = stats::median(after),
    # One sample before the peak is required only where there is an ascending
    # limb to sample. An intravenous bolus peaks at the dose, and asking for a
    # sample before it would refuse every such study a distribution phase.
    rich = stats::median(distinct) >= 4 && stats::median(after) >= 2 &&
      (peak_time <= min(profile$time) + 1e-8 || stats::median(before) >= 1)
  )
}

.model_detect_design <- function(source, roles, observations, pk_endpoint) {
  interval <- .model_richest_interval(observations)
  route <- .model_detect_route(source, roles, observations, pk_endpoint,
                               interval)
  richness <- .model_sampling_richness(observations, pk_endpoint, interval)

  # One compartment, and that is the whole default candidate set.
  #
  # A two-compartment model is not what this generator is for. It exists to make
  # simulated profiles resemble the source study, and a distribution phase is a
  # refinement of a shape the one-compartment model already has -- while costing
  # a fit that takes five times as long and, on a study a one-compartment model
  # describes, spends that time against a flat likelihood. Measured on a
  # thirty-subject oral study: 12 s against 49 s, for a worse AIC.
  #
  # It remains available. `pk = "2cmt_oral"` or `pk = "2cmt_iv"` forces it and
  # skips the search, and `richness$rich` says whether the sampling would
  # support one, so a caller who wants it is told when it is worth asking for.
  candidates <- switch(route$route,
                       infusion = "1cmt_infusion",
                       iv = "1cmt_iv",
                       oral = "1cmt_oral",
                       both = c("1cmt_iv", "1cmt_oral"))

  list(route = route$route, rising = route$rising,
       reason = route$reason, interval = interval, richness = richness,
       candidates = candidates)
}
