# Estimating the candidates --------------------------------------------------
#
# The only stage that reads patient data and the only one that needs `nlmixr2`.
# The candidate set is exactly the models `.pk_single_dose()` can evaluate in
# closed form: a candidate the fitter could estimate and the generator could not
# simulate would be a model that fits and then generates nothing, so the two
# lists are one list.

# Non-compartmental analysis, per subject, within one dose interval.
#
# One interval, and that restriction is the whole of what makes this correct on
# a multiple-dose study. Sorting a subject's samples by time after dose pools
# every cycle they were in: on a twelve-week daily regimen the resulting
# sequence interleaves a sample drawn 0.1 h into cycle 85 with one drawn 0.2 h
# into cycle 1, which is not a concentration-time profile and has no area and no
# terminal slope. Read that way, `case1_pkpd` reports a half-life of 371 hours
# and the fit is started from a number that describes nothing.
#
# The terminal slope is the quantity everything else leans on. It is fitted by
# log-linear regression over the points after the peak -- at least three of
# them, all positive -- which is the standard reading.
#
# `extrapolate` says whether anything follows this interval. Where nothing does,
# the tail is real and the area runs to infinity, so dose over it is clearance.
# Where another dose follows, the area is taken over the interval as it stands:
# at steady state dose over the interval area is clearance too, and before
# steady state it is an overestimate that is still the right order of magnitude.
.model_nca_subject <- function(part, extrapolate = TRUE) {
  part <- part[order(part$actual_tad), , drop = FALSE]
  part <- part[is.finite(part$actual_tad) & is.finite(part$dv), , drop = FALSE]
  empty <- list(cmax = NA_real_, tmax = NA_real_, auc = NA_real_,
                aumc = NA_real_, lambda_z = NA_real_)
  if (nrow(part) < 3L) return(empty)
  time <- part$actual_tad
  value <- part$dv
  peak <- which.max(value)

  width <- diff(time)
  auc <- sum(width * (utils::head(value, -1L) + utils::tail(value, -1L)) / 2)
  aumc <- sum(width * (utils::head(time * value, -1L) +
                         utils::tail(time * value, -1L)) / 2)

  # The terminal phase: everything from the peak on, positive, and at least
  # three points. A slope that comes out non-negative is not a terminal phase
  # and is discarded rather than used.
  terminal <- seq(peak, length(value))
  terminal <- terminal[value[terminal] > 0]
  lambda_z <- NA_real_
  if (length(terminal) >= 3L) {
    slope <- stats::coef(stats::lm(log(value[terminal]) ~ time[terminal]))[2L]
    if (is.finite(slope) && slope < 0) lambda_z <- -as.numeric(slope)
  }
  if (extrapolate && is.finite(lambda_z)) {
    last_value <- value[length(value)]
    last_time <- time[length(time)]
    auc <- auc + last_value / lambda_z
    aumc <- aumc + last_time * last_value / lambda_z +
      last_value / lambda_z^2
  }
  list(cmax = value[peak], tmax = time[peak], auc = auc, aumc = aumc,
       lambda_z = lambda_z)
}

# Absorption from the peak position rather than from a rule of thumb. For a
# one-compartment oral model the peak sits at log(ka/ke)/(ka - ke), which is one
# equation in one unknown once the terminal slope has given `ke`. Solving it is
# worth the ten lines: `4 / tmax` was the previous answer and it is off by
# whatever the elimination rate happens to be.
.model_ka_from_tmax <- function(tmax, ke) {
  if (!is.finite(tmax) || tmax <= 0 || !is.finite(ke) || ke <= 0) return(1)
  gap <- function(ka) log(ka / ke) / (ka - ke) - tmax
  lower <- ke * 1.01
  upper <- ke * 1000
  if (!is.finite(gap(lower)) || !is.finite(gap(upper)) ||
      gap(lower) * gap(upper) > 0) {
    return(max(4 / tmax, ke * 1.5))
  }
  as.numeric(stats::uniroot(gap, c(lower, upper))$root)
}

# Starting values read off the data rather than guessed. A population fit
# started far from the answer either converges slowly, reports the starting
# values back, or -- as the two-compartment candidates did before this was
# written -- grinds against a flat likelihood for a minute where the same fit
# from good starts takes seconds.
#
# Every value here is a textbook non-compartmental reading. Clearance is dose
# over the area extrapolated to infinity, not over the trapezoid alone, which on
# a study sampled to four half-lives understates the area by about a tenth and
# overstates clearance by the same. Volume is clearance over the terminal slope,
# which is the quantity the terminal phase actually identifies -- dose over the
# peak, the previous answer, is not a volume of any kind for an oral dose.
# Clearance and volume from single samples, for a design that cannot be read
# non-compartmentally. `given` is how much drug the subject had received by the
# time the sample was drawn and `first_dose_at` when their first dose went in,
# both on recorded times.
.model_crude_estimates <- function(rows) {
  empty <- list(cl = NA_real_, v = NA_real_)
  if (!nrow(rows) || is.null(rows$given)) return(empty)
  usable <- rows[is.finite(rows$dv) & rows$dv > 0 &
                   is.finite(rows$given) & rows$given > 0, , drop = FALSE]
  if (!nrow(usable)) return(empty)
  volume <- stats::median(usable$given / usable$dv)
  elapsed <- usable$time - usable$first_dose_at
  at_rate <- is.finite(elapsed) & elapsed > 0
  clearance <- if (any(at_rate)) {
    stats::median(usable$given[at_rate] / elapsed[at_rate] / usable$dv[at_rate])
  } else NA_real_
  list(cl = clearance, v = volume)
}

# A route needs this many subjects before it is read on its own. Below it the
# pooled read is the worse of two bad options rather than the wrong one.
.model_route_minimum <- 3L

# SIM-076. The starting values for a mixed study are read one route at a time.
#
# Pooled, the cohort's median profile is an intravenous decline and an
# extravascular rise averaged together, and no non-compartmental quantity read
# off it describes either route: on a simulated study with true `cl` 4, `v` 40
# and `ka` 0.4 the pooled read starts at 6.17, 55.6 and 1.37, and from there
# `focei` does not move at all -- which is `SIM-075`, and is how this was found.
#
# Split, each half answers the question it can. The intravenous records give
# clearance and volume with no bioavailability in the way, which is the whole
# reason `f` is identifiable in a study like this. The extravascular records
# give the absorption rate, since `tmax` only means anything where there is an
# absorption phase to peak. And the two clearances together give `f` a starting
# value that is read from the data rather than assumed: an extravascular read
# returns `cl/f`, so their ratio is `f`.
.model_initial_estimates_mixed <- function(rows, structural, pk_endpoint) {
  if (is.null(rows$route)) return(NULL)
  subjects_in <- function(route) {
    length(unique(rows$subject[!is.na(rows$route) & rows$route == route]))
  }
  if (subjects_in("iv") < .model_route_minimum ||
      subjects_in("extravascular") < .model_route_minimum) {
    return(NULL)
  }
  iv <- rows[!is.na(rows$route) & rows$route == "iv", , drop = FALSE]
  ev <- rows[!is.na(rows$route) & rows$route == "extravascular", , drop = FALSE]
  disposition <- .model_initial_estimates(
    iv, if (grepl("^2cmt", structural)) "2cmt_iv" else "1cmt_iv", pk_endpoint)
  absorption <- .model_initial_estimates(ev, "1cmt_oral", pk_endpoint)

  # `cl` from the extravascular half is `cl/f`, so the ratio is bioavailability
  # -- and it is also a check on the intravenous read. Bioavailability cannot
  # exceed one, so an implied `f` above one says the two reads disagree in a way
  # `f` cannot express, and the intravenous one is the suspect: it rests on a
  # terminal slope and an area, and a study sampled at troughs across an
  # accumulating regimen gives it neither.
  #
  # Measured on a real mixed study 2026-09-09: the intravenous read returned a
  # clearance of 4.87 where fitting the extravascular half alone found `cl/f`
  # 0.37, so the intravenous number was high by more than an order of magnitude
  # -- and, taken as the start, held the whole fit there. Where the two disagree
  # this way the extravascular read is used for disposition instead, scaled by a
  # neutral `f`, which keeps its own apparent clearance exactly where it read it.
  ratio <- disposition[["cl"]] / absorption[["cl"]]
  if (!is.finite(ratio) || ratio <= 0 || ratio > 1) {
    f <- 0.7
    disposition <- absorption[names(disposition)] * f
  } else {
    f <- max(ratio, 0.05)
  }

  out <- c(disposition[c("cl", "v")], ka = unname(absorption[["ka"]]), f = f)
  if (grepl("^2cmt", structural)) {
    out <- c(out, disposition[c("q", "v2")])
  }
  out[.required_pk_params[[structural]]]
}

# Starting values the caller supplied, over the ones read off the curve.
#
# The automatic read is non-compartmental, and on a study it cannot read -- one
# sampled only at troughs, one whose two routes disagree, one whose units are
# not what they look like -- it can start the search somewhere the optimizer
# cannot leave. `SIM-075` reports that after the fact; this is how a caller who
# knows the compound fixes it in advance.
#
# Keyed by parameter and nothing else, because the population model is fitted to
# exactly one endpoint: the concentration `endpoint_roles = c(pk = )` names. The
# PD endpoints are least-squares time courses with their own parameters and are
# not reached from here.
#
# Partial by design. A caller who knows the clearance and not the absorption
# says so, and the rest is read off the curve as before.
# With several PK endpoints the flat form is ambiguous, so it is keyed by
# endpoint: `start_param = list(parent = c(cl = 4), metabolite = c(cl = 9))`.
# With one it stays flat, because naming the only endpoint there is would be
# ceremony. Returns a list keyed by endpoint either way.
.model_split_start_param <- function(start_param, pk_endpoints) {
  if (is.null(start_param)) {
    return(stats::setNames(vector("list", length(pk_endpoints)), pk_endpoints))
  }
  if (is.list(start_param)) {
    unknown <- setdiff(names(start_param), pk_endpoints)
    if (is.null(names(start_param)) || length(unknown)) {
      stop(.condition_text(
        "`start_param` is a list, so it is read as one set of starting values ",
        "per concentration endpoint, and it names endpoint(s) that are not ",
        "fitted as one:",
        items = if (is.null(names(start_param))) "<unnamed>" else unknown,
        why = paste("The concentration endpoint(s) here are:",
                    paste(pk_endpoints, collapse = ", "))), call. = FALSE)
    }
    out <- stats::setNames(vector("list", length(pk_endpoints)), pk_endpoints)
    out[names(start_param)] <- start_param
    return(out)
  }
  if (length(pk_endpoints) > 1L) {
    stop(.condition_text(
      "`start_param` is a flat vector but ", length(pk_endpoints),
      " endpoints are fitted as concentrations, so it is not clear which one ",
      "it describes.",
      items = pk_endpoints,
      fix = paste0("Key it by endpoint: `start_param = list(`",
                   pk_endpoints[[1L]], "` = c(cl = ...))`.")), call. = FALSE)
  }
  stats::setNames(list(start_param), pk_endpoints)
}

.model_validate_start_param <- function(start_param, candidates) {
  if (is.null(start_param)) return(NULL)
  if (!is.numeric(start_param) || is.null(names(start_param)) ||
      anyNA(start_param) || any(!is.finite(start_param)) ||
      any(start_param <= 0) || any(!nzchar(names(start_param))) ||
      anyDuplicated(names(start_param))) {
    stop(.condition_text(
      "`start_param` must be a named vector of distinct, finite, positive ",
      "starting values, as `start_param = c(cl = 4, v = 40)`.",
      why = paste("Every parameter here is estimated on the log scale, so a",
                  "zero or negative starting value has no logarithm.")),
      call. = FALSE)
  }
  usable <- unique(unlist(.required_pk_params[candidates], use.names = FALSE))
  unknown <- setdiff(names(start_param), usable)
  if (length(unknown)) {
    stop(.condition_text(
      "`start_param` names parameter(s) that no candidate model has: ",
      paste(unknown, collapse = ", "), ".",
      why = paste0("The candidate(s) here are ",
                   paste(candidates, collapse = ", "), ", which take ",
                   paste(usable, collapse = ", "), "."),
      fix = paste("Name the model with `pk` if you meant a different one, or",
                  "drop the parameter.")), call. = FALSE)
  }
  start_param
}

# Applied per candidate, because a two-compartment candidate takes parameters a
# one-compartment candidate does not. A value that does not apply to this
# candidate is left out rather than an error: the caller named the parameters
# of the model they have in mind, and the search may be over several.
.model_apply_start_param <- function(start, start_param) {
  if (is.null(start_param)) return(start)
  named <- intersect(names(start_param), names(start))
  start[named] <- start_param[named]
  start
}

.model_initial_estimates <- function(observations, structural, pk_endpoint) {
  rows <- observations[observations$endpoint == pk_endpoint, , drop = FALSE]
  if (structural %in% names(.pk_mixed_forms)) {
    split_read <- .model_initial_estimates_mixed(rows, structural, pk_endpoint)
    if (!is.null(split_read)) return(split_read)
  }
  by_subject <- split(rows, rows$subject)
  dose <- stats::median(vapply(by_subject, function(part) {
    part$first_dose_amt[1L]
  }, numeric(1)), na.rm = TRUE)
  if (!is.finite(dose) || dose <= 0) dose <- 1

  # Each subject read inside the cohort's richest dose interval, which is the
  # same interval every other signal is read in.
  interval <- .model_richest_interval(rows)
  nca <- lapply(by_subject, function(part) {
    within <- part[!is.na(part$interval) & part$interval == interval, ,
                   drop = FALSE]
    if (nrow(within) < 3L) return(.model_nca_subject(part[0L, ], TRUE))
    # Nothing follows the subject's last interval, so its tail is real.
    last <- max(part$interval, na.rm = TRUE)
    .model_nca_subject(within, extrapolate = identical(interval, last))
  })
  middle <- function(name) {
    stats::median(vapply(nca, function(x) x[[name]], numeric(1)), na.rm = TRUE)
  }
  auc <- middle("auc")
  aumc <- middle("aumc")
  cmax <- middle("cmax")
  tmax <- middle("tmax")
  lambda_z <- middle("lambda_z")

  # A study with fewer than three samples inside any one dose interval has no
  # non-compartmental reading at all, and what used to stand in for one was a
  # constant: clearance 1 and volume ten. `pheno_sd` is the study that shows
  # what that costs -- one concentration per neonate per interval, so `focei`
  # starts a phenobarbital fit two orders of magnitude from the answer and lands
  # somewhere that generates a cohort at the assay floor.
  #
  # Two textbook identities need one sample each and are read here instead.
  # Amount given over concentration is a volume of distribution. Dose rate over
  # average concentration is clearance -- exactly so at steady state, and an
  # overestimate before it, which is the same caveat the interval area carries
  # above. Both are medians over the cohort's observations, and both are only
  # consulted where the non-compartmental quantity is missing.
  crude <- .model_crude_estimates(rows)

  cl <- if (is.finite(auc) && auc > 0) dose / auc else
    if (is.finite(crude$cl) && crude$cl > 0) crude$cl else 1
  volume_terminal <- if (is.finite(lambda_z) && lambda_z > 0) cl / lambda_z else
    if (is.finite(cmax) && cmax > 0) dose / cmax else
      if (is.finite(crude$v) && crude$v > 0) crude$v else 10 * cl
  ka <- .model_ka_from_tmax(tmax, if (is.finite(lambda_z)) lambda_z else
    cl / volume_terminal)

  out <- c(cl = cl, v = volume_terminal)
  if (grepl("oral|mixed", structural)) out <- c(out, ka = ka)
  # Bioavailability starts at 0.7, which is neither of the two values that
  # would make the search start on a boundary of what it can mean: 1 says the
  # extravascular dose is fully absorbed and 0 says none of it is. It is a
  # starting value rather than an assumption -- the two routes' contrast is
  # what moves it -- and nothing here can read it off the curve, since a
  # non-compartmental area is an area under whatever reached the blood.
  if (grepl("mixed", structural)) out <- c(out, f = 0.7)
  if (grepl("^2cmt", structural)) {
    # Steady-state volume from the mean residence time, which is what the first
    # moment of the curve is for. For an oral dose the residence time includes
    # the time spent absorbing, so the absorption mean is taken back off.
    mrt <- if (is.finite(aumc) && is.finite(auc) && auc > 0) aumc / auc else NA
    if (grepl("oral", structural) && is.finite(mrt) && is.finite(ka) && ka > 0) {
      mrt <- mrt - 1 / ka
    }
    vss <- if (is.finite(mrt) && mrt > 0) cl * mrt else volume_terminal / 2
    # The central compartment is the smaller half and the peripheral the rest,
    # with intercompartmental clearance started at elimination clearance. These
    # are starting values for a search, not a decomposition of the curve: what
    # matters is that they are the right order of magnitude and that the central
    # volume is below the terminal one, which the previous `v * 2` was not.
    central <- max(min(vss / 2, volume_terminal * 0.75), 1e-6)
    out[["v"]] <- central
    out <- c(out, q = cl, v2 = max(vss - central, central))
  }
  out
}

# The estimation method is `"focei"` rather than `"saem"`, which reverses what
# the design assumed, and the reason is the selection criterion. Choosing among
# candidates on AIC needs every candidate to have one. SAEM's log-likelihood is
# a Gaussian-quadrature step run after the fit, and on cohorts the size of a
# phase 1 study it returns a non-finite value: `theo_sd` fits perfectly well
# under SAEM -- clearance 2.75, volume 32.3, absorption 1.51, which are the
# textbook values -- and reports `AIC = Inf`, so a search over two candidates
# has nothing to compare. Under FOCEi the same fit reports AIC -42.1. SAEM
# remains available through `estimation` for a study large enough to give it a
# likelihood, and a candidate whose AIC is not finite is recorded as not
# converged whichever method produced it.

# One `nlmixr2` model function per candidate, written as text because that is
# what the shape of these functions is: a fixed block of parameter declarations
# and a fixed block of assignments, keyed by which parameters the structural
# model needs. `linCmt()` reads the parameter names and picks the same solution
# `.pk_single_dose()` evaluates, which is what keeps the two lists one list.
#
# `weight` folds allometric scaling into the same model rather than fitting a
# second one to compare against. The exponents are the standard 0.75 and 1 and
# are not estimated: this is a shape that makes a synthetic cohort's spread look
# right, not a covariate analysis, and testing it against AIC would double the
# cost of the only fit this function performs.
.model_allometric_exponents <- c(cl = 0.75, v = 1, q = 0.75, v2 = 1)

# The compartments a mixed model is written with. Their order is the `CMT`
# numbering the estimation data doses into: depot 1, central 2, peripheral 3,
# which is what `.model_estimation_data()` assigns each dose record by route.
.model_mixed_states <- function(structural) {
  central_out <- if (identical(structural, "2cmt_mixed")) {
    paste("    d/dt(central) <- ka * depot - (cl / v) * central -",
          "(q / v) * central + (q / vp) * peripheral")
  } else {
    "    d/dt(central) <- ka * depot - (cl / v) * central"
  }
  c("    d/dt(depot) <- -ka * depot",
    central_out,
    if (identical(structural, "2cmt_mixed"))
      "    d/dt(peripheral) <- (q / v) * central - (q / vp) * peripheral",
    "    cp <- central / v")
}

.model_nlmixr_function <- function(structural, start, error, error_start,
                                   weight = NULL) {
  parameters <- names(start)
  # Bioavailability carries no between-subject term. It is identified here only
  # by the contrast between the two routes, and asking a study to place a
  # per-subject distribution on that contrast as well is asking more than a
  # design with one dose each way can answer.
  random <- setdiff(parameters, "f")
  ini <- c(
    sprintf("    t%s <- log(%.10g)", parameters, start),
    sprintf("    eta.%s ~ %.10g", random, .model_eta_init),
    sprintf("    %s.err <- %.10g", error, error_start)
  )
  assignments <- vapply(parameters, function(parameter) {
    scaling <- if (!is.null(weight) &&
                   parameter %in% names(.model_allometric_exponents)) {
      sprintf(" * (%s / %.10g)^%.2f", weight$covariate, weight$reference,
              .model_allometric_exponents[[parameter]])
    } else ""
    if (identical(parameter, "f")) {
      return(sprintf("    f <- exp(tf)%s", scaling))
    }
    sprintf("    %s <- exp(t%s + eta.%s)%s", parameter, parameter, parameter,
            scaling)
  }, character(1))
  # `linCmt()` names the central volume `v` and the peripheral one `vp`.
  assignments <- sub("^    v2 <- ", "    vp <- ", assignments)
  # SIM-077. A mixed study is written as explicit compartments rather than
  # `linCmt()`, and the reason is bioavailability. `f(depot)` under `linCmt()`
  # does not reach the doses that `CMT` 1 sends to the depot: evaluated at known
  # parameters it scales the intravenous doses instead, so `f` came back as its
  # own reciprocal and `cl` and `v` came back as `cl/f` and `v/f` -- a fit whose
  # own `PRED` tracks the data, reporting parameters that mean something else.
  # With states the compartments are real, `f()` binds to the one it names, and
  # the predictions match `.pk_profile()` to four figures on both routes.
  #
  # Only the mixed models pay for this. A single-route study has no `f` to place
  # and keeps the closed-form solution, which needs no solver.
  if ("f" %in% parameters) {
    assignments <- c(assignments, "    f(depot) <- f",
                     .model_mixed_states(structural))
    predicted <- switch(error,
                        prop = "    cp ~ prop(prop.err)",
                        add = "    cp ~ add(add.err)")
  } else {
    predicted <- switch(error,
                        prop = "    linCmt() ~ prop(prop.err)",
                        add = "    linCmt() ~ add(add.err)")
  }
  text <- paste(c(
    "function() {", "  ini({", ini, "  })", "  model({", assignments,
    predicted, "  })", "}"
  ), collapse = "\n")
  eval(parse(text = text))
}

# The dataset `nlmixr2` reads: the concentration endpoint and the dosing
# records, on recorded times.
#
# Recorded, not planned. Fitting against the planned schedule where a patient's
# dose was reduced would push the drop in concentration that followed into
# clearance, and the model would report a population that eliminates the drug
# faster than the real one.
.model_estimation_data <- function(source, roles, pk_endpoint) {
  time <- suppressWarnings(as.numeric(source[[roles$time]]))
  dv <- suppressWarnings(as.numeric(source[[roles$dv]]))
  endpoint <- .endpoint(source, roles)
  dosed <- .dose_rows(source, roles)
  observed <- .observation_rows(source, roles, require_present = TRUE) &
    endpoint == pk_endpoint
  keep <- (dosed | observed) & is.finite(time)

  amount <- if (is.null(roles$amt)) rep(0, nrow(source)) else
    suppressWarnings(as.numeric(source[[roles$amt]]))
  amount[!is.finite(amount)] <- 0

  out <- data.frame(
    ID = as.character(source[[roles$id]])[keep],
    TIME = time[keep],
    DV = ifelse(observed[keep], dv[keep], NA_real_),
    AMT = ifelse(dosed[keep], amount[keep], 0),
    EVID = ifelse(dosed[keep], 1L, 0L),
    stringsAsFactors = FALSE
  )
  if (!is.null(roles$rate)) {
    rate <- suppressWarnings(as.numeric(source[[roles$rate]]))
    rate[!is.finite(rate)] <- 0
    out$RATE <- ifelse(out$EVID != 0L, rate[keep], 0)
  }
  # Where the study doses both ways, the solver is told which compartment each
  # dose enters: `linCmt()` with an absorption rate constant is a depot model,
  # and `CMT` 1 is that depot while `CMT` 2 is the central compartment every
  # observation is drawn from. This is the whole of what makes one patient's
  # intravenous and subcutaneous doses different events rather than one average
  # of the two.
  routes <- .dose_routes(source, roles)
  if (length(.study_routes(source, roles)) > 1L) {
    out$CMT <- ifelse(out$EVID == 0L, 2L,
                      ifelse(routes[keep] == "extravascular", 1L, 2L))
    out$CMT[is.na(out$CMT)] <- 2L
  }
  # Censoring is handed to the fitter rather than imputed away, on the same
  # convention `pmx_roles(cens=, limit=)` already carries and `nlmixr2` already
  # reads: `CENS` 1 is left-censored with `DV` holding the limit, -1 is
  # right-censored, and `LIMIT` closes the interval where the study reports the
  # other bound. A censored row then enters the likelihood as the probability
  # of falling below the limit, which is what it is, instead of as a number
  # nobody measured.
  if (!is.null(roles$cens)) {
    cens <- suppressWarnings(as.numeric(as.character(source[[roles$cens]])))
    cens[!is.finite(cens)] <- 0
    out$CENS <- ifelse(out$EVID == 0L, cens[keep], 0)
    if (!is.null(roles$limit)) {
      limit <- suppressWarnings(as.numeric(source[[roles$limit]]))
      out$LIMIT <- ifelse(out$EVID == 0L & out$CENS != 0, limit[keep],
                          NA_real_)
    }
  }
  # A subject dosed but never sampled for this endpoint adds no term to the
  # likelihood. Its dosing and its visits are real and still reach the arm
  # models, which read the source rather than this table, but a dose-only
  # record here is at best ignored by the solver and at worst an error in it.
  out <- out[out$ID %in% unique(out$ID[!is.na(out$DV)]), , drop = FALSE]
  # Subjects in the order the study lists them, not in the order their
  # identifiers sort as text. `ID` is written to the fit table as character
  # because that is what the solver wants, and sorting on it puts subject 10
  # before subject 2 whenever the study numbers its patients -- a different
  # order for the same study, which `focei` answers with a slightly different
  # optimum. Ordering on first appearance makes the fit the same fit whether
  # the identifiers are numbers or text (`SIM-071`).
  out$ID <- factor(out$ID, levels = .unique_in_order(out$ID))
  out <- out[order(out$ID, out$TIME, out$EVID == 0L), , drop = FALSE]
  out$ID <- as.character(out$ID)
  rownames(out) <- NULL
  out
}

# Repeated dosing, written as one record plus `ADDL`/`II`.
#
# `nlmixr2` sums a contribution per dose per subject at every likelihood
# evaluation, so a twelve-week daily regimen costs two orders of magnitude more
# per subject than a single dose: `case1_pkpd` hands the solver 12,750 dose
# records for 3,600 observations. `ADDL`/`II` says "and 84 more like it, every
# 24 hours", which is the representation the solver is fastest on.
#
# The catch, and the reason this is not a plain row-collapse. Real studies
# record actual dose times -- 0, 24.22, 48.28 -- and `ADDL` can only place
# doses on an exact interval, so the schedule has to be regularised first.
# Regularising the doses alone would move every sample's time after dose by
# whatever that dose had drifted, and with an absorption peak an hour after the
# dose that is the part of the profile the fit is reading. So the observations
# move with their own dose: each sample keeps the exact interval between it and
# the dose it followed, and only the absolute clock is redrawn. What is
# approximated is the spacing of the OLDER doses, whose contribution at that
# sample has been decaying for at least one interval.
#
# A subject is compressed only where the schedule is genuinely regular: the
# same amount (and infusion rate) throughout, at least three doses, and every
# interval within `tolerance` of the median. A dose reduction, a skipped cycle
# or an early stop fails all three tests and is left alone -- which is the
# behaviour that matters, since those are the studies whose dosing carries the
# information.
#
# `tolerance` is 0 by default, which is to say the schedule has to be exactly
# regular already and no sample moves at all. Measured on `case1_pkpd`
# 2026-09-02, whose doses drift by up to 0.4 h: allowing the drift compresses
# 12,750 dose records to 150 and the fit still takes 485 s against 647 s -- a
# factor of 1.3, because the solver's cost there is not only the dose count --
# while landing on a different and much worse optimum (`cl` 28.4 against 8.17,
# AIC 742,777 against -8,633). The transformation is not what moved it: solving
# both tables at one set of parameters agrees to 0.6% sample for sample, and
# `theo_md`, whose doses are exactly 24 h apart, compresses 84 records to 12
# and returns the same fit to four significant digits. A study whose dosing is
# already on a grid gets the win for free; one that records actuals is left
# alone rather than have its fit moved for a third off the run time.
.compress_dose_schedule <- function(data, tolerance = 0) {
  pieces <- lapply(split(data, factor(data$ID, levels = unique(data$ID))),
                   function(part) {
    doses <- part[part$EVID != 0L, , drop = FALSE]
    if (nrow(doses) < 3L || length(unique(doses$AMT)) != 1L) return(part)
    if (!is.null(part$RATE) && length(unique(doses$RATE)) != 1L) return(part)
    gaps <- diff(doses$TIME)
    interval <- stats::median(gaps)
    if (!is.finite(interval) || interval <= 0) return(part)
    if (max(abs(gaps - interval)) > tolerance * interval) return(part)

    nominal <- doses$TIME[[1L]] + interval * (seq_len(nrow(doses)) - 1L)
    observations <- part[part$EVID == 0L, , drop = FALSE]
    # The dose each sample followed, or the first dose for a baseline sample
    # taken before any dosing: its offset is negative and carried as such.
    after <- pmax(findInterval(observations$TIME, doses$TIME), 1L)
    offset <- observations$TIME - doses$TIME[after]
    # A trough drawn just before a dose that ran late is more than one interval
    # after the dose it belongs to, and placing it at that offset on the
    # regular grid would carry it past the next nominal dose -- turning the
    # lowest sample of an interval into a near-peak one. It is held at the end
    # of its own interval instead, which is where the protocol put it.
    offset <- pmin(offset, interval * (1 - 1e-6))
    observations$TIME <- nominal[after] + offset

    first <- doses[1L, , drop = FALSE]
    first$ADDL <- nrow(doses) - 1L
    first$II <- interval
    # `rep()` rather than the scalar, because a subject can arrive here dosed
    # and never sampled -- and a one-row value assigned to a zero-row frame is
    # an error rather than an empty column.
    observations$ADDL <- rep(0L, nrow(observations))
    observations$II <- rep(0, nrow(observations))
    rows <- rbind(first, observations)
    rows[order(rows$TIME, rows$EVID == 0L), , drop = FALSE]
  })
  # A subject left alone keeps the columns it arrived with, so a study that
  # mixes a regular regimen with a single dose, a dose reduction or an early
  # stop reaches here holding both shapes -- and `rbind()` on frames whose
  # columns differ is an error rather than a fill (`SIM-081`). The pieces that
  # did not compress are given the columns back, saying what is true of them:
  # this record is the whole dose, and no more follow it.
  compressed <- vapply(pieces, function(part) "ADDL" %in% names(part),
                       logical(1))
  if (any(compressed)) {
    pieces[!compressed] <- lapply(pieces[!compressed], function(part) {
      part$ADDL <- rep(0L, nrow(part))
      part$II <- rep(0, nrow(part))
      part
    })
  }
  # A study where nothing compressed keeps the columns off the data entirely,
  # so that the fitted dataset is the one this function was handed.
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

# What the compression did, as one line before a fit that may take minutes.
# `synpmx_model_estimate()` knows the dose-record count before it fits, and a
# user who is told "12,750 dose records" understands the wait; a user who is
# told nothing suspects a hang.
.dose_record_message <- function(before, after) {
  doses_before <- sum(before$EVID != 0L)
  doses_after <- sum(after$EVID != 0L)
  if (doses_after < doses_before) {
    return(paste0(doses_before, " dose records compressed to ", doses_after,
                  " with `ADDL`/`II`, keeping each sample's time after dose."))
  }
  if (doses_before < 1000L) return(NULL)
  paste0(doses_before, " dose records, and every likelihood evaluation sums a ",
         "contribution per dose: expect minutes rather than seconds. The ",
         "schedule is not exactly regular -- it varies in amount, or in ",
         "interval, or its dose times are recorded actuals -- so it cannot be ",
         "written as one record plus `ADDL`/`II`.")
}

# The search. Every candidate is fitted, the ones that converge are compared on
# AIC, and the ones that do not stay in the table carrying their reason -- a
# search that came down to one survivor should not look like a search that had
# one candidate.
.model_fit_candidates <- function(data, candidates, observations, pk_endpoint,
                                  error, estimation, quiet, weight = NULL,
                                  start_param = NULL) {
  fits <- list()
  rows <- list()
  starts <- list()
  for (candidate in candidates) {
    start <- .model_apply_start_param(
      .model_initial_estimates(observations, candidate, pk_endpoint),
      start_param)
    starts[[candidate]] <- start
    error_start <- if (identical(error, "prop")) 0.2 else
      stats::sd(data$DV, na.rm = TRUE) / 5
    spec <- .model_nlmixr_function(candidate, start, error, error_start,
                                   weight)
    # Wall clock rather than CPU time: the fitter threads, and what the caller
    # waited through is the number they are comparing against.
    started <- proc.time()[["elapsed"]]
    # Neither the covariance step nor the residual tables is read by anything
    # downstream: the generator takes the fixed effects, the omega matrix, the
    # residual error and the empirical Bayes estimates, and every one of those
    # is on the fit without them. Measured on a 40-patient cut of `onc_sim`
    # (2026-09-08): 244 s with them, 98 s without, same estimates to the digit.
    # The dose count is not where the time goes -- cutting each subject's dose
    # history to 30 days saved 23% -- so this is the lever, not the schedule.
    control <- if (identical(estimation, "focei")) {
      nlmixr2est::foceiControl(print = 0, covMethod = "", calcTables = FALSE)
    } else list(print = 0L)
    fit <- try(suppressWarnings(suppressMessages(
      nlmixr2est::nlmixr(spec, data, est = estimation, control = control)
    )), silent = TRUE)
    seconds <- proc.time()[["elapsed"]] - started
    converged <- !inherits(fit, "try-error") &&
      is.finite(suppressWarnings(stats::AIC(fit)))
    rows[[length(rows) + 1L]] <- data.frame(
      model = candidate, converged = converged,
      aic = if (converged) as.numeric(stats::AIC(fit)) else NA_real_,
      seconds = seconds,
      note = if (converged) "" else
        if (inherits(fit, "try-error")) .model_first_line(fit) else
          "no finite objective function",
      stringsAsFactors = FALSE
    )
    if (converged) fits[[candidate]] <- fit
    if (!quiet) {
      message(sprintf("  %-14s %-20s %s", candidate,
                      if (converged) sprintf("AIC %.1f", stats::AIC(fit))
                      else "did not converge",
                      .model_duration(seconds)))
    }
  }
  table <- do.call(rbind, rows)
  table <- table[order(table$aic, na.last = TRUE), , drop = FALSE]
  rownames(table) <- NULL
  if (!any(table$converged)) {
    stop(.condition_text(
      "No candidate model converged, so there is nothing to generate from. ",
      "Candidates tried:",
      items = table$model,
      fix = paste("`synpmx_avatar()` and `synpmx_pca()` need no identifiable",
                  "structure and will run on this study.")), call. = FALSE)
  }
  list(table = table, selected = table$model[which.min(table$aic)],
       fits = fits, starts = starts)
}

.model_first_line <- function(x) {
  message <- conditionMessage(attr(x, "condition"))
  trimws(strsplit(message, "\n", fixed = TRUE)[[1L]][1L])
}

# Every between-subject term starts here, so this is what "unchanged" is
# measured against as well as what the model text declares.
.model_eta_init <- 0.1

# A fit that reports an objective and an AIC has not necessarily estimated
# anything. Where the optimizer takes no effective step -- a bad starting point,
# a flat or badly scaled objective -- `nlmixr2` returns the starting values, and
# every number downstream is then a starting value wearing the costume of an
# estimate: the report prints them, the generator simulates from them, and the
# scorecard passes, because it asks whether the output copies anybody or changed
# the study's shape and a fit that never moved does neither.
#
# Found on a mixed-route study whose report read `f 0.7` -- the starting value
# exactly -- beside three between-subject terms all reading 0.316, which is
# `sqrt(0.1)`, the eta init. Three identical values is the tell, and nothing in
# the package said a word about it.
#
# Judged per parameter against its own start, because the scales differ by
# orders of magnitude. `tolerance` is deliberately loose: this is not asking
# whether the optimizer converged tightly, it is asking whether it moved at all.
.model_fit_movement <- function(fixed, omega, start, tolerance = 0.01) {
  names_fixed <- names(fixed)
  start <- start[names_fixed]
  moved_fixed <- abs(as.numeric(fixed) / as.numeric(start) - 1)
  eta <- if (is.null(omega) || !nrow(omega)) numeric() else diag(omega)
  moved_eta <- if (length(eta)) abs(eta / .model_eta_init - 1) else numeric()
  changes <- data.frame(
    parameter = c(names_fixed,
                  if (length(eta)) paste0("omega.", rownames(omega))),
    start = c(as.numeric(start), rep(.model_eta_init, length(eta))),
    estimate = c(as.numeric(fixed), as.numeric(eta)),
    relative_change = c(moved_fixed, moved_eta),
    stringsAsFactors = FALSE
  )
  # Reported separately, because the fixed effects moving while every
  # between-subject term sits on its init is a fit that found a population mean
  # and never looked for its spread. The stored `case1_pkpd` fit was that for
  # months: `cl`, `v` and `ka` estimated, all three omegas at 0.1, and the
  # synthetic spread 0.68 of the source's as a result. Not the banner -- the
  # fixed effects are estimates -- but a line the reader has to see.
  list(moved = any(changes$relative_change > tolerance, na.rm = TRUE),
       omega_moved = if (length(eta)) any(moved_eta > tolerance, na.rm = TRUE)
                     else NA,
       tolerance = tolerance, changes = changes)
}

# Which subjects the population model is fitted to, under `max_fit_subjects`.
# `NULL` where the cap does not bind. Proportional to the arms, so a study of
# eleven cohorts keeps eleven cohorts in the fit: each arm gets its share of the
# cap, rounded, every arm with a fitted subject keeps at least one, and the
# largest arms give up the rounding remainder. The draw uses the run's seed, so
# the same call fits the same people.
.model_fit_subset <- function(fitted_ids, subjects, subject_group, cap, seed) {
  fitted_ids <- unique(as.character(fitted_ids))
  if (length(fitted_ids) <= cap) return(NULL)
  arm_of <- stats::setNames(as.character(subject_group), as.character(subjects))
  arm <- arm_of[fitted_ids]
  arm[is.na(arm)] <- "<none>"
  sizes <- table(arm)
  share <- pmax(1L, as.integer(round(as.numeric(sizes) * cap / length(fitted_ids))))
  names(share) <- names(sizes)
  # Rounding can overshoot or undershoot the cap; settle it on the largest arms.
  while (sum(share) != cap) {
    largest <- names(share)[order(-as.numeric(sizes[names(share)]))]
    for (nm in largest) {
      if (sum(share) == cap) break
      if (sum(share) > cap && share[[nm]] > 1L) share[[nm]] <- share[[nm]] - 1L
      if (sum(share) < cap && share[[nm]] < sizes[[nm]]) share[[nm]] <- share[[nm]] + 1L
    }
    if (all(share >= sizes[names(share)]) && sum(share) < cap) break
  }
  draw <- function() {
    unlist(lapply(names(share), function(nm) {
      pool <- fitted_ids[arm == nm]
      if (length(pool) <= share[[nm]]) pool else sample(pool, share[[nm]])
    }), use.names = FALSE)
  }
  chosen <- if (is.null(seed)) draw() else .with_local_seed(seed, draw())
  fitted_ids[fitted_ids %in% chosen]
}

# Reading a converged fit back into the numbers the generator needs, and no
# others. Empirical Bayes estimates are read to report a correlation and are
# then discarded, because they are a description of each real patient.
.model_read_fit <- function(fit, structural, error) {
  fixed <- stats::setNames(
    exp(as.numeric(fit$parFixedDf[paste0("t", .required_pk_params[[structural]]),
                                  "Estimate"])),
    .required_pk_params[[structural]]
  )
  omega <- fit$omega
  rownames(omega) <- sub("^eta\\.", "", rownames(omega))
  colnames(omega) <- rownames(omega)
  keep <- rownames(omega) %in% names(fixed)
  omega <- omega[keep, keep, drop = FALSE]
  residual_value <- as.numeric(
    fit$parFixedDf[paste0(error, ".err"), "Estimate"]
  )
  list(
    fixed = fixed, omega = omega,
    residual = if (identical(error, "prop")) {
      list(kind = "proportional", cv = residual_value)
    } else {
      list(kind = "additive", sd = residual_value)
    },
    # `$eta` carries the subject identifier in its first column. Only the
    # random effects themselves are read, and only to report a correlation.
    etas = {
      frame <- as.data.frame(fit$eta)
      frame[, grepl("^eta\\.", names(frame)), drop = FALSE]
    }
  )
}

# The PD shapes --------------------------------------------------------------
#
# Three time courses with no exposure dependence, fitted by least squares on the
# pooled observations and selected on AIC. Least squares rather than `nlmixr2`
# deliberately: these are three-parameter curves on one endpoint, the fit is
# well conditioned, and routing them through a population fitter would put a
# compiler in the path of every PD endpoint for no gain in what the generator
# then draws.
#
# A PD endpoint driven by concentration is reproduced as a time course that
# happens to resemble the average subject's response. A dataset whose point is
# the exposure-response relationship is not served by this generator, and this
# is where that is true.
# Start values for the exponential PD shape, tried in order.
#
# The first set is the one this function always used: the median of the late
# half as the plateau and of the early half as the baseline. It works on a
# monotone response and fails on one that falls and then recovers, because the
# two medians land close together and the rate has nothing to fit. The rest
# split the time range into thirds, which separates the two ends of a V, and
# offer three time constants across the range.
.pd_exponential_starts <- function(time, value) {
  span <- max(time) - min(time)
  if (!is.finite(span) || span <= 0) span <- 1
  midpoint <- stats::median(time)
  thirds <- stats::quantile(time, c(1 / 3, 2 / 3), names = FALSE, na.rm = TRUE)
  early <- stats::median(value[time <= thirds[1L]])
  late <- stats::median(value[time >= thirds[2L]])
  starts <- list(list(plateau = stats::median(value[time > midpoint]),
                      baseline = stats::median(value[time <= midpoint]),
                      rate = 1 / max(midpoint, 1e-6)))
  for (fraction in c(4, 10, 2)) {
    starts[[length(starts) + 1L]] <- list(plateau = late, baseline = early,
                                          rate = log(2) / (span / fraction))
  }
  Filter(function(start) all(vapply(start, is.finite, logical(1))), starts)
}

# The first line of a condition, for a table cell.
.trimmed_condition <- function(x) {
  text <- conditionMessage(attr(x, "condition"))
  sub("\n.*$", "", trimws(text))
}

# Each subject's own baseline under the chosen shape, and what is left over.
#
# `.pd_profile()` is linear in `baseline`: evaluating it at 0 and at 1 gives the
# intercept and the per-point derivative, so a subject's baseline is
# `sum((y - intercept) * derivative) / sum(derivative^2)` and no optimizer is
# needed. A subject whose derivative is zero everywhere carries no information
# about the baseline and keeps the typical one.
.pd_subject_baselines <- function(pd, typical, time, value, subject) {
  shape <- list(pd = pd)
  at <- function(b) .pd_profile(shape, time, numeric(), numeric(),
                                replace(typical, "baseline", b))
  intercept <- at(0)
  derivative <- at(1) - intercept

  baseline <- rep(NA_real_, length(unique(subject)))
  names(baseline) <- unique(subject)
  residual <- numeric(length(value))
  for (id in names(baseline)) {
    rows <- which(subject == id)
    denominator <- sum(derivative[rows]^2)
    own <- if (denominator > 0) {
      sum((value[rows] - intercept[rows]) * derivative[rows]) / denominator
    } else {
      typical[["baseline"]]
    }
    baseline[[id]] <- own
    residual[rows] <- value[rows] - (intercept[rows] + own * derivative[rows])
  }
  list(baseline = unname(baseline), residual = residual)
}

# A shape is a candidate where it has something left over to be wrong about.
#
# `AIC` on a saturated fit is `-Inf`: two points fit a line exactly, so a line
# would win every two-point comparison and generate a curve through both with no
# residual at all. That is an interpolation presented as a selection. One
# residual degree of freedom is the condition that rules it out, and it is the
# same condition for all three shapes rather than a count chosen per shape.
.pd_estimable <- function(model) {
  if (is.null(model) || inherits(model, "try-error")) return(FALSE)
  aic <- suppressWarnings(try(stats::AIC(model), silent = TRUE))
  df <- suppressWarnings(try(stats::df.residual(model), silent = TRUE))
  !inherits(aic, "try-error") && !inherits(df, "try-error") &&
    is.finite(aic) && is.finite(df) && df >= 1
}

.model_fit_pd <- function(observations, endpoint, shapes = NULL) {
  rows <- observations[observations$endpoint == endpoint &
                         is.finite(observations$dv) &
                         is.finite(observations$aligned), , drop = FALSE]
  if (!nrow(rows)) return(NULL)
  # Study time from the first dose, which is the axis generation evaluates the
  # shape on. Fitting against time after dose instead makes every sample in a
  # daily regimen land at nearly the same place.
  time <- rows$aligned
  value <- rows$dv

  candidates <- list()
  candidate_notes <- list()
  failed <- list()
  constant <- stats::lm(value ~ 1)
  if (.pd_estimable(constant)) {
    candidates$constant <- list(
      pd = "constant",
      typical = c(baseline = unname(stats::coef(constant)[1L])),
      aic = stats::AIC(constant)
    )
  }
  # A slope needs two distinct times to be a slope at all, whatever the number
  # of observations sitting on them.
  linear <- if (length(unique(time)) >= 2L) stats::lm(value ~ time) else NULL
  if (.pd_estimable(linear)) {
    coefficients <- stats::coef(linear)
    candidates$linear <- list(
      pd = "linear",
      typical = c(baseline = unname(coefficients[1L]),
                  slope = unname(coefficients[2L])),
      aic = stats::AIC(linear)
    )
  }
  # Several start sets rather than one. The single median-based set fails on any
  # response that falls and then recovers, which is the shape a turnover
  # endpoint has, and the failure was silent: the candidate simply vanished from
  # the table and a flat line won on AIC against two other flat lines.
  exponential <- NULL
  note <- ""
  starts <- if (length(unique(time)) >= 2L && nrow(rows) >= 4L) {
    .pd_exponential_starts(time, value)
  } else {
    note <- paste0(nrow(rows), " observation(s) at ", length(unique(time)),
                   " distinct time(s): three parameters and nothing left over")
    list()
  }
  for (start in starts) {
    attempt <- try(stats::nls(
      value ~ plateau + (baseline - plateau) * exp(-rate * pmax(time, 0)),
      start = start
    ), silent = TRUE)
    if (.pd_estimable(attempt)) {
      exponential <- attempt
      break
    }
    note <- if (inherits(attempt, "try-error")) .trimmed_condition(attempt) else
      "fitted exactly, with nothing left over to select it on"
  }
  if (!is.null(exponential)) {
    candidates$exponential <- list(
      pd = "exponential", typical = stats::coef(exponential),
      aic = stats::AIC(exponential)
    )
  } else {
    failed <- list(exponential = note)
  }
  if (!is.null(shapes) && endpoint %in% names(shapes)) {
    chosen <- candidates[[shapes[[endpoint]]]]
    if (is.null(chosen)) {
      stop("`pd` names the shape `", shapes[[endpoint]], "` for endpoint `",
           endpoint, "`, which could not be fitted to it.", call. = FALSE)
    }
  } else if (length(candidates)) {
    chosen <- candidates[[which.min(vapply(candidates, function(c) c$aic,
                                           numeric(1)))]]
  } else {
    # One observation, or several at one time: nothing about a time course can
    # be read, and the endpoint is a level. That is still worth generating --
    # the alternative is a synthetic study missing an endpoint the source has --
    # so it is a constant at the mean, with whatever spread the observations
    # show and none where they show none. Named in the candidate table as the
    # only thing that could be fitted, so the report says so rather than
    # implying a search happened.
    chosen <- list(pd = "constant", typical = c(baseline = mean(value)),
                   aic = NA_real_)
    candidates$constant <- chosen
    candidate_notes$constant <- paste0(
      nrow(rows), " observation(s) at ", length(unique(time)),
      " distinct time(s): a level, with no time course to select"
    )
  }

  # Between-subject variability and residual error, both read around each
  # subject's own curve rather than around the population one.
  #
  # Taking the residual as the spread of every point about the typical curve
  # puts all of the between-subject variation into it, and generation then emits
  # that as independent noise per observation. On a study where subjects share a
  # shape but sit at different levels, almost the whole of the structure comes
  # back as scatter and no synthetic subject has a profile at all. Splitting the
  # two is what makes the generated profiles profiles.
  #
  # The split is exact because `.pd_profile()` is linear in `baseline` for every
  # shape, so a subject's own baseline is a least-squares projection and needs
  # no second optimizer. `generation` varies `baseline` and nothing else, so
  # this is the same quantity the generator draws.
  levels <- .pd_subject_baselines(chosen$pd, chosen$typical, time, value,
                                  rows$subject)
  positive <- levels$baseline[is.finite(levels$baseline) &
                                levels$baseline > 0]
  chosen$baseline_cv <- if (length(positive) > 1L) stats::sd(log(positive)) else 0
  # `sd()` of one number is `NA`, and an `NA` residual reaches generation as an
  # `NA` observation. One observation shows no scatter, which is a measurement
  # this study cannot make rather than a claim that the endpoint is noiseless.
  spread <- if (length(levels$residual) > 1L) stats::sd(levels$residual) else 0
  chosen$residual <- list(kind = "additive",
                          sd = if (is.finite(spread)) spread else 0)
  chosen$candidates <- data.frame(
    shape = c(names(candidates), names(failed)),
    converged = c(rep(TRUE, length(candidates)), rep(FALSE, length(failed))),
    aic = c(vapply(candidates, function(c) c$aic, numeric(1)),
            rep(NA_real_, length(failed))),
    note = c(vapply(names(candidates), function(name) {
      candidate_notes[[name]] %||% ""
    }, character(1)), unlist(failed) %||% character()),
    row.names = NULL, stringsAsFactors = FALSE
  )
  chosen
}

# One shape per arm, when the caller asks for it (`SIM-060` option (a)).
#
# The pooled fit predicts one number for every arm, so a synthetic patient's
# dose does not reach their response: on `case1_pkpd` the source's mean PD after
# 1500 h runs 68, 70, 126, 237, 341 across placebo and five ascending doses, and
# the pooled shape answers 149 to all six. Fitting the shape per arm is the same
# per-arm summary the dosing and visit models already are, and it asserts no
# dose-response form -- an arm is fitted on its own observations or not at all.
#
# It is not the default, because it is a different claim about the study: the
# pooled shape says these endpoints are a time course, and the per-arm shape
# says each arm has its own. An arm can also carry too little of an endpoint to
# fit, and where it does it keeps the pooled shape rather than losing the
# endpoint, so no arm is left with nothing to generate.
.model_fit_pd_arms <- function(observations, endpoint, shapes, subject_group,
                               pooled) {
  arms <- .unique_in_order(unname(subject_group))
  fits <- lapply(arms, function(arm) {
    subjects <- names(subject_group)[subject_group == arm]
    rows <- observations[observations$subject %in% subjects, , drop = FALSE]
    .model_fit_pd(rows, endpoint, shapes) %||% pooled
  })
  stats::setNames(fits, arms)
}

# The covariate the allometric scaling is applied to, or nothing. Weight-like by
# name, positive, and numeric -- there is no way to recognise a body weight from
# its values alone, and guessing from a distribution would be worse than asking.
#
# Applied where it is found rather than tested for improvement. Testing costs a
# second fit, which is the whole budget of the default path, and allometry on a
# weight is a shape this generator asserts rather than a hypothesis it examines.
# `covariate_effects = "none"` switches it off.
#
# The cost of the default is explicit. A covariate that influences the real
# profiles and is not in the model is generated independently of them, so the
# synthetic data carries no relationship between the two. `synpmx_avatar()`
# preserves those relationships without modelling them, because a blended
# subject's covariates and profile come from the same donors. `model_report()`
# reports the correlation between each declared covariate and the individual
# random effects, which is where an unmodelled relationship shows up.
.model_weight_covariate <- function(source, roles) {
  for (covariate in roles$covariates) {
    values <- suppressWarnings(as.numeric(source[[covariate]]))
    values <- values[is.finite(values)]
    if (!length(values) || any(values <= 0)) next
    # A baseline weight is rarely called `WT`. `WEIGHTB`, `WTBL` and `BWT` are
    # all the same column with a suffix or a prefix saying it was measured at
    # baseline, and refusing them means silently fitting no covariate at all on
    # a study that declared one. Anything not on this list is not scaled, which
    # is the safe direction: `HEIGHT` and `AGE` must not match.
    if (grepl("^(b|base|baseline)?_?(wt|wgt|weight|bw|bodywt|body_?weight)_?(b|bl|base|baseline|0)?$",
              covariate, ignore.case = TRUE)) {
      return(list(covariate = covariate, reference = stats::median(values)))
    }
  }
  NULL
}

# The per-subject weight, on the estimation dataset.
.model_attach_weight <- function(data, source, roles, weight) {
  if (is.null(weight)) return(data)
  by_subject <- vapply(split(source[[weight$covariate]],
                             as.character(source[[roles$id]])),
                       function(x) suppressWarnings(as.numeric(x[1L])),
                       numeric(1))
  data[[weight$covariate]] <- unname(by_subject[data$ID])
  if (any(!is.finite(data[[weight$covariate]]))) return(NULL)
  data
}

# What the assay limit cost, per endpoint.
#
# The concentration's censored rows go to `nlmixr2` as censored (`SIM-068`):
# the population fit has a likelihood that can express "below this limit", and
# using it is both more nearly right and materially different -- measured on 30
# subjects of `case1_pkpd`, 45% of them below the limit, the residual error
# falls from 0.517 to 0.290 against a uniform draw inside the censoring region,
# and residual error is what sets the scatter of everything generated. It costs
# a factor of 1.7 in fitting time there, and less where less of the endpoint is
# censored.
#
# Everything else still reads the imputed value, because nothing else has a
# likelihood: the apparatus, the covariate model and the PD shapes are `lm()`
# and `nls()`, and a stack of identical boundary substitutions would bend each
# of them toward the limit. So the share below the limit is still worth
# reporting, and this measures it.
# The floor below which the generator will not emit a value, per endpoint.
#
# A study that declares a censoring column says where its assay stopped, and
# `.censor_latent()` puts that boundary back at emit. A study that declares none
# still had an assay, and its smallest reported value is the only evidence of
# where the limit sat. Without a floor the residual draw can put a synthetic
# concentration orders of magnitude below anything the study could have
# measured, which on a log axis is the whole of what makes a figure look wrong.
#
# Half the smallest positive value reported, which is where a below-the-limit
# value is conventionally substituted, and which cannot sit above anything the
# study actually reported.
#
# Only for an endpoint that lives on a positive scale, because an endpoint
# recording a zero is recording something a floor would contradict and a PD
# score with a true zero is the ordinary case of that. What decides that is the
# *share* of non-positive values rather than whether any exists: a concentration
# assay reports a reading near its limit as a small negative number now and
# again, and one such row in three hundred is noise around the limit rather than
# a statement that this endpoint reaches zero. An endpoint whose scale really
# includes zero says so in many rows, not one.
.model_nonpositive_tolerance <- 0.05

.model_assay_floor <- function(values) {
  values <- values[is.finite(values)]
  if (!length(values)) return(NULL)
  if (mean(values <= 0) >= .model_nonpositive_tolerance) return(NULL)
  positive <- values[values > 0]
  if (!length(positive)) return(NULL)
  min(positive) / 2
}

# The floor exists only where the study declared no limit of its own for that
# endpoint. Where it declared one, `.censor_latent()` puts the real boundary
# back at emit and a second floor underneath it would be counted and warned
# about while changing nothing. A `cens` column is per-endpoint evidence, not
# per-study: `case1_pkpd` declares one and sets it on the concentration only, so
# its PD endpoint is a study that declared no limit and does get a floor.
.model_quantification_floor <- function(source, roles, endpoints) {
  observed <- .observation_rows(source, roles, require_present = TRUE)
  endpoint <- .endpoint(source, roles)
  dv <- suppressWarnings(as.numeric(source[[roles$dv]]))
  censored <- if (is.null(roles$cens)) rep(FALSE, nrow(source)) else {
    flag <- suppressWarnings(as.numeric(as.character(source[[roles$cens]])))
    is.finite(flag) & flag != 0
  }
  floors <- lapply(endpoints, function(name) {
    at <- observed & endpoint == name
    if (any(at & censored)) return(NULL)
    .model_assay_floor(dv[at])
  })
  names(floors) <- endpoints
  floors <- floors[!vapply(floors, is.null, logical(1))]
  if (!length(floors)) NULL else floors
}

.model_censoring_summary <- function(source, roles, endpoints) {
  if (is.null(roles$cens)) return(NULL)
  observed <- .observation_rows(source, roles, require_present = TRUE)
  endpoint <- .endpoint(source, roles)
  cens <- suppressWarnings(as.numeric(as.character(source[[roles$cens]])))
  dv <- suppressWarnings(as.numeric(source[[roles$dv]]))
  rows <- lapply(endpoints, function(name) {
    at <- observed & endpoint == name
    if (!any(at)) return(NULL)
    censored <- at & is.finite(cens) & cens != 0
    data.frame(
      endpoint = name, observations = sum(at), censored = sum(censored),
      fraction = sum(censored) / sum(at),
      limit = if (any(censored)) stats::median(dv[censored]) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Estimate a population model from a trial
#'
#' The only stage that reads patient data, and the only one that needs
#' `nlmixr2`. It works out which endpoint is the drug concentration and what
#' design produced it, fits the candidate models that design admits, picks one
#' on AIC, and returns that fit alongside the dosing and visit models the
#' generated subjects are built from. No patient row survives it.
#'
#' **The fitted parameters are not estimates to report.** They exist to make
#' simulated profiles look like the source study. The candidate set is five
#' linear models and the covariate model is allometric scaling or nothing, which
#' is too little to answer a scientific question, and the object prints that
#' warning with itself because its contents look exactly like the output of a
#' real population analysis.
#'
#' `nominal_time` is required, for two reasons. The dosing and visit models sit
#' on the nominal grid, and a grid inferred from recorded times is a statement
#' about the protocol only the caller can make. Estimation, separately, reads
#' the recorded times and the recorded dosing history, because a population fit
#' is a statement about the dose that was actually given.
#'
#' No formal privacy guarantee is offered. No patient's measured value reaches
#' the output, which is the claim `synpmx_pca()` makes and is stronger than
#' `synpmx_avatar()`'s, but the fixed effects and the covariance matrix are
#' functions of the individuals in the source and neither is noised. The cohort
#' floor is the whole defence and it is a threshold rather than an accounting.
#'
#' @param data Source PMX event data.
#' @param roles Explicit column roles from [pmx_roles()], including
#'   `nominal_time`.
#' @param pk One of the five built-in structural models, forcing it and skipping
#'   the search. `NULL` searches the candidates the design admits.
#' @param pd Named character vector of PD shapes per endpoint, skipping that
#'   search. One of `"constant"`, `"linear"` or `"exponential"` each.
#' @param pd_by_arm Fit each PD endpoint's shape per arm rather than once over
#'   the pooled cohort. `FALSE`, the default, gives every arm the same time
#'   course, so a synthetic patient's dose does not reach their response.
#'   `TRUE` gives each arm its own shape, selected on AIC within that arm and
#'   asserting no dose-response form; an arm holding too little of an endpoint
#'   to fit keeps the pooled shape. Costs nothing measurable, since these are
#'   least-squares fits.
#' @param endpoint_roles Which endpoint is the drug concentration, as
#'   `c(pk = "cp")`, overriding the inference.
#'
#'   **More than one may be named**, for a study that measures two
#'   concentrations — two drugs, or a parent and its metabolite. Each gets its
#'   own structural model, its own parameters and its own residual error,
#'   evaluated against the one dose schedule they share; no correlation between
#'   their random effects is estimated. Write it as
#'   `list(pk = c("parent", "metabolite"))`, or `c(pk = c("parent",
#'   "metabolite"))`, which R renames to `pk1`/`pk2` and which is read the same
#'   way.
#'
#'   Naming several is a declaration, never an inference. Left to itself the
#'   classification picks one concentration and treats every other continuous
#'   endpoint as a pharmacodynamic time course, because a second endpoint that
#'   passes the concentration signals is at least as often a biomarker as a
#'   metabolite — `onc_sim`'s tumour size passes them. Where a demoted endpoint
#'   really is a concentration, that time course has no dose term in it and the
#'   generated values lose their dose ordering, which is why naming both is
#'   worth doing.
#' @param covariate_effects `"auto"` fits allometric scaling on clearance and
#'   volume where a weight-like covariate is declared and keeps it where it
#'   improves AIC. `"none"` fits nothing.
#' @param min_subjects Cohort size the fit should have. Below it the covariance
#'   matrix describes the subjects it was fitted to rather than a population,
#'   which warns rather than refuses: the fit runs on whatever the study has.
#' @param min_arm_patients Minimum patients in every arm, as
#'   [synpmx_pca_summarize()] uses. Patients in a shorter arm are dropped with a
#'   warning before anything is fitted, so that arm is absent from the fitted
#'   model and from the data generated from it.
#' @param start_param Starting values for the population PK parameters, as
#'   `start_param = c(cl = 4, v = 40, ka = 0.5)`. Anything not named is read
#'   off the curve as usual, so a caller who knows the clearance and not the
#'   absorption says only the clearance.
#'
#'   Where the study fits **more than one** concentration endpoint, key it by
#'   endpoint — `start_param = list(parent = c(cl = 4), metabolite =
#'   c(cl = 9))` — because a flat vector would not say which one it describes,
#'   and is refused with that message. With one concentration the flat form is
#'   what to write. The PD endpoints are least-squares time courses that this
#'   does not reach.
#'
#'   The automatic starting values are a non-compartmental read of the cohort's
#'   median profile. On a study that read cannot describe — sampled only at
#'   troughs, dosed by two routes whose reads disagree, or recorded in units
#'   that are not what they appear — it can start the search somewhere the
#'   optimizer cannot leave, which `model_report()` then reports as a fit that
#'   did not move. This is how a caller who knows the compound fixes that in
#'   advance, and `model_report()` says which values were declared.
#' @param max_fit_subjects Most subjects the population model is fitted to,
#'   60 by default. A study above it has its PK parameters estimated on a subset
#'   drawn in proportion to the arms, using `seed`; the dosing, visit and
#'   covariate models still read every subject. Fit time is linear in subjects
#'   and the parameters a synthetic study needs are settled well before the
#'   sixtieth, so this is where a twenty-minute fit becomes a five-minute one.
#'   `model_report()` states the count. Cannot be below `min_subjects`.
#' @param min_time_bins Distinct nominal times after a dose the cohort should
#'   hold. Below it a one-compartment model is not identifiable, which warns
#'   rather than refuses: the fit runs and its parameters sit close to their
#'   starting values. No post-dose observation at all is an error.
#' @param estimation Passed to `nlmixr2`. `"focei"` by default because the
#'   selection criterion is AIC and `"saem"` does not reliably produce one at
#'   these cohort sizes.
#' @param seed Seed for the one random step, which is imputing censored values
#'   before the fit.
#' @param quiet Suppress the per-candidate progress messages.
#'
#' @return A `pmx_fitted_model`.
#' @seealso [synpmx_model_generate()], [synpmx_model()], [model_report()],
#'   [model_candidates()], [model_parameters()].
#' @export
synpmx_model_estimate <- function(data, roles, pk = NULL, pd = NULL,
                                  pd_by_arm = FALSE,
                                  endpoint_roles = NULL, start_param = NULL,
                                  covariate_effects = "auto",
                                  min_subjects = 20L, min_arm_patients = 3L,
                                  min_time_bins = 6L, max_fit_subjects = 60L,
                                  estimation = "focei",
                                  seed = NULL, quiet = FALSE) {
  if (!inherits(roles, "pmx_roles")) {
    stop("`roles` must come from `pmx_roles()`.", call. = FALSE)
  }
  max_fit_subjects <- .positive_integer(max_fit_subjects, "max_fit_subjects")
  if (max_fit_subjects < .positive_integer(min_subjects, "min_subjects")) {
    stop("`max_fit_subjects` (", max_fit_subjects, ") is below `min_subjects` ",
         "(", min_subjects, "): the fit cannot be capped under the floor it ",
         "is told to hold.", call. = FALSE)
  }
  # Arguments are checked before the fitter is required, so that a
  # mistyped model name reads as a mistyped model name rather than as a
  # missing suggested package.
  if (!is.null(pk)) {
    unknown <- setdiff(pk, .pk_models)
    if (length(unknown)) {
      stop("`pk` names model(s) outside the closed-form set: ",
           paste(unknown, collapse = ", "), ". Available: ",
           paste(.pk_models, collapse = ", "), ".", call. = FALSE)
    }
  }
  if (!identical(covariate_effects, "auto") &&
      !identical(covariate_effects, "none")) {
    stop("`covariate_effects` must be \"auto\" or \"none\".", call. = FALSE)
  }
  if (!is.logical(pd_by_arm) || length(pd_by_arm) != 1L || is.na(pd_by_arm)) {
    stop("`pd_by_arm` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!requireNamespace("nlmixr2est", quietly = TRUE)) {
    stop(.condition_text(
      "`synpmx_model_estimate()` needs the nlmixr2 package, which is in ",
      "Suggests.",
      fix = paste("Install it, or use `synpmx_avatar()` or `synpmx_pca()`,",
                  "which fit no structural model.")), call. = FALSE)
  }
  started <- proc.time()[["elapsed"]]
  data <- as.data.frame(data)
  source <- data[, intersect(.retained_role_columns(roles), names(data)),
                 drop = FALSE]
  source <- pmx_expand_doses(source, roles)

  # Arms too small to summarize leave first, so every gate, count and fit below
  # sees the cohort that will actually be modelled.
  subjects <- .unique_in_order(source[[roles$id]])
  subject_group <- .model_subject_arms(source, roles)
  keep <- .drop_short_arms(subject_group, min_arm_patients,
                           "synpmx_model_estimate()")
  if (!all(keep)) {
    subject_group <- subject_group[keep]
    source <- source[source[[roles$id]] %in% subjects[keep], , drop = FALSE]
  }

  n_source <- length(.unique_in_order(source[[roles$id]]))
  .model_note_subjects(n_source, min_subjects)
  .model_require_nominal_time(source, roles)
  .model_require_time_coverage(source, roles, min_time_bins)

  censoring_source <- source
  source <- if (is.null(seed)) .impute_censored(source, roles) else
    .with_local_seed(seed, .impute_censored(source, roles))

  observations <- .model_observations(source, roles)
  classified <- .model_classify_endpoints(source, roles, observations,
                                          endpoint_roles)
  # One design read per concentration endpoint: the route is a property of the
  # study, but the sampling richness that prunes the candidate set is a property
  # of the endpoint, and a metabolite is not always sampled like its parent.
  designs <- stats::setNames(lapply(classified$pk, function(endpoint) {
    one <- .model_detect_design(source, roles, observations, endpoint)
    if (!is.null(pk)) {
      one$candidates <- pk
      one$reason <- if (length(pk) == 1L) {
        paste0("declared through `pk = \"", pk, "\"`")
      } else {
        paste0("searched over the ", length(pk), " models named in `pk`")
      }
    }
    one
  }), classified$pk)
  design <- designs[[1L]]

  # The apparatus, on the nominal grid, exactly as `synpmx_pca_summarize()`
  # builds it. `.model_cells()` is this generator's adapter over the same
  # `.arm_models()` the other one calls.
  fittable <- c(classified$pk, classified$pd)
  planned <- source
  planned[[roles$time]] <-
    suppressWarnings(as.numeric(source[[roles$nominal_time]]))
  cells <- .model_cells(source, roles, c(fittable, classified$discrete),
                        min_arm_patients)
  arm_models <- .arm_models(planned, roles, cells, subject_group,
                            min_arm_patients)

  # `censoring_source` rather than `source`: the fitter is given the study's own
  # values with their censoring flags, while the summaries below keep the
  # imputed ones. The two frames differ only in the censored `DV` values.
  recorded_data <- .model_estimation_data(censoring_source, roles,
                                          classified$pk)
  # Who that left out. Said plainly, because the fit is then a statement about
  # fewer people than the study has, and a warning where it takes the fitted
  # cohort under the floor the run was told to hold.
  fitted_subjects <- length(unique(recorded_data$ID))
  # The cap. Fit time is linear in subjects and the inner per-subject loop is
  # where it goes -- measured 2026-09-08 on `onc_sim`, 40 patients fit in four
  # minutes and 200 in twenty-five -- while what the generator needs from the
  # fit, population parameters that put synthetic values where the real ones
  # are, is settled long before the two-hundredth patient. So the population
  # model is fitted to a subset drawn in proportion to the arms, and every
  # other model -- dosing, visits, covariates, cells -- still reads the whole
  # study. The report says so, because a reader who sees 200 patients and a
  # fit on 60 is owed the number.
  fit_draw <- .model_fit_subset(recorded_data$ID, subjects, subject_group,
                                max_fit_subjects, seed)
  if (!is.null(fit_draw)) {
    recorded_data <- recorded_data[recorded_data$ID %in% fit_draw, ,
                                   drop = FALSE]
  }
  fit_subjects <- list(fitted = length(unique(recorded_data$ID)),
                       of = fitted_subjects, cap = max_fit_subjects)
  if (!quiet && !is.null(fit_draw)) {
    message(.wrap_plain(paste0(
      "Fitting the population model to ", fit_subjects$fitted, " of ",
      fitted_subjects, " patients (`max_fit_subjects` = ", max_fit_subjects,
      "), drawn in proportion to the arms; the dosing, visit and covariate ",
      "models read all ", n_source, ".")))
  }
  if (fitted_subjects < n_source) {
    note <- paste0(n_source - fitted_subjects, " of ", n_source,
                   " subjects have no `", classified$pk,
                   "` observation and are not fitted; their dosing and visits ",
                   "still reach the arm models.")
    # A warning only where it is news: a cohort already under the floor was
    # warned about before the endpoints were classified.
    if (fitted_subjects < min_subjects && n_source >= min_subjects) {
      warning(note, " That leaves ", fitted_subjects,
              " subjects under `min_subjects` = ", min_subjects,
              ", so the covariance describes those subjects rather than a ",
              "population.", call. = FALSE)
    } else if (!quiet) {
      message(note)
    }
  }
  # A source that arrived compressed is folded back by runs, which is exact and
  # keeps a dose change as its own record; one that wrote every dose out is
  # compressed only where the whole schedule is regular.
  estimation_data <- if (is.null(roles$addl)) {
    .compress_dose_schedule(recorded_data)
  } else {
    .compress_dose_runs(recorded_data)
  }
  # Proportional error unless the concentration lives on a scale that includes
  # zero. A handful of non-positive readings does not make that scale: they are
  # what the assay returns near its limit, and treating them as evidence costs
  # the whole error model. `nimoData` reports one negative concentration in 321,
  # and reading that one row as "this endpoint reaches zero" fitted an additive
  # residual of 1.46 to values whose median is 3 -- which then generated a
  # cohort scattered from zero upward. Below the limit they are substituted at
  # the same floor the generator will not emit below, which is the LLOQ/2
  # convention applied to a value the assay reported as if it were a reading.
  values <- estimation_data$DV[!is.na(estimation_data$DV)]
  assay_floor <- .model_assay_floor(values)
  error <- if (is.null(assay_floor)) "add" else "prop"
  substituted <- if (is.null(assay_floor)) integer(0) else
    which(!is.na(estimation_data$DV) & estimation_data$DV <= 0)
  if (length(substituted)) {
    estimation_data$DV[substituted] <- assay_floor
    if (!quiet) {
      message(length(substituted), " of ", length(values), " `",
              classified$pk, "` observations are not positive and were fitted ",
              "at ", signif(assay_floor, 4),
              ", half the smallest positive value the study reports.")
    }
  }

  # Allometric scaling is folded into the fit rather than compared against one
  # without it, so the default path performs exactly one fit.
  weight <- if (identical(covariate_effects, "auto")) {
    .model_weight_covariate(source, roles)
  } else NULL
  with_weight <- .model_attach_weight(estimation_data, source, roles, weight)
  if (is.null(with_weight)) weight <- NULL else estimation_data <- with_weight

  if (!quiet) {
    note <- .dose_record_message(recorded_data, estimation_data)
    if (!is.null(note)) message(note)
  }
  starts_by_endpoint <- .model_split_start_param(start_param, classified$pk)

  # One population model per concentration endpoint. They share the dosing
  # records -- the same doses drive a parent and its metabolite -- and nothing
  # else: each has its own structural model, its own parameters and its own
  # residual error, and no correlation between their random effects is
  # estimated, which is a limitation worth knowing rather than a claim.
  pk_models <- stats::setNames(lapply(classified$pk, function(endpoint) {
    own_design <- designs[[endpoint]]
    own_data <- if (identical(endpoint, classified$pk[[1L]])) estimation_data else {
      one <- .model_estimation_data(censoring_source, roles, endpoint)
      one <- if (is.null(roles$addl)) .compress_dose_schedule(one) else
        .compress_dose_runs(one)
      attached <- .model_attach_weight(one, source, roles, weight)
      if (is.null(attached)) one else attached
    }
    if (!quiet) {
      message("Fitting ", length(own_design$candidates), " model(s) for `",
              endpoint, "` (", own_design$route, ")",
              if (!is.null(weight)) paste0(", allometric on ", weight$covariate),
              ":")
    }
    own_start <- .model_validate_start_param(starts_by_endpoint[[endpoint]],
                                             own_design$candidates)
    search <- .model_fit_candidates(own_data, own_design$candidates,
                                    observations, endpoint, error,
                                    estimation, quiet, weight, own_start)
    selected <- search$selected
    parameters <- .model_read_fit(search$fits[[selected]], selected, error)
    movement <- .model_fit_movement(parameters$fixed, parameters$omega,
                                    search$starts[[selected]])
    if (!movement$moved) .model_warn_unmoved(selected, movement)
    list(endpoint = endpoint, structural = selected, parameters = parameters,
         candidates = search$table, movement = movement, design = own_design,
         start_param = own_start,
         seconds = sum(search$table$seconds, na.rm = TRUE),
         effects = if (is.null(weight)) list() else stats::setNames(
           lapply(intersect(names(parameters$fixed),
                            names(.model_allometric_exponents)),
                  function(parameter) {
                    list(covariate = weight$covariate,
                         reference = weight$reference,
                         exponent = unname(
                           .model_allometric_exponents[[parameter]]))
                  }),
           intersect(names(parameters$fixed),
                     names(.model_allometric_exponents))))
  }), classified$pk)

  primary <- pk_models[[1L]]
  selected <- primary$structural
  parameters <- primary$parameters
  movement <- primary$movement
  effects <- primary$effects
  fit_seconds <- sum(vapply(pk_models, function(m) m$seconds, numeric(1)))

  pd_seconds <- stats::setNames(numeric(length(classified$pd)), classified$pd)
  pd_fits <- stats::setNames(lapply(classified$pd, function(endpoint) {
    started_pd <- proc.time()[["elapsed"]]
    pooled <- .model_fit_pd(observations, endpoint, pd)
    if (!is.null(pooled) && pd_by_arm) {
      pooled$arms <- .model_fit_pd_arms(observations, endpoint, pd,
                                        subject_group, pooled)
    }
    pd_seconds[[endpoint]] <<- proc.time()[["elapsed"]] - started_pd
    pooled
  }), classified$pd)
  pd_fits <- pd_fits[!vapply(pd_fits, is.null, logical(1))]
  # An endpoint nothing could be fitted to (`SIM-069`). There is no minimum
  # number of observations -- one is a level, and a level is generated -- so
  # this is the endpoint with no usable row at all: no value, or no time to fit
  # against. The gates make that hard to reach, since `nominal_time` is required
  # on every observation before this runs, so treat it as the guard it is. What
  # it guards against is silence: such an endpoint stays in the source, in the
  # visit grid and in the schema, and generation reaches its cells, finds no
  # shape and no discrete marginal, and emits nothing, so the synthetic study
  # would be missing an endpoint the source has with nothing said.
  unfitted <- setdiff(classified$pd, names(pd_fits))
  if (length(unfitted)) {
    counts <- vapply(unfitted, function(endpoint) {
      sum(observations$endpoint == endpoint & is.finite(observations$dv) &
            is.finite(observations$aligned))
    }, integer(1))
    warning(.condition_text(
      "`synpmx_model_estimate()` fitted no shape to ", length(unfitted),
      " endpoint(s):",
      items = sprintf("`%s` (%d usable observation(s))", unfitted, counts),
      why = paste("Every row of the endpoint is missing a value or a time on",
                  "the nominal grid. Those endpoints are absent from the",
                  "generated data entirely, rather than generated badly.")),
      call. = FALSE)
  }

  # The empirical Bayes estimates are one row per FITTED subject, in source
  # order, so the frame they are set beside has to be the same subjects in the
  # same order.
  fitted_ids <- unique(recorded_data$ID)
  in_fit <- as.character(source[[roles$id]]) %in% fitted_ids
  correlations <- .model_covariate_correlations(
    source[in_fit, , drop = FALSE], roles,
    subject_group[as.character(subjects) %in% fitted_ids], parameters$etas)
  parameters$etas <- NULL

  fitted <- .pmx_fitted_model(
    structural = selected, candidates = primary$candidates,
    parameters = parameters, pk_models = pk_models,
    movement = movement, fit_subjects = fit_subjects,
    start_param = start_param,
    endpoints = list(pk = classified$pk, pd = names(pd_fits),
                     discrete = classified$discrete, signals = classified$signals,
                     decided_by = classified$decided_by),
    arms = list(arms = arm_models$arms, sizes = arm_models$sizes),
    dosing = arm_models$dosing, visits = arm_models$visits,
    schema = .source_schema(censoring_source, roles, fittable, subject_group),
    roles = roles,
    settings = list(min_subjects = min_subjects,
                    min_arm_patients = min_arm_patients,
                    min_time_bins = min_time_bins,
                    max_fit_subjects = max_fit_subjects,
                    estimation = estimation,
                    covariate_effects = covariate_effects, error = error),
    n_source = n_source,
    cells = cells, pd = pd_fits, covariate_effects = effects,
    covariates = .covariate_model(source, roles, subject_group),
    discrete = .discrete_model(source, roles, cells, subject_group),
    design = design, correlations = correlations,
    censoring = .model_censoring_summary(censoring_source, roles, fittable),
    quantification_floor = .model_quantification_floor(censoring_source, roles,
                                                       fittable),
    # Reported with the fit because it is the number a caller weighs a rerun
    # against: `fit` is what the fitter took, `total` what the call took.
    # Broken down rather than totalled. A fit is the slow thing this package
    # does, and "ten minutes" is not actionable while "nine of them in the
    # second candidate" is: the answer to that is `pk = "1cmt_oral"`.
    timing = list(fit = fit_seconds,
                  total = proc.time()[["elapsed"]] - started,
                  # Every candidate of every concentration endpoint, so a
                  # study fitting two of them accounts for both waits.
                  candidates = do.call(rbind, lapply(pk_models, function(m) {
                    row <- m$candidates[, c("model", "converged", "seconds")]
                    if (length(pk_models) > 1L) {
                      row$model <- paste0(m$endpoint, ": ", row$model)
                    }
                    row
                  })),
                  pd = pd_seconds)
  )
  if (!quiet) {
    message("Fitted in ", .model_duration(fitted$timing$fit), "; ",
            .model_duration(fitted$timing$total), " in total.")
  }
  fitted
}

# Named by subject, always. `vapply()` names its result from a character input
# and leaves a numeric one unnamed, so a study whose `ID` is a number produced
# an unnamed vector and anything reading "which subjects are in this arm" off
# the names got nothing -- silently, since every arm then looks empty rather
# than wrong.
.model_subject_arms <- function(source, roles) {
  subjects <- .unique_in_order(source[[roles$id]])
  strata_key <- as.character(.subject_strata(source, roles))
  arms <- vapply(subjects, function(subject) {
    rows <- which(!is.na(source[[roles$id]]) & source[[roles$id]] == subject)
    strata_key[rows[1L]]
  }, character(1))
  stats::setNames(arms, as.character(subjects))
}

# Where an unmodelled covariate relationship shows up. The random effects are
# per-subject quantities and are not stored; only these correlations are, which
# is a number per covariate and parameter rather than a description of anybody.
.model_covariate_correlations <- function(source, roles, subject_group, etas) {
  if (is.null(etas) || !length(roles$covariates)) return(NULL)
  subjects <- .unique_in_order(source[[roles$id]])
  first_row <- vapply(subjects, function(subject) {
    which(!is.na(source[[roles$id]]) & source[[roles$id]] == subject)[1L]
  }, integer(1))
  rows <- list()
  for (covariate in roles$covariates) {
    values <- suppressWarnings(as.numeric(source[[covariate]][first_row]))
    if (sum(is.finite(values)) < 3L) next
    for (parameter in names(etas)) {
      eta <- etas[[parameter]]
      if (length(eta) != length(values)) next
      correlation <- suppressWarnings(
        stats::cor(values, eta, use = "complete.obs")
      )
      if (!is.finite(correlation)) next
      rows[[length(rows) + 1L]] <- data.frame(
        covariate = covariate, parameter = sub("^eta\\.", "", parameter),
        correlation = correlation, stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Estimate a population model and generate a synthetic dataset from it
#'
#' A single call for [synpmx_model_estimate()] followed by
#' [synpmx_model_generate()]. Use the two separately to look at the fit before
#' generating from it; it is on the result either way, as the `pmx_fitted_model`
#' attribute.
#'
#' @param data Source PMX event data.
#' @param roles Explicit column roles from [pmx_roles()], including
#'   `nominal_time`.
#' @param n_subjects Number of synthetic subjects. Defaults to the source count.
#' @param seed Seed, used for both stages.
#' @param ... Passed to [synpmx_model_estimate()].
#'
#' @return A data frame in the source's shape, carrying the fitted model as an
#'   attribute.
#' @seealso [synpmx_model_estimate()], [synpmx_model_generate()],
#'   [synpmx_pca()], [synpmx_avatar()].
#' @export
synpmx_model <- function(data, roles, n_subjects = NULL, seed = NULL, ...) {
  synpmx_model_generate(
    synpmx_model_estimate(data, roles, seed = seed, ...),
    n_subjects = n_subjects, seed = seed
  )
}
