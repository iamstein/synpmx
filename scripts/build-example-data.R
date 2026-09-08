# Build the two simulated studies that ship as package data, `mixroute_sim` and
# `onc_sim`, and write them to `data/`.
#
#   Rscript scripts/build-example-data.R
#
# Both exist to exercise shapes the public-data survey has no example of. Both
# are simulated from a written-down model, so unlike the eight public datasets
# the truth is known and the survey can say how close a generator came to it.
#
# Neither is a substitute for real data, and neither is fitted here: this script
# writes the source studies, and `scripts/build-model-fits.R` fits them.

set.seed(20260908)

# ---------------------------------------------------------------- mixroute_sim
#
# Three arms, one drug, two routes: intravenous only, subcutaneous only, and an
# intravenous loading dose followed by subcutaneous maintenance. That third arm
# is the point -- a patient who receives both routes is what a mixed model is
# for, and no public dataset in the survey has an administration column at all.
#
# One compartment with first-order absorption. Bioavailability is 0.7 and is
# identifiable here precisely because the same study doses both ways.

mix_pop <- c(cl = 2, v = 10, ka = 0.5)   # L/day, L, 1/day
MIX_F <- 0.7
mix_bsv <- c(cl = 0.30, v = 0.25, ka = 0.40)
MIX_PROP_ERR <- 0.15
MIX_DOSE <- 100                           # mg
MIX_LLOQ <- 0.05                          # mg/L

# One dose's contribution, by route. Superposition is exact because the model is
# linear in dose, which is the same reason `synpmx` can sum a mixed regimen.
mix_one_dose <- function(t, dose, route, p) {
  ke <- p[["cl"]] / p[["v"]]
  out <- numeric(length(t))
  live <- t >= 0
  if (!any(live)) return(out)
  tt <- t[live]
  out[live] <- if (identical(route, "iv")) {
    dose / p[["v"]] * exp(-ke * tt)
  } else {
    MIX_F * dose * p[["ka"]] / (p[["v"]] * (p[["ka"]] - ke)) *
      (exp(-ke * tt) - exp(-p[["ka"]] * tt))
  }
  out
}

mix_arms <- list(
  `IV only`   = list(routes = c("iv", "iv", "iv")),
  `SC only`   = list(routes = c("extravascular", "extravascular",
                               "extravascular")),
  `IV then SC` = list(routes = c("iv", "extravascular", "extravascular"))
)
mix_dose_days <- c(0, 7, 14)
# Rich after the first and last dose, troughs in between. A real study of this
# shape samples the first and last cycle and takes pre-dose levels otherwise.
mix_rich <- c(0.042, 0.25, 0.5, 1, 2, 3, 5, 7)
mix_sparse <- c(0, 7)

mix_rows <- list()
mix_id <- 0L
for (arm in names(mix_arms)) {
  for (k in seq_len(30L)) {
    mix_id <- mix_id + 1L
    wt <- round(exp(stats::rnorm(1, log(72), 0.18)), 1)
    p <- mix_pop * exp(stats::rnorm(3, 0, mix_bsv))
    p[["cl"]] <- p[["cl"]] * (wt / 72)^0.75
    p[["v"]] <- p[["v"]] * (wt / 72)
    routes <- mix_arms[[arm]]$routes

    ntime <- sort(unique(c(
      mix_dose_days[1L] + mix_rich,
      mix_dose_days[2L] + mix_sparse,
      mix_dose_days[3L] + mix_rich)))
    # Recorded times drift off the plan the way real visits do.
    time <- ntime + stats::rnorm(length(ntime), 0, 0.02)
    time[ntime == 0] <- 0
    time <- pmax(time, 0)

    conc <- rep(0, length(time))
    for (d in seq_along(mix_dose_days)) {
      conc <- conc + mix_one_dose(time - mix_dose_days[d], MIX_DOSE,
                                  routes[d], p)
    }
    conc <- conc * (1 + stats::rnorm(length(conc), 0, MIX_PROP_ERR))
    cens <- as.integer(conc < MIX_LLOQ)
    conc[cens == 1L] <- MIX_LLOQ

    doses <- data.frame(
      ID = mix_id, TIME = mix_dose_days, NTIME = mix_dose_days,
      DV = NA_real_, AMT = MIX_DOSE, EVID = 1L,
      CMT = ifelse(routes == "iv", 2L, 1L),
      ADM = ifelse(routes == "iv", 1L, 2L),
      CENS = 0L, ARM = arm, WT = wt)
    obs <- data.frame(
      ID = mix_id, TIME = time, NTIME = ntime, DV = round(conc, 4),
      AMT = 0, EVID = 0L, CMT = 2L,
      ADM = ifelse(routes[pmax(findInterval(ntime, mix_dose_days), 1L)] == "iv",
                   1L, 2L),
      CENS = cens, ARM = arm, WT = wt)
    mix_rows[[length(mix_rows) + 1L]] <- rbind(doses, obs)
  }
}
mixroute_sim <- do.call(rbind, mix_rows)
mixroute_sim <- mixroute_sim[order(mixroute_sim$ID, mixroute_sim$TIME,
                                   -mixroute_sim$EVID), ]
rownames(mixroute_sim) <- NULL

# --------------------------------------------------------------------- onc_sim
#
# A simulated oncology study shaped like RECORD-1 (Stein et al., BMC Cancer
# 2012;12:311), which modelled the sum of longest tumour diameters against the
# everolimus dose each patient actually received.
#
# Three shapes the survey has no other example of, and one it says outright it
# has none of:
#
#  * a slow endpoint anchored on a per-patient baseline, months long, where
#    every other endpoint in the survey is a concentration or a fast PD signal;
#  * a dose that changes within a patient for reasons the patient's own data
#    explains -- crossover on progression, reduction for toxicity;
#  * `ADDL`/`II`. Daily dosing for a year is 365 rows a patient written out, and
#    the survey's own text says "none of these datasets uses steady state,
#    `ADDL`, or `II`", so compressed dose records have never been exercised
#    end to end.
#
# The tumour model is the paper's model 2, with its published parameters:
#
#   dy/dt = r_i - E_dose,i * y_i,  E_dose = E10 (10 mg), E5 (5 mg), 0 (none)
#   r_i    = r   * (y0_i / y0_hat)^theta1 + N(0, eta_r)
#   E10_i  = E10 * (y0_i / y0_hat)^theta2 + N(0, eta_E10)
#   E5_i   = E5  * (y0_i / y0_hat)^theta2 + N(0, eta_E5)
#
# Pharmacokinetics are NOT what drives the tumour here, exactly as in the paper:
# the effect is indexed by the dose in force, and the trough concentrations ride
# alongside as a second endpoint. They are a one-compartment oral model at the
# published everolimus disposition, so a PK fit of them means something.

ONC_Y0_HAT <- 14.4           # cm, the paper's median baseline SLD
# The paper reports the median baseline but not its spread, and the spread is
# what sets the placebo arm's percentage change: change is proportional to
# y0^(theta1 - 1), so a small tumour grows by a large percentage. 0.35 on the log
# scale reproduces the paper's reported one-year placebo change of +142.1 +/-
# 98.3% (simulated: +139.8 +/- 86.2) and puts the 5th to 95th centile of baseline
# SLD at 8.1 to 25.5 cm.
ONC_Y0_SD <- 0.35
onc_r      <- 46.0e-3        # cm/day
onc_E10    <- 3.9e-3         # 1/day
onc_E5     <- 2.3e-3         # 1/day
onc_theta1 <- 0.4
onc_theta2 <- -0.7
onc_eta    <- c(r = 3.5e-3, E10 = 0.3e-3, E5 = 0.2e-3)  # Table 2
ONC_SIGMA  <- 1.1            # cm, additive residual

# Everolimus disposition: CL/F 19.6 L/h and a 30 h terminal half-life, which put
# the 10 mg steady-state trough near the 15-20 ng/mL the label reports.
onc_pk  <- c(cl = 470, v = 850, ka = 5)   # L/day, L, 1/day
onc_pk_bsv <- c(cl = 0.40, v = 0.30, ka = 0.50)
ONC_PK_PROP_ERR <- 0.20
ONC_PK_LLOQ <- 1.0                        # ng/mL

# SLD between two scans under a dose held constant across the gap.
onc_sld_step <- function(y0, r, E, dt) {
  if (E <= 0) return(y0 + r * dt)
  ss <- r / E
  ss + (y0 - ss) * exp(-E * dt)
}

# Trough on day `t` from every daily dose given before it. Superposition, which
# is exact for a linear model and is what `ADDL`/`II` is shorthand for.
onc_trough <- function(t, dose_days, dose_mg, p) {
  ke <- p[["cl"]] / p[["v"]]
  live <- dose_days < t & dose_mg > 0
  if (!any(live)) return(0)
  el <- t - dose_days[live]
  sum(dose_mg[live] * 1000 * p[["ka"]] / (p[["v"]] * (p[["ka"]] - ke)) *
        (exp(-ke * el) - exp(-p[["ka"]] * el)))
}

ONC_SCANS <- c(0, seq(42, 378, by = 42))       # baseline then ~6-weekly
ONC_PK_DAYS <- c(14, 28, 56, 112, 168, 252)
ONC_N <- 200L

onc_rows <- list()
for (id in seq_len(ONC_N)) {
  y0 <- exp(stats::rnorm(1, log(ONC_Y0_HAT), ONC_Y0_SD))
  scale1 <- (y0 / ONC_Y0_HAT)^onc_theta1
  scale2 <- (y0 / ONC_Y0_HAT)^onc_theta2
  r_i   <- max(onc_r   * scale1 + stats::rnorm(1, 0, onc_eta[["r"]]),   1e-4)
  E10_i <- max(onc_E10 * scale2 + stats::rnorm(1, 0, onc_eta[["E10"]]), 1e-5)
  E5_i  <- max(onc_E5  * scale2 + stats::rnorm(1, 0, onc_eta[["E5"]]),  1e-5)
  p <- onc_pk * exp(stats::rnorm(3, 0, onc_pk_bsv))

  arm <- if (id %% 3L == 0L) "Placebo" else "Everolimus 10 mg"
  age <- round(stats::rnorm(1, 61, 10))
  sex <- sample(c("M", "F"), 1, prob = c(0.74, 0.26))

  # The dose in force, day by day, is the whole input to the tumour model. Built
  # as a daily vector first because that is what the model integrates against,
  # then compressed into `ADDL`/`II` blocks for the table.
  horizon <- max(ONC_SCANS)
  daily <- rep(if (identical(arm, "Placebo")) 0 else 10, horizon + 1L)

  if (!identical(arm, "Placebo")) {
    # A quarter of the everolimus arm reduce to 5 mg for toxicity, and a few of
    # those interrupt entirely for a couple of weeks first.
    if (stats::runif(1) < 0.25) {
      at <- sample(28:200, 1)
      if (stats::runif(1) < 0.4) {
        gap <- sample(7:21, 1)
        daily[(at + 1L):min(at + gap, horizon + 1L)] <- 0
        at <- at + gap
      }
      if (at < horizon) daily[(at + 1L):(horizon + 1L)] <- 5
    }
  }

  # Walk the scan schedule, integrating the tumour and testing for progression.
  sld <- numeric(0); sld_at <- numeric(0)
  y <- y0; nadir <- y0; crossed_at <- NA_real_; last_day <- horizon
  for (k in seq_along(ONC_SCANS)) {
    day <- ONC_SCANS[k]
    if (k > 1L) {
      span <- (ONC_SCANS[k - 1L] + 1L):day
      # One step per distinct dose level inside the gap, which is exact because
      # the dose is constant within each run.
      levels_in <- rle(daily[span + 1L])
      at <- ONC_SCANS[k - 1L]
      for (j in seq_along(levels_in$lengths)) {
        E <- switch(as.character(levels_in$values[j]),
                    `10` = E10_i, `5` = E5_i, 0)
        y <- onc_sld_step(y, r_i, E, levels_in$lengths[j])
        at <- at + levels_in$lengths[j]
      }
    }
    sld <- c(sld, y); sld_at <- c(sld_at, day)
    nadir <- min(nadir, y)
    progressed <- y >= 1.2 * nadir && (y - nadir) >= 0.5 && day > 0
    if (progressed) {
      if (identical(arm, "Placebo") && is.na(crossed_at) &&
          stats::runif(1) < 0.8) {
        # Crossover: the placebo patient starts everolimus at this scan and is
        # followed on. This is the shape the paper's Figure 1 top row shows.
        crossed_at <- day
        daily[(day + 1L):(horizon + 1L)] <- 10
        nadir <- y
      } else {
        last_day <- day          # off study after progression
        break
      }
    }
  }
  keep <- sld_at <= last_day
  sld <- sld[keep]; sld_at <- sld_at[keep]
  sld_obs <- round(pmax(sld + stats::rnorm(length(sld), 0, ONC_SIGMA), 0.1), 2)

  dose_days <- which(daily > 0) - 1L
  dose_amt <- daily[dose_days + 1L]
  pk_at <- ONC_PK_DAYS[ONC_PK_DAYS <= last_day]
  pk <- vapply(pk_at, onc_trough, numeric(1), dose_days = dose_days,
               dose_mg = dose_amt, p = p)
  pk <- pk * (1 + stats::rnorm(length(pk), 0, ONC_PK_PROP_ERR))
  pk_cens <- as.integer(pk < ONC_PK_LLOQ)
  pk[pk_cens == 1L] <- ONC_PK_LLOQ

  # Compress the daily vector into one row per constant-dose run.
  runs <- rle(daily[seq_len(last_day + 1L)])
  starts <- c(0L, cumsum(runs$lengths)[-length(runs$lengths)])
  dose_rows <- data.frame(
    ID = id, TIME = starts, NTIME = starts, DV = NA_real_,
    AMT = runs$values, EVID = 1L, CMT = 1L,
    ADDL = runs$lengths - 1L, II = 1,
    NAME = NA_character_, CENS = 0L)

  onc_rows[[length(onc_rows) + 1L]] <- rbind(
    dose_rows,
    data.frame(ID = id, TIME = sld_at, NTIME = sld_at, DV = sld_obs,
               AMT = 0, EVID = 0L, CMT = 3L, ADDL = 0L, II = 0,
               NAME = "SLD", CENS = 0L),
    data.frame(ID = id, TIME = pk_at, NTIME = pk_at, DV = round(pk, 3),
               AMT = 0, EVID = 0L, CMT = 2L, ADDL = 0L, II = 0,
               NAME = "Everolimus trough", CENS = pk_cens))
  n_new <- nrow(onc_rows[[length(onc_rows)]])
  onc_rows[[length(onc_rows)]]$ARM <- arm
  onc_rows[[length(onc_rows)]]$CROSSOVER <- !is.na(crossed_at)
  onc_rows[[length(onc_rows)]]$BSLD <- round(y0, 2)
  onc_rows[[length(onc_rows)]]$AGE <- age
  onc_rows[[length(onc_rows)]]$SEX <- sex
}
onc_sim <- do.call(rbind, onc_rows)
onc_sim <- onc_sim[order(onc_sim$ID, onc_sim$TIME, -onc_sim$EVID), ]
rownames(onc_sim) <- NULL

# ------------------------------------------------------------------- write out
dir.create("data", showWarnings = FALSE)
save(mixroute_sim, file = "data/mixroute_sim.rda", compress = "xz", version = 2)
save(onc_sim, file = "data/onc_sim.rda", compress = "xz", version = 2)
message("wrote data/mixroute_sim.rda: ", nrow(mixroute_sim), " rows, ",
        length(unique(mixroute_sim$ID)), " patients")
message("wrote data/onc_sim.rda: ", nrow(onc_sim), " rows, ",
        length(unique(onc_sim$ID)), " patients")
