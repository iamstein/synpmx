# Exercise the generation screen against known population models across seeds.
# This evaluates the screening thresholds, not optimizer convergence or recovery.
# Time is hours, dose mg, volumes L, clearances L/h, absorption 1/h;
# concentrations are mg/L. Population log-variance is 0.09 on each parameter.
# Run from the repository root: Rscript scripts/evaluate-model-acceptance.R
devtools::load_all(".", quiet = TRUE)

models <- list(
  `1cmt_oral` = c(cl = 4, v = 40, ka = 1),
  `2cmt_oral` = c(cl = 4, v = 20, q = 6, v2 = 40, ka = 1),
  `2cmt_infusion` = c(cl = 4, v = 20, q = 6, v2 = 40),
  `2cmt_mixed` = c(cl = 4, v = 20, q = 6, v2 = 40, ka = 1, f = 0.7)
)

simulate_design <- function(structural, parameters, seed) {
  .with_local_seed(seed, {
    effects <- .draw_random_effects(parameters$omega, 30L)
    do.call(rbind, lapply(seq_len(30L), function(id) {
      times <- c(0.25, 0.5, 1, 2, 4, 8, 12, 24, 24.5, 28, 36, 48)
      dose_times <- c(0, 24)
      dose <- rep(if (id <= 15) 100 else 300, 2)
      routes <- if (grepl("mixed$", structural)) c("iv", "extravascular") else NULL
      duration <- if (grepl("infusion$", structural)) c(2, 2) else
        if (grepl("mixed$", structural)) c(2, 0) else c(0, 0)
      p <- .subject_parameters(parameters$fixed, list(), effects[id, ], list())
      value <- .pk_profile(list(pk = structural), times, dose, dose_times,
                            p, duration, routes)
      observations <- data.frame(ID = id, TIME = times,
        DV = .add_residual_error(value, parameters$residual, floor = 0),
        AMT = 0, EVID = 0L, RATE = 0, CMT = 2L)
      doses <- data.frame(ID = id, TIME = dose_times, DV = NA_real_,
        AMT = dose, EVID = 1L, RATE = ifelse(duration > 0, dose / duration, 0),
        CMT = if (is.null(routes)) 1L else ifelse(routes == "iv", 2L, 1L))
      rbind(doses, observations)
    }))
  })
}

rows <- list()
for (structural in names(models)) {
  fixed <- models[[structural]]
  omega <- diag(0.09, length(fixed))
  dimnames(omega) <- list(names(fixed), names(fixed))
  parameters <- list(fixed = fixed, omega = omega,
                       residual = list(kind = "proportional", cv = 0.2))
  for (seed in 1:20) {
    source <- simulate_design(structural, parameters, seed)
    for (scenario in c("known_population", "gross_level_error", "explosive_covariance")) {
      given <- parameters
      if (scenario == "gross_level_error") given$fixed <- given$fixed *
        ifelse(names(fixed) %in% c("cl", "v", "q", "v2"), 1e-4, 1)
      if (scenario == "explosive_covariance") given$omega <- given$omega * 1e8
      check <- .model_check_generation(given, structural, source)
      rows[[length(rows) + 1L]] <- data.frame(
        model = structural, seed = seed, scenario = scenario,
        accepted = is.null(check$failure), note = check$note,
        failure = check$failure %||% "", stringsAsFactors = FALSE)
    }
  }
}
result <- do.call(rbind, rows)
dir.create("output/model-acceptance", recursive = TRUE, showWarnings = FALSE)
utils::write.csv(result, "output/model-acceptance/generation-screen.csv", row.names = FALSE)
print(aggregate(accepted ~ model + scenario, result, sum), row.names = FALSE)
stopifnot(all(result$accepted[result$scenario == "known_population"]),
          !any(result$accepted[result$scenario != "known_population"]))
