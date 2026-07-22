private_fixture <- function(n = 8L) {
  pieces <- lapply(seq_len(n), function(subject) {
    time <- c(0, 0, 0.25, 1, 2, 4, 6, 8, 12, 12, 12.25, 13, 14, 16, 18, 20)
    endpoint <- c(
      "cp", "pd", "cp", "cp", "cp", "pd", "cp", "pd",
      "cp", "pd", "cp", "cp", "cp", "pd", "cp", "pd"
    )
    evid <- c(1L, 0L, rep(0L, 6L), 1L, rep(0L, 7L))
    cp_base <- c(NA, NA, 1, 8, 4, NA, 1.5, NA,
                 NA, NA, 1.2, 8.5, 4.2, NA, 1.6, NA)
    pd_base <- c(NA, 80, NA, NA, NA, 70, NA, 55,
                 NA, 40, NA, NA, NA, 50, NA, 65)
    dv <- ifelse(endpoint == "cp", cp_base, pd_base)
    dv[evid != 0] <- 0
    dv[evid == 0] <- dv[evid == 0] * (1 + (subject - 4.5) * 0.015)
    dose_times <- c(0, 12)
    occasion <- pmax(1L, findInterval(time, dose_times))
    tad <- pmax(0, time - dose_times[pmin(occasion, 2L)])
    data.frame(
      ID = as.integer(subject), TIME = time, NTIME = time, TAD = tad,
      OCC = occasion, DV = dv,
      AMT = ifelse(evid != 0, 100 + subject, 0), RATE = 0,
      EVID = evid,
      CMT = ifelse(evid != 0, 1L, ifelse(endpoint == "cp", 2L, 3L)),
      DVID = factor(endpoint, levels = c("cp", "pd")),
      MDV = ifelse(evid == 0, 0L, 1L), CENS = 0L, LIMIT = NA_real_,
      WT = rep(55 + 4 * subject, length(time)),
      AGE = rep(as.integer(25 + subject), length(time)),
      SEX = factor(rep(ifelse(subject %% 2L, "female", "male"), length(time)),
                   levels = c("female", "male")),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

private_roles <- function() {
  pmx_roles(
    id = "ID", time = "TIME", nominal_time = "NTIME", tad = "TAD",
    occasion = "OCC", dv = "DV", amt = "AMT", evid = "EVID",
    cmt = "CMT", dvid = "DVID", mdv = "MDV", rate = "RATE",
    cens = "CENS", limit = "LIMIT", covariates = c("WT", "AGE", "SEX")
  )
}

private_endpoints <- function() {
  list(
    cp = pmx_endpoint(
      dvid = "cp", alignment = "dose_relative", transform = "log",
      shape = "occasion", grid = c(0.25, 1, 2, 6), cmt = 2
    ),
    pd = pmx_endpoint(
      dvid = "pd", alignment = "study_time", transform = "identity",
      shape = "global", grid = c(0, 4, 8, 12, 16, 20), cmt = 3
    )
  )
}

private_bounds <- function() {
  pmx_bounds(
    time = c(0, 24), endpoints = list(cp = c(0, 20), pd = c(0, 120)),
    amt = c(0, 200), rate = c(-200, 200),
    covariates = list(WT = c(40, 120), AGE = c(18, 90)),
    limit = list(cp = c(0, 20), pd = c(0, 120))
  )
}

private_design <- function(data = private_fixture()) {
  pmx_public_design(
    schema = pmx_schema(data), dose_times = c(0, 12), dose_interval = 12,
    n_doses = 2, dose_amount = 100, dose_rate = 0,
    dose_evid = 1, dose_cmt = 1,
    endpoint_grids = lapply(private_endpoints(), `[[`, "grid"),
    endpoint_cmt = list(cp = 2, pd = 3), time_jitter_sd = 0.01
  )
}

private_limits <- function(max_rows = 40L, max_cells = 8L) {
  pmx_contribution_limits(
    max_rows = max_rows, max_doses = 4, max_occasions = 4,
    max_observations_per_endpoint = c(cp = 8, pd = 8),
    max_timing_cells = max_cells
  )
}

private_budget <- function() {
  pmx_budget_allocation(
    subject_count = 0.10, event = 0.15, timing = 0.15,
    covariates = 0.10, endpoints = 0.45, censoring = 0.05
  )
}

fit_public_fixture <- function(data = private_fixture(), ...) {
  suppressWarnings(fit_private_pmx(
    data = data, roles = private_roles(), endpoints = private_endpoints(),
    epsilon = 5, delta = 0, bounds = private_bounds(),
    public_design = private_design(data),
    contribution_limits = private_limits(),
    budget_allocation = private_budget(), backend = "public",
    public_source = TRUE, ...
  ))
}

# Minimal normalized endpoint plus trajectory field map for decode-level tests
# (SIM-020 / REV-001). Built through the public constructors so the working
# bounds, transform, and grid match what `fit_private_pmx()` would produce.
.threshold_fixture <- function(cells = 5L) {
  bounds <- pmx_bounds(c(0, 24), list(cp = c(0, 200)), amt = c(0, 500))
  roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                     evid = "EVID")
  design <- pmx_public_design(
    pmx_schema(data.frame(
      ID = 1L, TIME = 0, DV = 0, AMT = 0, EVID = 0L
    )),
    dose_evid = 1
  )
  limits <- pmx_contribution_limits(40, 8, 8, 20, cells)
  endpoints <- .normalize_endpoints(
    list(cp = pmx_endpoint(alignment = "dose_relative", transform = "log",
                           shape = "occasion")),
    roles, bounds, design, limits
  )
  list(
    endpoints = endpoints,
    map = list(cp = list(
      presence = paste("cp", "presence", seq_len(cells), sep = "__"),
      value = paste("cp", "value", seq_len(cells), sep = "__"),
      local_presence = character(), local_value = character()
    ))
  )
}
