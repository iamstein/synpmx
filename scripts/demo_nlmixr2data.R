# End-to-end demonstrations for the three nlmixr2data integration datasets.
# Install pmxmock first with: R CMD INSTALL .

if (!requireNamespace("pmxmock", quietly = TRUE)) {
  stop("Install pmxmock before running this script: R CMD INSTALL .")
}
if (!requireNamespace("nlmixr2data", quietly = TRUE)) {
  stop("Install nlmixr2data to run the integration demonstrations.")
}

load_dataset <- function(name) {
  environment <- new.env(parent = emptyenv())
  utils::data(list = name, package = "nlmixr2data", envir = environment)
  get(name, envir = environment, inherits = FALSE)
}

run_demo <- function(name, roles, seed) {
  source <- load_dataset(name)
  mock <- pmxmock::mock_pmx(source, roles, seed = seed)
  validation <- pmxmock::validate_pmx(mock, roles)
  comparison <- pmxmock::compare_pmx(source, mock, roles)
  stopifnot(validation$valid)
  print(comparison)
  if (length(comparison$plots) && interactive()) {
    print(comparison$plots$faceted)
  }
  invisible(list(source = source, mock = mock, comparison = comparison))
}

theophylline <- run_demo(
  "theo_md",
  pmxmock::pmx_roles(
    id = "ID", time = "TIME", dv = "DV", amt = "AMT",
    evid = "EVID", cmt = "CMT", covariates = "WT"
  ),
  seed = 101
)

warfarin_result <- run_demo(
  "warfarin",
  pmxmock::pmx_roles(
    id = "id", time = "time", dv = "dv", amt = "amt",
    evid = "evid", dvid = "dvid",
    covariates = c("wt", "age", "sex")
  ),
  seed = 202
)

wbc <- run_demo(
  "wbcSim",
  pmxmock::pmx_roles(
    id = "ID", time = "TIME", dv = "DV", amt = "AMT",
    evid = "EVID", cmt = "CMT", rate = "RATE",
    covariates = c("V2I", "V1I", "CLI")
  ),
  seed = 303
)
