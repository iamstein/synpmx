# End-to-end demonstrations for the three nlmixr2data integration datasets.
# Install pmxSynthData first with: R CMD INSTALL .

if (!requireNamespace("pmxSynthData", quietly = TRUE)) {
  stop("Install pmxSynthData before running this script: R CMD INSTALL .")
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
  mock <- pmxSynthData::mock_pmx(source, roles, seed = seed)
  validation <- pmxSynthData::validate_pmx(mock, roles)
  comparison <- pmxSynthData::compare_pmx(source, mock, roles)
  stopifnot(validation$valid)

  # Identify the concrete source object and generated data directly on every
  # plot. The existing legends/facets continue to distinguish source and mock.
  if (length(comparison$plots)) {
    dataset_caption <- paste0(
      "Source: nlmixr2data::", name,
      " | Mock: pmxSynthData::mock_pmx()"
    )
    comparison$plots <- lapply(comparison$plots, function(plot) {
      plot + ggplot2::labs(caption = dataset_caption)
    })
  }

  print(comparison)
  if (length(comparison$plots) && interactive()) {
    print(comparison$plots$faceted)
  }
  invisible(list(source = source, mock = mock, comparison = comparison))
}

theophylline <- run_demo(
  "theo_md",
  pmxSynthData::pmx_roles(
    id = "ID", time = "TIME", dv = "DV", amt = "AMT",
    evid = "EVID", cmt = "CMT", covariates = "WT"
  ),
  seed = 101
)

warfarin_result <- run_demo(
  "warfarin",
  pmxSynthData::pmx_roles(
    id = "id", time = "time", dv = "dv", amt = "amt",
    evid = "evid", dvid = "dvid",
    covariates = c("wt", "age", "sex")
  ),
  seed = 202
)

wbc <- run_demo(
  "wbcSim",
  pmxSynthData::pmx_roles(
    id = "ID", time = "TIME", dv = "DV", amt = "AMT",
    evid = "EVID", cmt = "CMT", rate = "RATE",
    covariates = c("V2I", "V1I", "CLI")
  ),
  seed = 303
)
