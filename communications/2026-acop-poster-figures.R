# Generate the public mad example figures for the ACOP poster.
# From the repository root:
#   Rscript communications/2026-acop-poster-figures.R
# Optional output directory and generation seed:
#   Rscript communications/2026-acop-poster-figures.R /tmp/acop-figures 909
# Uses the working-tree package and the stored fit; never estimates a new fit.
# Inputs are fixed to public xgxr data. No internal-study input is accepted.
# Outputs: PNG (200 dpi), vector PDF, aggregate evaluation CSV and run provenance.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2L) stop("Expected at most an output directory and a seed.")
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) != 1L) stop("Run this script with Rscript.")
script <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
root <- dirname(dirname(script))
out <- if (length(args)) args[[1L]] else
  file.path(root, "communications", "2026-acop-poster-figures")
seed <- if (length(args) == 2L) suppressWarnings(as.integer(args[[2L]])) else 909L
if (is.na(seed) || seed < 0L) stop("Seed must be a non-negative integer.")
packages <- c("devtools", "xgxr", "ggplot2", "patchwork")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))
devtools::load_all(root, quiet = TRUE)
library(ggplot2)
library(patchwork)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out, mustWork = TRUE)

source_data <- as.data.frame(get(utils::data(list = "mad", package = "xgxr")))
roles <- pmx_roles(
  id = "ID", time = "TIME", dv = "LIDV", amt = "AMT", evid = "EVID",
  cmt = "CMT", dvid = "NAME", mdv = "MDV", nominal_time = "NOMTIME",
  strata = c("TRTACT", "DOSE"), covariates = c("WEIGHTB", "SEX")
)
fit_file <- file.path(root, "inst", "extdata", "mad-model-fit.rds")
fit <- readRDS(fit_file)
synthetic <- synpmx_model_generate(fit, seed = seed)
card <- synpmx_scorecard(source_data, synthetic, roles)
utils::write.csv(as.data.frame(card), file.path(out, "mad-scorecard.csv"),
                 row.names = FALSE)
summaries <- compare_pmx_distributions(source_data, synthetic, roles, output = "tables")
for (name in names(summaries)) {
  if (!is.null(summaries[[name]])) {
    utils::write.csv(summaries[[name]], file.path(out, paste0("mad-", name, ".csv")),
                     row.names = FALSE)
  }
}

arm_table <- unique(source_data[c("TRTACT", "DOSE")])
arms <- as.character(arm_table$TRTACT[order(arm_table$DOSE)])
colours <- c(Source = "#1B6CA8", Synthetic = "#D95F02")
poster_theme <- theme_minimal(base_size = 24) + theme(
  panel.grid.minor = element_blank(),
  strip.text = element_text(size = 22, face = "bold"),
  strip.clip = "off",
  axis.text = element_text(size = 20),
  axis.title = element_text(size = 24),
  plot.title = element_text(size = 27, face = "bold"),
  plot.subtitle = element_text(size = 20),
  plot.margin = margin(15, 65, 15, 15),
  legend.position = "none"
)

# Exactly the same plotting transformation reads synthetic and source tables.
# Select the first PK interval on the declared nominal clock, then plot actual
# recorded times. No source-derived cutoffs or per-source plotting patches.
profile_rows <- function(data, label, endpoint, first_interval = FALSE) {
  keep <- data$EVID == 0 & data$MDV == 0 & data$NAME == endpoint &
    is.finite(data$TIME) & is.finite(data$LIDV)
  if (first_interval) keep <- keep & data$NOMTIME >= 0 & data$NOMTIME < 24
  keep[is.na(keep)] <- FALSE
  rows <- data[keep, c("ID", "TIME", "LIDV", "TRTACT")]
  rows$Dataset <- factor(label, levels = c("Source", "Synthetic"))
  rows$TRTACT <- factor(rows$TRTACT, levels = arms)
  rows
}

profile_plot <- function(endpoint, title, y_label, first_interval = FALSE,
                         log_y = FALSE) {
  # Exercise the helper on the synthetic table first, then on the public source.
  generated <- profile_rows(synthetic, "Synthetic", endpoint, first_interval)
  original <- profile_rows(source_data, "Source", endpoint, first_interval)
  stopifnot(nrow(generated) > 0L, nrow(original) > 0L)
  rows <- rbind(original, generated)
  if (log_y && any(rows$LIDV <= 0)) stop("Non-positive PK values need review before log plotting.")
  plot <- ggplot(rows, aes(TIME, LIDV, group = interaction(Dataset, ID),
                           colour = Dataset)) +
    geom_line(alpha = 0.45, linewidth = 0.55) +
    geom_point(alpha = 0.55, size = 0.8) +
    facet_grid(Dataset ~ TRTACT, drop = FALSE) +
    scale_colour_manual(values = colours) +
    labs(title = title, x = "Time (hours)", y = y_label) + poster_theme
  if (first_interval) plot <- plot + scale_x_continuous(breaks = c(0, 12, 24))
  else plot <- plot + scale_x_continuous(breaks = c(0, 96, 192))
  if (log_y) plot <- plot + scale_y_log10()
  plot
}

save_figure <- function(plot, stem, width, height) {
  png_file <- file.path(out, paste0(stem, ".png"))
  pdf_file <- file.path(out, paste0(stem, ".pdf"))
  ggsave(png_file, plot, width = width, height = height, units = "in",
         dpi = 200, bg = "white")
  ggsave(pdf_file, plot, device = grDevices::pdf, width = width, height = height,
         units = "in", bg = "white", useDingbats = FALSE)
  stopifnot(all(file.info(c(png_file, pdf_file))$size > 0))
  message("Wrote ", png_file, " and PDF")
}

save_figure(profile_plot("PK Concentration", "PK profiles: first dosing interval",
                          "Concentration (ng/mL)", TRUE, TRUE),
            "mad-pk-profiles", 18, 7)
save_figure(profile_plot("PD - Continuous", "Continuous PD profiles",
                          "Response (IU/L)"),
            "mad-pd-profiles", 18, 7)

# Use the public distribution API for both comparisons. The compact figure
# selects PK, continuous PD and baseline weight; the full figure retains every
# endpoint and baseline covariate. Density curves use relative peak height.
compact_roles <- pmx_roles(
  id = "ID", time = "TIME", dv = "LIDV", amt = "AMT", evid = "EVID",
  cmt = "CMT", dvid = "NAME", mdv = "MDV", nominal_time = "NOMTIME",
  strata = c("TRTACT", "DOSE"), covariates = "WEIGHTB"
)
compact_data <- function(data) {
  data[data$EVID != 0 | data$NAME %in% c("PK Concentration", "PD - Continuous"), ]
}
distribution_theme <- theme(
  text = element_text(size = 22),
  plot.title = element_text(size = 24, face = "bold"),
  axis.title = element_text(size = 22),
  axis.title.y = element_text(size = 22),
  axis.text = element_text(size = 18),
  legend.text = element_text(size = 22),
  plot.margin = margin(12, 12, 12, 12)
)
compact <- compare_pmx_distributions(compact_data(source_data),
                                      compact_data(synthetic), compact_roles)
save_figure(compact & distribution_theme, "mad-distributions", 18, 4.5)
full <- compare_pmx_distributions(source_data, synthetic, roles)
save_figure(full & distribution_theme, "mad-distributions-full", 18, 12)

provenance <- c(
  "Input: publicly available xgxr::mad (not an internal clinical dataset).",
  paste("Generated:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("Generation seed:", seed),
  "Fit: stored mad-model-fit.rds; no model re-estimation in this run.",
  paste("Fit MD5:", unname(tools::md5sum(fit_file))),
  paste("Script MD5:", unname(tools::md5sum(script))),
  "Working-tree R source MD5s:",
  capture.output(print(tools::md5sum(list.files(file.path(root, "R"), full.names = TRUE)))),
  "Profile helper executed on synthetic and source data without modification.",
  "This checks plotting compatibility, not a full AI-assisted analysis workflow.",
  "Profile axes: actual time in hours; PK ng/mL; continuous PD IU/L (xgxr metadata).",
  "Distribution panels pool arms and visits; they do not establish arm-specific agreement.",
  "PK has no placebo observations: the empty placebo panels are intentional.",
  capture.output(print(card)), capture.output(sessionInfo())
)
writeLines(provenance, file.path(out, "run-info.txt"))
message("Figures and supporting evaluation saved to ", out)
