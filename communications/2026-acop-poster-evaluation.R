# Rebuild the public-data evaluation table for the poster.
# Run from the repository root:
#   Rscript communications/2026-acop-poster-evaluation.R
# Reads the dataset preparation and stored-fit generation chunks in the model
# survey, preserving its roles and seeds. Estimation calls are never executed.
# Writes aggregate CSV, a Markdown table and run information beside this script.
# Paste the Markdown table into the outline; internal rows are maintained there.

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) != 1L) stop("Run with Rscript.")
script <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
root <- dirname(dirname(script))
out <- file.path(root, "communications", "2026-acop-poster-evaluation")
required <- c("devtools", "xgxr", "nlmixr2data")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))
devtools::load_all(root, quiet = TRUE)

article_file <- file.path(root, "vignettes", "articles", "pmxmodel-public-data-examples.Rmd")
article <- readLines(article_file)
headers <- grep("^```\\{r ", article)
labels <- sub("^```\\{r ([^,}]+).*", "\\1", article[headers])
chunks <- setNames(lapply(headers, function(i) {
  end <- which(article == "```" & seq_along(article) > i)[1L]
  article[seq.int(i + 1L, end - 1L)]
}), labels)
selected <- c("helpers", "case1", "mad", "warfarin", "wbc", "mavo",
              "theo-grid", "theo", "nimo", "nimo-run", "pheno",
              "pheno-grid", "pheno-run", "mixroute", "onc")
if (!all(selected %in% names(chunks))) stop("Survey chunk names changed; review this script.")

survey <- new.env(parent = globalenv())
warnings <- character()
withCallingHandlers({
  for (label in selected) {
    message("Evaluating survey chunk: ", label)
    for (expr in parse(text = chunks[[label]])) {
      # The nimo and pheno preparation chunks end in deliberate fitting
      # refusals. Skip those expressions; only stored fits drive this survey.
      if ("synpmx_model_estimate" %in% all.names(expr)) next
      invisible(eval(expr, envir = survey))
    }
    if (label == "helpers") {
      survey$stored_fit <- function(file) readRDS(file.path(root, "inst", "extdata", file))
    }
  }
}, warning = function(w) warnings <<- c(warnings, conditionMessage(w)))
runs <- survey$runs
datasets <- c("case1_pkpd", "mad", "warfarin", "wbcSim", "mavoglurant",
              "theo_md", "nimoData", "pheno_sd", "mixroute_sim", "onc_sim")
if (!identical(names(runs), datasets)) stop("Survey dataset list changed; review the table.")

# Editorial findings describe the stored-fit configuration in the survey.
# Review these against the outputs when fits or generation behaviour change.
design <- c(
  "Multiple arms; censored concentrations",
  "Multiple doses; mixed endpoint types",
  "Single oral dose; delayed response",
  "Infusions; white-cell nadir and recovery",
  "Repeated occasions; resetting clock",
  "Repeated oral doses; small cohort",
  "Weekly infusions; small cohort",
  "Neonatal care; sparse, irregular sampling",
  "Intravenous and subcutaneous dosing",
  "Trough sampling; tumour response; dose changes"
)
findings <- c(
  "Censoring represented; pooled PD loses arm differences.",
  "Endpoint types retained; early PK peak underrepresented.",
  "PD decline represented; recovery absent.",
  "Response fitted as PK; nadir and recovery missed.",
  "Second occasions lost; fewer observations generated.",
  "Repeated dosing retained; fit rests on a small cohort.",
  "Infusion schedule retained; concentration spread inflated.",
  "Nominal grid required; fewer doses generated.",
  "Both routes retained; bioavailability close to simulation truth.",
  "Pooled tumour curve loses arm differences; fewer doses generated."
)
table <- data.frame(
  Dataset = datasets,
  Patients = vapply(runs, function(x) length(unique(x$source[[x$roles$id]])), integer(1)),
  `Design feature` = design, `Main finding` = findings,
  row.names = NULL, check.names = FALSE
)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(table, file.path(out, "public-data-table.csv"), row.names = FALSE)
rows <- vapply(seq_len(nrow(table)), function(i) {
  paste0("| `", table$Dataset[i], "` | ", table$Patients[i], " | ",
         table$`Design feature`[i], " | ", table$`Main finding`[i], " |")
}, character(1))
writeLines(c("| Dataset | Patients | Design feature | Main finding |",
             "|---|---:|---|---|", rows), file.path(out, "public-data-table.md"))
cards <- do.call(rbind, lapply(runs, function(x) {
  cbind(Dataset = x$label, as.data.frame(x$card))
}))
utils::write.csv(cards, file.path(out, "supporting-scorecards.csv"), row.names = FALSE)
fit_files <- list.files(file.path(root, "inst", "extdata"),
                       pattern = "model-fit[.]rds$", full.names = TRUE)
writeLines(c(
  paste("Generated:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  "Public data only. Generation and evaluation rerun; population fits not re-estimated.",
  "Patient counts computed from source subjects; qualitative findings reviewed separately.",
  "Seeds and source preparation follow the named survey chunks.",
  "Warnings:", warnings,
  "Survey, script, stored-fit and working-tree R source MD5s:",
  capture.output(print(tools::md5sum(c(article_file, script, fit_files,
    list.files(file.path(root, "R"), full.names = TRUE))))),
  capture.output(sessionInfo())
), file.path(out, "run-info.txt"))
message("Wrote public-data table to ", out)
