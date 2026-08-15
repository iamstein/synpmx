# Developer conveniences. In the R console:
#
#     source("dev.R")
#
# Then use the functions below. Nothing here is part of the package --- this
# file is excluded from the build by .Rbuildignore.
#
# Two speeds, and the difference matters:
#
#   dev_preview("avatar-alg")   fast. Renders ONE document against the working
#                               tree via devtools::load_all(). Use while writing.
#
#   dev_articles()              slow, honest. Builds a tarball, installs it into
#                               a throwaway library, and renders against THAT.
#                               Use before pushing. Only this catches "it works
#                               because my session has something loaded".

.dev_root <- normalizePath(dirname(sys.frame(1)$ofile %||% "."), mustWork = FALSE)
if (!nzchar(.dev_root) || is.na(.dev_root)) .dev_root <- getwd()

`%||%` <- function(x, y) if (is.null(x)) y else x

.dev_sh <- function(...) {
  args <- paste(c(...), collapse = " ")
  cat("$ ./build.sh", args, "\n\n")
  status <- system2(file.path(.dev_root, "build.sh"), args)
  invisible(status)
}

#' Full verification: roxygen, tarball, R CMD check. What CI runs.
dev_check <- function() .dev_sh("check")

#' Render the three shipped vignettes to output/vignettes/, clean library.
dev_vignettes <- function() .dev_sh("vignettes")

#' Render every pkgdown article to output/articles/, clean library.
#' R CMD check never executes articles, so this is the only local check on them.
dev_articles <- function() .dev_sh("articles")

#' Both of the above.
dev_docs <- function() .dev_sh("docs")

#' Every document that can be rendered, with where its HTML lands.
dev_list <- function() {
  v <- Sys.glob(file.path(.dev_root, "vignettes", "*.Rmd"))
  a <- Sys.glob(file.path(.dev_root, "vignettes", "articles", "*.Rmd"))
  out <- data.frame(
    name = c(tools::file_path_sans_ext(basename(v)),
             tools::file_path_sans_ext(basename(a))),
    kind = c(rep("vignette", length(v)), rep("article", length(a))),
    stringsAsFactors = FALSE
  )
  out$rendered <- file.exists(file.path(
    .dev_root, "output",
    ifelse(out$kind == "vignette", "vignettes", "articles"),
    paste0(out$name, ".html")
  ))
  # dev_preview() writes through pkgdown into docs/articles/, so a document can
  # be rendered without output/ knowing about it.
  out$previewed <- file.exists(file.path(
    .dev_root, "docs", "articles", paste0(out$name, ".html")
  ))
  out[order(out$kind, out$name), ]
}

.dev_find <- function(pattern) {
  files <- c(Sys.glob(file.path(.dev_root, "vignettes", "*.Rmd")),
             Sys.glob(file.path(.dev_root, "vignettes", "articles", "*.Rmd")))
  hit <- files[grepl(pattern, basename(files), ignore.case = TRUE, fixed = FALSE)]
  if (length(hit) == 0L) {
    stop("no document matches ", sQuote(pattern), ". Try dev_list().", call. = FALSE)
  }
  if (length(hit) > 1L) {
    stop("several documents match ", sQuote(pattern), ":\n  ",
         paste(basename(hit), collapse = "\n  "), call. = FALSE)
  }
  hit
}

#' Fast single-document preview against the working tree, then open it.
#'
#' Renders with devtools::load_all(), so edits to R/ take effect without an
#' install. That speed is the whole point, and also the catch: it can succeed
#' on something a clean install would fail. Run dev_articles() before pushing.
#'
#' Rendered through pkgdown, so the page looks the way readers will see it:
#' contents in the right-hand sidebar, site CSS, real navbar. The .Rmd's own
#' `output:` format puts its table of contents inline at the top instead, which
#' is what makes a long page hard to judge.
#'
#' `new_process = FALSE` is load-bearing. pkgdown defaults to rendering in a
#' fresh R session, where `library(synpmx)` would pick up whatever is installed
#' and the working tree would go untested. In-process, the load_all() namespace
#' below is already attached and that call is a no-op.
#'
#' @param pattern Any substring of the file name, e.g. "pyraz" or "avatar-alg".
#' @param open Open the rendered HTML when done.
#' @param site Render through the site template. FALSE falls back to the .Rmd's
#'   own format in output/preview/, which is a little faster and does not need
#'   the site assets.
dev_preview <- function(pattern, open = TRUE, site = TRUE) {
  rmd <- .dev_find(pattern)
  # Attach the working tree BEFORE rendering. Each document calls
  # library(synpmx) in its setup chunk, which would otherwise pick up whatever
  # is installed; with the namespace already loaded from source, that call is a
  # no-op and the render genuinely exercises R/.
  devtools::load_all(.dev_root, quiet = TRUE)
  name <- tools::file_path_sans_ext(basename(rmd))

  if (site) {
    # The site assets are copied once and then reused; docs/ is gitignored.
    if (!file.exists(file.path(.dev_root, "docs", "pkgdown.yml"))) {
      message("initialising the site (once)")
      pkgdown::init_site(.dev_root)
    }
    # pkgdown names an article by its path under vignettes/, without extension.
    article <- sub("\\.Rmd$", "",
                   sub(file.path(.dev_root, "vignettes", ""), "", rmd,
                       fixed = TRUE))
    message("rendering ", basename(rmd), " through pkgdown (working tree)")
    elapsed <- system.time(
      pkgdown::build_article(article, pkg = .dev_root, new_process = FALSE,
                             quiet = TRUE)
    )[["elapsed"]]
    html <- file.path(.dev_root, "docs", "articles", paste0(name, ".html"))
  } else {
    dest <- file.path(.dev_root, "output", "preview")
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    message("rendering ", basename(rmd), " (working tree via load_all)")
    elapsed <- system.time(
      rmarkdown::render(rmd, output_dir = dest, quiet = TRUE,
                        envir = new.env(parent = globalenv()))
    )[["elapsed"]]
    html <- file.path(dest, paste0(name, ".html"))
  }

  message(sprintf("%.1fs -> %s", elapsed, html))
  if (open) utils::browseURL(html)
  invisible(html)
}

#' Open an already-rendered document without re-rendering it.
dev_open <- function(pattern) {
  name <- tools::file_path_sans_ext(basename(.dev_find(pattern)))
  candidates <- c(
    file.path(.dev_root, "docs", "articles", paste0(name, ".html")),
    file.path(.dev_root, "output", c("preview", "articles", "vignettes"),
              paste0(name, ".html"))
  )
  found <- candidates[file.exists(candidates)]
  if (length(found) == 0L) {
    stop("no rendered HTML for ", sQuote(name),
         ". Render it first with dev_preview(\"", name, "\").", call. = FALSE)
  }
  newest <- found[which.max(file.mtime(found))]
  message("opening ", newest)
  utils::browseURL(newest)
  invisible(newest)
}

message(
  "dev helpers loaded:\n",
  "  dev_preview(\"name\")  render one doc as the site shows it + open\n",
  "  dev_open(\"name\")     open what is already rendered\n",
  "  dev_list()           what exists and what is rendered\n",
  "  dev_articles()       clean-library render of all articles (before pushing)\n",
  "  dev_vignettes()      clean-library render of all vignettes\n",
  "  dev_docs()           both\n",
  "  dev_check()          full R CMD check (what CI runs)"
)
