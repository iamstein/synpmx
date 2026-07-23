#!/usr/bin/env Rscript
#
# Remove agent-facing repository conventions from the built pkgdown site.
#
# pkgdown's package_mds() renders every top-level .md into the site. As of
# pkgdown 2.2.1 its exclusion list is hardcoded (README, LICENCE, NEWS, 404,
# and three template names) with no configuration hook, and .Rbuildignore is
# not consulted. AGENTS.md and CLAUDE.md are instructions for coding agents,
# not user documentation, so they are pruned here instead.
#
# Run after pkgdown::build_site(). Safe to run repeatedly.

exclude <- c("AGENTS", "CLAUDE")
docs <- "docs"

if (!dir.exists(docs)) {
  stop("No 'docs' directory; build the site before pruning it.", call. = FALSE)
}

pages <- paste0(exclude, ".html")

# 1. The rendered pages and their copied markdown sources.
targets <- file.path(docs, c(pages, paste0(exclude, ".md")))
removed <- file.remove(targets[file.exists(targets)])
message("pruned ", length(removed), " file(s) from ", docs, "/")

# 2. Search index entries, or the pages stay reachable from the search box.
search_json <- file.path(docs, "search.json")
if (file.exists(search_json)) {
  index <- jsonlite::fromJSON(search_json, simplifyVector = FALSE)
  keep <- !vapply(index, function(entry) {
    path <- entry$path
    if (is.null(path) || length(path) != 1L || is.na(path)) return(FALSE)
    any(vapply(pages, function(p) grepl(p, path, fixed = TRUE), logical(1)))
  }, logical(1))
  jsonlite::write_json(index[keep], search_json, auto_unbox = TRUE)
  message("pruned ", sum(!keep), " search index entr(ies)")
}

# 3. Sitemap entries, so the pages are not advertised to crawlers.
sitemap <- file.path(docs, "sitemap.xml")
if (file.exists(sitemap)) {
  lines <- readLines(sitemap, warn = FALSE)
  drop <- Reduce(`|`, lapply(pages, function(p) grepl(p, lines, fixed = TRUE)))
  writeLines(lines[!drop], sitemap)
  message("pruned ", sum(drop), " sitemap entr(ies)")
}
