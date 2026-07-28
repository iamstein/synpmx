# How to build things, and where the output goes

Internal note. Nothing here is cited from shipped documentation.

## Where the HTML is

**Everything rendered lands under `output/`, which is gitignored.** That is why
it never appears in `git status`, and why most editors hide it by default. It
is disposable: delete the whole directory any time and rebuild.

| What you ran | Where the HTML lands |
|---|---|
| `./build.sh vignettes` | `output/vignettes/` |
| `./build.sh articles` | `output/articles/` |
| `./build.sh docs` | both of the above |
| `dev_preview("name")` in R | `output/preview/` |
| `./build.sh` (plain check) | no HTML — see below |

So for a pkgdown article such as the AVATAR relationships example:

```
output/articles/example-avatar-PKPD-covariate-treatment-effect.html
```

Open it from the shell with `open output/articles/<name>.html`, or from R with
`dev_open("example-avatar")`, which finds the newest copy for you.

**A copy can exist in more than one of those directories** if you have used
several commands — `dev_preview()` writes to `output/preview/` while
`./build.sh articles` writes to `output/articles/`, and neither cleans the
other. `dev_open()` resolves this by opening whichever is newest. If they ever
disagree, trust `output/articles/`: that one was rendered against a clean
install rather than your working tree.

## Three ways to build, and what each one actually proves

They differ in what they would catch, which is the only reason to have three.

| | Speed | Renders against | Catches |
|---|---|---|---|
| `dev_preview("name")` | ~7 s | your working tree, via `load_all()` | broken code or prose in one document |
| `./build.sh articles` / `vignettes` / `docs` | ~3 min | a freshly built tarball installed into a throwaway library | the above, plus "it only worked because something was loaded in my session" |
| `./build.sh` | ~5 min | same clean library | `R CMD check`: tests, examples, `NAMESPACE`/`man` drift, shipped vignettes |
| pkgdown job on GitHub | ~10 min | a clean CI machine | everything above, plus the site actually assembling |

The fast path can succeed where a clean install fails. That is the trade, and
it is worth making while writing — just not the last thing you run.

## The gotcha that motivated all of this

**`R CMD check` never executes anything in `vignettes/articles/`.**
`.Rbuildignore` contains `^vignettes/articles$`, so articles are not part of
the package at all. They are executed only by:

- `./build.sh articles` (locally), and
- the pkgdown job on GitHub (on every push to `main`).

Before this was set up, a broken article passed every local check and failed in
CI. **Run `./build.sh articles` before pushing anything that touches an
article.** The plain `./build.sh` run now prints a reminder saying exactly this.

Shipped vignettes (`vignettes/*.Rmd`) do not have this problem — `R CMD check`
rebuilds them every run.

## `build.sh` reference

```bash
./build.sh              # roxygen, tarball, R CMD check. What CI runs.
./build.sh check        # same, explicit
./build.sh vignettes    # install, render vignettes/*.Rmd        -> output/vignettes/
./build.sh articles     # install, render vignettes/articles/*   -> output/articles/
./build.sh docs         # both render sets
./build.sh --keep-lib   # leave the throwaway library for inspection
./build.sh --help
```

Every mode builds a source tarball first and installs it into a fresh temporary
library placed ahead of the user library, so nothing is ever validated against
a stale installed `synpmx` or an already-loaded namespace. Suggested packages
(ggplot2, nlmixr2data, opendp, …) still resolve from the user library.

The render modes deliberately skip `R CMD check` — they are for reading output,
not for verification. Run `./build.sh` for that.

Logs for every step, including per-document render logs, are in `output/logs/`.
A failed render prints the path to its log.

## `dev.R` reference

For the R console, and the fast loop while writing:

```r
source("dev.R")

dev_preview("example-avatar")  # render one doc against the working tree + open
dev_open("avatar-alg")         # open the newest already-rendered copy
dev_list()                     # every doc, its kind, whether it is rendered

dev_articles()                 # ./build.sh articles
dev_vignettes()                # ./build.sh vignettes
dev_docs()                     # ./build.sh docs
dev_check()                    # ./build.sh
```

`dev_preview()` and `dev_open()` take any substring of the file name and error
if it is ambiguous. `dev.R` is excluded from the package by `.Rbuildignore`.

One subtlety worth knowing: every document calls `library(synpmx)` in its setup
chunk, which would load the *installed* package. `dev_preview()` therefore calls
`devtools::load_all()` first, so that `library()` call becomes a no-op and the
render genuinely exercises the working tree. Without that step the "fast
preview" would silently be testing whatever was last installed.

## The published site

`_pkgdown.yml` defines the site. The `pkgdown` GitHub Actions workflow builds it
on every push to `main` and deploys to the `gh-pages` branch, which serves
<https://iamstein.github.io/synpmx/>. A build takes roughly ten minutes, plus a
minute for Pages to publish.

There is no need — and no supported path — to build the full site locally:
`pkgdown::build_site()` would write into `docs/`, which is not how this repo
deploys. Use `./build.sh articles` to check that articles execute, and let CI
assemble the site.

Two things about that workflow worth remembering:

- **`rxode2` is deliberately not installed in the docs job.** It pulls in
  `qs` → `stringfish`, whose TBB symbols do not resolve against RcppParallel
  6.0.0, and that broke every site build from 2026-07-23. The one `rxode2`
  block in the documentation is a plain fence that never executed anyway.
- **Concurrency cancels superseded runs.** Pushing twice in quick succession
  shows the earlier run as `cancelled`, which is normal — the later run carries
  both commits.

To confirm something is live, fetch the page and look for text you added rather
than trusting the workflow's green tick:

```bash
curl -s https://iamstein.github.io/synpmx/articles/<name>.html | grep -c "some phrase"
```

Beware line wrapping when grepping rendered HTML: pandoc breaks lines mid
sentence, so a phrase that exists can still fail a line-based `grep`. Strip
tags and collapse whitespace before matching if a check comes back zero.
