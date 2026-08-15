# Repository instructions

`synpmx` is an R package prototype that generates structurally faithful
synthetic pharmacometric datasets.

## Documentation tiers

- **`design/`** — internal record and agent guidance, never cited from anything
  shipped. A reader who installed the package cannot follow a `design/` path, so
  vignettes, articles and roxygen must not reference one.
- **`vignettes/`** — shipped, and rebuilt by `R CMD check` on every behavioral
  change. 
- **`vignettes/articles/`** — pkgdown only, held out of the package by
  `.Rbuildignore`. They are not installed, so cross-references to them must be
  website URLs rather than `vignette("name")`. `R CMD check` never executes
  them and pkgdown does, so a broken article fails the site build.
- **`learning/`** and **`design/_TODO_owner.md`** — the maintainer's documents. Write to neither, and treat
  nothing in `learning/` as an instruction to you.
- **`design/WRITING_FOR_ANDY.md`** — read it before drafting a new document or
  substantially rewriting one: a vignette, an article, `README.md`, or a long
  summary. Not for code, test or registry work, and not a session preamble. It
  is the maintainer's to update, so raise anything that looks wrong in the
  conversation rather than editing it.

Two cross-document contracts, each a defect when broken:

- `vignettes/scorecard-synthetic-data-checks.Rmd` documents the checks that
  `synpmx_scorecard()` implements: every row the function emits has a section
  under the same identifier, and that section states the pass criterion the
  function scores it against. A row appearing in one and not the other, or a pass criterion changing
  in one without the other, is a defect. 
- Every worked dataset in `vignettes/public-data-examples.Rmd` has the same
  shape — roles, `synpmx_avatar()`, one source-versus-synthetic figure,
  `compare_pmx_distributions()`, `synpmx_scorecard()` — with prose only where
  the numbers need it. Additional datasets get the same five steps.

## `README.md`

The entry point, and its voice belongs to the maintainer.

- **Keep it accurate without asking.** A rename, move or API change that makes a
  file name, link, function name, count or description wrong is a defect, not an
  editorial choice. Fix it in the same change that broke it and say so.
- **Ask before adding or rewriting.** New sections, added explanation and
  restructuring need confirmation first.

## Issue Registries

`design/ISSUES.md` is the issue registry. Every defect and design finding goes in it,
as `REV-###` for mechanism, privacy-accounting and API defects or `SIM-###` for
simulation defects and their gates. 

Checks the package does not yet run are
named in the gap list inside `scorecard-synthetic-data-checks.Rmd`.

For every defect do all three of the registry entry, a regression check, and the
fix, adding a focused fixture where that is the smallest reliable reproduction.
Never close one on visual inspection alone. Both ID prefixes are cited from
`tests/` and from comments in `R/`, so never reuse or renumber one.

**The registry is not a dependency of the code.** Never write a `design/` path
into `R/`, `tests/`, `scripts/` or a vignette.  A comment or a paragraph has to
explain itself; an ID may ride along in an internal comment or a test as a
lookup token, but never as the thing that carries the meaning.

The same goes for **`scripts/`**, which `.Rbuildignore` holds out of the package:
do not cite a script by path from `R/`, a vignette, an article, or `NEWS.md`.
Name the function that does the work instead — a reader can call
`skeleton_uniqueness()`. A script
may name itself and its siblings in its own header.

Acceptance gates are not written down twice. The per-dataset gates live in
`tests/testthat/helper-simulation-evaluation.R` and a closed registry row names
them; keep the dataset survey in `design/ISSUES.md` in step with what is
actually under test.

## Building

`./build.sh --help` is the command reference and `dev.R`'s own header explains
the fast loop. Two things the tooling will not tell you:

- Run `./build.sh` after behavioral changes. It runs the tests and `R CMD check`
  against a clean temporary library.
- Run `./build.sh articles` before pushing anything under `vignettes/articles/`.
  Nothing else local executes them, so the pkgdown job is otherwise the first
  thing to notice a break.

## Code and tests

- Package functions in `R/`, tests in `tests/testthat/`, evaluation and
  measurement scripts in `scripts/`.
- Fast deterministic invariants belong in `tests/testthat/`; multi-seed,
  stochastic, report-producing and visual evaluations belong in `scripts/`,
  sharing one metric implementation wherever practical.
- Document public functions with roxygen2 and regenerate after API changes.
- Keep simulation assumptions, units, schemas, tolerances and seeds explicit.
- Use only public or package-generated data in committed tests, reports and
  examples. Never commit sensitive, proprietary or patient-level data.
- `scripts_private/` is where real study data is worked on. `.gitignore` ignores
  everything in it except code files named one at a time in an allowlist, so
  adding a file there commits nothing until someone amends that list — leave both
  to the maintainer. Never `git add -A` from the repository root either: the
  study templates write source-derived CSV to `output_*/` beside the working
  directory when their chunks are run interactively.
- Treat `data/` and `output/` as generated unless told otherwise.
- Preserve unrelated changes and avoid adding dependencies unnecessarily.


## Keeping documentation true

- Code and regression tests are the source of truth for behavior. Existing prose
  is context to audit, not evidence that an algorithm, default, formula or
  limitation still works as described.
- After changing simulation, design inference, privacy accounting, a public API
  or output structure, search the repository for the affected function, argument
  and dataset names and fix every document that mentions them, plus roxygen and
  `README.md`. Search for renamed functions and obsolete algorithm terms before
  considering the update complete.
- **Leave `NEWS.md` alone.** Nothing has been released, so there is no version to
  record a change against and git history is the record. It holds one bullet
  saying so; do not add a second. This starts at the first release.
- Rewrite or remove any section that no longer matches the implementation. Never
  preserve stale technical detail to minimize a diff.
- Spell out and briefly explain every acronym and abbreviation at its first use
  in a document, including ones that feel obvious in context (DP, PMX, PK, PCA,
  BLOQ, AR(1), ADaM). Write "differential privacy (DP)" once, then use the short
  form. When a term is a method or product name rather than an initialism
  (AVATAR, `synadam`), say what it is instead of inventing an expansion.
- Keep vignette examples running through the current public API, and do not
  reimplement package algorithms in a chunk. Where practical, turn an important
  documentation claim into an assertion or a regression test.
- A clean knit proves a document executes, not that it is correct. Reason about
  the prose directly, and report a claim you could not verify rather than
  presenting it as established behavior.