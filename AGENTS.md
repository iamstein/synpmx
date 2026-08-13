# Repository instructions

This repository contains the `synpmx` R package prototype for simulating
structurally faithful synthetic pharmacometric datasets.

## Where documentation lives

Three tiers, and the tier decides both the audience and the maintenance cost.

**`design/` — internal record. Never cited from anything shipped.** A reader who
installed the package cannot follow a `design/` path, so vignettes, articles,
and roxygen comments must not reference one.

- `design/TODO.md` — what to do next. Read this first; it is the working queue.
- `design/REVIEW_BACKLOG.md` — defects and design findings (`REV-###`).
- `design/TEST_SIM.md` — simulation defects and their regression gates (`SIM-###`).
- `design/METHOD_DISCUSSION.md` — AVATAR-versus-DP tradeoffs, `synadam` parity.
- `design/PROTOTYPE_SPEC.md` — the specification being implemented.
- `design/SYNTHETIC_DATA_CHECKS.md` — the taxonomy of checks to run on
  generated data, and the specification for the vignette that will present it.
- `design/BUILD_DOCUMENTATION.md` — every build command, where the rendered
  HTML lands, and which check catches what.
- `design/WRITING_STYLE.md` — prose rules for vignettes, articles and README,
  derived from the maintainer's own edits. Read before drafting or revising any
  document. Its companion `design/_THEORY_OF_MIND.md` records how he reads.
- `design/learning/` — the maintainer's own learning record, not a package
  document. `EFFECTIVE_LEARNING.md` is the plan, `GLOSSARY.md` holds definitions
  written from memory, and `QUESTIONS.md` logs questions verbatim and dated.
  Add to `QUESTIONS.md` when a question is answered; never rewrite his wording.

**`vignettes/` — shipped, and rebuilt by `R CMD check` on every behavioral
change.** Keep this set small; each one is a recurring cost, not just a
document.

- `demo.Rmd` — one dataset (`xgxr::case1_pkpd`) end to end, and deliberately
  short: roles, generation, source-versus-synthetic plots and distributions,
  then `synpmx_scorecard()` in one line. Written 2026-08-04, modelled on the
  private per-study try files; cut to this shape on 2026-08-13 because the
  checking half belongs in `scorecard-synthetic-data-checks.Rmd`, which the demo now
  points at rather than restating. Keep it that way: anything that explains a
  check rather than running it goes in that vignette. Chosen for `case1_pkpd`
  because all 180 patients have a unique observation schedule as recorded and
  coarsening onto `NOMTIME` takes that to zero, which is the clearest
  demonstration of the visit grid available in public data.
- `synpmx-4-methods.Rmd` — the four generation modes at a high level.
- `avatar-algorithm.Rmd` — the default generator step by step, and the seven
  masking mechanisms M1-M6. Promoted from `articles/` on 2026-07-31 because it
  is the method reference, not supporting evidence.
- `avatar-evaluation-public-data.Rmd` — `synpmx_avatar()` run over every public
  dataset in the evaluation set, then the exposure measurements across all
  eight. Renamed from `synpmx-demo.Rmd` on 2026-08-04: it is an evaluation of
  algorithm performance, not a tutorial. It carries a table describing every
  dataset it uses; keep that table in step with the datasets actually run.
  Reworked on 2026-08-13 so all eight worked examples have the same shape as
  `demo.Rmd` — roles, `synpmx_avatar()`, one source-versus-synthetic figure,
  `compare_pmx_distributions()`, `synpmx_scorecard()` — followed by a subsection
  per dataset only where something fails or reads oddly, which is where
  `pmx_masking_report()` and any dataset-specific digging belong. Keep it that
  way: a new dataset gets the same five steps, and prose only where the numbers
  need it.
- `scorecard-synthetic-data-checks.Rmd` — the six categories of check to run on generated
  data (A validity, B who is singled out, C same study, D distributions, E what
  the checks cannot tell you, F check the output not the algorithm), worked on
  `xgxr::case1_pkpd` with `nlmixr2data::pheno_sd` as the case where they fail.
  Written 2026-08-04 from the specification in `design/SYNTHETIC_DATA_CHECKS.md`;
  it names its own gaps, and that list is the queue for what to build next.
  Renamed from `synthetic-data-checks.Rmd` on 2026-08-13: it is the reference
  for `synpmx_scorecard()`, and the name should say so. Its scorecard table is
  the contract that function implements — a row's pass criterion changing in one
  without the other is a defect.

**`vignettes/articles/` — pkgdown only.** Excluded from the build by
`.Rbuildignore`, so `R CMD check` never touches these and they are not shipped
in the tarball. Use this tier for teaching and evidence that does not need to be
rebuilt on every change. Note that pkgdown *executes* article code, so a broken
article fails the site build.

- `synpmx-privacy.Rmd` — the trust-boundary decision rule and choosing
  epsilon. Moved here from `vignettes/` on 2026-08-04: it belongs with the
  differential-privacy material rather than with the AVATAR guides, and it
  is not needed offline by someone running `synpmx_avatar()`. Cross-
  references to it must be website URLs, not `vignette("synpmx-privacy")`,
  because an article is not installed.
- `privacy-background.Rmd` — `d`, `f`, sensitivity, the error law.
- `privacy-argument.Rmd` — the formal mechanism-level argument, for a reviewer.
- `feasibility.Rmd` — what is achievable at which cohort size.
- `model-elicitation.Rmd` / `data-elicitation.Rmd` — producing the public
  structural model, priors, and design without reading data.

`README.md` is the entry point. It is critical that it be human readable and
understandable, and its voice belongs to the maintainer rather than to an agent.
Two different rules apply to it, and the difference is between keeping it true
and changing what it says.

- **Keep it accurate without asking.** When a rename, move, or API change makes a
  file name, link, function name, count, or description in `README.md` wrong,
  fix it in the same change that broke it, and say so in the summary. A stale
  link or a wrong file name is a defect, not an editorial choice. The same goes
  for a documentation table that no longer matches the vignettes that exist.
- **Ask before adding or rewriting.** New sections, added explanation,
  restructuring, or reworking existing prose are substantial changes to a
  document the maintainer owns. Describe what you would change and why, and get
  confirmation first.

Keep `design/TODO.md` current: tick items as they close, add newly discovered
work, and record the reasoning in the registry that owns it rather than in the
task list itself.

- Put package functions in `R/`, tests in `tests/testthat/`, and runnable
  demonstrations in `scripts/`.
- Document public functions with roxygen2 and regenerate documentation after
  API changes.
- Keep simulation assumptions, units, schemas, and seeds explicit.
- Do not commit sensitive, proprietary, or patient-level data.
- Treat `data/` and `output/` as local/generated unless told otherwise.
- Preserve unrelated changes and avoid adding dependencies unnecessarily.
- Run the full tests and `R CMD check` after behavioral changes. `./build.sh`
- `./build.sh articles` before pushing anything that touches
  `vignettes/articles/`. `.Rbuildignore` keeps articles out of the package, so
  `R CMD check` never executes them and the pkgdown job on GitHub is otherwise
  the first thing to notice a break. `./build.sh docs` does vignettes and
  articles together.
- `source("dev.R")` in the R console for the fast loop: `dev_preview("name")`
  renders one document against the working tree in seconds and opens it,
  `dev_list()` shows what exists. Rendered HTML lands in `output/`, which is
  gitignored --- that is why it does not show up in `git status`.
  does both against a clean temporary library.

## Simulation testing and evaluation

- Treat `design/TEST_SIM.md` as the living evaluation specification. Keep its
  dataset registry, issue registry, metrics, and acceptance gates synchronized
  with the implemented evaluator.
- Evaluate every public dataset used by the evaluation vignette. Add a focused
  regression fixture when that is the smallest reliable way to reproduce a
  defect.
- For every newly discovered simulator defect, add or update all three of: the
  issue entry in `design/TEST_SIM.md`, an automated regression check, and the
  implementation fix. Do not close an issue based only on visual inspection.
- Put fast, deterministic invariants in `tests/testthat/`. Put multi-seed,
  stochastic, report-producing, and visual evaluations in `scripts/`, while
  sharing one metric implementation wherever practical.
- Use only public or package-generated data in committed tests, reports, and
  examples. Keep dataset-specific assumptions, tolerances, units, and seeds
  explicit.

## Keeping documentation synchronized

- Treat current package code and regression tests as the source of truth for
  behavior. Existing vignette prose is context to audit, not evidence that an
  algorithm, default, formula, or limitation still works as described.
- In vignettes, design documents, and `README.md`, spell out and briefly explain
  every acronym and abbreviation at its first use in that document — including
  ones that feel obvious in context (DP, PMX, PK, PCA, BLOQ, AR(1), ADaM). Write
  "differential privacy (DP)" once, then use the short form. When a term is a
  method or product name rather than an initialism (AVATAR, `synadam`), say what
  it is instead of inventing an expansion.
- Preserve each vignette's audience, purpose, and broad information structure
  by default, but rewrite or remove any section that no longer matches the
  implementation. Never preserve stale technical detail merely to minimize a
  documentation diff.
- After changes to simulation, design inference, privacy accounting, public
  APIs, or output structure, make an explicit documentation-impact pass:
  update the evaluation vignette for workflow/output changes, the
  simulation-method vignette for generator changes, the privacy-introduction
  vignette for privacy changes,
  and the epsilon-exploration vignette for privacy--utility behavior. Also
  check roxygen documentation, `README.md`, `NEWS.md`, and design specifications
  where relevant.
- Verify every implementation-specific statement against the exact functions
  and tests that establish it. This includes defaults, constants, equations,
  grids, budget allocation, randomness, inferred dosing/sampling behavior,
  schemas, and known limitations. Prefer references to public functions over
  descriptions of internal call sequences that are likely to drift.
- Keep vignette examples executable through the current public API. Avoid
  reimplementing package algorithms in vignette chunks; when practical, turn
  important documentation claims into assertions or regression tests.
- After documentation-affecting changes, update the vignette code and prose to
  match the new behavior. A full clean-library re-render of every vignette is
  not required for each change; `./build.sh` rebuilds them once inside
  `R CMD check`, which is enough to prove they still execute. Use
  `./build.sh vignettes` when you want inspectable HTML to read the tables and
  plots. A successful knit is necessary but does not prove that the explanation
  is semantically correct, so reason about the prose directly rather than
  treating a clean render as verification. When you do render, never validate
  against a previously installed package, an already loaded namespace, or stale
  rendered HTML.
- Search the repository for renamed functions, old vignette names, removed
  arguments, and obsolete algorithm terms before considering the update
  complete. Report any claim that cannot be verified instead of presenting it
  as established behavior.

## Theory of Mind

Consult the design/_THEORY_OF_MIND.md and update.