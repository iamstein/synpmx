# Repository instructions

This repository contains the `synpmx` R package prototype for simulating
structurally faithful synthetic pharmacometric datasets.

## Design documents

Read `design/TODO.md` first. It is the working task queue and it links to the
document that owns each kind of decision:

- `design/TODO.md` — what to do next.
- `design/REVIEW_BACKLOG.md` — defects and design findings (`REV-###`).
- `design/TEST_SIM.md` — simulation defects and their regression gates (`SIM-###`).
- `design/FEASIBILITY.md` — what is achievable at which cohort size.
- `design/PRIVACY_BACKGROUND.md` — how the privacy arithmetic works: `d`, `f`,
  sensitivity, and the error law, with worked examples.
- `design/PRIVACY_ARGUMENT.md` — the formal mechanism-level argument, for a reviewer.
- `design/MODEL_ELICITATION.md` — the interview producing the public structural
  model and priors, and the rules that keep those inputs data-independent.
- `design/DATA_ELICITATION.md` — the trial-design complexity ladder, and which
  parts of a protocol are genuinely public.
- `design/PROTOTYPE_SPEC.md` — the specification being implemented.

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
  does both against a clean temporary library.

## Simulation testing and evaluation

- Treat `design/TEST_SIM.md` as the living evaluation specification. Keep its
  dataset registry, issue registry, metrics, and acceptance gates synchronized
  with the implemented evaluator.
- Evaluate every public dataset used by the demo. Add a focused
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
  update the demo for workflow/output changes, the simulation-method vignette
  for generator changes, the privacy-introduction vignette for privacy changes,
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
