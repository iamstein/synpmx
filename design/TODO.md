# TODO

Living task list. One line per actionable item, newest scope at the top.

How this relates to the other documents. `AGENTS.md` explains the three tiers
and which one new writing belongs in.

Internal design record (`design/`, cited by nothing shipped):

- `design/TODO.md` (this file) — **what to do next.** The working queue.
- `design/REVIEW_BACKLOG.md` — **why**, for defects and design findings. `REV-###`.
- `design/TEST_SIM.md` — **evidence**, for simulation defects and their gates. `SIM-###`.
- `design/METHOD_DISCUSSION.md` — **tradeoffs**, AVATAR blending vs formal DP,
  and why AVATAR is the trajectory-level analogue of synadam.
- `design/PROTOTYPE_SPEC.md` — **contract**, the specification being implemented.


## Next Steps: adoption, external validation, and the ACoP poster

The package, its documentation, and the live site are all in place. The work now
is to get the method used by others, understood well enough to defend, and
presented. The ACoP (American Conference on Pharmacometrics) poster is the
culminating deliverable, and it paces the rest.

Try it on real data — no mode has run on a real study yet, and that is the test
that matters: role declaration against a real schema, event grammar the template
sampler has not seen, and whether the output is actually useful for workflow
development. Keep all of it in `scripts_private/`.

- [ ] Run AVATAR and the calibrated path on the internal INTERNAL_STUDY data.
- [ ] Ask Alex to try it on his own dataset.
- [ ] Consider asking Anwesha or Bambang to try it as well.
- [ ] Turn tester friction into `REV-###` entries in `REVIEW_BACKLOG.md` as it
      surfaces, so external use feeds back into the package.

Understand it well enough to defend it.

- [ ] Read the original AVATAR paper (`references/Guillaudeux23.pdf`). Feeds the
      novelty positioning owed under "Verification owed" below.
- [ ] Write my own mental map of how AVATAR works, checked against
      `vignettes/articles/avatar-mathematics.Rmd`, to confirm my understanding.

Documentation the owner asked for (2026-07-25). Both are writing tasks that
depend on inputs not yet in hand, so they are queued, not started:

- [ ] A simplified, made-basic-and-clear explanation of how AVATAR works, from
      the 2023 paper (`references/Guillaudeux23.pdf`; now readable — poppler
      installed 2026-07-25). Likely lives in the method vignette or a new short
      article. Note the paper's k = 20 headline / k = 4 floor when writing the
      k-choice justification.
- [ ] A "what synpmx adds to AVATAR" section: the BLOQ/censoring handling, the
      full event-table (dosing + observation grammar) rather than a flat matrix,
      and — once designed and built — the minimum-donor pooling and rare-event
      pooling from `REV-025`/`REV-026`. Do not document the pooling as a feature
      until it exists.

Present it.

- [ ] Follow up with David Zhang: a presentation following their template. Then
      potentially with Greg Pinhault.
- [ ] Confirm the ACoP abstract status and the poster deadline. The abstract is
      drafted (`communications/2026-acop-abstract.md`) — is it submitted and
      accepted, and by when is the poster due? This gates everything above.
- [ ] Make the ACoP poster (`communications/2026-acop-poster-notes.md`). Public
      or fixture data only, and internal review before submission — see
      `communications/README.md`; any internal-study figure stays in
      `scripts_private/` and never on the poster.

Deferred.

- [ ] Decide how date/datetime columns should be handled. Low priority; dates
      are rarely analysis-relevant. Today `time` must be numeric elapsed time,
      and a raw `RFSTDTC`-style datetime column is either converted by the user
      beforehand or dropped as undeclared. Not needed for INTERNAL_STUDY.

## Done: documentation reorganization (2026-07-23)

Decided and executed. The reasoning, the audience analysis, and the rationale
for each call are in `design/DOCUMENTATION_SCOPE.md`; delete that file once
this section is stale. `AGENTS.md` now records the resulting three-tier rule.

- [x] Adopt pkgdown. `_pkgdown.yml` with a grouped reference index over all 31
      exports, plus a GitHub Actions workflow deploying to `gh-pages`.
- [x] `README.Rmd` → `README.md` as the entry point: pitch, one runnable
      example, the four-mode table, and the documentation map. 247 lines of
      API reference and limitations came out.
- [x] Vignette set cut from five to three: `synpmx-method` (all four modes,
      high level), `synpmx-demo`, `synpmx-privacy`. `synpmx-intro` and
      `synpmx-epsilon-exploration` were merged away.
- [x] Deep AVATAR mathematics moved out of the method vignette into
      `vignettes/articles/avatar-mathematics.Rmd`.
- [x] Five design documents moved to `vignettes/articles/`: feasibility,
      privacy background, privacy argument, model and data elicitation.
- [x] Every citation into `design/` removed from shipped documentation and
      roxygen comments; roxygen now links to the website.
- [x] Deleted `design/METHODS_VIGNETTE_SPEC.md` (stale) and `scripts/README.md`
      (stale). `NEWS.md` reduced to a stub.

### Earlier scope: one entry point, four named modes

- [x] An entry point covering all four modes applied to `theo_md`, a properties
      table, and a table mapping environments (trusted / partner / published)
      to acceptable modes. Now the four-mode tour in `synpmx-method`.
- [x] Privacy vignette: explain what AVATAR *is* and what DP *is* before
      comparing them, with the formal `(epsilon, delta)` definition and the
      kind-not-degree table.
- [x] Method vignette: all options at the top, then the default AVATAR
      algorithm, then the model-based alternatives at the end.
- [x] Demo vignette: state the four modes up front and run the model-based path
      on theophylline, as `scripts/demo_nlmixr2data.R` does.
- [x] `README.md`: table of contents explaining how the documentation set is
      organized and which document answers which question. Superseded by the
      full `README.Rmd` rewrite above.
- [x] House style: spell out every acronym on first use in a document
      (`AGENTS.md`). The word "mock" is gone: prose says "synthetic data", the
      `compare_pmx()` argument and outputs are `synthetic`, and generated
      character/factor IDs are `syn_001` rather than `mock_001`.
- [x] `design/DOCUMENTATION_SCOPE.md` — inventory of all 23 documents with
      guessed audiences, since rewritten as the decision record.

## Version 4 — return to AVATAR blending as the primary method

Scope decision (2026-07-22): after comparing to Novartis's `synadam` (which
resamples each column marginally from the data with no formal guarantee), AVATAR
is the trajectory-level analogue of the same governance-based approach, and it is
the right default when the synthetic data reaches no one the source data could
not. The DP (v2) and structural
(v3) engines are kept as superseded alternatives for when a formal guarantee is
required, not removed.

- [x] `design/METHOD_DISCUSSION.md`: the AVATAR vs DP tradeoff essay, the
      trajectory-is-a-fingerprint asymmetry, and the synadam parity argument.
- [x] `PROTOTYPE_SPEC.md`: Version 4 section (history) plus a section 0 banner
      making AVATAR the default and the trust-boundary the decision rule.
- [x] Restore the AVATAR engine as `synpmx_avatar()` (renamed from the v1
      `mock_pmx`), `synthesis.R`, `profiles.R`, plus the ported `utils.R`
      helpers. Exported and working; no name collisions with v2/v3.
- [x] AVATAR tests, including all five nlmixr2data datasets. `test-avatar.R`
- [x] Rebuilt the demo vignette around AVATAR, keeping all five nlmixr2data
      datasets, and the method vignette from the Version 1 "How synpmx
      Works".
- [x] Method vignette explains the (epsilon, delta) vs AVATAR distinction and
      the synadam parity argument.
- [x] Slimmed the two DP vignettes to short "formal-privacy alternative" asides.
- [x] `./build.sh` clean, all four vignettes knit, 395 tests pass.

Version 4 is complete.

Hardening completed after review: `REV-018` and `REV-019` are fixed with
regression tests in `test-structural-v3.R` and `test-avatar.R`. The remaining
open findings below apply to the superseded formal-DP v2 path.

## Superseded: v3 low-dimensional structural generator (kept as an alternative)

The scope decision (2026-07-22): prefer small trials over pooled corpora. That
requires releasing a handful of parameters against public structural priors
instead of a dense grid. See `vignettes/articles/feasibility.Rmd` section 8 and
`design/PROTOTYPE_SPEC.md` "Version 3 scope".

### Core: model in, correction out

- [x] `pmx_structural_model()` — public structural model with built-in analytic
      1-cmt IV/oral/infusion PK and direct-effect or indirect-response PD.
      Built-ins need no compiler; `rxode2` is accepted and validated but not yet
      wired through. `R/structural.R`
- [x] `pmx_trial_design()` — dose levels, cohort sizes, protocol sampling
      schedule, dosing interval, infusion duration, visit windows.
- [x] Multiplicative correction release rather than absolute parameters.
      `synpmx_calibrated()`, `d = 2-3`.
- [x] `pmx_prior(range, source)` / `pmx_priors()` with mandatory provenance.
- [x] Per-subject NCA estimator. The correction is the ratio of predicted to
      observed AUC on the subject's own grid, which avoids needing `F` or an
      extrapolation to infinity.
- [x] Sampling schedule from the protocol. The timing release group is gone.
- [x] `at_prior_boundary` diagnostic, warned on and shown by `print()`.
- [x] Prior mode is first-class: `pmx_generate(model, design)` takes no data and
      no epsilon.
- [ ] Support ladder levels 3-4 (intra-patient escalation, titration):
      occasion-varying dose applied to *generated* subjects from the public
      rule. Never replay a source subject's sequence.
- [ ] `REV-017` Record realized trial-design quantities in `proof_assumptions`
      with a source field.
- [x] Two-compartment PK (`2cmt_iv`, `2cmt_oral`), analytic and verified
      against `Dose/CL`.
- [ ] Transit-absorption option.
- [x] Covariates: **out of scope**, stated in `PROTOTYPE_SPEC.md` section 6.
      Users join their own; covariate-handling code is not exercised by this
      package's output.
- [ ] Schema flexibility. Output is a fixed column set; real datasets are not.

### Guardrails

- [ ] Enforce public model selection ergonomically: the model is built and
      passed in before `synpmx_empirical()` sees data, so there is no easy route
      to fit, look, and revise.
- [ ] Guard against the prior range being set from the data without budget.
      Mean +/- 2 SD is the right target for a private stage 1 to *discover*;
      computing it directly voids the guarantee undetectably.
- [ ] Document that confidential data must never be sent to an external service,
      including an LLM. Outside the accounting entirely.

### Optional, if priors turn out too wide

- [ ] Two-stage range-finding: coarse log-scale histogram at `epsilon_0` (~20%
      of budget, L1 sensitivity 1 regardless of bin count), then clip and
      release means at `epsilon_1`. Measured worth with an absolute prior:
      2.44-fold -> 1.52-fold at N = 20, epsilon 1. Likely unnecessary once the
      release is a correction factor, since that prior is already tight.
- [ ] Small budget for realized messiness rates (dropout, missed dose, BLQ) if
      public assumptions prove too crude.

### Verification owed

- [ ] Verify the `rxode2` templates in `vignettes/articles/model-elicitation.Rmd` compile and
      produce sensible profiles. Authored without a working C compiler and never
      executed.
- [ ] Measure the correction-factor parameterization. The ~1.37-fold estimate at
      N = 20, epsilon 1 is arithmetic from the error law, not a measurement.
- [x] Measure the PD correction. Exact without residual error; biased low by
      about a third with 15% residual on a small deviation. Documented in
      `vignettes/articles/feasibility.Rmd`; PD is experimental.
- [x] Improve the PD estimator. Solved by changing the endpoint rather than the
      statistic: simple time-course shapes take a level correction (ratio of
      means), which is unbiased under residual error. Exposure-driven PD keeps
      the signed-area statistic and stays experimental.
- [x] Retire the exposure-driven PD shapes. Removed entirely; PD is now a
      simple exposure-independent time course.
- [ ] Literature check before claiming novelty: DP + non-compartmental analysis,
      DP + popPK, DP synthetic data under informative structural priors.

## Next: correctness and privacy findings

- [ ] `REV-003` Data-dependent `stop()` on confidential rows before any noise is
      applied. An unaccounted output channel in the core DP claim.
- [ ] `REV-016` Pooling breaks the stated adjacency when one person appears more
      than once (rollover, extension, crossover, re-enrollment). Group privacy
      silently degrades the guarantee to `k * epsilon`.
- [ ] `REV-004` `delta` is validated and requires justification but is never
      spent. Either restrict to `delta = 0` or wire it to a real mechanism.
- [ ] `REV-002` Pre-flight feasibility check, so an infeasible `(N, epsilon, d)`
      is refused before budget is spent rather than discovered in a plot.
- [ ] `REV-020` `pmx_structural_model(rx = )` is stored but never used; it now
      warns. Either wire it through `rxode2::rxSolve()` with a regression test
      against the analytic solution, or reject it outright. (Was stranded in the
      completed documentation-reorganization list.)
- [~] `REV-025` **Owner-flagged, 2026-07-25. Core fix landed.** AVATAR now
      borrows the nearest donors across dose/schedule groups to blend `k` (=5)
      real patients into every avatar, with a loud red alert only when the
      source has fewer than `k + 1` subjects (`.select_donors()`,
      `.loud_warn()`; `test-avatar-pooling.R`). Discovery: weight-based dosing
      made nearly every subject its own singleton, so this was the common case,
      not an edge case. Outlier detector done too:
      `flag_identifiable_subjects()` screens follow-up time, dose count, dose
      magnitude, and DV. Still open: whether the headline floor should exceed 5.
- [ ] `REV-026` **Consider, not committed (2026-07-25 decision).** "Skeleton
      sampling." Today each avatar copies one real anchor's *event skeleton* —
      the number of doses, their sizes, and the observation times — verbatim,
      and only the DV values are blended. So a source subject with a unique
      structure (the long-followed `wbcSim` subject, a one-off dose) yields an
      avatar with that same unique structure, which is identifying even though
      its values are blended. Skeleton sampling would instead *draw* each
      avatar's schedule (follow-up length, dose count, observation times) from
      the cohort distribution, so no avatar carries any one real subject's
      unique structure and there is nothing to flag. It is the thorough,
      preventive cure.
      **Why it is paused:** `flag_identifiable_subjects()` +
      `remediate_identifiable_subjects()` (with replacement) already remove
      structural outliers after generation and refill the cohort, which covers
      the same privacy goal from the other direction. The owner chose the
      detect-and-remediate path for now; skeleton sampling stays a possible
      future direction, to weigh against the diversity it would add versus its
      cost. Design notes in `design/METHOD_DISCUSSION.md` §6a; would build on
      the `REV-025` relaxed pooling.

## Then: utility headroom in the existing dense-grid path

Keep this path for pooled corpora; it is not superseded by v3.

- [ ] `REV-005` Derive sensitivity from the declared contribution limits instead
      of matrix width.
- [ ] `REV-006` Gaussian/zCDP accounting. **Note:** only worthwhile at high
      dimension. At `d` around 6 pure-DP Laplace is the better mechanism, so
      this belongs to the dense-grid path, not to v3.
- [ ] `REV-007` Low-dimensional trajectory representation. Largely subsumed by
      the v3 work above.
- [ ] `REV-008` Name and document the post-processing rescue constants; add a
      diagnostic for how much output shape comes from the release versus the
      fallbacks.

## Hygiene

- [ ] `REV-009` `.public_input_manifest()` ignores all six arguments and returns
      a constant. Derive it or delete it.
- [ ] `REV-010` Replace the nine-name payload denylist with an allowlist.
- [ ] `REV-011` Add `print.pmx_private_validation` and `print.pmx_backend_tests`.
- [ ] `REV-012` `run_dp_backend_tests()` advertises a `privacy_map` test it does
      not actually assert.
- [ ] `REV-013` `enable_features("contrib")` mutates global OpenDP state on every
      backend resolve.
- [ ] `REV-015` Convert load-bearing vignette prose into assertions or tests.

## Done

- [x] `REV-001` / `SIM-020` Scale-aware support threshold in decoding. `44db89f`
- [x] `REV-014` Land the working tree in reviewable commits.
- [x] `./build.sh` for `R CMD check` and clean-library vignette rendering. `778848b`
- [x] `vignettes/articles/feasibility.Rmd` scoping assessment. `d5b0e30`
- [x] `REV-023` Session-level `synpmx_enable_dp_engines()` gate on the DP
      engines' unaudited status, so it is enforced rather than only documented.
- [x] `REV-024` Trust boundary reframed as organizational rather than
      geographic, so taking synthetic data to a local machine for code
      development is a supported use rather than a boundary crossing.
- [x] Publish the pkgdown site. Live and verified 2026-07-24: landing page,
      articles, and reference all render, and `AGENTS.html` / `CLAUDE.html` are
      pruned. The stringfish/RcppParallel ABI break that failed every build is
      fixed by keeping `rxode2` out of the docs job. `f1326fe`
- [x] Delete `.github/workflows/r.yml`, which could never pass (R 3.6.3 leg
      against `DESCRIPTION` R >= 4.1.0); `./build.sh` checks locally. `2521f8e`
- [x] One README. Deleted `README.Rmd`; `README.md` is the only entry point,
      its example pinned by `tests/testthat/test-readme.R`. `e09011d` `215cf50`
