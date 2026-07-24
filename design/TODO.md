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

Published as pkgdown articles (`vignettes/articles/`, no `R CMD check` cost):

- `feasibility.Rmd` — **scope**, what is achievable at which cohort size.
- `privacy-background.Rmd` — **intuition**, where `d`, `f`, and the error law
  come from. Start here if the arithmetic is unfamiliar.
- `privacy-argument.Rmd` — **proof**, the formal mechanism-level argument.
- `avatar-mathematics.Rmd` — **algorithm**, the default generator step by step.
- `model-elicitation.Rmd` — **inputs**, the interview that produces a public
  structural model and priors before any data is read.
- `data-elicitation.Rmd` — **structure**, the trial-design ladder and which
  parts of a protocol are actually public.

Keep items here short and link out. When an item closes, tick it and update the
registry entry it points at rather than deleting the history.

---

The `synpmx` rename is complete (2026-07-23), in all three places: the package
(`DESCRIPTION`, `R/synpmx-package.R`, `man/`, `build.sh`, `synpmx.Rproj`, all
five vignettes by filename, title, and `\%\VignetteIndexEntry`, and every
mention in prose), the GitHub repository `iamstein/synpmx` with `origin`
repointed, and the local clone directory `~/git/synpmx`. The exported API kept
its `pmx_*` / `synpmx_avatar()` names deliberately: `pmx_` says what the data
is, `syn` says what the package does, and `synadam` does not prefix its own
functions either.

## Next

1. **Finish publishing the pkgdown site.** `usethis::use_github_pages()` has
   been run (2026-07-23): it created the `gh-pages` branch and activated Pages
   through the API, so no manual Settings step is needed to turn it *on*. Two
   things still block the site:

   - **Pages is publishing from the wrong source.** The API reports
     `source: main /`, but the workflow deploys to `gh-pages`, so nothing the
     workflow produces is ever served. Repoint it, either in Settings → Pages
     or with `usethis::use_github_pages(branch = "gh-pages")`.
   - **The push token needs `workflow` scope.** The credential git uses has
     only `repo`. GitHub refuses any push that creates or modifies a file
     under `.github/workflows/`, so pushing `pkgdown.yaml` fails until the
     scope is added at <https://github.com/settings/tokens>.

   The site serves at <https://iamstein.github.io/synpmx/>, the URL already
   declared in `DESCRIPTION`, `_pkgdown.yml`, and every roxygen and vignette
   cross-link — those links stay broken until this is done. Confirm the first
   deploy rendered all nine documents and that `AGENTS.html` and `CLAUDE.html`
   are absent (`pkgdown/prune-site.R` removes them).

2. **Decide what to do with `.github/workflows/r.yml`.** Added through the
   GitHub UI on 2026-07-23 from GitHub's default R template. It will fail on
   every push as written: it checks against R 3.6.3 and 4.1.1, but
   `DESCRIPTION` requires R >= 4.1.0, so the 3.6.3 leg cannot even install the
   package. Either delete it — `./build.sh` already does a stricter check
   locally — or replace it with `usethis::use_github_action("check-standard")`,
   which uses the maintained r-lib matrix.

3. **Try the approach on the internal PIT565 data.** The methods have only been
   exercised on public `nlmixr2data` sources and package fixtures. Running
   AVATAR and the calibrated structural path on a real internal study is the
   test that matters: role declaration against a real schema, event grammar
   that the template sampler has not seen, and whether the generated data is
   actually useful for workflow development. Keep it in `scripts_private/`.

4. **Decide how date/datetime columns should be handled.** Low priority; dates
   are rarely analysis-relevant. Today `time` must be numeric elapsed time, and
   a raw `RFSTDTC`-style datetime column is either converted by the user
   beforehand or dropped as undeclared. Options if this ever matters: accept a
   datetime `time` and derive elapsed time from the first event per subject; or
   offer a `keep`-like path that shifts dates to a synthetic origin so they stay
   internally consistent without carrying a real calendar date. Not needed for
   PIT565.

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
- [ ] `REV-020` — `pmx_structural_model(rx = )` is stored but never used. It now
      warns; either wire it through `rxode2::rxSolve()` with a regression test
      against the analytic solution, or reject it outright.

## Version 4 — return to AVATAR blending as the primary method

Scope decision (2026-07-22): after comparing to Novartis's `synadam` (which
resamples each column marginally from the data with no formal guarantee), AVATAR
is the trajectory-level analogue of the same governance-based approach, and it is
the right default for trusted-environment synthetic data. The DP (v2) and structural
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
- [ ] `REV-024` Reframe the trust boundary as organizational rather than
      geographic. As written, the docs say AVATAR output must stay in the
      validated environment, which forbids the use case the package exists
      for: taking synthetic data out to a local machine for code development.
      `REV-024` carries the approved wording, the full site inventory, and the
      three sites that must *not* change. Own commit, own branch.

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
