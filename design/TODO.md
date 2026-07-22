# TODO

Living task list. One line per actionable item, newest scope at the top.

How this relates to the other design documents:

- `design/TODO.md` (this file) — **what to do next.** The working queue.
- `design/REVIEW_BACKLOG.md` — **why**, for defects and design findings. `REV-###`.
- `design/TEST_SIM.md` — **evidence**, for simulation defects and their gates. `SIM-###`.
- `design/FEASIBILITY.md` — **scope**, what is achievable at which cohort size.
- `design/PRIVACY_BACKGROUND.md` — **intuition**, where `d`, `f`, and the error
  law come from. Start here if the arithmetic is unfamiliar.
- `design/PRIVACY_ARGUMENT.md` — **proof**, the formal mechanism-level argument.
- `design/MODEL_ELICITATION.md` — **inputs**, the interview that produces a
  public structural model and priors before any data is read.
- `design/DATA_ELICITATION.md` — **structure**, the trial-design ladder and
  which parts of a protocol are actually public.
- `design/PROTOTYPE_SPEC.md` — **contract**, the specification being implemented.

Keep items here short and link out. When an item closes, tick it and update the
registry entry it points at rather than deleting the history.

---

## Now: v3 low-dimensional structural generator

The scope decision (2026-07-22): prefer small trials over pooled corpora. That
requires releasing a handful of parameters against public structural priors
instead of a dense grid. See `design/FEASIBILITY.md` section 8 and
`design/PROTOTYPE_SPEC.md` "Version 3 scope".

### Core: model in, correction out

- [x] `pmx_structural_model()` — public structural model with built-in analytic
      1-cmt IV/oral/infusion PK and direct-effect or indirect-response PD.
      Built-ins need no compiler; `rxode2` is accepted and validated but not yet
      wired through. `R/structural.R`
- [x] `pmx_trial_design()` — dose levels, cohort sizes, protocol sampling
      schedule, dosing interval, infusion duration, visit windows.
- [x] Multiplicative correction release rather than absolute parameters.
      `fit_calibrated_pmx()`, `d = 2-3`.
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
      passed in before `fit_private_pmx()` sees data, so there is no easy route
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

- [ ] Verify the `rxode2` templates in `design/MODEL_ELICITATION.md` compile and
      produce sensible profiles. Authored without a working C compiler and never
      executed.
- [ ] Measure the correction-factor parameterization. The ~1.37-fold estimate at
      N = 20, epsilon 1 is arithmetic from the error law, not a measurement.
- [x] Measure the PD correction. Exact without residual error; biased low by
      about a third with 15% residual on a small deviation. Documented in
      `design/FEASIBILITY.md`; PD is experimental.
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
- [x] `design/FEASIBILITY.md` scoping assessment. `d5b0e30`
