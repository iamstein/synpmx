# TODO

Living task list. One line per actionable item, newest scope at the top.

How this relates to the other design documents:

- `design/TODO.md` (this file) — **what to do next.** The working queue.
- `design/REVIEW_BACKLOG.md` — **why**, for defects and design findings. `REV-###`.
- `design/TEST_SIM.md` — **evidence**, for simulation defects and their gates. `SIM-###`.
- `design/FEASIBILITY.md` — **scope**, what is achievable at which cohort size.
- `design/MODEL_ELICITATION.md` — **inputs**, the interview that produces a
  public structural model and priors before any data is read.
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

- [ ] `pmx_structural_model(rx, typical, source)` accepting an `rxode2` model as
      a public input. `rxode2` is a Suggests dependency needed only at
      generation time. Fall back to the v2 grid or Mode A when absent.
- [ ] Release a **multiplicative correction** to the model's prediction, not an
      absolute parameter. Draft vector `d = 3`: cohort size, PK correction, PD
      correction. The prior on "how wrong is the prediction" is ~8-fold and
      free; a prior on absolute CL is ~100-fold.
- [ ] `pmx_prior(range, source)` with required data-independent provenance,
      recorded in the release ledger. Replaces `pmx_bounds()` as the dominant
      sensitivity term.
- [ ] Per-subject NCA (trapezoidal AUC, terminal slope) as the estimator. Each
      subject's value must depend only on that subject's own rows — this is what
      makes the sensitivity argument work and why NLME fitting cannot be used.
- [ ] Sampling schedule from the protocol as a public input. Removes the
      `endpoint_timing` release group entirely (36 dimensions in the fixture).
- [ ] Report the released correction as a diagnostic. Already public and
      accounted, so free, and a correction pressed against the clipping boundary
      is the one signal that the prior was wrong.
- [ ] Public-design-only (Mode A) as a first-class entry point taking no data
      and no epsilon — not reachable only via the `backend = "public"` backdoor.

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
- [ ] Measure PD correction, PD baseline, and t-half. Only CL has been confirmed
      against the error law; terminal-slope estimates are noisier per subject.
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
