# TODO

Living task list. One line per actionable item, newest scope at the top.

How this relates to the other design documents:

- `design/TODO.md` (this file) — **what to do next.** The working queue.
- `design/REVIEW_BACKLOG.md` — **why**, for defects and design findings. `REV-###`.
- `design/TEST_SIM.md` — **evidence**, for simulation defects and their gates. `SIM-###`.
- `design/FEASIBILITY.md` — **scope**, what is achievable at which cohort size.
- `design/PROTOTYPE_SPEC.md` — **contract**, the specification being implemented.

Keep items here short and link out. When an item closes, tick it and update the
registry entry it points at rather than deleting the history.

---

## Now: v3 low-dimensional structural generator

The scope decision (2026-07-22): prefer small trials over pooled corpora. That
requires releasing a handful of parameters against public structural priors
instead of a dense grid. See `design/FEASIBILITY.md` section 8 and
`design/PROTOTYPE_SPEC.md` "Version 3 scope".

- [ ] Decide and document the released parameter vector. Target `d <= 8`. Draft
      in `PROTOTYPE_SPEC.md` section 6: cohort size, CL, t-half, PD baseline, PD
      effect magnitude, PD onset rate.
- [ ] Add public structural priors to the configuration: `pmx_prior(range,
      source)`, where `source` records data-independent provenance. This range
      replaces `pmx_bounds()` as the dominant sensitivity term.
- [ ] Implement two-stage range-finding: a coarse log-scale histogram at
      `epsilon_0` (~20% of budget, L1 sensitivity 1 regardless of bin count),
      then clip and release means at `epsilon_1`. Measured worth: 2.44-fold ->
      1.52-fold at N = 20, epsilon 1.
- [ ] Guard against the range being set from the data without budget. Mean +/-
      2 SD is the right target for stage 1 to *discover*, but computing it
      directly voids the guarantee and nothing downstream detects it.
- [ ] Measure `t-half`, PD baseline, and PD magnitude. Only CL has been
      confirmed against the error law; terminal-slope estimates are noisier per
      subject and may behave worse.
- [ ] Implement per-subject non-compartmental summaries (trapezoidal AUC,
      terminal slope) as the estimator. Each subject's value must depend only on
      that subject's own rows — this is what makes the sensitivity argument
      work, and it is why NLME/popPK fitting cannot be used here.
- [ ] Take the sampling schedule from the protocol as a public input rather than
      learning it. This removes the `endpoint_timing` release group entirely
      (36 dimensions in the current fixture) and frees its budget.
- [ ] Spend a small budget on realistic messiness instead: dropout rate, missed
      dose rate, BLQ fraction. Mock data that is too clean tests pipelines
      poorly, which defeats the package's stated purpose.
- [ ] Add a declared-assumption check for linear PK. Dose-normalization is only
      valid when it holds; the assumption must be explicit and testable.
- [ ] Measure the realized frontier and replace the arithmetic in
      `design/FEASIBILITY.md` section 8 with measurements.
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
