# Review backlog

Prioritized findings from a code review on 2026-07-22, covering `R/`,
`design/`, `tests/`, and the vignettes at commit `bf1ebd2` plus the
uncommitted working tree.

This is a living backlog in the same spirit as `design/TEST_SIM.md`. `TEST_SIM.md`
tracks *simulation* defects found by evaluating output; this file tracks
*mechanism, privacy-accounting, and API* defects found by reading the code.
When an item here produces a reproducible output failure, mirror it into the
`TEST_SIM.md` issue registry with a gate and close it here.

Status values: `open`, `in-progress`, `closed`, `wontfix`.

---

## Part 1 — Why epsilon has to be so large right now

This section exists because the headline symptom ("utility only appears at
epsilon 50-500") has *two* independent causes that look identical in a plot.
Separating them changes what you should fix first.

### The noise arithmetic

Every source-dependent number this package releases is a **column sum over
subjects**, computed from a per-subject row clipped to `[0, 1]`. The Laplace
mechanism adds noise with scale

```
b = sensitivity / epsilon_group
```

`sensitivity` is the L1 sensitivity: the largest total change one subject can
cause across the whole released vector. `.release_matrix_sum()`
(`R/representation.R:496`) charges `sensitivity = ncol(matrix)` — the worst case
where a single subject moves every coordinate by the full 1.0.

Measured on `pmx_simulated_fixture(60)` at `epsilon = 5` with the default
budget split, using the real OpenDP backend:

| query | epsilon | sensitivity | dims | Laplace scale `b` |
|---|---:|---:|---:|---:|
| `subject_count` | 0.50 | 1 | 1 | 2 |
| `event_and_regimen` | 0.75 | 9 | 9 | 12 |
| `endpoint_timing` | 0.75 | 36 | 36 | 48 |
| `endpoint_trajectories` | 2.00 | 40 | 40 | 20 |
| `baseline_covariates` | 0.50 | 8 | 8 | 16 |
| `censoring` | 0.50 | 12 | 12 | 24 |

So the noise added to each trajectory coordinate has scale 20, while the true
value of that coordinate is a count bounded by N = 60. That is signal-to-noise
below 1 on essentially the entire release.

### The part that is not actually a problem: N

Almost everything decoded from these sums is a **ratio** — a mean or a
probability — of the form `released_sum / released_count`. The noise on that
ratio is approximately

```
error on a decoded unit-scale mean  ~=  b / N  =  sensitivity / (epsilon_group * N)
```

The `N` in the denominator is the whole story. For the trajectory group above
(`b = 20`):

| N | approx. error on a `[0,1]` unit-scale mean |
|---:|---:|
| 60 | 0.33 — unusable, the domain is only 1.0 wide |
| 200 | 0.10 |
| 600 | 0.033 |
| 2000 | 0.010 — good |

This is confirmed empirically. At N = 2000 and `epsilon = 5`, the released
per-cell presence counts came back as 1988.7, 1997.5, 2058.9, 1994.5 against a
true value of 2000 — accurate to well under 1%.

**Conclusion: the mechanism is already fine at realistic clinical N.** The
demos and fixtures are 8-60 subjects, which is the regime where *any* honest DP
mechanism fails. A pooled dataset of several hundred subjects at `epsilon = 5`
is already in usable territory with today's code. What is missing is not
primarily a better mechanism — it is telling the user, before they spend
budget, whether their (N, epsilon, dimension) combination can work at all. See
`REV-002`.

### The part that *is* a bug, and hurts at every N

Separately from noise magnitude, the decode layer compares **raw noisy counts**
against the literal constant `0.25`, as if they were probabilities in `[0, 1]`.
They are not — they are unnormalized sums of order N.

`.decode_trajectories()` (`R/representation.R:624`):

```r
value_unit <- ifelse(presence > 0.25, value_sum / pmax(presence, 1e-8), NA_real_)
```

`presence > 0.25` is meant to ask "did this grid cell have any released
support?". At N = 2000 with `b = 20`, a cell with *zero* real support routinely
draws noise of +43 and sails through the gate. Its `value_sum` is separately
noised, clamped at 0 by `pmax(..., 0)`, so the cell decodes to `0 / 43 = 0`,
which `.from_unit()` maps to the **bottom of the working domain**.

Observed directly at N = 2000, `epsilon = 5`, log-transformed `cp` endpoint
(working scale spans about -8.5 to 5.3):

```
true  presence : 2000 2000 2000    0 2000    0    0    0   0   0
noisy presence : 1989 1998 2059    0 1995 43.3 10.7    0 0.7 0.7
true  curve    : 0.088 2.11 1.41 1.03 0.434 0.434 0.434 0.434 0.434 0.434
noisy curve    : 0.14  2.43 1.06 0.873 0.584 -2.09 -8.52 -8.52 -8.52 -8.52
```

The first five cells — the ones with real support — are *accurate*: errors of
0.2-0.35 on a domain 13.8 wide, exactly the `b/N` prediction. The last five are
pinned to the domain floor. That is not noise; that is a broken gate. And
because `.fill_unoccupied_curve()` only interpolates cells that fail the same
threshold, it cannot rescue them — the bad cells are labelled "occupied".

So a meaningful share of the visual damage in the epsilon vignette is a
scale-inconsistency defect, not an inherent privacy cost. Fix `REV-001` before
concluding anything about how much epsilon this design needs.

---

## Part 2 — Prioritized issues

### P0 — Fix before drawing further conclusions about privacy/utility

| ID | Area | Issue | Suggested direction | Status |
|---|---|---|---|---|
| `REV-001` | Decoding | Raw noisy counts compared against the literal `0.25` as if they were probabilities. Occurred at `R/representation.R:592, 605, 624, 640, 656, 658, 685`. Unoccupied cells decoded to a domain endpoint instead of "no support", at every N. | **Fixed.** Added `.support_threshold(count, noise_scale)` and a backend-supplied `noise_scale()`; decoders now gate on `max(0.25, min(3 * b, count / 2))`. Tracked as `SIM-020` with regression tests in `tests/testthat/test-decode-support-threshold.R`. Measured effect at epsilon 5: empty cells pinned to the working-domain floor drop from 83% to 0% at N >= 600, and curve error at N = 2000 halves (0.099 -> 0.045). No effect at N <= 100, and none on the noiseless public-fixture path by construction. | closed |
| `REV-002` | API / guidance | Nothing warns a user that their N is too small for their epsilon and release dimension until after the budget is irreversibly spent. The existing check (`R/fit.R:252`) fires on `dimensions > 6 * count`, which is dimension-vs-N only and ignores epsilon entirely. | Add a public pre-flight helper that takes the configuration (no data) and reports the predicted `sensitivity / (epsilon_group * N)` error per release group against a stated usability threshold. Document the `N >~ sensitivity / (epsilon_group * target_error)` envelope in the epsilon vignette. This converts "unusable" into "usable within a stated envelope", which is the honest and far more useful claim. | open |

### P1 — Privacy correctness

| ID | Area | Issue | Suggested direction | Status |
|---|---|---|---|---|
| `REV-003` | DP guarantee | Data-dependent `stop()` on confidential rows *before* any noise is applied: `validate_pmx(..., strict = TRUE)` at `R/fit.R:206`, plus `R/representation.R:90` ("An observed DVID is not covered...") and `R/representation.R:215` ("A subject-property value is outside the declared public category levels"). Whether the call throws — and its message — is a function of one individual's record. This is an unaccounted output channel in exactly the place the package claims formal DP. | Drop or clip out-of-domain records rather than erroring, matching how numeric bounds are already handled by `.clip()`. If a hard failure must be retained, route it through a private test and account for it. Either way, disclose it in `proof_assumptions` (`R/fit.R:271`) and `vignettes/articles/privacy-argument.Rmd`. | open |
| `REV-004` | Accounting | `delta` is validated and requires a written justification when positive (`R/fit.R:14-24`), but every mechanism is Laplace and every accounting entry hardcodes `delta = 0` (`R/privacy.R:158`). `unspent_delta` always equals the request. Users are asked to justify a parameter that is never spent. | Either restrict the API to `delta = 0` and say so, or wire delta to a real approximate-DP mechanism — which is the same work as `REV-006`. Do not leave it as decorative. | open |

### P2 — Utility headroom (real work, do after P0/P1)

| ID | Area | Issue | Suggested direction | Status |
|---|---|---|---|---|
| `REV-005` | Sensitivity | `sensitivity = ncol(matrix)` is sound but pessimistic. A subject occupies a handful of trajectory/timing cells, not all of them; the true row L1 is bounded by `2 * min(cells, observation_limit)` per endpoint, and the timing matrix's true bound is `n_cells + 2 * max_occasions` rather than `ncol`. | Derive the bound from the contribution limits already declared in `pmx_contribution_limits()` rather than from matrix width. Cheap, no change to the proof structure, and it makes `max_timing_cells` do real work — today it only caps grid size (`R/endpoints.R:157, 182`) and never enters a sensitivity argument. | open |
| `REV-006` | Mechanism | Pure-DP Laplace with basic sequential composition costs linearly in dimension. Gaussian noise under zCDP/RDP accounting costs roughly `sqrt(d)` instead of `d`, and composes better across the six query groups. Rough estimate for the trajectory group at epsilon 5: Laplace sd ~28 vs Gaussian sd ~12 under zCDP. | Add a zCDP accountant alongside the existing basic-composition one, and a Gaussian measurement in the OpenDP adapter. This is the natural consumer of `delta` (`REV-004`). Keep the pure-DP path as the default and make the approximate-DP path opt-in. | open |
| `REV-007` | Representation | The trajectory group spends 40 numbers on per-cell presence/value pairs to describe what is qualitatively a 3-4 parameter curve. Dimension is the dominant cost term under any accountant. | Investigate releasing a low-dimensional shape (spline/basis coefficients, or a small set of quantiles) instead of a dense grid. Highest ceiling of anything on this list, and the most work. | open |
| `REV-008` | Post-processing | The decoders carry a lot of unlabelled rescue logic: `.ensure_grid_presence()` forces at least 3 cells regardless of the release (`R/generate.R:119`), variance floors `0.05^2`/`0.01^2`/`0.25` (`R/representation.R:658-660`), the `0.95` rate cap (`:686`). All are privacy-safe post-processing, but collectively they mean that at high noise the generator is mostly sampling from hardcoded priors — which undercuts "source-calibrated". | Promote these to named, documented constants. Add a diagnostic reporting what fraction of the output shape is attributable to the release versus to the fallbacks, so the honest claim is measurable rather than asserted. | open |

### P3 — Integrity and hygiene

| ID | Area | Issue | Suggested direction | Status |
|---|---|---|---|---|
| `REV-009` | Governance | `.public_input_manifest()` (`R/fit.R:103-124`) accepts six arguments, uses none of them, and returns a hardcoded 8-row constant that ships into the release ledger as `public_inputs`. It is formatted like an audit record and contains no information about the fit it documents. | Derive it from the actual inputs, or delete it. In a governance context a convincing-looking empty artifact is worse than no artifact. | open |
| `REV-010` | Leakage guard | The prohibited-payload check (`R/privacy.R:247-252`) is a nine-name denylist. It would not catch `subjects`, `piece`, or `rows`. | Invert to an allowlist over the expected released structure. A denylist is a smoke detector, not a guard, and it is currently presented as the latter. | open |
| `REV-011` | API | No `print.pmx_private_validation` or `print.pmx_backend_tests` method exists, so both exported classed returns dump as raw lists. Compare the five `print` methods that are registered in `NAMESPACE`. | Add both methods. | open |
| `REV-012` | Testing | `run_dp_backend_tests()` advertises a `"privacy_map"` entry in its `tests` vector (`R/privacy.R:135`) but `passed` only checks lengths and finiteness. The map check happens incidentally inside `release()`. | Either assert the privacy relation explicitly or stop naming a test that is not run. | open |
| `REV-013` | Side effects | `.opendp_backend()` calls `enable_features("contrib")` on every resolve (`R/privacy.R:5`), mutating global OpenDP state as a side effect of ordinary package use. | Do it once, guarded, and document it. | open |
| `REV-014` | Process | 8 commits total against a ~1,200-line uncommitted working tree that contains a whole feature (`R/properties.R`, the evaluation harness, `synpmx-epsilon-exploration.Rmd`). Large unreviewed surface for a package making formal privacy claims. | Land the working tree in reviewable commits before starting on this backlog. | open |
| `REV-015` | Docs | ~2,100 lines of vignette against ~3,500 lines of R, with `AGENTS.md` imposing a manual synchronization burden on every change. | Convert load-bearing prose claims into assertions or regression tests, per `AGENTS.md`'s own guidance, so they cannot silently go stale. | open |

---

## Suggested order of work

1. `REV-014` — land the working tree so the rest is reviewable.
2. `REV-001` — the decode threshold bug. Cheap, and it changes what the
   epsilon vignette shows, so everything downstream depends on it.
3. `REV-003` — the data-dependent `stop()` paths. Cheap, and it is a
   correctness issue in the package's central claim.
4. Re-run `scripts/evaluate_simulations.R` and the epsilon sweep. **Re-measure
   before planning further.** The remaining utility gap may be much smaller
   than it looks today.
5. `REV-002` + `REV-005` — planning helper and tighter sensitivity. Together
   these likely make `epsilon = 5` defensible at realistic N without touching
   the mechanism.
6. `REV-004` + `REV-006` — Gaussian/zCDP, if step 4 shows it is still needed.
7. `REV-007` — representation redesign, only if steps 5-6 are insufficient.
8. P3 items opportunistically.

---

## Part 3 — Is this feasible at N = 8, 20, 40, 100?

Measured after the `REV-001` fix, on `pmx_simulated_fixture(N)`, dose-relative
log `cp` endpoint, 8 repetitions per cell. `curve_err` is the median absolute
error of the decoded population curve on the working scale, **divided by the
true curve's own dynamic range**. So `curve_err >= 1` means the error is larger
than the entire signal being estimated — the output carries no information
about the source.

| N | epsilon 1 | epsilon 5 | epsilon 50 |
|---:|---:|---:|---:|
| 8 | 1.76 | 4.34 | 1.03 |
| 20 | 3.27 | 2.68 | 0.58 |
| 40 | 4.35 | 2.17 | 0.36 |
| 100 | 2.55 | 1.46 | 0.16 |
| 600 | 0.85 | 0.21 | 0.03 |
| 2000 | 0.32 | 0.09 | 0.01 |

At defensible epsilon (1-5), every cohort at or below N = 100 produces
`curve_err > 1`. The generated data at those sizes is a draw from the
package's post-processing priors with a privacy-noise decoration on top. It is
not source-calibrated in any meaningful sense, and no amount of decoder
engineering changes that.

### This is a lower bound, not an implementation defect

Known lower bounds for releasing `d` attribute means over `N` records under
`(epsilon, delta)`-DP put the per-coordinate error at roughly

```
error  >=  sqrt(d) / (epsilon * N)     (times a modest sqrt(log(1/delta)) factor)
```

Order-of-magnitude, with `d = 40` and a delta factor of ~5:

| N | floor at epsilon 5, d = 40 | floor at epsilon 5, d = 6 |
|---:|---:|---:|
| 8 | 0.79 | 0.31 |
| 20 | 0.32 | 0.12 |
| 40 | 0.16 | 0.06 |
| 100 | 0.06 | 0.02 |
| 600 | 0.01 | 0.004 |

Two things follow. First, **N < 20 is not achievable under any formal guarantee
at a defensible epsilon**, by any implementation — stop trying. Second, the gap
between today's measured 1.46 at N = 100 and the ~0.06 floor is roughly 25x,
which is approximately what `REV-005` + `REV-006` + `REV-007` claim to recover
together. The small-to-mid range is an engineering problem; the very small
range is not.

### Recommended tiering

| Cohort | Verdict | Recommended mode |
|---|---|---|
| N >= ~500 | Works today at epsilon 5, after `REV-001` | Current DP path unchanged |
| N ~100-500 | Not viable today; plausibly viable after P2 | Do `REV-005`, `REV-006`, `REV-007`, then re-measure |
| N ~20-100 | At the edge of theory even with a perfect implementation | Only reachable by cutting `d` to single digits (`REV-007`). Treat as research, not roadmap |
| N < 20 | Not achievable, provably | Do not fit. See below |

### What to do instead at small N

The package already contains almost everything needed for the honest answer,
but it is currently reachable only as a testing backdoor.

1. **Make "public design only" a first-class generation mode.** Users with 8 or
   40 patients overwhelmingly want structurally realistic PMX tables to
   exercise cleaning, joins, control-file plumbing, and censoring logic. None
   of that requires touching patient data. `pmx_public_design()`,
   `pmx_bounds()`, `pmx_endpoint()`, and `.resolved_regimen()`'s fallbacks
   already support generating entirely from declared public inputs. Expose a
   `synpmx_generate()` path that consumes **no** source data and **no** privacy
   budget, and therefore needs no DP claim. Today the only way to get this is
   `backend = "public"`, which is gated behind `public_source = TRUE` and is
   framed as a fixture hack rather than a supported answer.
2. **Separate the fitting cohort from the generated cohort in the
   documentation.** The privacy unit is one subject in the *fit*; `n_subjects`
   in `synpmx_generate()` is unrelated. A user can fit on 2,000 pooled subjects
   and generate 20. This already works and is the single most useful thing to
   tell a small-study user, but neither the README nor the demo vignette says
   it. Pooling across studies is the real answer to small N.
3. **Make release dimension adapt to (N, epsilon).** Rather than always
   releasing a 40-number curve, choose a tier from the planning calculation in
   `REV-002`: a large cohort gets the dense grid; a small cohort gets a handful
   of scalars (subject count, mean dose, mean observation count) with
   everything else supplied by the public design. A few scalars at N = 40 have
   error `~1 / (epsilon * 40)`, which is perfectly acceptable — it is only the
   40-dimensional curve that is hopeless.
4. **Refuse, or warn hard, rather than silently producing noise.** The existing
   `population$private_subject_count < 6` warning (`R/fit.R:246`) fires far too
   late and far too rarely. `REV-002`'s pre-flight check should make the
   infeasible case loud *before* budget is spent.

The net position worth documenting: this package is a **pooled-data tool**. It
turns a large confidential corpus into a reusable generator. It is not, and
cannot be, a way to make a 20-patient study shareable.

## Caveats on this review

- The `b/N` convergence claim and the `REV-001` diagnosis were verified
  empirically against `pmx_simulated_fixture()` at N = 60, 600, and 2000. The
  *proposed fix* for `REV-001` was not implemented or tested.
- The Gaussian/zCDP noise estimate in `REV-006` is an order-of-magnitude
  calculation, not a derived bound. Derive it properly before committing to it.
- Sensitivity bounds in the current code were checked for soundness and appear
  correct everywhere reviewed. `REV-005` is about tightness, not correctness.

---

## Addendum: findings added after the initial review

| ID | Area | Issue | Suggested direction | Status |
|---|---|---|---|---|
| `REV-016` | DP guarantee | The declared adjacency is "add-or-remove one complete subject" and `.bound_subject_contributions()` groups rows by the ID column, which is correct only if one person appears exactly once. Rollover and extension studies, crossovers pooled with parallel-group studies, and re-enrolled subjects all violate it. A person contributing `k` records receives roughly `k * epsilon` by group privacy, and nothing in the accounting, the ledger, or `validate_private_model()` reveals it. This matters more now that pooling is a recommended path. | Require an explicit assertion that IDs are unique persons, or accept a person-level grouping column and bound contributions on that rather than on the study subject ID. | open |
| `REV-017` | Public inputs | Realized trial-design quantities (dose levels reached, cohort sizes after expansion) are formally functions of the source dataset, so treating them as fixed public inputs is an assumption rather than a derivation. Practically this is low risk: such facts are usually already published in registries and presentations, the inference from a stopping point to any cause is many-to-one and confounded, and a committee's decision about a program is not one individual's record. The genuine per-subject case is separate and is a generation rule, not an accounting problem. | Record realized design in `proof_assumptions` alongside the other public-input assertions, with a source field, so a reviewer can see it. Separately, enforce that per-subject escalation and titration sequences are generated from the public rule and never copied from source subjects. See `vignettes/articles/data-elicitation.Rmd` section 3. | open |
| `REV-018` | Calibrated DP diagnostics | Correction diagnostics exposed the exact number of usable source subjects and continuous covariate means used the raw source count as a denominator. These were source-dependent output channels outside the released mechanism. | Use the private subject-count release for all released means and preflight accounting; keep exact counts internal only. | fixed |
| `REV-020` | Structural model API | `pmx_structural_model(rx = )` accepts and stores an `rxode2` model and errors when `rxode2` is missing, but nothing ever reads `model$rx`: `.pk_profile()`/`.pd_profile()` always use the built-in analytic solutions, and there is no test. A user supplying an ODE model silently gets the analytic 1-/2-compartment curve instead, with no warning. Found on 2026-07-23 while documenting the model-based engines. | Either wire `rx` into profile evaluation through `rxode2::rxSolve()` with a regression test against the analytic solution for a shared model, or reject `rx` with a "not yet supported" error until it is implemented. Until then, no vignette or roxygen text may claim ODE support. | open |
| `REV-019` | AVATAR schema and donors | AVATAR copied columns declared in `roles$exclude`, failed for factor-valued IDs, and matched dose events only by sign, allowing incompatible dose magnitudes to exchange trajectories. | Remove excluded columns before synthesis, extend factor levels for fresh IDs, and include rounded dose/rate magnitudes in event signatures. | fixed |
| `REV-021` | Validation coverage | `validate_pmx()` accepted a `CENS` flag with any `DV`: it checked the flag's value set, that event rows are uncensored, that boundaries are finite, and interval direction, but never that the flag agrees with the value. On the pre-fix AVATAR output (copied `CENS`, independently blended `DV`) it returned `valid = TRUE` on data with 36 rows flagged left-censored whose DV was above the limit and 12 unflagged rows below it. Found on 2026-07-23 after adding AVATAR censoring; the AVATAR defect that produced such data is fixed, but the validation blind spot is independent and applies to any caller-supplied dataset. | Fixed. Added a `cens_dv_coherence` check: per endpoint, a point left-censored DV (at the limit) may not exceed an uncensored DV, and a right-censored DV may not fall below one, or the flag and value were produced independently. **Scope:** the check assumes a single censoring threshold per endpoint, which is what AVATAR reconstructs; a genuine multi-assay-limit source where an uncensored value legitimately sits below a higher limit could be flagged. Revisit if a multi-limit dataset is added. Regression tests in `test-avatar-censoring.R`. | fixed |
| `REV-022` | AVATAR coherence | A baseline covariate and a longitudinal endpoint can be two views of the same quantity, and AVATAR keeps them consistent only by accident. In `CPIT565A1` a `B0` covariate holds baseline B-cell count while a separate DVID carries the B-cell kinetic trajectory. AVATAR blends covariates (`.synthesize_covariates()`) and endpoint trajectories (`.synthesize_trajectories()`) in independent passes, with independent donor weighting after the shared neighbour set and independent subject/residual perturbation, and copies the anchor's event template besides. So a synthetic subject's `B0` need not equal the baseline (t = 0, or the pre-dose value) of its own generated B-cell trajectory. The same applies to any baseline-covariate/endpoint pair — baseline biomarker vs its time course, screening creatinine vs a renal endpoint. Found 2026-07-24 on a real dataset. **Judged acceptable for now**: AVATAR targets workflow-development realism, not analysis-grade internal consistency, and the pieces are individually plausible. | No automatic detection: nothing in the schema says `B0` *is* the baseline of a given DVID; it is semantic knowledge only the modeller has. Options when it matters: (a) document the limitation so a user knows not to rely on covariate/endpoint-baseline consistency; (b) let a role declare that a covariate is an endpoint's baseline, and derive it from the generated trajectory instead of blending it independently; (c) offer a post-generation reconciliation that overwrites such a covariate with the generated baseline. Start with (a). | open |
