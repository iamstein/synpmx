# Feasibility assessment: what can be released, and from how many patients

Written 2026-07-22, after the `SIM-020` decode fix. Companion to
`design/REVIEW_BACKLOG.md` (defect-level findings) and
`design/PRIVACY_ARGUMENT.md` (the formal argument). This document asks a
scoping question rather than an implementation question:

> Is the current approach viable for the datasets we actually care about? If
> formal privacy cannot be delivered at those cohort sizes, what *can* be?

The short answer: privacy is entirely achievable at any cohort size — the thing
that becomes impossible at small N is **useful source-calibrated release**. That
is a statement about utility, not about privacy, and it constrains the avatar
approach in v1 exactly as much as the differential-privacy approach in v2. The
two designs fail in opposite, equally informative directions.

---

## 1. The constraint that binds both designs

A synthetic-data generator has exactly one job: carry information from the
source cohort into the generated cohort. Privacy engineering is the business of
limiting how much of that information is attributable to any single person.

At N = 12, those two goals are the same quantity. There is no "population
structure" in twelve patients that is separable from "these twelve patients."
The mean is a twelfth of each person. A PK curve shape is the average of twelve
individual curves, each contributing 8% of it. Any statistic rich enough to make
generated data look like the source is, at that cohort size, a lightly
obfuscated transcription of individuals.

This is not a limitation of Laplace noise, or of OpenDP, or of this codebase. It
is the reason both of the following are true:

- **v2 (DP) tells you honestly that it cannot help.** It adds noise calibrated
  to one person's worst-case influence. At N = 12 that influence is ~8% of every
  statistic, so the noise swamps the signal. The output is visibly useless.
- **v1 (avatar) did not tell you.** It produced attractive output at N = 12
  precisely because it carried individual-level information through. The
  plausibility *was* the leak.

The instinct that moved you off v1 was correct. But it is worth being precise
about why, because the same reasoning bounds what v2 can do.

---

## 2. What v1 actually did, and why it does not escape the constraint

From `21eb6e2` (`R/synthesis.R`, `R/profiles.R`), v1 was a `template` +
`avatar_blend` design:

1. Group subjects by an `.event_signature()` — a token built from EVID/CMT/DVID
   patterns, the *number of dose starts*, and the *rounded inter-dose gaps*.
2. Copy a source subject's event skeleton as the template for a generated
   subject, with `time_jitter` applied to unique times.
3. Project subject profiles with PCA, retaining `pca_variance` of variance.
4. Take the `k` nearest neighbors in that space and blend their covariates and
   DV trajectories with distance weights.
5. Add `subject_noise_sd`, `residual_noise_sd`, and AR(1) `residual_phi` noise.

Three properties of this design are load-bearing, and all three get *worse* as N
shrinks:

- **The event skeleton is copied, not generated.** A visit schedule is a
  quasi-identifier. Jittering unique times preserves the pattern — the number of
  doses, the gaps, the number and rough placement of samples. In a 12-patient
  study where one patient missed a visit, that patient's skeleton is unique and
  survives into the output.
- **k nearest neighbors is a large fraction of a small cohort.** With N = 12 and
  k = 5, each generated subject is a weighted average of 42% of the study. As N
  falls, neighbors are further away and the weights concentrate on the closest
  donor, so the blend degenerates toward copying.
- **The signature grouping partitions before blending.** A patient with an
  unusual regimen lands in a small — possibly singleton — signature group, and
  then gets "blended" with themselves.

None of this is fixable by tuning `k` or raising the noise. Push the noise high
enough to defeat a nearest-neighbor attack and you have destroyed the same
signal DP would have destroyed, without the accounting to prove it. **v1 is not
a way around the constraint in section 1; it is the same constraint with the
failure mode moved from "visibly useless" to "invisibly unsafe."**

That said, v1 is not worthless. It has a legitimate niche described in
section 6 — it just cannot be called private.

---

## 3. What differential privacy does and does not promise

Worth stating plainly, because "is true privacy possible" depends on what is
meant by privacy.

DP promises: **the released output is nearly the same whether or not you
participated.** Formally, for adjacent datasets differing in one subject, any
output's probability changes by at most a factor `exp(epsilon)`.

DP does **not** promise:

- that nothing can be inferred about you. If a study establishes that a
  biomarker predicts toxicity, that conclusion applies to you whether or not you
  enrolled. DP explicitly permits this and does not consider it a violation.
- anonymity, or any legal release authorization.
- protection against a bad public-input assertion. If a "public" bound was
  actually chosen by looking at confidential extrema, the guarantee is void, and
  no amount of noise detects that.

The practical reading of epsilon:

| epsilon | `exp(epsilon)` | Informal reading |
|---:|---:|---|
| 0.1 | 1.1 | Very strong |
| 1 | 2.7 | Strong; a common target |
| 5 | 148 | Weak worst case; often used anyway |
| 50 | 5.2e21 | Not a guarantee |
| 500 | — | Decorative |

This matters for reading the epsilon-exploration vignette: its "Large" column at
epsilon 500 is not a privacy setting, it is a demonstration of what the
generator does when the mechanism is effectively switched off.

**So: is true privacy impossible?** No. At N = 12 you can have perfect privacy
trivially — release nothing derived from the data. What is impossible is
*simultaneously* having privacy and source-calibrated utility. The impossibility
lives on the utility side of the tradeoff, and no framework moves it.

---

## 4. The measured feasibility frontier

From `pmx_simulated_fixture(N)`, dose-relative log `cp` endpoint, 8 repetitions,
after the `SIM-020` fix. Error is the median absolute deviation of the decoded
population curve divided by the true curve's own dynamic range, so **>= 1 means
the error exceeds the entire signal**.

| N | epsilon 1 | epsilon 5 | epsilon 50 |
|---:|---:|---:|---:|
| 8 | 1.76 | 4.34 | 1.03 |
| 20 | 3.27 | 2.68 | 0.58 |
| 40 | 4.35 | 2.17 | 0.36 |
| 100 | 2.55 | 1.46 | 0.16 |
| 600 | 0.85 | 0.21 | 0.03 |
| 2000 | 0.32 | 0.09 | 0.01 |

Error scales as `sensitivity / (epsilon * N)`. Extrapolating that law, with
today's 40-dimensional trajectory release:

| Target | Required N at epsilon 1 | Required N at epsilon 5 |
|---|---:|---:|
| error 0.10 (good) | ~6,000 | ~1,300 |
| error 0.25 (usable) | ~2,500 | ~500 |

Two conclusions:

- **At a defensible epsilon of 1, today's design needs thousands of subjects.**
  The epsilon 5 column is where the package currently looks acceptable, and
  epsilon 5 is a weak guarantee.
- **Every dataset in the demos and the epsilon sweep is far below the frontier**:
  `theo_md` 12, `nimoData` 12, `wbcSim` 45, `warfarin` 32, `mavoglurant` 120.
  Not one of them is in feasible territory at any defensible epsilon.

### Where the budget currently goes

The default allocation spends **half the epsilon on trajectory shape**
(`endpoints = 0.50` in the epsilon vignette; 0.40 in the review measurements).
That is the largest single line item, and it is spent learning the shape of a PK
or PD curve.

This is the most questionable allocation in the design, because **curve shape is
usually the least confidential thing in the dataset.** That theophylline
concentrations rise to a peak in one to two hours and decline log-linearly is
public knowledge, published, and in some cases already encoded in
`nlmixr2data`. The package is spending most of its privacy budget rediscovering
literature from twelve patients.

---

## 5. What can actually be achieved

Distinguish three things that get conflated as "privacy":

1. **Formal guarantee** — a worst-case bound on one person's influence (DP).
   Needs N. Provable.
2. **Empirical risk assessment** — measured attack success: membership
   inference, attribute inference, nearest-neighbor distance ratios. No
   guarantee, but evidence. This is what avatar/synthpop-style tools and most
   regulatory de-identification frameworks actually rely on.
3. **Governance** — controlled environment, contractual limits, access logging.
   Privacy from the process rather than from the artifact.

Most real pharmaceutical workflows run on (2) + (3). This package attempts (1),
which is a stronger and rarer claim, and the honest cost of that claim is the
cohort-size frontier in section 4.

Given all of the above, here is the achievable menu, ordered by how much source
information each carries:

### Tier A — Public-design generation (no source data)
Generate entirely from declared public inputs: schema, roles, endpoint clocks,
bounds, contribution limits, nominal regimen. **Zero privacy cost, no DP claim
needed, works at N = 0.** Fully sufficient for the package's stated purpose —
exercising cleaning, joins, reshaping, control-file plumbing, censoring
conventions, repeated-dose code. Already implemented; currently reachable only
via `backend = "public"` behind a `public_source = TRUE` gate, framed as a
fixture hack rather than the answer.

### Tier B — Literature-informed shape plus a few private scalars
Take curve shape from a published popPK/popPD model. Spend the entire epsilon
budget on the handful of quantities that are genuinely study-specific and not
public: cohort size, mean dose, doses per subject, observations per subject,
censoring rate. That is `d` around 5-8 instead of 40, and it changes the
arithmetic completely:

| Design | `d` | Error at N = 40, epsilon 1 |
|---|---:|---:|
| Today (dense grid, Laplace) | 40 | ~1.0 (useless) |
| Scalars only, Laplace | 6 | ~0.15 |
| Scalars only, Gaussian/zCDP | 6 | ~0.06 (usable) |

**This is the design that makes N = 40 feasible at a defensible epsilon**, and
it is a much better fit for pharmacometrics than the current one, because the
field already has strong public structural priors. It is not on the backlog yet;
it should be.

### Tier C — Full DP population generator (current design, improved)
Viable at N in the high hundreds after `REV-005`/`REV-006`/`REV-007`; needs
thousands at epsilon 1 without them. Appropriate for **pooled** corpora: many
studies of a compound, a legacy database, a consortium dataset.

### Tier D — Avatar/blend synthesis (v1)
No formal guarantee. Defensible only inside Tier-3 governance, with Tier-2
empirical attack testing, and never as an external release. Its legitimate niche
is *inside* a trusted environment where the alternative is passing the real data
around — but that is a different product from the one this package claims to be.

---

## 6. Recommendation

> **Superseded in part by section 8.** This section was written before the
> target cohort sizes were fixed. The scope decision of 2026-07-22 prioritizes
> small trials, which rules out point 2 below as the primary direction and
> promotes point 3 (Tier B) from research direction to main architecture. Points
> 1 and 4 stand unchanged.

**Reposition the package around two supported modes, and be explicit that the
middle ground does not exist.**

1. **Make Tier A a first-class, prominently documented mode.** For studies below
   roughly 100 patients, this is the answer, and it is a *good* answer: it
   delivers exactly what the README says the package is for. It needs a public
   entry point that reads no source data and spends no budget — not a fixture
   backdoor.
2. **Scope the DP path explicitly to pooled data.** State a minimum viable
   cohort in the documentation. On today's implementation that is roughly N >=
   500 at epsilon 5, or N >= 2,500 at epsilon 1. Implement `REV-002`'s
   pre-flight check so an infeasible configuration is refused *before* budget is
   spent, not discovered afterward in a plot.
3. **Investigate Tier B as the main research direction.** It is the only path
   that makes small-cohort DP genuinely feasible, and it exploits the one
   structural advantage this domain has over generic tabular synthesis: strong,
   public, well-validated priors for curve shape. Concretely: stop spending half
   the budget learning what a PK curve looks like.
4. **Do not revive v1 as a privacy tool.** It may be worth keeping as an
   explicitly-labelled in-environment convenience, but the label matters more
   than the code.

The one-sentence scope statement worth putting in the README:

> This package turns a **large pooled** confidential corpus into a reusable
> generator with a formal privacy guarantee. It is not a way to make a small
> study shareable; for that, generate from public design inputs and spend no
> privacy budget at all.

---

## 7. What would change this assessment

- **Tier B prototyping.** The `d = 6` numbers in section 5 are arithmetic from
  the error law, not measurements. Build it and measure before trusting them.
- **A tighter sensitivity analysis** (`REV-005`). If the true per-subject L1 is
  much smaller than `ncol`, every N threshold here drops proportionally.
- **Empirical attack testing on v1** would replace my structural argument in
  section 2 with evidence. A nearest-neighbor distance-ratio test and a
  membership-inference test against the `21eb6e2` code would settle how bad it
  actually was, which is worth knowing before reusing any of it.
- **A real target dataset.** The frontier in section 4 assumes the fixture's
  release dimension. Section 8 maps it onto typical development-program cohort
  sizes; confirm which of those rows are the actual targets.

---

## 8. Assessment against the target cohort sizes

Target sizes in scope (2026-07-22): **6, 20, 60, 100, 500, 1000, 10000**.
Scope preference: prioritize the small end.

That preference rules out the section 6 recommendation, which pointed at pooled
corpora. It forces the Tier B design, because Tier B is the only column in which
the small sizes appear at all.

### Today's dense-grid design against these sizes

**Measured**, not extrapolated. Same protocol as section 4; 8 repetitions below
N = 500, 5 at N = 500-1000, 3 at N = 10000.

| N | epsilon 1 | epsilon 5 | Verdict |
|---:|---:|---:|---|
| 6 | 2.43 | 2.71 | Tier A only |
| 20 | 3.74 | 3.42 | Tier A only |
| 60 | 3.49 | 1.46 | Tier A only |
| 100 | 2.69 | 1.06 | no |
| 500 | 1.08 | 0.26 | marginal, and only at a weak epsilon |
| 1000 | 0.55 | 0.09 | usable at a weak epsilon |
| 10000 | 0.08 | 0.02 | good; the first size that works at epsilon 1 |

Values below N = 100 are non-monotone in both N and epsilon. That is expected
and is itself the finding: in that regime the output is pure noise, so the
metric is measuring the normalization rather than any signal.

Only the bottom two rows of the target list are served, and N = 1000 only at an
epsilon that is not a real guarantee. **The current architecture does not serve
the stated scope.**

### The Tier B arithmetic

Release `d` bounded per-subject scalars instead of a dense grid. Each subject's
value is clipped to a **public prior range**, and the released quantity is a
mean. Under Laplace with basic composition the error on each released mean, as a
fraction of the prior range, is

$$f \;=\; \frac{d}{\varepsilon N}.$$

With a draft `d = 6` (cohort size, CL, t-half, PD baseline, PD magnitude, PD
onset rate):

| N | eps 1 | eps 3 | eps 5 |
|---:|---:|---:|---:|
| 6 | 1.00 | 0.33 | 0.20 |
| 20 | 0.30 | 0.10 | 0.06 |
| 60 | 0.10 | 0.033 | 0.020 |
| 100 | 0.060 | 0.020 | 0.012 |
| 500 | 0.012 | 0.004 | 0.002 |
| 1000 | 0.006 | 0.002 | 0.001 |
| 10000 | 0.0006 | 0.0002 | 0.0001 |

`f` is a fraction of the *prior range*, so it converts to a parameter error only
once that range is fixed. For a two-fold public prior on CL (0.69 log units),
relative error is approximately `exp(0.69 f) - 1`:

| `f` | Relative error on CL |
|---:|---:|
| 0.30 | 23% |
| 0.20 | 15% |
| 0.10 | 7.1% |
| 0.06 | 4.2% |
| 0.033 | 2.3% |
| 0.012 | 0.8% |

Taking `f <= 0.10` as the usability bar gives the frontier `epsilon * N >= 60`:

| N | Minimum epsilon for `f <= 0.10` | Verdict |
|---:|---:|---|
| 6 | 10 | **No.** Not a guarantee. Tier A |
| 20 | 3 | Marginal; defensible only with governance |
| 60 | 1 | **Works at a real epsilon** |
| 100 | 0.6 | Works at a strong epsilon |
| 500 | 0.12 | Very strong |
| 1000 | 0.06 | Very strong |
| 10000 | 0.006 | Essentially free |

**This is the innovation.** Today epsilon 1 needs roughly 6,000 subjects. Under
Tier B it needs about 60. That is a hundredfold reduction in the cohort size at
which a defensible guarantee becomes possible, and it moves Phase 1 from
impossible into range.

### Measured confirmation

The arithmetic above was validated against real OpenDP releases. Per-subject CL
was computed by non-compartmental analysis (`CL = Dose / AUC`, trapezoidal over
the first occasion) from `pmx_simulated_fixture`, clipped to a public range on
the log scale, and released as a mean with `d = 5` scalars sharing the budget.
200 replicate releases per cell. The table reports **median fold-error on CL**;
1.5 means "within 50%".

| N | 100-fold prior, eps 0.5 | eps 1 | eps 5 | 5-fold prior, eps 0.5 | eps 1 | eps 5 |
|---:|---:|---:|---:|---:|---:|---:|
| 6 | 9.90 | 9.60 | 1.71 | 2.18 | 2.15 | 1.27 |
| 20 | 6.05 | 2.44 | 1.26 | 1.98 | 1.52 | 1.07 |
| 60 | 2.00 | 1.46 | 1.06 | 1.23 | 1.10 | 1.03 |
| 100 | 1.48 | 1.21 | 1.04 | 1.13 | 1.07 | 1.02 |
| 300 | 1.14 | 1.07 | 1.01 | 1.05 | 1.03 | 1.00 |
| 500 | 1.09 | 1.04 | 1.01 | 1.03 | 1.01 | 1.00 |
| 1000 | 1.03 | 1.02 | 1.00 | 1.01 | 1.01 | 1.00 |

**The error law is confirmed.** For N = 20, epsilon 1, 5-fold prior:
`f = 5/(1 * 20) = 0.25`, and a 5-fold prior spans `log(5) = 1.61` units, so the
predicted error is `exp(0.25 * 1.61) = 1.49`-fold. Measured: 1.52.

The formula is mildly conservative at wide priors, because it uses the Laplace
*scale* while the table reports the *median* absolute error, which is
`ln(2) = 0.69` times the scale. Treat `f = d/(epsilon N)` as a safe planning
bound rather than a point prediction.

### The two prior columns are the argument for range-finding

Comparing the two halves of that table at epsilon 1: N = 20 improves from
2.44-fold to 1.52-fold, and N = 60 from 1.46 to 1.10, purely by narrowing the
clipping range from 100-fold to 5-fold. That is roughly a 3x effective privacy
gain — the same accuracy at a third of the epsilon.

A coarse private histogram can buy that narrowing for about 20% of the budget,
because **a histogram has L1 sensitivity 1 regardless of bin count**: each
subject falls in exactly one bin. Locating the data is nearly free; estimating
its mean precisely is what costs `d/(epsilon N)`. Spending a fifth of the budget
to triple the value of the remainder is clearly correct.

This also answers the objection that a new drug has no published model. See
`design/PROTOTYPE_SPEC.md` section 5: "public" in the DP sense means
*independent of this dataset*, not *published*. Preclinical allometric scaling,
the protocol, the assay validation report, and physiological ceilings are all
available before the data exist, and a preclinical human PK prediction is
typically within 3-5 fold — already enough for stage 2 even before
range-finding.

### Measured against the implementation

The tables above measure the *arithmetic*. The v3 generator now exists, so the
frontier can be measured end to end: simulate a cohort whose true clearance is
25 against a public prediction of 10, calibrate, and compare. Forty replicate
OpenDP releases per cell, `d = 2`, an 8-fold correction-factor prior.

Three error sources are separated. `dp` is the noisy release against a noiseless
release of the same cohort, isolating the privacy cost. `bias` is the noiseless
release against the population truth. `total` is what a user actually gets.

| N | eps | `f` | predicted | dp | bias | total |
|---:|---:|---:|---:|---:|---:|---:|
| 6 | 0.25 | 1.33 | 40.3 | 2.40 | 1.12 | 2.32 |
| 20 | 0.25 | 0.40 | 3.03 | 1.71 | 1.07 | 1.60 |
| 60 | 0.25 | 0.13 | 1.45 | 1.40 | 1.08 | 1.40 |
| 100 | 0.25 | 0.08 | 1.25 | 1.21 | 1.07 | 1.23 |
| 300 | 0.25 | 0.03 | 1.08 | 1.07 | 1.07 | 1.08 |
| 20 | 0.5 | 0.20 | 1.74 | 1.63 | 1.07 | 1.60 |
| 60 | 0.5 | 0.07 | 1.20 | 1.17 | 1.08 | 1.21 |
| 300 | 0.5 | 0.01 | 1.04 | 1.04 | 1.07 | 1.08 |
| 20 | 1 | 0.10 | 1.32 | 1.24 | 1.07 | 1.26 |
| 60 | 1 | 0.03 | 1.10 | 1.12 | 1.08 | 1.13 |
| 300 | 1 | 0.01 | 1.02 | 1.02 | 1.07 | 1.06 |
| 20 | 2 | 0.05 | 1.15 | 1.16 | 1.07 | 1.19 |
| 300 | 2 | 0.00 | 1.01 | 1.01 | 1.07 | 1.07 |

**The error law holds for `f` below about 0.25**, typically within 10%. Above
that it is increasingly pessimistic, and at `f = 1.33` it predicts 40-fold
against a measured 2.4-fold.

The reason is worth recording: **clipping bounds the damage.** A release cannot
land outside the prior range, so however large `f` grows the error saturates
near the prior's half-width rather than diverging. `pmx_preflight()` now caps
its reported fold-error accordingly. The `f >= 1` verdict is unaffected and
remains correct — at that point the release is a draw from the prior and
conveys nothing — but the number attached to it was nonsense.

### A systematic floor the error law does not describe

`bias` sits near 1.07 and **does not shrink with N or epsilon**. It is not
sampling error; it is estimator bias, and it sets a floor on total accuracy that
no privacy budget can buy through.

Measuring it found a real defect, now fixed. The generator drew between-subject
variability centred on the *arithmetic* mean while the released correction, being
a mean on the log scale, targets the *geometric* mean. The two differ by
`exp(sigma^2 / 2)`, which for 30% CV is about 4.4%, and nothing about that gap
shrinks with more subjects or more budget. `typical` is now the median
throughout, which is also the usual population-PK convention. That change took
the floor from 1.12 to 1.07.

The residual ~7% is AUC truncation: the correction compares trapezoidal AUC over
a finite window, and a prediction 2.5-fold wrong captures a different fraction
of the true curve than the observation does. It shrinks as the prediction
improves, and at 7% it is far inside the accuracy bar. It is recorded rather
than chased.

### Structural misspecification does not break exposure

The measurements above fit the same model family that generated the fixture. To
test the assumption directly, a two-compartment truth was fitted with a
one-compartment public model at epsilon 0.5:

| N | Recovered CL (true 25) | Generated/source AUC | Generated/source Cmax |
|---:|---:|---:|---:|
| 60 | 29.5 | 0.99 | 0.89 |
| 300 | 22.9-25.5 | 1.02-1.10 | 0.93 |

**Getting the compartment count wrong does not break exposure magnitude.** The
correction is an AUC ratio, so it aligns total exposure regardless of how the
curve gets there; Cmax lands within about 10%. The generated profile is
mono-exponential where the truth is biphasic, which is visible in a plot and
irrelevant to the accuracy bar. This is the design working as intended: shape is
asserted, magnitude is calibrated.

### PD calibrates cleanly with an exposure-independent endpoint

PD is a simple time course — constant, linear, or exponential approach to a
plateau — with no exposure coupling. That is a modeling choice, but it is also
what makes the correction well conditioned. The released quantity is a **level
correction**, the ratio of mean observed to mean predicted response, where both
terms are the response itself rather than a small deviation from a large
baseline.

Measured on an exponential-decay endpoint, true baseline 100 against a predicted
40, so a true correction of 2.5x:

| Residual error | Recovered correction | Recovered baseline |
|---|---:|---:|
| 0% | 2.51 | 100.4 |
| **15%** | **2.53** | **101.3** |
| 30% | 2.61 | 104.4 |

The estimator is essentially unbiased. At realistic budgets, with 15% residual
and a 100-fold PD prior at epsilon 0.5, the recovered baseline is 155 at N = 20,
106 at N = 60, and 102 at N = 300 — the error at small N being privacy noise
rather than bias, exactly as the error law predicts. A tighter baseline prior,
usually available from healthy-volunteer literature, would shrink all three.

An earlier design used an exposure-driven Emax endpoint. Its correction had to
estimate a deviation from baseline and was biased low by about a third under the
same residual error, because a geometric mean of noisy per-subject ratios sits
below the ratio of their means. Rather than chase that bias with a better
estimator, the endpoint itself was simplified. Exposure-driven PD has been
removed entirely.

### Caveat on scope of the measurement### Caveat on scope of the measurement

Only the PK correction was measured. The PD correction follows the same error
law but is unconfirmed, and should be expected to behave worse: its prior is
wider and its per-subject estimate noisier. Measure it before relying on the
parameter vector as a whole.

Misspecification is now covered for PK, above, but only for a two-compartment
truth fitted by a one-compartment model. Absorption misspecification, nonlinear
elimination, and time-varying clearance are untested.

### Why the prior range is the real lever

The reason this works is worth stating precisely, because it determines where
the effort should go.

The DP release is not "the PK parameters". It is a **correction to a public
prior**. Noise matters only relative to the width of that prior. A prior saying
CL is somewhere in 0.1-100 L/h needs far more budget to locate the truth than
one saying 3-6 L/h. Halving the prior width halves the required epsilon at fixed
N.

Pharmacometrics is unusually well placed to exploit this. The field already
builds informative structural priors as a matter of course — published popPK
models, allometric scaling, known elimination pathways, established PD
mechanisms. **The tightest defensible public prior is worth more than any
mechanism improvement**, and tightening priors is cheaper than any of
`REV-005`/`REV-006`/`REV-007`.

The corresponding hazard is exact: the prior now does most of the work, so the
temptation to tune it against the data is much stronger, and the consequence of
doing so is much worse. A prior fitted to the confidential data voids the
guarantee entirely and nothing downstream detects it. This needs the same
provenance discipline as `pmx_bounds()`, applied harder.

### Why non-compartmental analysis, and not popPK fitting

The estimator choice is a privacy constraint, not a modeling preference.

- **NCA is DP-compatible.** Trapezoidal AUC and a terminal slope are computed
  from one subject's own rows. One subject's data influences only that
  subject's value, so clipping to a public range bounds the sensitivity
  directly.
- **NLME/popPK fitting is not.** Shrinkage couples subjects: every individual
  post-hoc estimate depends on the population fit, which depends on everyone.
  One subject perturbs all N estimates by an amount with no simple bound.
  Making that private requires DP-SGD-style machinery, far more budget, and a
  much harder proof.

This constraint aligns favorably with the stated scope. NCA needs rich sampling,
and rich sampling is exactly what small early-phase studies have. Sampling
richness and cohort size are inversely related across development, so the
small-N target and the DP-compatible estimator want the same designs.

### Mechanism note: do not use Gaussian here

`REV-006` recommends Gaussian noise under zCDP. That advice is **dimension
dependent and does not apply to Tier B**. At `d = 6` and `epsilon = 1`, Laplace
gives a noise scale of 6 on the sum, while a Gaussian at `delta = 1e-6` gives
roughly 13. Gaussian wins only once `d` is large enough for the `sqrt(d)`
sensitivity to beat the `sqrt(2 ln(1.25/delta))` constant, which is well above
`d = 6`. Keep pure-DP Laplace for the low-dimensional path.

### Residual risks in the Tier B design

1. **Prior provenance**, as above. The dominant threat to the claim.
2. **Confidently wrong output.** Today's noise looks like noise. A tight but
   incorrect prior produces clean, plausible, wrong data. Failure becomes
   invisible again — the same property that made v1 dangerous, arriving through
   a different door.
3. **Too clean to be useful.** Mock data generated from a smooth structural
   model lacks the outliers, BLQ runs, missed doses, and protocol deviations
   that actually break analysis pipelines. This directly undercuts the stated
   purpose. Mitigate by injecting messiness from public and protocol knowledge
   rather than learning it — dropout rates, BLQ rules, and visit windows are
   design facts, not patient data.
4. **Linear PK is an assumption, not a fact.** Dose-normalization is invalid
   under TMDD or saturable elimination. Superseded by the structural-model
   design in `design/PROTOTYPE_SPEC.md` section 6, where the model states its
   own assumptions and a nonlinear structure can simply be supplied.

### Very small cohorts

The DP guarantee at a given epsilon is independent of `N` — that is the
definition, and resisting an attacker who already knows five of six subjects is
precisely what DP provides and k-anonymity does not. Cohort size affects only
how much signal survives the noise.

Under the correction-factor design (`d = 3`, 8-fold prior), N = 6 at epsilon 1
gives `f = 0.5`, halving the prior uncertainty and taking 8-fold to about
2.8-fold. Against the accuracy bar in `design/PROTOTYPE_SPEC.md` section 1 that
is useful. The margin is thin, though — at epsilon 0.5, `f = 1.0` and the
release conveys nothing beyond the prior — so run the pre-flight check rather
than assuming. Governance may separately decline any release derived from six
patients, which is a reasonable policy position but not a mathematical one.

### On novelty

The combination — public structural priors from the pharmacometric literature,
a handful of DP-released scalars, and per-subject NCA as the estimator — is a
credible methodological contribution and I am not aware of published work doing
exactly this. But DP under informative priors is not itself new, and I cannot
verify novelty from inside this repository. Before making a novelty claim in a
paper or a package abstract, check the literature for DP with non-compartmental
analysis, DP applied to popPK, and DP synthetic data under strong structural
priors. The engineering position stands regardless of what that search finds.

### A pooling hazard that matters more now

The recommendation to pool makes an assumption the current code does not check.
The declared adjacency is **add-or-remove one complete subject**, and
`.bound_subject_contributions()` groups rows by the ID column. That is correct
only when one person appears exactly once in the pooled dataset.

Common situations that violate it:

- **Rollover and extension studies.** A patient completing a Phase 2 study and
  enrolling in its open-label extension appears twice, usually under different
  subject numbers.
- **Crossover designs** pooled with the parallel-group studies they accompany.
- **Re-screened or re-enrolled subjects** within a program.
- **Any pooling that concatenates studies without a subject-level crosswalk.**

If one person contributes `k` records, the guarantee they actually receive
degrades to roughly `k * epsilon` by group privacy. At epsilon 5 with a single
rollover, that person has epsilon 10 — and nothing in the accounting, the
ledger, or `validate_private_model()` would reveal it.

This is not hypothetical for the recommended use case: pooled submission
datasets frequently contain rollover subjects. Before pooling is recommended in
user-facing documentation, the package should either require an assertion that
IDs are unique persons, or accept a person-level grouping column and bound
contributions on that instead of on the study subject ID. Tracked as a
prerequisite for the section 6 recommendation.
