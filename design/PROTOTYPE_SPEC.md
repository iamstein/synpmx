# synpmx implementation specification

Working specification for the package. Read `design/TODO.md` for the current
task queue and `design/FEASIBILITY.md` for the evidence behind the scope
decisions here.

**Document order:** objective, scope, architecture, then the detailed contract,
then version history newest-first at the bottom.

---

# 0. Current method (Version 4): AVATAR is the default

**Read this first.** As of Version 4 (2026-07-22), the primary and default
method is **AVATAR-style blending**, exposed as `synthesize_pmx()`. It resamples
and blends whole source subject trajectories, works at any cohort size, needs no
elicitation, and makes **no formal privacy guarantee**. It is the right tool for
**synthetic data that stays inside a trusted computing environment**, which is the
common case, and it is the trajectory-level analogue of what Novartis's
`synadam` already does column by column. See `design/METHOD_DISCUSSION.md`.

The differentially private engines described in most of this document —
`fit_calibrated_pmx()` (structural correction, Version 3) and `fit_private_pmx()`
(dense grid, Version 2) — are **retained as alternatives** for the case that
actually needs a formal guarantee: when the generated data **crosses a trust
boundary** (shared externally, published, or moved to a less-controlled system).
They are not removed and not deprecated; they are the right tool for a different
question.

The rest of this specification (sections 1-10) documents the DP engines in
detail, because they carry the formal obligations that need specifying. AVATAR's
method is documented in the method vignette and `design/METHOD_DISCUSSION.md`.
The decision rule between them is a single question: **does the output cross a
trust boundary?** No → AVATAR. Yes → a DP engine.

---

# 1. Objective

Build a generator of **synthetic pharmacometric data for model-workflow
exploration**. It reads a real dataset inside a restricted environment, learns a
small number of differentially private summaries, and generates structurally
coherent synthetic PMX event datasets outside that environment.

The generated data must be good enough to exercise:

- data cleaning, joins, reshaping, filtering, and derivation code;
- exploratory visualizations;
- control-file and model-run plumbing;
- PK, PD, biomarker, and repeated-dose workflows;
- preliminary model-code debugging before returning to the restricted
  environment.

It must capture:

- the correct order of magnitude for each endpoint;
- whether an endpoint rises and falls after each dose;
- whether an endpoint changes over a longer study-level timescale;
- **monotone dose-exposure behavior** — higher doses give higher exposures;
- approximate dosing and sampling-time scales;
- coarse subject-to-subject and observation-level variability;
- endpoint, censoring, and event-table conventions.

It explicitly does **not** preserve exact distributions, covariance or
covariate-response relationships, rare individual trajectories, parameter
estimates, scientific conclusions, inferential validity, or model-selection
results.

The target is "vaguely right, structurally exact". A user should be able to
develop and debug an analysis pipeline against this data and have it run
unchanged against the real data. They should not be able to draw a scientific
conclusion from it. The deliberate lack of statistical fidelity is an advantage
for privacy, not a limitation to be engineered away.

## The accuracy bar, stated explicitly

This is a design input, not a consequence, and it drives every other decision.

**Structure must be right. Numbers need only be in the right region.**

| Aspect | Requirement |
|---|---|
| Sampling times | **Exact.** From the protocol. Nominal times, visit windows, and the relationship between doses and samples must be correct |
| Event-table structure | **Exact.** EVID, CMT, DVID, censoring conventions, ordering, schema |
| Dose-exposure ordering | **Correct.** Higher dose gives higher exposure |
| Curve shape | **Approximately correct.** Rises and falls when it should; nadir and recovery where they should be |
| Parameter magnitudes | **Within a few fold is acceptable.** A 2-3 fold error on CL is fine |
| Distributions, covariance, parameter estimates | **Not required at all** |

The consequence is important and easy to miss: **the quantities that matter most
cost no privacy budget.** Sampling times come from the protocol. Structure comes
from the schema and the event conventions. Shape comes from the structural
model. All are public inputs.

Privacy budget is spent only on the part the user cares least about — roughly
where the exposure magnitude sits. That is a fortunate alignment and the design
should exploit it: **spend the surplus on a smaller epsilon, not on accuracy
nobody needs.**

---

# 2. Scope

## Target cohort sizes

**6, 20, 60, 100, 500, 1000, 10000 subjects, with the small end prioritized.**

Small cohorts are the hard case and the main design driver. A method that only
works at N = 1000 does not serve this package's users.

## What is achievable

Measured, not assumed. See `design/FEASIBILITY.md` sections 4 and 8.

The controlling relationship is that the error on any released per-subject mean,
expressed as a fraction of the range it was clipped to, is approximately

$$f \;=\; \frac{d}{\varepsilon N}$$

for `d` released scalars under Laplace with basic composition. Three levers:
fewer released quantities, a narrower clipping range, more subjects. Only the
first two are under our control at a fixed study size, and **the clipping range
is the cheaper of the two**.

## Background: `d`, `f`, and the error law

Full derivation with worked examples in `design/PRIVACY_BACKGROUND.md`. The
short version, because these two symbols appear throughout:

**`d` is how many separate quantities you release.** Privacy budget is spent,
not shared, so releasing `d` things splits epsilon `d` ways under basic
composition. Each release gets `epsilon/d`, and since a clipped per-subject
value has sensitivity 1, the Laplace noise on each released *sum* has scale
`b = d/epsilon`. **Releasing twice as many numbers doubles the noise on every
one of them.**

We want means, not sums, so we divide by the released count — post-processing,
which is free. That divides the error by `N`:

$$f \;=\; \frac{b}{N} \;=\; \frac{d}{\varepsilon N}$$

**`f` is the error expressed as a fraction of the prior range's width.** It is
dimensionless because every subject's value was scaled into `[0, 1]` by clipping
to that range, so 0 is the bottom of the prior and 1 is the top.

That reading is what makes `f` useful. `f = 1` means the noise is as wide as the
whole prior, so the release said nothing you had not already assumed. `f = 0.25`
means it quartered your uncertainty.

To get a fold-error, multiply by the prior's log span `S` and exponentiate: a
prior spanning a factor of `k` has `S = ln k`, and the fold-error is
`exp(f * S)`. Note that `f` says nothing about `S` — a wide prior and a narrow
one at identical `f` give very different answers, which is exactly why the
correction-factor parameterization in section 6 is worth so much.

## The decision rule: does the release beat the prior?

The right question is **not** "is the error small enough?" It is "is the
released value more informative than the prior I already had?" Below a
threshold, the release adds nothing over the prior and paying budget for it is
strictly worse than not paying — Prior mode gives the same output for free.

With `d` releases, epsilon `e`, `N` subjects and a prior spanning `S` log units,
the released mean has error about `f * S` where `f = d/(e N)`. So `f` *is* the
fraction of the prior's width that survives as noise:

| `f` | Interpretation | Verdict |
|---:|---|---|
| >= 1.0 | Noise equals or exceeds the prior width | **Worthless.** Use Prior mode |
| ~0.5 | Roughly halves the uncertainty | Marginal but real |
| ~0.25 | Quarters it | Clearly worthwhile |
| <= 0.1 | Prior contributes almost nothing | Consider lowering epsilon instead |

**When `f` falls below about 0.1, do not bank the accuracy — reduce epsilon.**
Given the accuracy bar above, a smaller guarantee is worth more than precision
nobody needs.

## Very small cohorts

The differential-privacy guarantee does not weaken as `N` falls. At a given
epsilon it is exactly as strong at N = 6 as at N = 6000 — that is the
definition, and resisting an attacker who already knows five of six subjects is
precisely what DP provides and k-anonymity does not. Cohort size affects only
how much signal survives the noise, which is a utility question.

At N = 6 with `d = 3` and a correction-factor prior, `f = 0.5` at epsilon 1: the
release roughly halves the prior uncertainty, taking an 8-fold prior to about
2.8-fold. Against the accuracy bar above, that is useful.

Two cautions apply. The margin is thin — at epsilon 0.5, `f = 1.0` and the
release conveys nothing — so run the pre-flight check rather than assuming. And
governance may decline to authorize any release derived from six patients as a
policy matter. That is a reasonable position, but it is a policy judgment rather
than a mathematical one, and this specification should not dress it up as the
latter.

---

# 3. Overall shape of what will be developed

Three generation modes sharing one API, one validator, and one event-table
constructor. They form an increasing-data-use axis:

| Mode | Uses the data | One-line description |
|---|---|---|
| **Prior** | Not at all | Generate entirely from declared design and priors. No budget, no DP claim |
| **Calibrated** | A few numbers | A structural model corrected by a handful of privately released scalars |
| **Empirical** | The shape itself | Learn trajectories directly from a dense private grid (Version 2) |

## Choosing a mode

The two questions that decide it: **how many patients**, and **how many numbers
are you willing to learn from them**.

| | Prior | Calibrated | Empirical |
|---|---|---|---|
| **Parameters learned from data (`d`)** | **0** | **2-3** | **~100** |
| Works from | any N, including no dataset at all | N >= ~6 | N >= ~1000 |
| Privacy budget | none | small | large |
| DP claim | none needed, none made | yes | yes |
| Where magnitudes come from | the prior alone | prior corrected by the data | the data |
| Where structure comes from | protocol + structural model | protocol + structural model | protocol + released grid |

`d` is the whole story. Error goes as `d/(epsilon N)`, so learning 100 numbers
instead of 3 costs roughly 30x the noise or 30x the epsilon. The Empirical mode
buys empirical trajectory shape with that; the Calibrated mode decides the shape
is not worth it and takes it from a structural model instead.

### Quick guide by cohort size

| Patients | Mode | `d` | Starting epsilon |
|---|---|---:|---:|
| No data, or 1-5 | Prior | 0 | n/a |
| 6-20 | Calibrated, or Prior if the margin is too thin | 3 | 0.5-2 |
| 20-100 | Calibrated | 3 | 0.2-0.5 |
| 100-1000 | Calibrated | 3 | 0.05-0.2 |
| 1000+ | Calibrated | 3 | 0.05 or lower |
| 1000+ **and** empirical shape is specifically needed | Empirical | ~100 | 5 or higher |

A useful rule of thumb for the Calibrated mode is **`epsilon ~ 10/N`**, which
holds `f` near 0.3 — roughly a 2-fold spread on exposure against an 8-fold
prior. Treat it as a starting point, then check it with the pre-flight rule in
section 2 rather than adopting it blindly. Section 4 works through three cohort
sizes in detail.

Note what the last row costs: the Empirical mode needs an epsilon two orders of
magnitude weaker than the Calibrated mode at the same cohort size. That is the
price of learning shape from data instead of asserting it.


```
                 public inputs                     confidential data
                       |                                   |
                       v                                   v
  +------------------------------+     +------------------------------+
  |  PRIOR                       |     |  optional range-finding      |
  |  schema, roles, protocol,    |     |  coarse histogram,           |
  |  structural model, priors,   |     |  sensitivity 1, cheap        |
  |  assay, dose levels          |     +--------------+---------------+
  |  no budget, no DP claim      |                    |
  +--------------+---------------+     +--------------v---------------+
                 |                     |  CALIBRATED                  |
                 |                     |  d <= 3 correction factors   |
                 |                     |  from per-subject NCA,       |
                 |                     |  clipped to public priors    |
                 |                     +--------------+---------------+
                 |                                    |
                 |                     +--------------v---------------+
                 |                     |  EMPIRICAL (v2)              |
                 |                     |  dense per-cell grid release |
                 |                     |  large pooled corpora only   |
                 |                     +--------------+---------------+
                 |            +---------------------- +
                 v            v
        +------------------------------------------+
        |  Structural generator                    |
        |  rxode2 model + protocol schedule +      |
        |  public variability + injected messiness |
        +--------------------+---------------------+
                             |
                             v
                  PMX event table + validation
```

- **Prior.** Reads no confidential data, spends no budget, needs no DP claim.
  The answer for very small studies, and the default for anyone who just wants a
  realistic table. Everything in sections 1 and 5 still applies.
- **Calibrated.** The main development work. A public structural model plus a
  small multiplicative correction learned from the data. Section 6.
- **Empirical.** The implemented Version 2 engine. Retained for large pooled
  corpora where its cost is affordable. Not the default.

The structural generator, event-table construction, censoring, schema
restoration, and validation are **shared by all three modes**. Only the source of
the population parameters differs.

---

# 4. Scenario guide

What to use, at what epsilon, and what to expect. Fold-error figures are measured
median error on CL from 200 replicate OpenDP releases; see
`design/FEASIBILITY.md` section 8.

Figures below are **measured against the implementation** — 40 replicate OpenDP
releases per cell, `d = 2`, an 8-fold correction-factor prior, reported as the
total fold-error on clearance including estimator bias. See
`design/FEASIBILITY.md` section 8.

**In every scenario, sampling times and event structure are exact**, because
they come from the protocol and the schema. What varies below is only how well
the exposure magnitude is pinned down.

## Smallest: N = 6

**Calibrated mode at epsilon 1-2, or Prior mode. Run the pre-flight check and decide.**

| Configuration | `f` | Measured total error |
|---|---:|---:|
| epsilon 0.25 | 1.33 | 2.3-fold |
| epsilon 1 | 0.33 | 1.7-fold |
| epsilon 2 | 0.17 | 1.6-fold |

Note how little epsilon buys here: eight times the budget moves the error from
2.3-fold to 1.6-fold, because at N = 6 the release is dominated by the prior and
by estimator bias rather than by noise. That is the signature of a configuration
where **Prior mode is the honest choice** — it gives roughly the same answer for
no privacy cost at all. Use the calibrated path at N = 6 only if the prior is
genuinely poor and the correction is expected to be large.

## Mid: N = 20

**Calibrated mode at epsilon 0.5. This is the design centre.**

| Configuration | `f` | Measured total error |
|---|---:|---:|
| epsilon 0.25 | 0.40 | 1.60-fold |
| **epsilon 0.5** | **0.20** | **1.60-fold** |
| epsilon 1 | 0.10 | 1.26-fold |
| epsilon 2 | 0.05 | 1.19-fold |

At epsilon 0.5 exposures land within about 1.6-fold, comfortably inside the
accuracy bar. Note that epsilon 0.25 measures the same 1.6-fold, so **the
cheaper setting is the better buy** — the extra budget at 0.5 is absorbed by
bias rather than converted into accuracy. Epsilon 1 does reach 1.26-fold, but
that is precision the use case does not need at twice the cost.

Stage-1 range-finding is not needed here. The correction-factor prior is already
tight enough.

## Large for this context: N = 300

**Calibrated mode at epsilon 0.1, or lower.**

| Configuration | `f` | Measured total error |
|---|---:|---:|
| epsilon 0.25 | 0.03 | 1.08-fold |
| epsilon 0.5 | 0.01 | 1.08-fold |
| epsilon 1 | 0.01 | 1.06-fold |

At N = 300 the constraint has stopped binding entirely, and the residual error
is estimator bias rather than privacy noise — epsilon 0.25 and epsilon 1 are
indistinguishable in the output. **Take the smaller epsilon and bank the
guarantee.** Anything below about 0.25 is free here.

Empirical mode becomes technically viable around N = 1000 at epsilon 5 and N = 10000 at
epsilon 1. Given the accuracy bar, it is hard to justify: it buys empirical
shape detail that is not required, at an epsilon an order of magnitude weaker.
Prefer Calibrated mode.

## Summary

| N | Mode | Recommended epsilon | Expected spread | Notes |
|---:|---|---:|---|---|
| 6 | Calibrated, or Prior | 1-2 | ~1.7-2.8-fold | Thin margin; check before spending |
| 20 | Calibrated | 0.5 | ~1.9-fold | Design centre |
| 60 | Calibrated | 0.25 | ~1.5-fold | Comfortable |
| 100 | Calibrated | 0.2 | ~1.4-fold | Comfortable |
| 300 | Calibrated | 0.1 | ~1.2-fold | Constraint no longer binds |
| 1000+ | Calibrated | 0.05 or lower | ~1.1-fold | Very strong guarantee |

These epsilons are roughly **an order of magnitude smaller** than a design
targeting parameter-estimate accuracy would require. That is the direct payoff
of the accuracy bar in section 1: a modest utility requirement converts into a
strong privacy guarantee.

---

# 5. Where the prior comes from

The clipping range is the dominant sensitivity term, so its provenance is the
most important governance question in the design.

## "Public" means independent of this dataset, not published

This is the key point, and it resolves the objection that new drugs have no
published models. Differential privacy requires that the prior be chosen
**without inspecting the confidential data**. It does not require that the prior
be in a journal.

For a first-in-human compound, all of the following are available before the
data exist, and all qualify:

| Source | Constrains |
|---|---|
| Preclinical allometric scaling / predicted human PK | CL, Vd, t-half — this is the model that picked the FIH dose |
| Protocol | Dose levels, sampling schedule, visit windows, study duration |
| Bioanalytical method validation | LLOQ, ULOQ, assay precision |
| Physiology | CL cannot exceed hepatic blood flow (~90 L/h); renal CL capped by GFR; Vd bounded below by plasma volume |
| Sampling window | t-half cannot be estimated beyond the observation period |
| Compound class | Typical ranges for a small molecule, a mAb, an ADC |
| Prior compounds in the same program or target class | Sponsor-internal but data-independent |

A predicted popPK model from preclinical scaling is often accurate to within
3-5 fold in humans. That is already a usable stage-2 prior and costs nothing.

## Stage 1: buy a tighter range privately

When the design-based prior is wide, narrow it with a small budget slice before
spending the rest on means.

**This is cheap because of an asymmetry that matters a great deal:** a histogram
has L1 sensitivity **1**, regardless of how many bins it has. Each subject falls
in exactly one bin. Learning *where* the data lies costs almost nothing; learning
a precise *mean* costs `d/(epsilon N)`.

So:

1. Build coarse log-scale bins over the wide design-based range.
2. Release the histogram at `epsilon_0`, typically 15-25% of the total budget.
3. Take the interval covering the bulk of the released mass, widened for noise.
4. Clip subject values to that interval and release the parameter means at
   `epsilon_1`.

Total loss is `epsilon_0 + epsilon_1` by basic composition. Measured payoff at
epsilon 1: N = 20 improves from 2.44-fold to 1.52-fold, N = 60 from 1.46 to
1.10. Roughly a 3x effective privacy gain for a 20% budget slice.

The stage-1 range must be widened enough that clipping does not itself become
data-dependent in a way the accounting misses. Prefer a conservative widening.

## What must never happen

Setting the clipping range from the data's own mean and standard deviation
**without spending budget** voids the guarantee entirely, and nothing downstream
detects it. Mean +/- 2 SD is a perfectly good *target* for what stage 1 should
discover, but it must be discovered through a private mechanism, not computed
directly. This is the single easiest way to accidentally destroy the privacy
claim while believing it is intact.

---

# 6. The v3 release design

## The structural model is a public input

The package should not invent its own notion of curve shape. The user supplies a
**structural model** — preferably an `rxode2` model — together with predicted
typical parameter values. The model and its predictions are public inputs in
exactly the same sense as the schema and the dose levels: they exist
independently of the confidential dataset.

### Why this is the right shape

Four things follow, and together they are worth more than any mechanism
improvement:

1. **It collapses the release.** If the model supplies structure, `t-half` is
   not a released quantity — it falls out of CL and V. Neither does curve shape,
   accumulation, or the infusion profile. The release shrinks from "a handful of
   parameters" to "a correction to the prediction", and `d` drops from 6-8 to
   3-4. Since error is `d/(epsilon N)`, that is a direct proportional gain.

2. **It makes the linear-PK assumption swappable rather than baked in.** The
   earlier design assumed linear PK so that dose-normalization was valid. With a
   supplied model, a compound with target-mediated disposition uses a TMDD model
   and the assumption disappears. Nonlinearity, accumulation on repeat dosing,
   and infusion behavior all come out right because the model computes them
   rather than because the generator approximates them.

3. **It gives dose-exposure monotonicity for free.** Higher doses produce higher
   exposures because the model says so, at zero privacy cost. This was an
   explicit objective in section 1 and it is now structural rather than learned.

4. **It sources the prior from where the knowledge actually lives.** For a
   first-in-human study, a predicted human PK model already exists — it is what
   selected the starting dose. For a later study, an earlier study in the same
   program is data-independent *with respect to this dataset*, which is exactly
   what differential privacy requires.

### Interface

`rxode2` is a `Suggests` dependency, needed only at generation time. Generation
is post-processing, so a heavy dependency there costs nothing privacy-wise and
does not burden users of Prior mode.

```r
model <- pmx_structural_model(
  rx        = my_rxode2_model,
  typical   = c(cl = 10, v = 70, ka = 1),
  source    = "allometric scaling from rat and dog, FIH prediction memo v3",
  endpoints = c(cp = "cp", response = "R")
)
```

`source` is required and recorded in the release ledger. A structural model
without data-independent provenance is a governance failure, not a modeling
choice.

Where no model is supplied, fall back to the Version 2 grid behavior with its
worse constant, or to Prior mode.

## What is released: a correction, not a parameter

**Release the multiplicative correction between the model's prediction and the
data, not the absolute parameter.** This is the single most valuable design
decision in Version 3, and the reasoning is worth stating carefully.

The clipping range drives the error. A prior on absolute CL for a new compound
might span 100-fold, because CL genuinely could be almost anything. But a prior
on *how wrong the preclinical prediction is* is much tighter and far better
characterized: allometric scaling lands within roughly 2-3 fold for most
compounds. So `[1/4, 4]` on the correction factor is defensible from public
knowledge alone.

| Parameterization | Prior span (log units) |
|---|---:|
| Absolute CL, new compound | ~4.6 (100-fold) |
| Absolute CL after stage-1 range-finding | ~1.6 (5-fold), costs ~20% of budget |
| **Correction factor on a predicted CL** | **~2.1 (8-fold), costs nothing** |

The correction-factor prior is nearly as tight as a privately-purchased range,
and it is free. Combined with the smaller `d`, the estimated effect at
N = 20, epsilon 1:

| Approach | `d` | Span | Estimated fold-error |
|---|---:|---:|---:|
| Absolute prior + stage-1 range-finding | 5 | 1.61 | ~1.65 |
| Correction factor, no stage 1 | 3 | 2.08 | **~1.37** |

These are planning estimates from the error law confirmed in
`design/FEASIBILITY.md` section 8, not measurements of an implementation.

### Draft release vector

| Released | Meaning | Prior source |
|---|---|---|
| Cohort size | count, sensitivity 1 | n/a |
| PK correction | `log(CL_observed / CL_predicted)`, per subject by NCA | Allometric scaling accuracy |
| PD correction | `log(level_observed / level_predicted)`, a ratio of means | Baseline literature, inclusion criteria |

Add a PD baseline correction, or a between-subject variability term, only when a
workflow specifically needs it. Every addition raises `d` and costs accuracy
proportionally.

Between-subject variability is a **public assumption** (for example 50% CV) by
default. A second moment costs roughly twice the sensitivity of a mean, so
release it only deliberately.

### The correction factor is a free diagnostic

The released correction is already public and already accounted, so reporting it
costs nothing. It directly answers "how wrong was my prior?"

A correction near 1 confirms the prediction. A correction pressed against the
clipping boundary means the prior was wrong and the release is censored — the
generated data is then driven by the boundary, not by the study, and the user
must be told. This is a real mitigation for the "confidently wrong output" risk
in `design/FEASIBILITY.md` section 8, and it is the one self-check the design
gets for free.

## Baseline covariates

The generator emits declared baseline covariate columns so that
covariate-handling pipeline code has something to run against. Fidelity is
secondary: the goal is exercising joins, filters, and covariate-model plumbing,
not reproducing a covariate analysis.

There are two ways to declare them, differing in privacy model.

**DP covariates** (`pmx_covariate()` / `pmx_covariates()`) take a public range
(continuous) or level set (categorical), chosen without inspecting the data. In
prior mode the values are drawn from those public declarations at no cost. In
calibrated mode each is released privately and costs exactly one budget slice,
regardless of its number of levels: a continuous covariate releases its clipped
mean, a categorical one its level-count vector, whose L1 sensitivity is one
because adding or removing a subject changes exactly one count by one. This
keeps the whole release differentially private.

**Bootstrap covariates** (`pmx_covariates_auto()`) take only a list of column
names. Values are resampled from the data: a uniform draw over the observed
range (continuous) or a proportional resample (categorical), the approach used
by Novartis's `synadam`. This is **not** differentially private -- it exposes
the data-derived support of each column -- and a model that uses it is marked
accordingly, with a warning and a `covariates_private = FALSE` flag in the
privacy record. It spends no budget and does not enter `d`. By default the
continuous range is clipped to the 1st and 99th percentiles so the two extreme
subjects are not exposed; `clip = NULL` restores the raw min/max, matching
`synadam` exactly. Use it only inside a trusted environment.

What remains out of scope is any **covariate-response relationship**. Covariates
are generated independently of the PK and PD, so a workflow that checks, say, a
weight effect on clearance will not find one. That is a modeling relationship,
not pipeline plumbing, and it is excluded by the same logic that keeps
exposure-response PD out. Note also that adding covariates raises `d` and so
dilutes the PK correction; the pre-flight in section 2 accounts for the full
`d`, including covariates.

## PD is a simple, exposure-independent time course

The PD endpoint is a constant, a linear trend, or an exponential approach to a
plateau (covering both decay and growth). It has a plausible time course and a
plausible magnitude, and it is deliberately **not** linked to the generated
concentrations. There is no Emax, no indirect-response model, no exposure
coupling.

This is the right choice for two reasons.

**It is adequate.** The accuracy bar in section 1 asks for a plausible PD column
so that longitudinal analysis code can be exercised. A pipeline reading a
biomarker does not check that the biomarker is Emax-consistent with exposure.

**It calibrates cleanly.** An exposure-driven correction would have to estimate a
*deviation from baseline*, and a PD deviation rides on a large baseline, so its
per-subject signal-to-noise is poor and the resulting estimator is biased low in
a way that does not shrink with N. A simple shape instead takes a **level
correction**: the ratio of the mean observed response to the mean predicted
response, where both terms are the response itself. Measured with 15% residual
error, that recovers 2.53 against a true 2.5.

There is a second advantage. A prior on a PD *baseline* can be tight, because
baselines are usually well characterized from healthy-volunteer literature or
the inclusion criteria. A prior on Emax for a new mechanism cannot be. Since the
prior span drives the whole error law, that difference is worth more than the
mechanistic realism given up.

The cost, stated plainly: a workflow that specifically checks exposure-response
coupling will not find it in this data. That is out of scope, by the same logic
that puts covariates out of scope — it is a modeling relationship, not a
pipeline-plumbing feature.

## When PD is still weaker than PK## When PD is still weaker than PK

Preclinical-to-clinical translation is less reliable for pharmacodynamics than
for pharmacokinetics. The correction-factor prior for a PD effect is
correspondingly wider — perhaps `[1/10, 10]` rather than `[1/4, 4]` — which
gives back part of the gain.

Plan for PD needing more budget than PK for the same relative accuracy. The
released quantity is a PD *level* correction; the shape comes from the
structural model and is never learned.

## Model selection is itself a privacy-relevant choice

Choosing the structural model by looking at the data is a data-dependent
decision, and the framework does not care whether a human or a language model
does the looking. Selecting a two-compartment model because the concentration
plot looks biphasic leaks information about the source, and that leak is outside
the accounting.

Two acceptable procedures:

1. **Public selection.** Choose from mechanism, compound class, modality, and
   preclinical data only. A library such as `nlmixr2lib` is a reasonable source
   of candidate structures, and an LLM agent may help navigate it, provided the
   agent is reasoning from the compound's public description and never from the
   data. See `design/MODEL_ELICITATION.md` for the interview that keeps this
   honest.
2. **Budgeted selection.** Choose privately from a small prespecified candidate
   set with an exponential mechanism, and account for it.

The API must make procedure 1 the path of least resistance: the structural model
is constructed **before** `fit_private_pmx()` is called and is passed in as a
public input, so there is no ergonomic route to fitting the data, looking at it,
and revising the model.

**Never send the confidential data to an external service, including a language
model.** That is a disclosure entirely outside the privacy accounting, and no
downstream validation in this package would detect it. State this in the user
documentation, not only here.

## Non-compartmental analysis, not popPK fitting

The estimator choice is a privacy constraint, not a modeling preference.

- **NCA is DP-compatible.** Trapezoidal AUC and a terminal slope are computed
  from one subject's own rows. One subject influences only their own value, so
  clipping to the stage-1 range bounds sensitivity directly.
- **NLME/popPK fitting is not.** Shrinkage couples subjects; every individual
  post-hoc estimate depends on the population fit, which depends on everyone.
  One subject perturbs all N estimates with no simple bound. Making that private
  needs DP-SGD-style machinery and far more budget.

This aligns with the scope: NCA needs rich sampling, and rich sampling is what
small early-phase studies have. Sampling richness and cohort size are inversely
related across development, so the small-N target and the DP-compatible
estimator want the same designs.

When sampling is too sparse for NCA, fall back to Prior mode rather than to a
population fit.

## Mechanism: pure-DP Laplace

Keep Laplace. Gaussian under zCDP helps only at high dimension: at `d = 6` and
epsilon 1, Laplace gives a noise scale of 6 on the sum while a Gaussian at
delta 1e-6 gives roughly 13. The `sqrt(d)` sensitivity advantage does not
overcome the `sqrt(2 ln(1.25/delta))` constant until `d` is well above 8.

Gaussian/zCDP remains appropriate for Empirical mode, where `d` is in the tens.

## Structure that costs nothing

Take from public inputs, never from the data:

- **The sampling schedule**, from the protocol. This removes the entire
  `endpoint_timing` release group — 36 dimensions in the current fixture — and
  frees its budget.
- **Curve shape, accumulation, and dose-exposure monotonicity**, from the
  structural model. Zero cost.
- **Visit windows, assay limits, BLQ rules.**

## Injected messiness

Synthetic data that is too clean does not exercise the pipelines this package exists
to test. A perfectly smooth model-generated dataset will not surface the bugs
that a real dataset surfaces.

Inject from public and protocol knowledge, not from the data:

- dropout, from the protocol's expected completion rate;
- missed doses and out-of-window visits;
- BLQ runs, from the assay LLOQ and the predicted concentration-time profile;
- occasional outliers and duplicate records.

Spend a small budget slice on realized rates only if the public assumption is
too crude to be useful.

## Structural assumptions are declared, not implied

With a supplied structural model there is no hidden linearity assumption: the
model states its own. What must still be declared and checkable is whether the
**NCA estimator** used for the correction factor is appropriate. Trapezoidal AUC
with a terminal slope assumes adequate sampling through the elimination phase.

Where sampling is too sparse for NCA, fall back to Prior mode rather than to a
population fit — see the estimator constraint above.

---

# 7. Privacy contract

Unchanged from Version 2 and binding on all private modes.

## Formal guarantee

Subject-level, add-or-remove `(epsilon, delta)` differential privacy. One
subject's complete longitudinal contribution is the privacy unit: all rows and
visits, dosing and infusion records, actual and nominal timing, baseline
covariates, all endpoint observations, DVID/censoring/missingness, and rare
schedules or protocol deviations.

For neighboring datasets `D` and `D'` differing by one complete subject, and
every output set `S`:

```text
Pr[M(D) in S] <= exp(epsilon) * Pr[M(D') in S] + delta
```

The accurate claim is:

> Generated from a subject-level `(epsilon, delta)`-differentially private model.

Never claim zero privacy risk, impossibility of re-identification, or legal
anonymity. Never use `private`, `anonymous`, `safe`, or `de-identified` as
unqualified binary labels.

## One person, one record

The adjacency above assumes each person appears exactly once. Rollover and
extension studies, crossovers pooled with parallel-group studies, and re-enrolled
subjects all violate it; a person contributing `k` records receives roughly
`k * epsilon` by group privacy.

Either require an explicit assertion that IDs are unique persons, or accept a
person-level grouping column and bound contributions on that instead of on the
study subject ID. Tracked as `REV-016`.

## Epsilon and delta in plain language

Epsilon is the **one-person influence limit**. Zero means the output cannot
depend on the source patients at all; smaller is stronger; a mechanism can
satisfy DP at a very large epsilon while providing weak practical protection.
Require explicit epsilon and do not ship a universal default.

Delta is a small additive allowance in the probability bound. It is **not** the
probability of re-identification, the fraction of unprotected patients, or a
breach probability. Prefer `delta = 0`; require justification for any positive
delta. The low-dimensional path uses pure DP, so delta is zero there.

## Identity, attributes, participation

Document that a release can expose that an unusual record existed without
identifying its owner; that absence of direct identifiers does not prevent
linkage; that membership and attribute disclosure are distinct protected
questions; and that DP limits the additional information attributable to one
person's participation even against an attacker with auxiliary information.

---

# 8. Public API

## Roles

```r
roles <- pmx_roles(
  id = "ID", time = "TIME", nominal_time = NULL, tad = NULL, occasion = NULL,
  dv = "DV", amt = "AMT", evid = "EVID", cmt = "CMT", dvid = NULL, mdv = NULL,
  rate = NULL, cens = NULL, limit = NULL, covariates = "WT",
  subject_properties = "ACTARM", assigned_dose = NULL, exclude = NULL
)
```

Critical roles are explicit. Never infer column meaning from names. Subject
properties such as `ACTARM`, `TRT`, or a nominal dose group are categorical
subject-level assignments whose released strata condition the generated regimen.
An occasion-varying assigned-dose column is derived from the generated event
amount, not sampled independently.

## Endpoints

```r
endpoints <- list(
  cp = pmx_endpoint(dvid = "cp", alignment = "dose_relative",
                    transform = "log", shape = "occasion")
)
```

Alignments: `dose_relative`, `study_time`, `occasion`, `hybrid`. Alignment is
user-declared scientific metadata. Any inference of endpoint behavior from
confidential data must be private and budgeted.

## Priors (new in v3)

```r
priors <- pmx_priors(
  cl      = pmx_prior(range = c(2, 50), source = "allometric scaling from rat/dog"),
  t_half  = pmx_prior(range = c(2, 24), source = "preclinical, sampling window"),
  ...
)
```

Every prior carries a `source` string recording its provenance. The fit records
these in the release ledger. A prior without a data-independent provenance is a
governance failure, not a modeling choice.

## Fit once, generate many

```r
private_model <- fit_private_pmx(data, roles, endpoints, epsilon, delta,
                                 priors, public_design, contribution_limits,
                                 budget_allocation)
synthetic <- generate_pmx(private_model, n_subjects = NULL, seed = 123)
privacy_report(private_model)
```

- `fit_private_pmx()` is the only stage that reads confidential data.
- `generate_pmx()` reads only the fitted model; repeated generation is
  post-processing and costs nothing.
- Refitting against the source composes and costs more budget.
- Generation defaults to the fitted privacy-accounted subject count.
- Privacy noise is never user-seeded, logged, or exposed.

A pre-flight helper must report predicted error per release group for a given
`(N, epsilon, d)` **before** budget is spent, and refuse or warn loudly on
infeasible configurations (`REV-002`).

## Public-design-only generation

Prior mode needs a first-class entry point that takes no confidential data and no
epsilon. It must not be reachable only as a testing backdoor.

---

# 9. Shared machinery

Binding on all modes.

## Contribution bounding

Deterministically bound each subject's contribution before computing any private
statistic: rows, doses, occasions, observations per endpoint and occasion,
timing cells, and bounded dose/rate/time/covariate/DV/limit domains. All
sensitivity calculations protect the subject's complete bounded contribution,
not a single row. Never release which subjects were clipped.

## Event structure

Construct coherent PMX records from the population model and public structural
rules. Never copy source row blocks. Generate EVID, AMT, RATE or duration, CMT,
DVID, ADDL/II, observation rows, MDV, censoring fields, and tied-time ordering
together. Never numerically average event-control fields.

## Nominal and actual time

Generate chronological actual-like times around the nominal design, preserving
nondecreasing within-subject time (or within subject and occasion when the clock
resets), correct dose occasion, postdose observations after the qualifying dose,
tied collection blocks, dose/predose ordering, nonnegative time, and consistent
TAD/occasion derivations. Never copy a subject's complete time vector,
missing-visit pattern, or infusion duration.

## Censoring

Monolix convention: `CENS = 0` uncensored; `CENS = 1` left censored with DV
holding the boundary; `CENS = -1` right censored; optional `LIMIT` for the other
boundary. Generate a latent trajectory first, then derive DV, CENS, and LIMIT
coherently. Never blend, average, or independently perturb CENS.

Validate allowed values, boundary direction, finite limits, lower not exceeding
upper, the DV-equals-boundary convention, and consistency among EVID, MDV, DV,
DVID, CENS, and LIMIT.

## Identifiers and schema

Generate new IDs. Never retain source IDs. Reject direct identifiers. Never
export calendar dates or datetimes; stop on unmodeled Date/POSIXct columns
unless explicitly excluded. Preserve declared public schema, names, order,
classes, and factor levels.

The fitted model must not contain raw rows, subject profiles, source IDs, event
templates, raw residuals, unnoised aggregates, confidential category sets or
bounds, or data-dependent debugging caches.

## Accounting

Account for every source-dependent released computation, including the stage-1
range-finding histogram, the subject count, parameter releases, any
data-dependent bound, and any released diagnostic. Realized loss must not exceed
the request. Emit a machine-readable release-ledger entry per fit.

**Data-dependent failures are releases.** A `stop()` whose occurrence depends on
one subject's values is an unaccounted output channel. Prefer clipping or
dropping out-of-domain records over erroring (`REV-003`).

## Validated backend

Never hand-code production mechanisms with ordinary R RNG. Use OpenDP behind a
small adapter. Fail closed when unavailable; never fall back to ordinary noise.
Document backend, version, adjacency, mechanism, accountant, epsilon, delta,
contribution bounds, proof assumptions, and floating-point handling.

## Post-processing discipline

Generation, constraint repair, plotting, and validation may use the fitted model
freely provided they do not consult the source. `compare_pmx()` is not
releasable merely because `synthetic` is private; any source-derived diagnostic
leaving the restricted environment must itself be budgeted.

## Validation

```r
validate_pmx(data, roles); validate_private_model(private_model)
privacy_report(private_model); compare_pmx(source, synthetic, roles)
```

Cover schema and classes, event coherence, chronological and tied-row ordering,
endpoint alignment, repeated-dose and study-time behavior, censoring
consistency, finite domain values, new IDs, and absence of direct identifiers.

---

# 10. Demonstrations

Public `nlmixr2data` datasets remain the structural test set: `theo_md`
(repeated Q24H dosing, dose-relative PK), `warfarin` (lower-case schema, `cp`
PK plus `pca` PD, factor preservation), `wbcSim` (infusion start/stop, study-time
decline/nadir/recovery), `nimoData` (nominal-dose property strata, declared
TAD/OCC, washout), `mavoglurant` (TIME reset within OCC, occasion-varying
assigned dose).

**These are structural demonstrations, not utility demonstrations.** At 12 to
120 subjects they are below the feasibility frontier for any defensible epsilon.
They prove the event tables are coherent; they cannot show an achievable
privacy-utility tradeoff, and the vignettes must say so.

Add a larger fully simulated public fixture for privacy-utility evaluation at
N = 300 and above.

---

# 11. Version history

Newest first.

## Version 4 — return to AVATAR blending as the primary method (2026-07-22)

**Driver:** comparison to Novartis's `synadam`. `synadam` generates synthetic
ADaM data by resampling each column marginally from the real data — uniform over
the observed range for continuous columns, proportional resample for categorical
— with no differential privacy and no formal guarantee, relying on governance
and the fact that the data stays in a controlled environment. That is accepted
practice. AVATAR-style blending is the same governance-based idea applied at the
granularity of a whole subject trajectory rather than a single column, so if
`synadam`'s model is acceptable, AVATAR's is acceptable for the same
trusted-environment use.

**Change:** the Version 1 AVATAR engine is restored as `synthesize_pmx()`
(renamed from `mock_pmx`) and made the primary, default method. It resamples a
source subject's event skeleton as a template and fills it with a
distance-weighted blend of similar subjects' covariates and trajectories, plus
subject and residual noise. It works at any cohort size, needs no elicitation,
and makes no formal privacy claim.

**Rationale.** For synthetic data that stays inside a trusted environment — the
common case, and this package's stated purpose — a resampling method is more
useful than differential privacy on every axis the user cares about, and DP's
formal guarantee buys nothing when no adversary can reach the output. See
`design/METHOD_DISCUSSION.md`.

**The asymmetry, recorded honestly.** A resampled covariate value is weakly
identifying; a resampled subject trajectory is strongly identifying, close to a
fingerprint. AVATAR therefore sits a notch higher on the risk ladder than
`synadam`'s column resampling, and its safety depends more heavily on the
governance context. This is the same concern that motivated moving away from
Version 1 originally; Version 4 accepts it deliberately, scoped to
trusted-environment use, with the DP engines available when the output crosses a
trust boundary.

**Retains** Versions 2 and 3 as alternatives, not deprecated. The decision rule
is the trust boundary: inside → AVATAR; crossing out → a DP engine.

## Version 3 — low-dimensional structural release (2026-07-22, in design)

**Driver:** the target cohort sizes were fixed at 6-10000 with the small end
prioritized, and the Version 2 dense-grid engine was measured to need N ~ 1000
at epsilon 5 and N ~ 10000 at epsilon 1. It does not serve the scope.

**Change:** replace the dense per-cell grid with `d <= 8` bounded scalars
released against public priors, structure taken from the protocol and preclinical
predictions, per-subject NCA as the estimator, and two-stage private
range-finding. Measured effect: a defensible epsilon 1 becomes viable at N ~ 60
rather than N ~ 6000, and N = 20 lands within about 1.5-fold on CL.

**Supersedes** these Version 2 sections: fixed-dimensional subject
representation, privately learned information, endpoint trajectory generation,
and sampling-design inference. Budget allocation is revised — no trajectory
group, and a new stage-1 range-finding group.

**Retains** the Version 2 privacy contract, roles, endpoint declarations,
contribution bounding, event and timing construction, censoring, identifiers and
schema, accounting, backend discipline, post-processing rules, and validation.

Version 2 remains implemented and supported as Empirical mode for large pooled corpora.

## Version 2 — subject-level DP population generator (implemented)

Replaced the Version 1 blending architecture with a fit-once/generate-many
subject-level differentially private population generator: explicit public
schema and design declarations, numeric bounds, contribution limits, budget
allocation, an OpenDP adapter that fails closed, fixed-dimensional bounded
subject summaries, composed accounting, release ledgers, and leakage validation.

Version 2 must not select a source subject as an event-template anchor, copy a
complete source event skeleton, choose raw subjects as nearest-neighbor donors,
blend raw subject trajectories, expose a `blend`/`avatar`/`template`/exact-timing
mode, or describe jitter, generalization, k-anonymity, or disclosure testing as a
substitute for differential privacy. These prohibitions carry forward to
Version 3.

Known limitations, all tracked in `design/REVIEW_BACKLOG.md`: sensitivity is
charged at `ncol` rather than derived from contribution limits (`REV-005`);
delta is validated but never spent (`REV-004`); data-dependent `stop()` paths
remain (`REV-003`).

## Version 1 — AVATAR-style blending (superseded, in git history)

Anchor-and-donor synthesis: subjects grouped by an event signature, event
skeletons copied as templates with time jitter, PCA over subject profiles,
k-nearest-neighbor donor blending of covariates and trajectories, plus subject,
residual, and AR(1) noise.

Abandoned because it offered no formal guarantee and its privacy degraded
precisely where it was most needed. Copied event skeletons are quasi-identifiers;
`k` neighbors is a large fraction of a small cohort; and signature grouping
partitions before blending, so an unusual regimen lands in a near-singleton group
and is blended with itself. See `design/FEASIBILITY.md` section 2.

The instructive comparison: Version 1 produced attractive output at N = 12
*because* it carried individual information through, while Version 3 tells you
honestly when it cannot help. Both face the same information-theoretic
constraint; only one of them reports it.
