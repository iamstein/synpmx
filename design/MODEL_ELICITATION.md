# Building the structural model: an elicitation guide

A question-by-question procedure for producing the **public structural model and
priors** that `pmxSynthData` needs before it touches any confidential data.

Written to be worked through by a pharmacometrician, or by an LLM agent
assisting one. The output is an `rxode2` model, a set of prior ranges, and a
provenance table.

Companion to `design/PROTOTYPE_SPEC.md` sections 5 and 6, which explain why the
design works this way. If the `d` and `f` arithmetic in Q7 is unfamiliar, read
`design/PRIVACY_BACKGROUND.md` first.

---

## The cardinal rule

> **Every answer must come from knowledge that exists independently of the
> dataset you are about to fit.**

Not "published" — *independent*. A preclinical scaling memo, a protocol, an
assay validation report, and an earlier study in the same program all qualify.
The dataset you are about to analyze does not.

This is not bureaucratic. The prior range is the dominant term in the privacy
arithmetic (`design/FEASIBILITY.md` section 8). A range set by peeking at the
data voids the differential-privacy guarantee entirely, and **nothing downstream
detects it** — the release will look completely normal, the accounting will
balance, and the guarantee will be false.

### Stop signals

If any of these occurs, stop and re-plan rather than continuing:

- "Let me check the data first." — Then it is not a prior.
- "I'll widen the range until it covers everything I saw." — Same.
- "I'll paste a few rows into the model to help pick a structure." — This is a
  disclosure outside the accounting entirely. Never send confidential data to an
  external service, including a language model.
- "The plot looks biphasic, so two compartments." — Data-dependent model
  selection. Either choose from mechanism and preclinical data, or use a
  budgeted private selection mechanism.

Answering "I don't know" is fine and often correct. A wide honest prior costs
budget; a narrow dishonest one costs the guarantee.

---

## Q0. Gate: how many subjects, and do you need a privacy claim?

| Answer | Route |
|---|---|
| Fewer than about 10 subjects | **Either.** Calibrated mode is viable at epsilon 1-2 with a thin margin; Prior mode is entirely reasonable. Run the Q7 pre-flight and decide |
| 10-300 subjects, claim needed | **Calibrated mode.** The main path. Continue |
| More than ~1000, and you want empirical shape rather than a structural model | **Empirical mode.** Version 2 dense grid. Different document |
| No privacy claim needed at all | **Prior mode.** Simpler, and honest |

Full selection guide, including how many parameters each mode learns from the
data, is in `design/PROTOTYPE_SPEC.md` section 3.

Prior mode still needs Q1-Q5. It just never runs a private fit. **Most of this
interview is about building a good simulator, not about privacy.**

---

## Q1. Compound and program context

1. **Is this first-in-human?**
   - Yes → PK priors come from preclinical allometric scaling. Expect the
     prediction to be within about 2-3 fold for CL.
   - No → priors come from earlier studies in the program. Much tighter;
     a 2-fold correction-factor prior may be defensible. Say which study.
2. **Modality?** Small molecule, mAb, ADC, peptide, oligonucleotide, cell
   therapy. This sets plausible CL and V scales and typical half-life ranges
   before any compound-specific information.
3. **Route and formulation?** IV bolus, IV infusion, oral, SC. Determines
   whether absorption parameters are needed at all.
4. **Are there prior compounds in the same class or against the same target?**
   Sponsor-internal knowledge counts as public here, provided it is not derived
   from the dataset being fitted.

---

## Q2. PK structure and predicted parameters

1. **How many compartments?** Choose from mechanism and preclinical profile
   shape, not from this study's data. When unsure, one compartment is usually
   adequate — the objective is "vaguely right", not model selection.
2. **Predicted CL**, with units and source.
3. **Predicted V** (central; plus peripheral V and Q if two-compartment).
4. **Predicted F** for extravascular routes.
5. **Predicted ka** or absorption half-life, for oral or SC.
6. **Is elimination linear over the dose range studied?** If not — TMDD,
   saturable metabolism — say so and use a nonlinear structure. This is exactly
   the case a supplied structural model handles and a dose-normalizing
   approximation does not.
7. **How confident are you in the CL prediction?** This becomes the
   correction-factor prior:

   | Confidence | Correction prior | Typical situation |
   |---|---|---|
   | High | `[1/2, 2]` | Later study, same population, same formulation |
   | Moderate | `[1/4, 4]` | FIH from good allometric scaling |
   | Low | `[1/10, 10]` | Novel modality, uncertain scaling, first in class |

   Choose honestly. This is the single number with the largest effect on how
   much privacy budget you need.

---

## Q3. Endpoints

For each endpoint:

1. **Name, units, and DVID coding.**
2. **Type:** PK concentration, PD biomarker, efficacy, safety lab.
3. **Scientific clock:** does it reset after each dose (`dose_relative`) or
   evolve over the study (`study_time`)?
4. **For PD endpoints:**
   - **Baseline** value and its plausible range. Often well characterized from
     healthy-volunteer literature or the inclusion criteria.
   - **Direction:** does it rise, fall, or both?
   - **Mechanism:** direct effect, indirect response (turnover), transit delay?
     Drives which `nlmixr2lib`-style structure to use.
   - **Predicted maximal effect** and the exposure producing it, with source.
   - **Onset timescale:** minutes, hours, days, weeks?
5. **Is the endpoint censored?** LLOQ, ULOQ, or both.

For PD, expect to need a wider correction prior than PK — preclinical-to-clinical
translation is less reliable. `[1/10, 10]` is a reasonable default.

---

## Q4. Protocol

All of this is a public document. None costs privacy budget.

1. **Dose levels and cohort assignment.**
2. **Dosing schedule:** single, QD, BID, Q2W, loading dose, titration.
3. **Number of doses and study duration.**
4. **Infusion duration**, if applicable.
5. **Sampling schedule** per endpoint — the nominal times from the protocol.
   Supplying this removes the entire timing-inference release group.
6. **Visit windows** — how far actual times may deviate from nominal.
7. **Planned N**, total and per cohort.
8. **Expected dropout rate**, from the protocol's assumptions.

---

## Q5. Assay and messiness

1. **LLOQ and ULOQ** per endpoint, from method validation.
2. **Assay precision**, as a CV. Becomes residual error.
3. **BLQ handling convention** in your datasets: dropped, set to zero, set to
   LLOQ/2, or `CENS` flagged.
4. **Expected messiness** you want reproduced, so the mock data exercises your
   pipeline: missed doses, out-of-window visits, missing covariates, duplicate
   records, occasional outliers.

Mock data that is too clean does not surface the bugs a real dataset surfaces.
Derive these rates from the protocol and from general experience, not from the
dataset.

---

## Q6. Variability

1. **Between-subject variability** on CL, V, and PD parameters. A public
   assumption of 30-50% CV is usually fine and costs nothing. Only release IIV
   privately if a workflow specifically needs realistic shrinkage behavior, and
   note that a second moment costs about twice the sensitivity of a mean.
2. **Residual variability**: proportional, additive, or combined, and roughly
   what magnitude. Usually the assay CV plus model misspecification.
3. **Covariate effects.** Prefer public structural relationships — allometric
   weight scaling on CL and V — over privately learned ones.

---

## Q7. Privacy budget (Calibrated mode only)

1. **What epsilon is approved?** This must come from governance, not from
   whichever value makes the plot look best.
2. **What will you release?** Target three: cohort size, PK correction, PD
   correction. Each addition costs accuracy proportionally.
3. **Run the pre-flight check** before spending anything. With `d` releases,
   epsilon `e`, and `N` subjects, compute `f = d / (e * N)`. This is the
   fraction of the prior's width that survives as noise.

   | `f` | Verdict |
   |---:|---|
   | >= 1.0 | The release tells you nothing the prior did not. **Use Prior mode** |
   | ~0.5 | Halves the uncertainty. Marginal but real |
   | ~0.25 | Clearly worthwhile |
   | <= 0.1 | **Lower epsilon instead.** You are buying accuracy you do not need |

   The question is not "is the error small" but "does this beat the prior".

4. **Start from the smallest epsilon that clears the bar, not the largest one
   approved.** Because the accuracy requirement is modest, most studies should
   land an order of magnitude below what a parameter-estimation exercise would
   need: epsilon 0.5 at N = 20, epsilon 0.1 at N = 300.

---

## Output

### 1. An rxode2 structural model

Template for one-compartment oral PK with an indirect-response PD endpoint.
**Not executed during authoring of this document** — the environment lacked a C
compiler. Verify it compiles and produces sensible profiles before use.

```r
library(rxode2)

pk_pd <- rxode2({
  cl <- tcl * exp(eta.cl)
  v  <- tv
  ka <- tka

  d/dt(depot) = -ka * depot
  d/dt(centr) =  ka * depot - (cl / v) * centr
  cp = centr / v

  # Indirect response: drug inhibits loss of R
  d/dt(R) = kin - kout * (1 - emax * cp / (ec50 + cp)) * R
})

typical <- c(tcl = 10, tv = 70, tka = 1,
             kin = 10, kout = 0.1, emax = 0.8, ec50 = 1)
```

Substitute a different structure where mechanism demands it. `nlmixr2lib` is a
reasonable source of candidates; select by mechanism, never by looking at data.

### 2. Prior ranges with provenance

```r
priors <- pmx_priors(
  pk_correction = pmx_prior(
    range  = c(1/4, 4),
    source = "allometric scaling, rat and dog; FIH prediction memo v3"
  ),
  pd_correction = pmx_prior(
    range  = c(1/10, 10),
    source = "target turnover from literature; no clinical PD precedent"
  )
)
```

### 3. A provenance table

One row per public input, to be recorded in the release ledger and reviewed by
whoever approves the release.

| Input | Value | Source | Data-independent? |
|---|---|---|---|
| Structural model | 1-cmt oral + IDR | Mechanism, preclinical | Yes |
| Predicted CL | 10 L/h | Allometric scaling memo v3 | Yes |
| PK correction prior | `[1/4, 4]` | Scaling accuracy literature | Yes |
| Sampling schedule | Protocol section 7.2 | Protocol v4 | Yes |
| LLOQ | 0.5 ng/mL | Method validation report | Yes |
| Expected dropout | 15% | Protocol assumptions | Yes |

Any row that cannot be answered "Yes" in the final column must either be removed
from the public inputs or moved inside the privacy budget.

---

## Worked example: FIH small molecule, 20 subjects

| Question | Answer |
|---|---|
| Q0 | 20 subjects, claim needed → **Calibrated mode** |
| Q1 | FIH, small molecule, oral |
| Q2 | 1-compartment; CL 10 L/h, V 70 L, F 0.5, ka 1/h; linear; moderate confidence → `[1/4, 4]` |
| Q3 | `cp` (ng/mL, dose-relative, LLOQ 0.5); `biomarker` (study-time, IDR, baseline 100, falls, `[1/10, 10]`) |
| Q4 | SAD, 5 cohorts, single dose; rich sampling to 48 h; 15% dropout |
| Q5 | Assay CV 15%; BLQ flagged with `CENS` |
| Q6 | 40% CV on CL, 30% on V; proportional residual 20% |
| Q7 | epsilon 1; release cohort size, PK correction, PD correction (`d = 3`) |

Pre-flight at the approved epsilon 1: `f = 3 / (1 * 20) = 0.15`. PK prior span
`log(4/(1/4)) = 2.08`, so PK error is `exp(0.15 * 2.08) = 1.37`-fold. PD prior
span `log(100) = 4.6`, so PD error is `exp(0.15 * 4.6) = 2.0`-fold.

**But epsilon 1 is more than this study needs.** At epsilon 0.5, `f = 0.30`,
giving PK 1.9-fold and PD 4.0-fold. PK still clears the accuracy bar
comfortably. PD does not.

**Recommendation:** take **epsilon 0.5** and spend it on PK only, dropping the
PD correction and generating PD from the public structural model alone. This
halves the privacy cost, keeps PK within about 2-fold, and avoids a PD release
that was never going to be informative given a 100-fold prior. `d` falls from 3
to 2, improving PK further to about 1.6-fold.

The general lesson: when a prior is too wide for the budget to narrow usefully,
**drop that release rather than funding it**. A release with `f` near 1 costs
privacy and returns nothing. Do not raise epsilon to make a number look
better.
