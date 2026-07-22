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
- **A real target dataset.** The whole frontier analysis assumes the fixture's
  release dimension. If the datasets of interest have more endpoints, more
  occasions, or more covariates, `d` grows and the required N grows with it.
  What are the actual cohort sizes and structures you need to support?
