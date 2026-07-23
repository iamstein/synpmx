# Method discussion: AVATAR blending vs formal differential privacy

Why `synpmx` uses AVATAR-style blending as its primary method, when it
also contains a fully differentially private engine. Written after building
both, measuring both, and comparing to Novartis's `synadam`.

The short version: for **synthetic data that stays in a trusted environment**, a
resampling method is the right tool, and AVATAR is the trajectory-level version
of exactly what `synadam` already does column by column. Formal differential
privacy is the right tool only when the output **crosses a trust boundary**, and
it is retained here for that case.

---

## 1. The two methods in one paragraph each

**AVATAR blending** (the `synthesize_pmx()` engine). For each generated subject,
sample a compatible source subject's event skeleton as a template, then fill its
covariates and endpoint trajectories with a distance-weighted blend of a handful
of similar source subjects, plus subject-level and within-trajectory noise. The
output looks like real data because it is built from real data. There is **no
formal privacy guarantee**.

**Differential privacy** (the `fit_calibrated_pmx()` and `fit_private_pmx()`
engines). Compute a small number of aggregate statistics, add mathematically
calibrated noise so that no single subject can move any released number by more
than a bounded amount, and generate from those noised aggregates against public
structural priors. The output carries a provable `(epsilon, delta)` guarantee,
at the cost of utility that degrades sharply as the cohort shrinks.

---

## 2. What synadam does, and why it matters here

`synadam` is Novartis's own ADaM synthetic-data package. Its `simulate_vec`
generates each column independently:

- **continuous:** `runif(n, min = observed_min, max = observed_max)` — a uniform
  draw over the observed range;
- **categorical:** `sample(unique_values, replace = TRUE)` — resample the
  observed categories;
- **flags:** sampled in the observed proportions.

There is no differential privacy, no epsilon. The method preserves each column's
marginal support and nothing else, and relies on the surrounding governance —
ADaM data is already analysis-level de-identified, and it stays inside a
controlled environment — rather than on a mathematical guarantee.

This is a completely standard and defensible position. Most real clinical-data
sharing runs on governance and heuristic de-identification (HIPAA Safe Harbor,
Expert Determination), not on differential privacy. `synadam` shipping this way
is evidence that the governance-based model is accepted practice.

**The key observation:** `synadam` resamples each *column* from the data.
AVATAR resamples each *subject trajectory* from the data. They are the same idea
at two different granularities. If `synadam`'s privacy model is acceptable for
its use, AVATAR's is acceptable for the same use — with one important caveat in
section 4.

---

## 3. Why a resampling method is the right default for this package

The package's stated purpose (`design/PROTOTYPE_SPEC.md` section 1) is **synthetic
data for model-workflow exploration**: exercising cleaning, joins, reshaping,
plotting, control-file plumbing, and repeated-dose or longitudinal analysis
code. The accuracy bar is "vaguely right, structurally exact." It is explicitly
*not* parameter estimation, inference, or scientific conclusions.

For that purpose, AVATAR is simply better than the DP engine on every axis the
user cares about:

- **It works at any cohort size.** The DP engine's utility collapses below a few
  hundred subjects (`design/FEASIBILITY.md`), which is most early-phase work.
  AVATAR produces plausible data from twelve subjects.
- **It preserves joint structure for free.** Real covariate correlations, real
  trajectory shapes, real timing patterns come through because whole real
  subjects come through (blended). The DP engine has to assert all of that from
  public structural models.
- **It needs no elicitation.** No public priors, no structural model, no
  protocol declaration. Point it at the data and go.
- **It spends no budget**, because there is no budget.

The DP engine, by contrast, spends most of its effort defending against an
adversary who wants to re-identify a patient from the output. If no such
adversary can reach the output — because it never leaves the safe environment —
that effort buys nothing.

---

## 4. The asymmetry that makes AVATAR riskier than synadam

Honesty requires stating where the `synadam` analogy strains.

A resampled **covariate value** is weakly identifying. A weight of 72 kg is
shared by thousands of people; releasing it reveals almost nothing about any
individual. This is why `synadam`'s per-column resampling is low-risk, and why
`synpmx`'s bootstrap covariates (`pmx_covariates_auto()`) are a reasonable
default.

A resampled **subject trajectory** is strongly identifying. A full
concentration-time profile — with its particular sampling times, missed visits,
noise pattern, and idiosyncratic shape — is close to a fingerprint. Blending
several donors and adding noise mitigates this, but the mitigation is not free
and not formal: push the noise high enough to defeat a nearest-neighbor linkage
attack and you have destroyed the same signal the DP engine would have
destroyed, without the accounting to prove it. This is the exact failure mode
that made the Version 1 design uncomfortable in the first place
(`design/FEASIBILITY.md` section 2).

So AVATAR sits a notch higher on the risk ladder than `synadam`'s column
resampling, for the same governance model. The conclusion is not that AVATAR is
unsafe; it is that **AVATAR's safety depends more heavily on the governance
context**, and the more identifying the trajectory, the more that dependence
matters.

---

## 5. The decision rule

The choice between the two engines is not about which is "more private" in the
abstract. It is about a single question:

> **Does the generated data cross a trust boundary?**

- **Stays inside** the safe computing environment — you are the only consumer,
  it never leaves, governance and access controls apply → **AVATAR.** It is more
  useful, works at any N, and its lack of a formal guarantee costs nothing
  because there is no adversary to guarantee against. This is `synadam`'s
  situation, and the reasoning is the same.
- **Crosses out** — shared with a partner, a vendor, a less-controlled system,
  or published → **the DP engine.** A formal guarantee is the only thing that
  survives a determined adversary, and here there might be one.

Differential privacy is expensive precisely because it defends against someone
who *wants* to break it. Buying that defense when no one can touch the output is
paying for a threat you do not have. Refusing to buy it when the output will be
handed to strangers is negligence. The engines exist side by side so the user
can match the tool to the boundary.

---

## 6. What the package keeps, and why

Both engines remain, deliberately.

- **`synthesize_pmx()`** — AVATAR, the primary and default method. The right
  answer for trusted-environment synthetic data, which is the common case.
- **`fit_calibrated_pmx()`** — the structural-correction DP engine. The right
  answer when a formal guarantee is needed and the cohort is small; it asserts
  shape from a public model and privately calibrates only the magnitude.
- **`fit_private_pmx()`** — the dense-grid DP engine. Retained for large pooled
  corpora where its cost is affordable.

Keeping all three is not indecision. It reflects that "synthetic clinical data"
is not one problem with one right method — it is a family of problems separated
by the trust boundary and the cohort size, and the honest package offers the
tool that fits each, clearly labeled with what it does and does not guarantee.

---

## 7. Reading guide

- `design/PROTOTYPE_SPEC.md` — the specification, with Version 4 (AVATAR) at the
  top and the DP versions retained below as alternatives.
- `design/FEASIBILITY.md` — the measurements behind "DP utility collapses at
  small N" and the Version 1 re-identification analysis.
- `design/PRIVACY_BACKGROUND.md` — how the DP arithmetic (`d`, `f`, epsilon)
  works, for when the DP engine is the right choice.
