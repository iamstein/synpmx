# Privacy in synpmx

## Scope

`synpmx` generates synthetic pharmacometric (PMX) datasets — dosing and
measurement event tables of the kind a population pharmacokinetic (PK)
or pharmacodynamic analysis consumes. It offers two families of method
that protect patients in fundamentally different ways. This vignette
explains what each one actually is, what differential privacy does and
does not promise, how to choose between them, and — once a private mode
is chosen — how to choose an epsilon.

The companion vignettes are
[`vignette("synpmx-method")`](https://iamstein.github.io/synpmx/articles/synpmx-method.md)
(the four generation modes) and
[`vignette("synpmx-demo")`](https://iamstein.github.io/synpmx/articles/synpmx-demo.md)
(the practical workflow).

## The two families in one page

### AVATAR blending: synthetic data built out of real subjects

The **default** method is
[`synthesize_pmx()`](https://iamstein.github.io/synpmx/reference/synthesize_pmx.md),
AVATAR-style blending. “AVATAR” is a method name rather than an
initialism: it comes from the patient-centric *avatarization*
literature, in which each synthetic record (“avatar”) is constructed
from the local neighborhood of real records rather than from a fitted
parametric model. This package implements an AVATAR-*inspired*
adaptation for longitudinal event tables, not published AVATAR software.

Mechanically, for each synthetic subject it:

1.  samples a real source subject as an **anchor** and keeps that
    subject’s event skeleton — dosing rows, observation times,
    missing-value pattern — so the generated subject is structurally a
    real trial subject;
2.  finds the **k nearest compatible neighbors** of the anchor in a
    standardized, principal-component (PCA) profile space (the default
    is `k = 5`);
3.  fills covariates and endpoint trajectories with a randomized,
    distance-weighted **blend** of those donors, with no single donor
    allowed more than 80% of the weight; and
4.  adds subject-level and within-trajectory random noise.

The output looks like trial data because it is assembled from trial
data. That is the source of both its utility and its risk: it works at
any cohort size, preserves real covariate correlations and real
trajectory shapes for free, needs no elicited model or priors — and
makes **no formal privacy guarantee**. Its safety rests on governance:
the data stays inside a trusted, access-controlled environment. This is
the same footing as Novartis’s `synadam`, an ADaM (Analysis Data Model)
synthetic-data package that resamples each column marginally from the
real data with no privacy accounting at all.

### Differential privacy: a bound on what any one subject can change

The **alternative** is the differentially private (DP) engines,
[`fit_calibrated_pmx()`](https://iamstein.github.io/synpmx/reference/fit_calibrated_pmx.md)
and
[`fit_private_pmx()`](https://iamstein.github.io/synpmx/reference/fit_private_pmx.md).
These never copy or blend a subject’s data into the output. Instead
they:

1.  compute a small number of **aggregate** statistics from the source
    data, with every per-subject contribution clipped into a publicly
    declared range;
2.  add **calibrated random noise** to each aggregate, with the noise
    scale derived from how much one subject could have moved that
    number; and
3.  generate the synthetic dataset from the noised aggregates and a
    public structural model.

Only the noised aggregates ever touch the generated data, so the release
carries a mathematical guarantee.

## What differential privacy actually is

Differential privacy is a property of the **release procedure**, not of
the released dataset. You cannot inspect a table and check whether it is
differentially private; you can only check the mechanism that produced
it.

A randomized mechanism `M` is `(epsilon, delta)`-differentially private
if, for every pair of **neighboring** datasets `D` and `D'` — identical
except that one complete subject is added or removed — and every set of
possible outputs `S`:

``` math
\Pr[M(D) \in S] \;\le\; e^{\varepsilon} \, \Pr[M(D') \in S] \;+\; \delta
```

Read it as a promise made to one patient: *whatever the analyst
concludes from this release, they would have concluded almost the same
thing had you never enrolled.* The bound holds against an adversary with
unlimited computing power and arbitrary side information, including one
who already knows every other subject in the study.

Making that true requires two ingredients, both visible in the package:

- **Sensitivity** — the most that adding or removing one subject can
  move the released quantity. It exists only because per-subject values
  are clipped into a range declared *without looking at the data*. A
  data-derived range is itself a leak, and the accounting will not catch
  it.
- **A noise scale tied to that sensitivity.** The engines use the
  Laplace mechanism: noise of scale `sensitivity / epsilon`. Smaller
  epsilon means a stronger promise and more noise. Budget is spent, not
  reused: releasing `d` separate quantities splits epsilon `d` ways.

The [privacy background
article](https://iamstein.github.io/synpmx/articles/privacy-background.html)
works the arithmetic through with examples.

## AVATAR and DP are not two points on one privacy scale

It is tempting to rank the methods as “less private” and “more private.”
That is the wrong axis. The difference is one of **kind**:

|  | AVATAR blending | Differentially private engines |
|----|----|----|
| What produces the output | Real subject trajectories, blended and perturbed | Noised aggregate statistics plus a public structural model |
| Privacy claim | None formal; heuristic mitigation by blending and noise | Proven `(epsilon, delta)` bound on one subject’s influence |
| Holds against a determined adversary | Not established | Yes, by construction |
| Rests on | Governance and access control | Mathematics, plus correct declared ranges |
| Utility at small N | Good — works from a dozen subjects | Degrades sharply below a few hundred |
| Elicitation required | None | Public structural model, priors, and clipping ranges |
| Cost | None | Epsilon budget, and accuracy paid for it |

The honest caveat about AVATAR is one of granularity. A resampled
covariate value — a body weight of 72 kg — is weakly identifying,
because thousands of people share it. A resampled subject **trajectory**
is much closer to a fingerprint: its particular sampling times, missed
visits, noise pattern, and curve shape can be nearly unique. Blending
several donors and adding noise mitigates this, but not formally, and
pushing the noise high enough to defeat a nearest-neighbor linkage
attack would destroy the same signal a DP mechanism would have destroyed
— without the accounting to prove it. AVATAR therefore leans on the
governance context more heavily than column-wise resampling like
`synadam` does. This package has undergone no attack-based privacy
validation.

## The decision rule

The choice is a single question:

> **Does the generated data cross a trust boundary?**

- **Stays inside** the safe environment — you are the only consumer,
  governance and access controls apply → use **AVATAR**
  ([`synthesize_pmx()`](https://iamstein.github.io/synpmx/reference/synthesize_pmx.md)).
  Its lack of a formal guarantee costs nothing, because there is no
  adversary to guarantee against.
- **Crosses out** — shared with a partner or vendor, published, or moved
  to a less-controlled system → use a **DP** engine. A formal guarantee
  is the only thing that survives a determined adversary.

Differential privacy is expensive precisely because it defends against
someone who wants to break it. Buy it when the output will be handed to
strangers; skip it when nothing can reach the output.

## The differentially private engines

Both make a subject-level guarantee: neighboring datasets differing by
one complete subject produce nearly indistinguishable output, so no one
person’s participation can be inferred from the release.

- **[`fit_calibrated_pmx()`](https://iamstein.github.io/synpmx/reference/fit_calibrated_pmx.md)**
  (structural correction). Asserts curve shape from a public structural
  model and privately calibrates only the exposure magnitude by a small
  correction factor. Because it releases very few quantities, it remains
  viable at small cohorts. This is the recommended DP path for
  early-phase studies.
- **[`fit_private_pmx()`](https://iamstein.github.io/synpmx/reference/fit_private_pmx.md)**
  (dense grid). Releases a larger set of noised population summaries.
  Retained for large pooled corpora where its cost is affordable.

Both are built on Laplace releases, which spend no delta, so the
realized accounting reports `delta = 0` — a pure epsilon guarantee —
even though the requested `delta` is carried as slack in the contract.
Both fail closed when the validated OpenDP backend is unavailable, and
neither ever substitutes ordinary random noise for a calibrated
mechanism:

``` r

dp_backend_status()
#>   backend available version production
#> 1  OpenDP      TRUE  0.15.1       TRUE
```

## What the guarantee does and does not mean

- **Epsilon** is the one-person influence limit: smaller is stronger. It
  is not a re-identification probability.
- **Delta** is a small additive slack in the probability bound. It is
  not the fraction of unprotected patients.
- Differential privacy bounds the information attributable to one
  person’s participation. It does **not** establish legal anonymity,
  authorize release, or secure the environment.
- The guarantee is only as good as the declared clipping ranges. Ranges
  chosen by looking at the data break it silently.

## Choosing an epsilon

Everything below applies to the **differentially private** modes only.
If your synthetic data stays inside a trusted environment, use
[`synthesize_pmx()`](https://iamstein.github.io/synpmx/reference/synthesize_pmx.md)
and ignore epsilon entirely.

Epsilon is the one-person influence limit: smaller means stronger
privacy and more noise, larger means weaker privacy and less noise.
There is no universal default; an approved value must come from
governance and threat modeling, not from whichever number makes a plot
look best.

For the structural-correction engine
([`fit_calibrated_pmx()`](https://iamstein.github.io/synpmx/reference/fit_calibrated_pmx.md)),
the usable accuracy is captured by one quantity:

``` math
f \;=\; \frac{d}{\varepsilon N}
```

where `d` is the number of released quantities, `N` the number of
subjects, and `f` the fraction of the prior’s width that survives as
noise. The
[`pmx_preflight()`](https://iamstein.github.io/synpmx/reference/pmx_preflight.md)
helper reports it before any budget is spent:

``` r

priors <- pmx_priors(pk = pmx_prior(c(1 / 4, 4), source = "example"))
pmx_preflight(priors, epsilon = 0.5, n_subjects = 60)
#> Pre-flight: d = 2, epsilon = 0.5, N = 60  ->  f = 0.067
#>  quantity prior_fold          f expected_fold_error
#>        pk         16 0.06666667            1.203025
#> 
#> Verdict: consider a smaller epsilon
#> The prior contributes almost nothing. Consider spending less epsilon rather than banking accuracy you do not need.
```

The decision rule is not “is the error small” but “does the release beat
the prior”: when `f` is near or above 1 the release conveys nothing the
prior did not, and generating from public inputs alone is strictly
better.

### The measured frontier

Utility degrades sharply in small cohorts, and this is a property of
differential privacy, not of the implementation. Measured total
fold-error on clearance, for the structural engine with a
correction-factor prior:

|   N | epsilon 0.25 | epsilon 1 |
|----:|-------------:|----------:|
|   6 |     2.3-fold |  1.7-fold |
|  20 |     1.6-fold |  1.3-fold |
|  60 |     1.4-fold |  1.1-fold |
| 300 |    1.08-fold | 1.06-fold |

From about 60 subjects upward the error is limited by estimator bias
rather than by the privacy mechanism, so a smaller epsilon is the better
buy: it strengthens the guarantee at no real cost in accuracy.

Epsilon and delta are governance decisions, not defaults. For anything
public facing they should be set and justified by whoever owns the data,
and recorded: every fit carries a release ledger, and
[`privacy_report()`](https://iamstein.github.io/synpmx/reference/privacy_report.md)
prints the realized accounting.

## Where to read more

- [Privacy
  background](https://iamstein.github.io/synpmx/articles/privacy-background.html)
  — how the arithmetic works: `d`, `f`, epsilon, and the error law, with
  worked examples.
- [Mechanism-level privacy
  argument](https://iamstein.github.io/synpmx/articles/privacy-argument.html)
  — the formal argument, for a reviewer.
- [Feasibility by cohort
  size](https://iamstein.github.io/synpmx/articles/feasibility.html) —
  the complete measured frontier, and why small cohorts are hard for any
  formal method.
- [`vignette("synpmx-method")`](https://iamstein.github.io/synpmx/articles/synpmx-method.md)
  — the four generation modes and why AVATAR is the default.
- [AVATAR
  mathematics](https://iamstein.github.io/synpmx/articles/avatar-mathematics.html)
  — the default generator in detail.
