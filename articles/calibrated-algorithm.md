# The synpmx_calibrated Algorithm

[`synpmx_calibrated()`](https://iamstein.github.io/synpmx/reference/synpmx_calibrated.md)
keeps a public structural model’s *shape* and spends a small privacy
budget correcting only its *magnitude*. Each subject is reduced to one
bounded number — how wrong the model’s prediction was for them — those
numbers are averaged, and the average is released with calibrated noise.
Two numbers leave the study: that correction and a noisy subject count.

Releasing two numbers rather than dozens is the whole design. Budget
buys accuracy in proportion to the number of subjects and in inverse
proportion to how many quantities are released, so a mode that releases
two stays usable at cohort sizes where a mode that releases fifty does
not.

**Out of scope, stated once and bluntly.** Everything not calibrated is
*asserted*: curve shape, between-subject variability, residual error,
and every covariate come from the public model, so the output is only as
realistic as that model. This is not an estimation procedure and the
released correction is not a parameter estimate. See [what prior-only
generation cannot
express](https://iamstein.github.io/synpmx/articles/prior-algorithm.html),
because all of it applies here too.

## Maintenance status

A secondary, provided-as-is path.
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
is the primary maintained method. The differentially private modes are
complete and tested but not under active development, carry known open
findings, and have not been independently privacy-audited. Use them to
demonstrate the privacy and utility tradeoff, not as a production
release mechanism; a real regulated release needs specialist review. A
session-level acknowledgment is required before the function runs:

``` r

synpmx_enable_dp_engines()
#> DP engines enabled for this session: the differentially private engines are complete and tested, but not under active development, carry known open findings, and have not been independently privacy-audited. See https://iamstein.github.io/synpmx/articles/synpmx-privacy.html for the trust-boundary decision rule and what a production release additionally needs.
```

`backend = "public"` is exempt from that gate because it makes no
privacy claim — and it adds **no noise at all**, so it shows the
correction mechanism and nothing about epsilon. Every number below that
is meant to be private comes from `backend = "opendp"`.

## Overview of Algorithm

| Step | Operation | Reads data? | Spends budget? |
|----|----|----|----|
| 1 | Validate the public model, design and priors | No | No |
| 2 | Resolve the privacy backend, failing closed | No | No |
| 3 | Reduce each subject to one bounded correction | **Yes** | No |
| 4 | Count the budget’s divisor and split it | No | No |
| 5 | Release the noisy subject count | **Yes** | **Yes** |
| 6 | Release each correction, clipped to its prior | **Yes** | **Yes** |
| 7 | Apply the correction to the public parameters | No | No |
| 8 | Warn where the release conveys nothing | No | No |
| 9 | Generate, which is post-processing | No | No |

Only steps 5 and 6 spend anything. Everything after step 6 reads the
release rather than the study, which is why any number of datasets can
be drawn from one payment with
[`synpmx_generate()`](https://iamstein.github.io/synpmx/reference/synpmx_generate.md).

## Step 1: Validate the Public Inputs

The model and design are the same public objects prior-only generation
takes, with the same required `source` on each.
[`pmx_priors()`](https://iamstein.github.io/synpmx/reference/pmx_priors.md)
adds one \[pmx_prior()\] per released correction, and a prior is **the
dominant sensitivity term in the whole design**, so its provenance
matters more than any other input:

``` r

priors <- pmx_priors(pk = pmx_prior(c(1 / 4, 4),
                                    source = "scaling literature"))
priors$pk
#> prior [0.25, 4] (16-fold, span 2.77 log units)
#>   source: scaling literature
```

The range brackets how wrong the prediction is believed to be. A prior
on “how wrong is this prediction” is far tighter than any prior on the
parameter itself, and it is free: nobody has to know the clearance to be
confident the model is within four-fold.

The `span` printed above is the width in log units, and it is what the
noise is measured against. A wider prior is safer against clipping and
costs accuracy proportionally.

## Step 2: Resolve the Backend, Failing Closed

`backend = "opendp"` resolves the validated OpenDP adapter and **errors
if it is unavailable** rather than falling back to ordinary R random
noise.

``` r

dp_backend_status()
#>   backend available version production
#> 1  OpenDP      TRUE  0.15.1       TRUE
```

## Step 3: Reduce Each Subject to One Bounded Correction

For each subject, the area under their observed concentration curve is
compared with the area under the public model’s prediction for their own
dosing history:

    correction_i = AUC(predicted) / AUC(observed)

Clearance is inversely proportional to area, so this ratio *is* the
clearance correction, with no fitting involved. A subject needs at least
three observations to contribute one; fewer and they contribute nothing.

This is the step that makes the sensitivity bound provable. Each
subject’s number is computed **only from that subject’s own rows**, so
adding or removing one subject changes one term in a sum. Clipping that
term to the public prior bounds how much it can change, and the bound is
what the noise is calibrated against.

A PD correction works the same way where the model has a PD component,
scaling the whole curve so its shape is kept.

## Step 4: Count the Divisor and Split the Budget

The budget is divided equally among the quantities released:

    d = 1 (subject count) + one per prior + one per DP-declared covariate
    per-query epsilon = epsilon / d

That `d` is the second term in every feasibility statement, and it is
why this mode releases as little as it can. A bootstrap-resampled
covariate is drawn from the data outside the budget and does not enter
`d` — and is not differentially private, which the release says of
itself.

[`pmx_preflight()`](https://iamstein.github.io/synpmx/reference/pmx_preflight.md)
reports what a given `d`, epsilon and cohort size will buy, reading no
data and spending nothing:

``` r

pmx_preflight(priors, epsilon = 0.1, n_subjects = 12)
#> Pre-flight: d = 2, epsilon = 0.1, N = 12  ->  f = 1.667
#>  quantity prior_fold        f expected_fold_error
#>        pk         16 1.666667                   4
#> 
#> Verdict: worthless
#> The noise is as wide as the prior. This release would tell you nothing you did not already assume.
#> Use prior-mode generation instead, or raise epsilon only if governance allows.
```

The quantity it reports is

    f = d / (epsilon * N)

the fraction of each prior’s width that survives as noise, and the
expected fold-error is `exp(f * span)`, capped at the prior’s half-width
because clipping cannot put a release outside its prior. Use it before
spending, not after. Above `f` of roughly 0.25 the uncapped form is
increasingly pessimistic, which the [public-data
evaluation](https://iamstein.github.io/synpmx/articles/calibrated-public-data-examples.html)
measures against realized draws.

## Step 5 and 6: Release the Count and the Corrections

Both releases go through the accountant with sensitivity one. The
subject count is a count, so removing a subject changes it by one. Each
correction is a sum of per-subject terms clipped to the unit interval,
so removing a subject changes it by at most one.

The noise is Laplace at scale `sensitivity / per-query epsilon`, `delta`
is zero, and this is where the privacy guarantee actually lives. **The
noise is never user-seeded.** A `seed` argument controls ordinary
generation only; seeding the privacy noise would make the release
reproducible, which would defeat it. A consequence for anyone reading
output: two runs of the same call give different releases, by design.

## Step 7: Apply the Correction

The released mean is mapped back out of the unit interval through the
prior, exponentiated, and multiplied onto the public parameters.
Clearance takes the PK correction. A PD correction multiplies
`baseline`, `plateau` and `slope` together so the curve’s shape survives
and only its level moves.

## Step 8: Warn Where the Release Conveys Nothing

Two warnings matter, and both describe a release that is technically
valid and practically empty.

**The noise is as wide as the prior.** When `f >= 1` the release carries
no information beyond the assumption, and the function says so:
prior-only generation would give the same output at no privacy cost.

**The correction is pressed against its prior boundary.** When the
released mean lands within 2% of either end of the unit interval, the
prior was probably wrong and the release is censored. The generated data
then reflects the boundary rather than the study. At small cohorts and
small epsilon this is common rather than exceptional, which the
public-data evaluation quantifies.

## Step 9: Generate, Which Costs Nothing

Generation from the release is post-processing: it reads no confidential
data and consumes no further budget. It is the prior-only algorithm with
corrected parameters, so every step of [that
algorithm](https://iamstein.github.io/synpmx/articles/prior-algorithm.html)
applies, including the `roles` declaration naming the output columns and
the refusal of a mixed-route model.

Draw further datasets with
[`synpmx_generate()`](https://iamstein.github.io/synpmx/reference/synpmx_generate.md)
rather than by calling
[`synpmx_calibrated()`](https://iamstein.github.io/synpmx/reference/synpmx_calibrated.md)
again. A second fit spends the budget a second time, and the function
warns when it notices one.

## What Leaves the Source Data

Two noisy numbers per fit, and nothing else. No patient’s record, no
parameter estimate, no summary of the observations beyond the single
clipped average. The release travels with the generated data so the
accounting is always at hand:

``` r

attr(synthetic, "synpmx_release")$privacy$accounting
```

What the guarantee does **not** cover is the honesty of the public
inputs. A model whose `typical` values were read off a fit to the same
study makes the prior data-dependent, and no check in the package can
detect that. The `source` fields are where that discipline is recorded.

## Where to go next

- [`vignette("calibrated-demo")`](https://iamstein.github.io/synpmx/articles/calibrated-demo.md)
  — one study end to end, with the preflight before the spend.
- [Evaluating calibration on public
  data](https://iamstein.github.io/synpmx/articles/calibrated-public-data-examples.html)
  — what a defensible epsilon actually produces at these cohort sizes.
- [Feasibility by cohort
  size](https://iamstein.github.io/synpmx/articles/feasibility.html) —
  where `f` comes from and how it behaves.
- [What differential privacy does and does not
  guarantee](https://iamstein.github.io/synpmx/articles/synpmx-privacy.html)
  — the trust-boundary decision rule.
