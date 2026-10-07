# Privacy checks on a fitted model

Five checks on what a fitted model releases, each with its pass
criterion. The first three ask whether the fingerprint holds anything
about one patient that it should not: a per-patient table, an
identifier, or a single patient's value. The fourth recounts the
smallest group of patients behind each kind of released frequency – arm
sizes, attendance, dose-change rates, categorical levels, carried
values, and the levels the schema keeps for factor columns and binary or
ordinal endpoints – on both sides of it, since "one patient missed this
visit" discloses as much as "one patient came". The fifth asks how far
any one patient moves a released estimate, from the individual random
effects for the population model, from each subject's baseline for the
PD shapes, and by leaving each patient out for the covariate summaries;
by default estimation has already left out any patient who moved one by
15 or more.

## Usage

``` r
model_privacy_checks(fitted_model)
```

## Arguments

- fitted_model:

  A `pmx_fitted_model` from
  [`synpmx_model_estimate()`](https://iamstein.github.io/synpmx/reference/synpmx_model_estimate.md),
  or a fingerprint from
  [`model_fingerprint()`](https://iamstein.github.io/synpmx/reference/model_fingerprint.md).
  A fingerprint carries only the verdict of the fifth check; the full
  fit carries the table behind it.

## Value

A `pmx_privacy_checks` data frame with columns `check`, `question`,
`result`, `criterion` and `verdict`. On a full fit, the reading behind
P5 is the `influence` attribute, one row per released estimate with the
released value, the largest move one patient causes, its unit and its
verdict; and the recount behind P4 is the `frequencies` attribute, one
row per kind of released frequency with the smallest group behind it,
its floor and how many values were changed to meet it. No row is about a
patient.

## Details

None of this is a formal privacy guarantee. The fifth check describes
this study only, and an adversary who knows every other patient in it
can detect a move smaller than any threshold. The article [Privacy
protections in the PMX model
generator](https://iamstein.github.io/synpmx/articles/pmxmodel-privacy.html)
states what each check establishes and what it does not.

## See also

[`model_fingerprint()`](https://iamstein.github.io/synpmx/reference/model_fingerprint.md),
[`synpmx_model_estimate()`](https://iamstein.github.io/synpmx/reference/synpmx_model_estimate.md),
[`synpmx_scorecard()`](https://iamstein.github.io/synpmx/reference/synpmx_scorecard.md).
