# Privacy checks on a fitted model

Five checks on what a fitted model releases, each with its pass
criterion. The first four ask whether the release holds anything about
one patient that it should not: a per-patient table, an identifier, a
single patient's value, or a frequency resting on too few patients. The
fifth asks how far any one patient moves a released estimate, from the
individual random effects for the population model, from each subject's
baseline for the PD shapes, and by leaving each patient out for the
covariate summaries.

## Usage

``` r
model_privacy_checks(fitted_model)
```

## Arguments

- fitted_model:

  A `pmx_fitted_model` from
  [`synpmx_model_estimate()`](https://iamstein.github.io/synpmx/reference/synpmx_model_estimate.md),
  or a release from
  [`model_release()`](https://iamstein.github.io/synpmx/reference/model_release.md).
  A release carries only the verdict of the fifth check; the full fit
  carries the table behind it.

## Value

A `pmx_privacy_checks` data frame with columns `check`, `question`,
`result`, `criterion` and `verdict`. On a full fit, the per-estimate
reading behind check P5 is the `influence` attribute: one row per
released estimate, with the released value, the largest move one patient
causes, its unit and its verdict. No row is about a patient.

## Details

None of this is a formal privacy guarantee. The fifth check describes
this study only, and an adversary who knows every other patient in it
can detect a move smaller than any threshold. The article [Privacy
protections in the PMX model
generator](https://iamstein.github.io/synpmx/articles/pmxmodel-privacy.html)
states what each check establishes and what it does not.

## See also

[`model_release()`](https://iamstein.github.io/synpmx/reference/model_release.md),
[`synpmx_model_estimate()`](https://iamstein.github.io/synpmx/reference/synpmx_model_estimate.md),
[`synpmx_scorecard()`](https://iamstein.github.io/synpmx/reference/synpmx_scorecard.md).
