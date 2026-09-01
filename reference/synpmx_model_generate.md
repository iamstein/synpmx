# Generate a synthetic PMX dataset from a fitted model

Draws new subjects from a
[`synpmx_model_estimate()`](https://iamstein.github.io/synpmx/reference/synpmx_model_estimate.md)
fit. This stage reads no patient data: its arguments are the model and a
subject count, so everything about the source that reaches the output
has already passed through the fit.

## Usage

``` r
synpmx_model_generate(fitted_model, n_subjects = NULL, seed = NULL)
```

## Arguments

- fitted_model:

  A `pmx_fitted_model` from
  [`synpmx_model_estimate()`](https://iamstein.github.io/synpmx/reference/synpmx_model_estimate.md).

- n_subjects:

  Number of synthetic subjects. Defaults to the source count.

- seed:

  Generation seed.

## Value

A data frame in the source's shape, carrying the fitted model as a
`pmx_fitted_model` attribute.

## Details

Per subject: an arm is assigned keeping the source arm shares,
covariates are drawn from the arm's covariate model, random effects from
the between-subject covariance matrix, and the dose schedule from the
arm's dosing model. The concentration is then evaluated at the visits
drawn from the arm's visit model, against the schedule that was drawn,
so a reduced or skipped dose reaches the concentrations rather than
appearing only in the dosing records.

## See also

[`synpmx_model_estimate()`](https://iamstein.github.io/synpmx/reference/synpmx_model_estimate.md),
[`synpmx_model()`](https://iamstein.github.io/synpmx/reference/synpmx_model.md),
[`model_report()`](https://iamstein.github.io/synpmx/reference/model_report.md).
