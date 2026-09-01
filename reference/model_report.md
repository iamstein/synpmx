# What a fitted model carries

An inventory of everything in a `pmx_fitted_model`, in two halves: what
`nlmixr2` estimated, and the dosing, visit and covariate models that are
summaries of the source rather than estimates. Nothing here is
per-subject.

## Usage

``` r
model_report(fitted_model)
```

## Arguments

- fitted_model:

  A `pmx_fitted_model` from
  [`synpmx_model_estimate()`](https://iamstein.github.io/synpmx/reference/synpmx_model_estimate.md).

## Value

A `pmx_model_report` list, printed as sections.

## See also

[`model_candidates()`](https://iamstein.github.io/synpmx/reference/model_candidates.md),
[`model_parameters()`](https://iamstein.github.io/synpmx/reference/model_parameters.md),
[`synpmx_model_estimate()`](https://iamstein.github.io/synpmx/reference/synpmx_model_estimate.md).
