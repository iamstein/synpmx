# The estimated parameters

Fixed effects, the between-subject covariance matrix and the residual
error. Not estimates to report: see
[`synpmx_model_estimate()`](https://iamstein.github.io/synpmx/reference/synpmx_model_estimate.md).

## Usage

``` r
model_parameters(fitted_model)
```

## Arguments

- fitted_model:

  A `pmx_fitted_model` from
  [`synpmx_model_estimate()`](https://iamstein.github.io/synpmx/reference/synpmx_model_estimate.md).

## Value

A list with `fixed`, `omega` and `residual`.

## See also

[`model_report()`](https://iamstein.github.io/synpmx/reference/model_report.md),
[`model_candidates()`](https://iamstein.github.io/synpmx/reference/model_candidates.md).
