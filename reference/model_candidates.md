# The candidate models the selection was made from

Every candidate the design admitted, whether or not it converged, with
the AIC it was compared on. A candidate that failed keeps its reason, so
a search that came down to one survivor does not look like a search that
had one candidate.

## Usage

``` r
model_candidates(fitted_model)
```

## Arguments

- fitted_model:

  A `pmx_fitted_model` from
  [`synpmx_model_estimate()`](https://iamstein.github.io/synpmx/reference/synpmx_model_estimate.md).

## Value

A data frame with columns `model`, `converged`, `aic` and `note`.

## See also

[`model_report()`](https://iamstein.github.io/synpmx/reference/model_report.md),
[`model_parameters()`](https://iamstein.github.io/synpmx/reference/model_parameters.md).
