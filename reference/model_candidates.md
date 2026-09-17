# The candidate models the selection was made from

Every candidate actually attempted, in attempt order, with numerical
convergence, acceptance after parameter and generation checks, and
reasons for rejection or warnings. An unneeded fallback has no row.

## Usage

``` r
model_candidates(fitted_model)
```

## Arguments

- fitted_model:

  A `pmx_fitted_model` from
  [`synpmx_model_estimate()`](https://iamstein.github.io/synpmx/reference/synpmx_model_estimate.md).

## Value

A data frame with columns `model`, `converged`, `accepted`, `aic`,
`seconds` (fitting and acceptance checks), and `note`. AIC may be
missing for an accepted stochastic approximation
expectation-maximization (SAEM) fit that was not compared with another
model. These are PK candidates; a PD-only fit returns an empty table. PD
candidates are in
`model_report(fitted_model)$pd[[endpoint]]$candidates`.

## See also

[`model_report()`](https://iamstein.github.io/synpmx/reference/model_report.md),
[`model_parameters()`](https://iamstein.github.io/synpmx/reference/model_parameters.md).
