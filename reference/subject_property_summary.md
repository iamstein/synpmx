# Summarize fitted subject properties and associated regimens

Subject properties are categorical treatment assignments or grouping
variables declared through
[`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md),
such as `ACTARM`, `TRT`, or a nominal dose group. Their released
probabilities and property-conditioned regimen summaries are
source-dependent and therefore part of the fitted model's privacy
accounting.

## Usage

``` r
subject_property_summary(private_model)
```

## Arguments

- private_model:

  A fitted model from
  [`fit_private_pmx()`](https://iamstein.github.io/synpmx/reference/fit_private_pmx.md).

## Value

A data frame with one row per declared property stratum. It has zero
rows when no subject properties were declared.
