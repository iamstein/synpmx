# Generate a new PMX event dataset from a fitted private model

This function never reads source data. Repeated calls are
post-processing and do not alter or consume the fitted model's privacy
accounting.

## Usage

``` r
generate_pmx(private_model, n_subjects = NULL, seed = 123)
```

## Arguments

- private_model:

  A fitted model from
  [`fit_private_pmx()`](https://iamstein.github.io/synpmx/reference/fit_private_pmx.md).

- n_subjects:

  Optional positive number of new synthetic subjects. By default,
  generation uses the privacy-accounted subject-count release stored in
  the fitted model (the exact source count for the explicitly public
  fixture backend).

- seed:

  Ordinary reproducibility seed for post-processing generation.

## Value

An ordinary data frame or tibble with the declared public schema and a
lightweight `pmx_privacy` attribute.
