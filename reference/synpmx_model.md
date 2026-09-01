# Estimate a population model and generate a synthetic dataset from it

A single call for
[`synpmx_model_estimate()`](https://iamstein.github.io/synpmx/reference/synpmx_model_estimate.md)
followed by
[`synpmx_model_generate()`](https://iamstein.github.io/synpmx/reference/synpmx_model_generate.md).
Use the two separately to look at the fit before generating from it; it
is on the result either way, as the `pmx_fitted_model` attribute.

## Usage

``` r
synpmx_model(data, roles, n_subjects = NULL, seed = NULL, ...)
```

## Arguments

- data:

  Source PMX event data.

- roles:

  Explicit column roles from
  [`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md),
  including `nominal_time`.

- n_subjects:

  Number of synthetic subjects. Defaults to the source count.

- seed:

  Seed, used for both stages.

- ...:

  Passed to
  [`synpmx_model_estimate()`](https://iamstein.github.io/synpmx/reference/synpmx_model_estimate.md).

## Value

A data frame in the source's shape, carrying the fitted model as an
attribute.

## See also

[`synpmx_model_estimate()`](https://iamstein.github.io/synpmx/reference/synpmx_model_estimate.md),
[`synpmx_model_generate()`](https://iamstein.github.io/synpmx/reference/synpmx_model_generate.md),
[`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md),
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md).
