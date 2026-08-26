# Summarize a PMX dataset and generate a synthetic one from the summary

A single call for
[`synpmx_pca_summarize()`](https://iamstein.github.io/synpmx/reference/synpmx_pca_summarize.md)
followed by
[`synpmx_pca_generate()`](https://iamstein.github.io/synpmx/reference/synpmx_pca_generate.md).
Use the two separately to look at what the summary contains before
generating from it; the model is on the result either way, as the
`pmx_pca_model` attribute.

## Usage

``` r
synpmx_pca(data, roles, n_subjects = NULL, seed = NULL, ...)
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

  Generation seed.

- ...:

  Passed to
  [`synpmx_pca_summarize()`](https://iamstein.github.io/synpmx/reference/synpmx_pca_summarize.md).

## Value

A data frame in the source's shape, carrying the model as an attribute.

## Details

Where
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
blends values from real neighbouring patients, this writes out no number
a patient measured. What it carries out of the source is a mean, a
scale, a set of principal-component loadings, one mean score vector per
arm, a residual covariance, and a dosing and visit model per arm.
[`pmx_pca_report()`](https://iamstein.github.io/synpmx/reference/pmx_pca_report.md)
inventories all of it.

`nominal_time` is required. No formal privacy claim is made.

## See also

[`synpmx_pca_summarize()`](https://iamstein.github.io/synpmx/reference/synpmx_pca_summarize.md),
[`synpmx_pca_generate()`](https://iamstein.github.io/synpmx/reference/synpmx_pca_generate.md),
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md),
[`pmx_pca_report()`](https://iamstein.github.io/synpmx/reference/pmx_pca_report.md).

## Examples

``` r
data <- pmx_simulated_fixture(60)
roles <- pmx_roles(
  id = "ID", time = "TIME", nominal_time = "NTIME", dv = "DV", amt = "AMT",
  evid = "EVID", cmt = "CMT", dvid = "DVID", mdv = "MDV"
)
synthetic <- synpmx_pca(data, roles, seed = 1)
nrow(synthetic) > 0
#> [1] TRUE
```
