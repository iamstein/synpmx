# Generate a synthetic PMX dataset from a principal-component model

Draws new subjects from a
[`synpmx_pca_summarize()`](https://iamstein.github.io/synpmx/reference/synpmx_pca_summarize.md)
model. This stage reads no patient data: its arguments are the model and
a subject count, so everything the synthetic dataset is built from is
visible in the model itself.

## Usage

``` r
synpmx_pca_generate(model, n_subjects = NULL, seed = NULL)
```

## Arguments

- model:

  A `pmx_pca_model` from
  [`synpmx_pca_summarize()`](https://iamstein.github.io/synpmx/reference/synpmx_pca_summarize.md),
  or a dataset generated from one.

- n_subjects:

  Number of synthetic subjects. Defaults to the number the model was
  fitted on.

- seed:

  Generation seed.

## Value

A data frame in the source's shape, carrying the model as an attribute.

## Details

Each generated subject is assigned an arm, keeping each arm's share of
the cohort. Its scores are that arm's mean plus a fresh residual, its
dose schedule is the one the arm holds in common, and it attends each
visit with the frequency the arm attended it. No individual's schedule
and no individual's visit set exists in the model to be copied.

## See also

[`synpmx_pca_summarize()`](https://iamstein.github.io/synpmx/reference/synpmx_pca_summarize.md),
[`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md).

## Examples

``` r
data <- pmx_simulated_fixture(60)
roles <- pmx_roles(
  id = "ID", time = "TIME", nominal_time = "NTIME", dv = "DV", amt = "AMT",
  evid = "EVID", cmt = "CMT", dvid = "DVID", mdv = "MDV"
)
model <- synpmx_pca_summarize(data, roles)
synthetic <- synpmx_pca_generate(model, seed = 1)
nrow(synthetic) > 0
#> [1] TRUE
```
