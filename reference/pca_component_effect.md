# What each component does to a profile, on the scale the data are reported in

A loading is a number on the standardized modelling scale, so a table of
them says which visits a component is concentrated on and not what
moving along it looks like. This runs the model's own inversion instead:
from the cohort centre, move one component by a given number of score
standard deviations, invert the standardization and the endpoint's log
transform, and report the feature values that come back.

## Usage

``` r
pca_component_effect(x, sds = 1)
```

## Arguments

- x:

  A dataset from
  [`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md),
  or its trial summary.

- sds:

  How far to move along each component, in score standard deviations.
  Defaults to 1. The spread used is the pooled residual standard
  deviation across arms, which is the scale
  [`pca_scores()`](https://iamstein.github.io/synpmx/reference/pca_scores.md)
  reports per arm.

## Value

A data frame with one row per component, feature and displacement:
`component`, `feature`, `kind`, `endpoint`, `time`, `covariate`,
`level`, `score_sd` (`-sds`, `0` or `+sds`), and `value` on the reported
scale. Marked `"restricted_not_releasable"`: the centre is a mean of
real data.

## Details

The result is the component in the units the study measured. A component
that raises the whole concentration profile shows as a curve shifted up
at every time; one that separates early visits from late shows as a
curve that crosses the centre. Read it as a description of the fitted
basis, not as pharmacology: no clearance was estimated and a score
standard deviation is a variance decomposition, so naming a mechanism
for a curve is a claim the fit does not support.

The reference is the cohort centre — score zero on every component —
which is the whole cohort's mean feature vector rather than any arm's.
Arms differ by their mean score vectors, which
[`pca_scores()`](https://iamstein.github.io/synpmx/reference/pca_scores.md)
reports.

## See also

[`pca_components()`](https://iamstein.github.io/synpmx/reference/pca_components.md)
for the loadings themselves,
[`pca_features()`](https://iamstein.github.io/synpmx/reference/pca_features.md)
for the grid,
[`pca_scores()`](https://iamstein.github.io/synpmx/reference/pca_scores.md)
for the per-arm score model.

## Examples

``` r
data <- pmx_simulated_fixture(60)
roles <- pmx_roles(
  id = "ID", time = "TIME", nominal_time = "NTIME", dv = "DV", amt = "AMT",
  evid = "EVID", cmt = "CMT", dvid = "DVID", mdv = "MDV"
)
effect <- pca_component_effect(synpmx_pca_summarize(data, roles))
head(effect)
#>    feature          kind endpoint  time covariate level component score_sd
#> 1 dv_cp__1 endpoint_cell       cp  0.25      <NA>  <NA>       PC1       -1
#> 2 dv_cp__2 endpoint_cell       cp  1.00      <NA>  <NA>       PC1       -1
#> 3 dv_cp__3 endpoint_cell       cp  2.00      <NA>  <NA>       PC1       -1
#> 4 dv_cp__4 endpoint_cell       cp  6.00      <NA>  <NA>       PC1       -1
#> 5 dv_cp__5 endpoint_cell       cp 12.25      <NA>  <NA>       PC1       -1
#> 6 dv_cp__6 endpoint_cell       cp 13.00      <NA>  <NA>       PC1       -1
#>      value
#> 1 1.085682
#> 2 8.686082
#> 3 4.342978
#> 4 1.628557
#> 5 1.302831
#> 6 9.228971
```
