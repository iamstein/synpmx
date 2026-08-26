# Component loadings over time, and the variance each component explains

The loadings are what makes a principal component readable. Plotted
against time rather than tabulated, a component that is flat and
positive is overall magnitude and one that crosses zero separates early
from late.

## Usage

``` r
pca_components(x)
```

## Arguments

- x:

  A dataset from
  [`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md),
  or the fit itself.

## Value

A data frame with one row per component and retained grid cell, carrying
`variance_explained` as an attribute.

## See also

[`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md),
[`pca_report()`](https://iamstein.github.io/synpmx/reference/pca_report.md).

## Examples

``` r
data <- pmx_simulated_fixture(60)
head(pca_components(synpmx_pca(data, pmx_generated_roles(), seed = 1)))
#>    feature          kind endpoint  time covariate patients component    loading
#> 1 dv_cp__1 endpoint_cell       cp  0.25      <NA>       60       PC1 -0.2672606
#> 2 dv_cp__2 endpoint_cell       cp  1.00      <NA>       60       PC1 -0.2672596
#> 3 dv_cp__3 endpoint_cell       cp  2.00      <NA>       60       PC1 -0.2672610
#> 4 dv_cp__4 endpoint_cell       cp  6.00      <NA>       60       PC1 -0.2672620
#> 5 dv_cp__5 endpoint_cell       cp 12.25      <NA>       60       PC1 -0.2672614
#> 6 dv_cp__6 endpoint_cell       cp 13.00      <NA>       60       PC1 -0.2672595
```
