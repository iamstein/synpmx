# What the PCA fit read out of the source data

One row per released quantity: what it is, how many numbers it holds,
and the smallest number of patients standing behind any one of them.
That last column is where disclosure risk sits. A grid cell is backed by
the patients who reached it, a covariate mean by the whole cohort, and a
covariate level by the patients who hold it, at least
`min_category_patients`.

## Usage

``` r
pca_report(x)
```

## Arguments

- x:

  A dataset from
  [`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md),
  or the fit itself.

## Value

A `pca_report` data frame.

## See also

[`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md),
[`pca_components()`](https://iamstein.github.io/synpmx/reference/pca_components.md).

## Examples

``` r
data <- pmx_simulated_fixture(60)
pca_report(synpmx_pca(data, pmx_generated_roles(), seed = 1))
#> What the PCA fit read out of the source data
#> 
#>   subjects: 60  components retained: 1 
#> 
#>             quantity
#>           visit grid
#>      feature centers
#>       feature scales
#>             loadings
#>          score means
#>     score covariance
#>  endpoint transforms
#>         assay limits
#>         dosing model
#>          visit model
#>        arm constants
#>                                                             what numbers
#>                             Nominal times modelled, per endpoint      14
#>                             Mean of each grid cell and covariate      14
#>                                   Standard deviation of the same      14
#>                               Component loadings on each feature      14
#>                                       Mean score vector, per arm       1
#>                           Residual covariance between components       1
#>                                    Log or identity, per endpoint       2
#>                                 Censoring boundary, per endpoint       0
#>  Planned cycles per arm; the dose ladder and three rates, pooled       8
#>     One attendance rate per endpoint, at the visits each arm has      14
#>                       Strata and kept columns, one value per arm       0
#>  min_patients
#>            60
#>            60
#>            60
#>            60
#>            60
#>            60
#>            60
#>            60
#>            60
#>            60
#>            60
```
