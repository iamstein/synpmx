# What the PCA fit read out of the source data

One row per released quantity: what it is, how many numbers it holds,
and the smallest number of patients standing behind any one of them.
That last column is where disclosure risk sits. A grid cell or a
covariate mean is backed by the whole cohort, while a rare covariate
level can be backed by a single patient.

## Usage

``` r
pmx_pca_report(x)
```

## Arguments

- x:

  A dataset from
  [`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md),
  or the fit itself.

## Value

A `pmx_pca_report` data frame.

## See also

[`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md),
[`pmx_pca_components()`](https://iamstein.github.io/synpmx/reference/pmx_pca_components.md).

## Examples

``` r
data <- pmx_simulated_fixture(60)
pmx_pca_report(synpmx_pca(data, pmx_generated_roles(), seed = 1))
#> What the PCA fit read out of the source data
#> 
#>   subjects: 60  components retained: 1 
#> 
#>             quantity                                               what numbers
#>           visit grid               Nominal times modelled, per endpoint      14
#>      feature centers               Mean of each grid cell and covariate      14
#>       feature scales                     Standard deviation of the same      14
#>             loadings                 Component loadings on each feature      14
#>          score means                         Mean score vector, per arm       1
#>     score covariance             Residual covariance between components       1
#>  endpoint transforms                      Log or identity, per endpoint       2
#>         assay limits                   Censoring boundary, per endpoint       0
#>         dosing model             Dose times and amounts each arm shares       4
#>          visit model Probability of a visit, per arm, endpoint and time      14
#>        arm constants         Strata and kept columns, one value per arm       0
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
