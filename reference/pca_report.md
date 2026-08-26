# What the PCA fit read out of the source data

One row per released quantity: what it is, how many numbers it holds,
and the smallest number of patients standing behind any one of them.
That last column is where disclosure risk sits. A grid cell or a
covariate mean is backed by the whole cohort, while a rare covariate
level can be backed by a single patient.

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
#>             quantity                                                      what
#>           visit grid                      Nominal times modelled, per endpoint
#>      feature centers                      Mean of each grid cell and covariate
#>       feature scales                            Standard deviation of the same
#>             loadings                        Component loadings on each feature
#>          score means                                Mean score vector, per arm
#>     score covariance                    Residual covariance between components
#>  endpoint transforms                             Log or identity, per endpoint
#>         assay limits                          Censoring boundary, per endpoint
#>         dosing model Planned cycles, the dose ladder, and three rates, per arm
#>          visit model        Probability of a visit, per arm, endpoint and time
#>        arm constants                Strata and kept columns, one value per arm
#>  numbers min_patients
#>       14           60
#>       14           60
#>       14           60
#>       14           60
#>        1           60
#>        1           60
#>        2           60
#>        0           60
#>        8           60
#>       14           60
#>        0           60
```
