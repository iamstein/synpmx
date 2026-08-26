# The dosing model each arm was generated from

One row per arm and dose, giving the time and the amount. This is the
schedule the arm holds in common: the modal set of times and amounts
among its patients, never an individual's. `share` is the fraction of
the arm holding it, so a schedule picked by a bare plurality is visible
as one, and `distinct` is how many schedules that arm actually
contained.

## Usage

``` r
pmx_pca_dosing(x)
```

## Arguments

- x:

  A dataset from
  [`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md),
  or its model.

## Value

A data frame with `arm`, `dose`, `time`, `amt`, `share`, `patients` and
`distinct`.

## Details

A study recording its dose times as actuals rather than as planned times
holds close to one schedule per patient. The generated data then holds
one per arm, which is the mechanism working and a real loss of variety.

## See also

[`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md),
[`pmx_pca_visits()`](https://iamstein.github.io/synpmx/reference/pmx_pca_visits.md),
[`pmx_pca_report()`](https://iamstein.github.io/synpmx/reference/pmx_pca_report.md).

## Examples

``` r
data <- pmx_simulated_fixture(60)
roles <- pmx_roles(
  id = "ID", time = "TIME", nominal_time = "NTIME", dv = "DV", amt = "AMT",
  evid = "EVID", cmt = "CMT", dvid = "DVID", mdv = "MDV"
)
pmx_pca_dosing(synpmx_pca(data, roles, seed = 1))
#>   arm dose time      amt      share patients distinct
#> 1 all    1    0 100.0885 0.01666667        1       60
#> 2 all    2   12 100.0885 0.01666667        1       60
```
