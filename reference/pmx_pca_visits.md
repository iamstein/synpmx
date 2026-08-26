# The visit model each arm was generated from

One row per arm, endpoint and modelled nominal time, giving the
probability that a generated subject in that arm has an observation
there. It is the fraction of the arm's patients who did, so attendance
is drawn per visit rather than a real patient's visit set being reused.

## Usage

``` r
pmx_pca_visits(x)
```

## Arguments

- x:

  A dataset from
  [`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md),
  or its model.

## Value

A data frame with `arm`, `endpoint`, `time`, `probability` and
`patients`, the last being how many patients across the study hold that
cell at all.

## See also

[`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md),
[`pmx_pca_dosing()`](https://iamstein.github.io/synpmx/reference/pmx_pca_dosing.md),
[`pmx_pca_report()`](https://iamstein.github.io/synpmx/reference/pmx_pca_report.md).

## Examples

``` r
data <- pmx_simulated_fixture(60)
roles <- pmx_roles(
  id = "ID", time = "TIME", nominal_time = "NTIME", dv = "DV", amt = "AMT",
  evid = "EVID", cmt = "CMT", dvid = "DVID", mdv = "MDV"
)
head(pmx_pca_visits(synpmx_pca(data, roles, seed = 1)))
#>   arm endpoint  time probability patients
#> 1 all       cp  0.25           1       60
#> 2 all       cp  1.00           1       60
#> 3 all       cp  2.00           1       60
#> 4 all       cp  6.00           1       60
#> 5 all       cp 12.25           1       60
#> 6 all       cp 13.00           1       60
```
