# The visit model each arm was generated from

One row per arm, endpoint and modelled nominal time, giving the
probability that a generated subject in that arm has an observation
there. It is the share of patients who did, pooled over the arms that
have the visit, so every such arm shows the same probability; an arm has
a visit where at least `min_arm_patients` of its patients were observed
there. Attendance is drawn per visit rather than a real patient's visit
set being reused.

## Usage

``` r
pca_visits(x)
```

## Arguments

- x:

  A dataset from
  [`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md),
  or its trial summary.

## Value

A data frame with `arm`, `endpoint`, `time`, `probability` and
`patients`, the last being how many patients across the study hold that
cell at all.

## See also

[`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md),
[`pca_dosing()`](https://iamstein.github.io/synpmx/reference/pca_dosing.md),
[`pca_report()`](https://iamstein.github.io/synpmx/reference/pca_report.md).

## Examples

``` r
data <- pmx_simulated_fixture(60)
roles <- pmx_roles(
  id = "ID", time = "TIME", nominal_time = "NTIME", dv = "DV", amt = "AMT",
  evid = "EVID", cmt = "CMT", dvid = "DVID", mdv = "MDV"
)
head(pca_visits(synpmx_pca(data, roles, seed = 1)))
#>   arm endpoint  time probability patients
#> 1 all       cp  0.25           1       60
#> 2 all       cp  1.00           1       60
#> 3 all       cp  2.00           1       60
#> 4 all       cp  6.00           1       60
#> 5 all       cp 12.25           1       60
#> 6 all       cp 13.00           1       60
```
