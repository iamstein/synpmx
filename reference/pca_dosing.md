# The planned dose schedule each arm was generated from

One row per arm and cycle, giving the nominal time and the amount the
arm was planned to receive there. The planned amount is read per cycle,
so a protocol-prescribed escalation is part of the plan rather than a
departure from it.

## Usage

``` r
pca_dosing(x)
```

## Arguments

- x:

  A dataset from
  [`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md),
  or its trial summary.

## Value

A data frame with `arm`, `cycle`, `time` and `planned_amt`.

## Details

This is what every generated patient in the arm starts from.
[`pca_dose_rates()`](https://iamstein.github.io/synpmx/reference/pca_dose_rates.md)
gives the three hazards that then move them off it: reductions,
interruptions and discontinuation. On a study with no dose modifications
those rates are zero and this schedule is what every patient receives.

The grid stops at the last cycle `min_arm_patients` of the arm's
patients reached, so a treatment duration only one patient had cannot be
generated.

## See also

[`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md),
[`pca_dose_rates()`](https://iamstein.github.io/synpmx/reference/pca_dose_rates.md),
[`pca_visits()`](https://iamstein.github.io/synpmx/reference/pca_visits.md).

## Examples

``` r
data <- pmx_simulated_fixture(60)
roles <- pmx_roles(
  id = "ID", time = "TIME", nominal_time = "NTIME", dv = "DV", amt = "AMT",
  evid = "EVID", cmt = "CMT", dvid = "DVID", mdv = "MDV"
)
head(pca_dosing(synpmx_pca(data, roles, seed = 1)))
#>   arm cycle time planned_amt
#> 1 all     1    0    100.0885
#> 2 all     2   12    100.0885
```
