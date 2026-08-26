# The dose-modification rates each arm was generated from

One row per arm, giving the three discrete-time hazards that turn a
planned schedule into the schedule a patient actually received, and the
dose ladder reductions move down.

## Usage

``` r
pca_dose_rates(x)
```

## Arguments

- x:

  A dataset from
  [`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md),
  or its trial summary.

## Value

A data frame with `arm`, `planned_cycles`, `levels`, `discontinuation`,
`interruption`, `reduction`, `patients`, `source_doses` and `distinct`.

## Details

A study where nobody reduces, skips or stops early has all three rates
at zero and a single level, and every generated patient then receives
the planned schedule exactly. A study with dose modifications — oncology
being the usual case — has non-zero rates, and generated patients differ
from one another in the same way and to the same degree the source
patients did.

`source_doses` and `distinct` describe what the arm actually contained,
so a generated dataset can be checked against them: `source_doses` is
the mean number of dosing events per patient, and `distinct` is how many
different schedules the arm held.

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
pca_dose_rates(synpmx_pca(data, roles, seed = 1))
#>   arm planned_cycles levels discontinuation interruption reduction patients
#> 1 all              2      1               0            0         0       60
#>   source_doses distinct
#> 1            2       60
```
