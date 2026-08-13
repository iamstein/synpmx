# The scorecard for one synthetic dataset

Runs the checks from
[`vignette("synthetic-data-checks")`](https://iamstein.github.io/synpmx/articles/synthetic-data-checks.md)
that can be run automatically, and returns them as one table with the
answer and whether that answer passes. The vignette is the reference for
what each check asks and why its pass criterion is what it is.

## Usage

``` r
pmx_scorecard(source, synthetic, roles, proximity = NULL)
```

## Arguments

- source:

  Source PMX data.

- synthetic:

  Generated synthetic PMX data from
  [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md),
  carrying its `"pmx_settings"` attribute.

- roles:

  Explicit roles from
  [`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md).

- proximity:

  An already-computed
  [`compare_pmx_proximity()`](https://iamstein.github.io/synpmx/reference/compare_pmx_proximity.md)
  result, to save recomputing it. Left `NULL` it is computed here, which
  is the slowest part of the scorecard.

## Value

A `pmx_scorecard` data frame with columns `check`, `question`, `reads`,
`result` and `verdict`, marked `"restricted_not_releasable"`.

## Details

The `reads` column decides where the table may go. Rows marked
`"source"` or `"both"` were computed from real patient data, so the
filled-in scorecard is itself restricted output and belongs in the
environment the source lives in. Only the `"synthetic"` and
`"run report"` rows can travel with the data.

A `"review"` verdict is not a soft `"pass"`. It marks a row where no
threshold would be honest, and it has to be read.

Three checks in the vignette are absent here because no function can
produce them: C2 (dose and observation ordering), the source-side
rare-level census behind B5, and the one that matters most – whether the
pipeline that will consume the real study runs unchanged against this
output.

## See also

[`compare_pmx()`](https://iamstein.github.io/synpmx/reference/compare_pmx.md),
[`pmx_masking_report()`](https://iamstein.github.io/synpmx/reference/pmx_masking_report.md),
[`vignette("synthetic-data-checks")`](https://iamstein.github.io/synpmx/articles/synthetic-data-checks.md).

## Examples

``` r
data <- pmx_simulated_fixture(30)
roles <- pmx_roles(
  id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
  cmt = "CMT", dvid = "DVID", covariates = "WT"
)
synthetic <- suppressWarnings(synpmx_avatar(data, roles, seed = 1))
#> synpmx_avatar(): dropped 9 undeclared column(s): NTIME, TAD, OCC, RATE, MDV, CENS, LIMIT, AGE, SEX.
#>   Declare a column in `keep` to carry it through verbatim.
pmx_scorecard(data, synthetic, roles)
#> Scorecard: see vignette("synthetic-data-checks") for what each asks
#> 
#>   A1    Synthetic table is a legal PMX dataset           synthetic   TRUE                     pass
#>   A2    Source is legal under the declared roles         source      TRUE                     pass
#>   A3    Every endpoint survived                          both        2 of 2                   pass
#>   A4    Cohort size survived                             both        30 -> 30                 pass
#>   A5    Observations per patient                         both        14 -> 14                 review
#>   A5    Doses per patient                                both        2 -> 2                   review
#>   B1a   Avatars wearing one real patient's visit set     run report  0                        pass
#>   B1b   Avatars wearing one real patient's dose schedule run report  0                        pass
#>   B2    Synthetic patients unusual within their stratum  synthetic   1 of 30                  review
#>   B3    Adversarial accuracy inside its null interval    both        0.767 in [0.248, 0.692]  FAIL
#>   B4a   Generated time vectors copying an exposed real one both        0                        pass
#>   B4b   Generated DV vectors copying an exposed real one both        0                        pass
#>   C4    Dose regimens represented                        both        1 of 1                   pass
#> 
#> 1 FAIL, 3 to review.
#> Rows reading `source` or `both` are restricted output.
```
