# The scorecard for one synthetic dataset

Runs the checks from
[`vignette("scorecard-synthetic-data-checks")`](https://iamstein.github.io/synpmx/articles/scorecard-synthetic-data-checks.md)
that can be run automatically, and returns them as one table with the
answer, whether that answer passes, and the call that explains it. The
vignette is the reference for what each check asks and why its pass
criterion is what it is.

## Usage

``` r
synpmx_scorecard(source, synthetic, roles, proximity = NULL)
```

## Arguments

- source:

  Source PMX data.

- synthetic:

  Generated synthetic PMX data. From
  [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
  it carries a `"pmx_settings"` attribute and the whole card can be
  filled in; without one, the three rows that need it read
  `"unavailable"`.

- roles:

  Explicit roles from
  [`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md).

- proximity:

  An already-computed
  [`compare_pmx_proximity()`](https://iamstein.github.io/synpmx/reference/compare_pmx_proximity.md)
  result, to save recomputing it. Left `NULL` it is computed here, which
  is the slowest part of the scorecard.

## Value

A `synpmx_scorecard` data frame with columns `check`, `question`,
`reads`, `result`, `verdict` and `explore`, marked
`"restricted_not_releasable"`.

## Details

The `reads` column decides where the table may go. Rows marked
`"source"` or `"both"` were computed from real patient data, so the
filled-in scorecard is itself restricted output and belongs in the
environment the source lives in. Only the `"synthetic"` and
`"run settings"` rows can travel with the data. `"run settings"` means
the value is the generation run's own record of what it did –
`attr(synthetic, "pmx_settings")` – rather than a measurement taken from
either table.

## Scoring a table this package did not generate

Everything measured from the two tables is measurable on any synthetic
dataset, whatever produced it, so a table carrying no `"pmx_settings"`
attribute is scored rather than refused: the three `"run settings"` rows
(B1a, B1b, C4) come back with the verdict `"unavailable"` and the rest
of the card is computed as usual. That covers another method's output,
and this package's own output read back from a file.

Those three rows cannot be recomputed from the finished table, and that
is a property of the generator rather than an omission here. Generated
times are the coarsened visit grid plus resampled deviations – applied
to dose rows too – so an avatar's schedule no longer matches any source
patient's key exactly, and matching it back by snapping to the grid
reports schedules that were never given. The run measures both
guarantees before the deviations are applied.
[`unmaskable_strata()`](https://iamstein.github.io/synpmx/reference/unmaskable_strata.md)
is the part of that question answerable without any run record: it reads
the source alone and names the arms whose patients no method could mask.

The `explore` column names the call to run when a row needs explaining.
The calls are written against this function's argument names (`source`,
`synthetic`, `roles`), so rename them to whatever the session calls
those objects. Every row carries one; printing lists them under the
table, for the rows that did not pass, rather than as a sixth column.

Printing and knitting differ on purpose.
[`print()`](https://rdrr.io/r/base/print.html) is a console layout: the
verdict table, then the calls to run, then the B5b levels. Knitting a
chunk that returns this object emits
[`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html) tables
instead, so a `.Rmd` or `.qmd` gets the whole card including the
`explore` column. Running a chunk interactively in an IDE shows the
console form, since nothing is knitting.

A `"review"` verdict is not a soft `"pass"`. It marks a row where no
threshold would be honest, and it has to be read. Nor is
`"unavailable"`: it marks a row nothing was measured for.

`"FAIL"` is reserved for the rows where the answer is always a defect:
the output is not a legal dataset (A1), it is not the study that went in
(A3, A6), or it reproduces one real patient's structure verbatim (B1a,
B1b, B4a, B4b). Everything else is `"review"`, including every row whose
answer can move for a legitimate reason – a subject dropped for want of
donors, a cohort statistic at a small sample size, a source a validator
objects to.

Two checks in the vignette are absent here because no function can
produce them: C2 (dose and observation ordering), and the one that
matters most – whether the pipeline that will consume the real study
runs unchanged against this output.

## See also

[`compare_pmx()`](https://iamstein.github.io/synpmx/reference/compare_pmx.md),
[`pmx_masking_report()`](https://iamstein.github.io/synpmx/reference/pmx_masking_report.md),
[`pmx_endpoint_types()`](https://iamstein.github.io/synpmx/reference/pmx_endpoint_types.md),
[`vignette("scorecard-synthetic-data-checks")`](https://iamstein.github.io/synpmx/articles/scorecard-synthetic-data-checks.md).

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
synpmx_scorecard(data, synthetic, roles)
#> Scorecard: see vignette("scorecard-synthetic-data-checks") for what each asks
#> 
#>   check question                                           reads        result                  verdict
#>   A1    Synthetic table is a legal PMX dataset             synthetic    TRUE                    pass
#>   A2    Source is legal under the declared roles           source       TRUE                    pass
#>   A3    Every endpoint survived                            both         2 of 2                  pass
#>   A4    Cohort size survived                               both         30 -> 30                pass
#>   A5    Observations per patient                           both         14 -> 14                review
#>   A5    Doses per patient                                  both         2 -> 2                  review
#>   B1a   Avatars with a visit set nobody else shares        run settings 0                       pass
#>   B1b   Avatars with a dose schedule nobody else shares    run settings 0                       pass
#>   B2    Synthetic patients unusual within their stratum    synthetic    1 of 30                 review
#>   B3    Adversarial accuracy inside its null interval      both         0.767 in [0.248, 0.692] review
#>   B4a   Generated time vectors copying an exposed real one both         0                       pass
#>   B4b   Generated DV vectors copying an exposed real one   both         0                       pass
#>   C4    Distinct dose-time schedules represented           run settings 1 of 1                  pass
#> 
#> To explore, with `source`, `synthetic` and `roles` named as you have them:
#>   A5    compare_pmx_distributions(source, synthetic, roles)
#>   A5    pmx_masking_report(synthetic, source, roles, section = "dose_schedules")
#>   B2    flag_identifiable_subjects(synthetic, roles)
#>   B3    compare_pmx_proximity(source, synthetic, roles)
#> 
#> no failures, 4 to review.
#> `run settings` rows come from the run's own record, `attr(synthetic, "pmx_settings")`.
#> Rows reading `source` or `both` are restricted output.
```
