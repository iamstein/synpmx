# A scorecard as a coloured HTML table

Displays a
[`synpmx_scorecard()`](https://iamstein.github.io/synpmx/reference/synpmx_scorecard.md)
as an interactive table with the verdicts coloured: `"FAIL"` in bold red
on a light red, `"review"` in bold orange on a light orange,
`"unavailable"` in muted grey, and `"pass"` left as ordinary text. A
card is five verdicts among thirty-odd rows of prose, and the rows that
need reading are the ones that have to be findable without reading all
of it.

## Usage

``` r
synpmx_scorecard_datatable(x, ...)
```

## Arguments

- x:

  A
  [`synpmx_scorecard()`](https://iamstein.github.io/synpmx/reference/synpmx_scorecard.md).

- ...:

  Passed to `DT::datatable()`. Paging is off and row numbers are
  suppressed by default, since the whole card is meant to be read at
  once and `check` already names each row.

## Value

An
[`htmltools::tagList`](https://rstudio.github.io/htmltools/reference/tagList.html)
holding the coloured card and the notes that knitting one carries.
Without `DT` installed, `x` invisibly, having printed it.

## Details

[`synpmx_scorecard()`](https://iamstein.github.io/synpmx/reference/synpmx_scorecard.md)
computes the card and this displays it, so the object is unchanged and
can be subset, saved or printed as usual. The colouring is the only
thing added.

What is emitted is what knitting the card itself emits, with the
colouring added: the card, then the B5b rare-level detail where a study
has any, then the reminder that D1 is a number and not a shape.

`DT` is a suggested package rather than a required one – `synpmx` has no
hard dependencies – so without it installed this says so and prints the
card in the console form instead. The verdicts are all still there; only
the colour is missing.

The same restriction applies as to the card itself. Rows reading
`"source"` or `"both"` were computed from real patient data, so an HTML
file holding the whole table belongs in the environment the source lives
in.

## See also

[`synpmx_scorecard()`](https://iamstein.github.io/synpmx/reference/synpmx_scorecard.md).

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
synpmx_scorecard_datatable(synpmx_scorecard(data, synthetic, roles))
#> DT is not installed, so the scorecard is printed uncoloured. Install DT for the coloured table.
#> Scorecard: see vignette("scorecard-synthetic-data-checks") for what each asks
#> 
#>   check question                                           reads        result                              verdict
#>   A1    Synthetic table is a legal PMX dataset             synthetic    TRUE                                pass
#>   A2    Source is legal under the declared roles           source       TRUE                                pass
#>   A3    Every endpoint survived                            both         2 of 2                              pass
#>   A4    Cohort size survived                               both         30 -> 30                            pass
#>   A5a   Observations per patient                           both         14 -> 14                            review
#>   A5b   Doses per patient                                  both         2 -> 2                              review
#>   A6    Discrete endpoints keeping their source scale      both         no discrete endpoint                pass
#>   B1a   Avatars with a visit set nobody else shares        run settings 0                                   pass
#>   B1b   Avatars with a dose schedule nobody else shares    run settings 0                                   pass
#>   B2    Synthetic patients unusual within their stratum    synthetic    1 of 30                             review
#>   B3    Adversarial accuracy inside its null interval      both         0.767 in [0.248, 0.692]             review
#>   B4a   Generated time vectors copying an exposed real one both         0                                   pass
#>   B4b   Generated DV vectors copying an exposed real one   both         0                                   pass
#>   B5a   Patients holding the least-held categorical level  synthetic    no categorical covariate or stratum pass
#>   B5b   Rare source levels copied into the output          both         no categorical covariate or stratum pass
#>   C1    Strata keeping their source size                   both         no strata declared                  pass
#>   C2    Distinct dose-time schedules represented           run settings 1 of 1                              pass
#>   D1    Values landing in the same range                   both         sd x1.4 on pd (furthest of 3)       review
#> 
#> To explore, with `source`, `synthetic` and `roles` named as you have them:
#>   A5a   compare_pmx_distributions(source, synthetic, roles)
#>   A5b   pmx_masking_report(synthetic, source, roles, section = "dose_schedules")
#>   B2    flag_identifiable_subjects(synthetic, roles)
#>   B3    compare_pmx_proximity(source, synthetic, roles)
#>   D1    compare_pmx_distributions(source, synthetic, roles)
#> 
#> D1 reports numbers, not shapes. Plot source and synthetic on the same axes
#> -- `DV` against time, and each covariate -- with whatever you normally use.
#> 
#> no failures, 5 to review.
#> `run settings` rows come from the run's own record, `attr(synthetic, "pmx_settings")`.
#> Rows reading `source` or `both` are restricted output.
```
