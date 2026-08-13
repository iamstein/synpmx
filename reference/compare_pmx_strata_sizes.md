# Stratum sizes, source against synthetic

One row per declared stratum level, with the source size, the synthetic
size, and the size the generator was aiming for. Reported for each
`strata` column on its own, and – when more than one is declared – for
the joint cells too, because the joint cell is what the generator
balances on.

## Usage

``` r
compare_pmx_strata_sizes(source, synthetic, roles)
```

## Arguments

- source:

  Source PMX data.

- synthetic:

  Generated synthetic PMX data.

- roles:

  Explicit roles from
  [`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md).

## Value

A `pmx_strata_sizes` data frame with `column`, `level`,
`source_patients`, `synthetic_patients`, `expected` and `balanced`. Zero
rows when the roles declare no strata.

## Details

`expected` is the source share of the cohort carried to the synthetic
cohort size, which is the generator's own target rule. It matters
whenever `n_subjects` differs from the source size: every stratum then
changes size by design, and comparing the raw counts says "everything
moved" when nothing is wrong.

`balanced` is whether the cell had enough source patients (3) to be held
exactly. Cells below that floor are deliberately left to vary by seed –
reproducing the size of a two-patient arm on every run discloses that
size. A cell above the floor landing far from `expected` is not that
mechanism, and is worth chasing.

With more than one `strata` column, `balanced` is `NA` on the
single-column rows: the floor applies to the joint cell, so a column
with 19 patients can still be split into joint cells that all fall below
it.

This reads real patient data on both sides and is marked
`"restricted_not_releasable"`.

## See also

[`synpmx_scorecard()`](https://iamstein.github.io/synpmx/reference/synpmx_scorecard.md),
which reports this as row C3.

## Examples

``` r
data <- pmx_simulated_fixture(20)
data$ARM <- ifelse(as.integer(data$ID) %% 2L == 0L, "A", "B")
roles <- pmx_roles(
  id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
  cmt = "CMT", dvid = "DVID", covariates = "WT", strata = "ARM"
)
synthetic <- suppressWarnings(synpmx_avatar(data, roles, seed = 1))
#> synpmx_avatar(): dropped 9 undeclared column(s): NTIME, TAD, OCC, RATE, MDV, CENS, LIMIT, AGE, SEX.
#>   Declare a column in `keep` to carry it through verbatim.
compare_pmx_strata_sizes(data, synthetic, roles)
#> Restricted PMX stratum sizes (source against synthetic)
#> 
#> 20 source patients -> 20 synthetic patients.
#> 
#> ARM:
#>   A
#>     10 source -> 10 synthetic (expected 10.0)
#>   B
#>     10 source -> 10 synthetic (expected 10.0)
#> 
#> Cells under the floor of 3 source patients are left to vary by seed on purpose.
#> Source-derived; not releasable.
```
