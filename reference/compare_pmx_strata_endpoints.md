# Endpoints held by each stratum, source against synthetic

One row per stratum level and endpoint, with how many source patients
and how many synthetic patients contributed an observation of that
endpoint. A zero on both sides is an endpoint the arm never held – a
placebo arm with no pharmacokinetic concentration – and is the normal
case rather than a finding. A nonzero source count against a zero
synthetic count is the row to act on: that arm lost an endpoint on the
way out.

## Usage

``` r
compare_pmx_strata_endpoints(source, synthetic, roles)
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

A `pmx_strata_endpoints` data frame with `level`, `endpoint`,
`source_patients` and `synthetic_patients`.

## Details

Reported for the first `strata` column only, which is the column
[`synpmx_scorecard()`](https://iamstein.github.io/synpmx/reference/synpmx_scorecard.md)
scores row C3 on. Zero rows when the roles declare no strata.

This reads real patient data on both sides and is marked
`"restricted_not_releasable"`.

## See also

[`synpmx_scorecard()`](https://iamstein.github.io/synpmx/reference/synpmx_scorecard.md),
which reports this as row C3, and
[`compare_pmx_strata_sizes()`](https://iamstein.github.io/synpmx/reference/compare_pmx_strata_sizes.md)
for the size half of the same question.

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
compare_pmx_strata_endpoints(data, synthetic, roles)
#> PMX endpoints by stratum (source against synthetic)
#> 
#> ARM, patients contributing each endpoint:
#> 
#>   A / cp
#>     10 source -> 10 synthetic
#>   A / pd
#>     10 source -> 10 synthetic
#>   B / cp
#>     10 source -> 10 synthetic
#>   B / pd
#>     10 source -> 10 synthetic
#> 
#> A zero on both sides is an endpoint the arm never held.
#> Source-derived; not releasable.
```
