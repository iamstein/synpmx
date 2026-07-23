# Public PMX censoring fixture

A tiny, fully simulated event table illustrating uncensored,
left-censored, right-censored, and interval-censored Monolix-style
records. It contains no patient or proprietary information.

## Usage

``` r
pmx_censoring_fixture()
```

## Value

A data frame with `CENS = 0`, `1`, and `-1`; interval censoring uses
`DV` as the reported upper boundary and `LIMIT` as the lower boundary.

## Examples

``` r
fixture <- pmx_censoring_fixture()
roles <- pmx_roles(
  id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
  dvid = "DVID", mdv = "MDV", cens = "CENS", limit = "LIMIT"
)
validate_pmx(fixture, roles)$valid
#> [1] TRUE
```
