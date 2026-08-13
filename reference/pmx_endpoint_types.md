# What kind of values each endpoint takes

Reports, per endpoint, whether its values are continuous or discrete,
and where that answer came from.
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
blends real trajectories, and a weighted mean of several patients' zeros
and ones is a number between them, so a discrete endpoint would come
back continuous unless the generated values are snapped back onto the
levels the source used. This is the function that decides which
endpoints that applies to.

## Usage

``` r
pmx_endpoint_types(data, roles)
```

## Arguments

- data:

  PMX data.

- roles:

  Explicit roles from
  [`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md).

## Value

A `pmx_endpoint_types` data frame with one row per endpoint and columns
`endpoint`, `type`, `levels`, `decided_by`, and `reason`. Marked
`"restricted_not_releasable"`: the level set is read from real data.

## Details

An endpoint is called:

- `"binary"` when every observed value is 0 or 1;

- `"ordinal"` when every observed value is a whole number and there are
  at most 12 distinct ones, which are then the scale;

- `"integer"` when every observed value is a whole number and there are
  more levels than that, so generated values are rounded rather than
  snapped onto a scale. Counts are the usual case, and `"count"` may be
  used to declare one, but the evidence is only that the values are
  whole numbers;

- `"continuous"` otherwise, including when there are fewer than 10
  observed values, which is too few to call.

Declare `endpoint_types` in
[`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md)
to override any of it. The `reason` column says what the data showed, so
a study whose endpoint was called wrongly can be corrected without
guessing at the rule.

Snapping a binary or ordinal endpoint means generated values are source
values: a 0/1 endpoint has no third value to emit. What protects a
patient on a discrete endpoint is the visit-set and dose-schedule
machinery, not the distinctness of any one number.

## See also

[`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md),
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md),
[`synpmx_scorecard()`](https://iamstein.github.io/synpmx/reference/synpmx_scorecard.md).

## Examples

``` r
data <- pmx_simulated_fixture(30)
roles <- pmx_roles(
  id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
  cmt = "CMT", dvid = "DVID", covariates = "WT"
)
pmx_endpoint_types(data, roles)
#> Endpoint value types
#> 
#>  endpoint       type levels decided_by
#>        cp continuous     --   inferred
#>        pd continuous     --   inferred
#>                                      reason
#>  not every observed value is a whole number
#>  not every observed value is a whole number
#> 
#> Generated values on a `binary` or `ordinal` endpoint are snapped to the
#> levels above, and on an `integer` endpoint rounded to whole numbers.
#> Override with `pmx_roles(endpoint_types = )`. Source-derived; not
#> releasable unless separately public or privately budgeted.
```
