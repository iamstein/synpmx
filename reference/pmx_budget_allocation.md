# Allocate an epsilon budget across private summary groups

Values are fractions of the requested epsilon. Active groups must
receive a positive fraction and the total may not exceed one.

## Usage

``` r
pmx_budget_allocation(
  subject_count,
  event,
  timing,
  covariates,
  endpoints,
  censoring
)
```

## Arguments

- subject_count, event, timing, covariates, endpoints, censoring:

  Nonnegative budget fractions.

## Value

A `pmx_budget_allocation` object.
