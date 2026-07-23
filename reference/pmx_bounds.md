# Declare public numeric domains for private PMX fitting

Bounds must be chosen without inspecting confidential patient values.
They define clipping domains and therefore enter the sensitivity
argument.

## Usage

``` r
pmx_bounds(
  time,
  endpoints,
  amt = NULL,
  rate = NULL,
  covariates = NULL,
  limit = NULL
)
```

## Arguments

- time:

  Bounds for actual study time.

- endpoints:

  Named list of DV bounds, one pair per endpoint declaration.

- amt, rate:

  Optional amount and rate bounds.

- covariates:

  Named list of numeric covariate bounds.

- limit:

  Named list of censoring-limit bounds by endpoint.

## Value

A public `pmx_bounds` configuration object.
