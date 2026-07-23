# Declare subject contribution limits

Declare subject contribution limits

## Usage

``` r
pmx_contribution_limits(
  max_rows,
  max_doses,
  max_occasions,
  max_observations_per_endpoint,
  max_timing_cells = 12L
)
```

## Arguments

- max_rows, max_doses, max_occasions, max_timing_cells:

  Positive public limits applied independently to every source subject.

- max_observations_per_endpoint:

  A positive scalar or a named vector giving per-endpoint limits.

## Value

A `pmx_contribution_limits` object.
