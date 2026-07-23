# Declare one public prior range

The range is the dominant sensitivity term in the whole design, so its
provenance matters more than any other input. It must be chosen without
inspecting the confidential data.

## Usage

``` r
pmx_prior(range, source)
```

## Arguments

- range:

  Two increasing positive multipliers bracketing the correction factor,
  for example `c(1/4, 4)` for a prediction believed accurate to about
  four-fold.

- source:

  Required provenance string.

## Value

A `pmx_prior`.
