# Declare one public baseline covariate

A covariate is either continuous, with a public plausible `range` used
for clipping, or categorical, with public `levels`. The range or level
set must be chosen without inspecting the confidential data.

## Usage

``` r
pmx_covariate(range = NULL, levels = NULL, source)
```

## Arguments

- range:

  Two increasing numbers bracketing a continuous covariate.

- levels:

  Character levels of a categorical covariate.

- source:

  Required provenance string.

## Value

A `pmx_covariate`.
