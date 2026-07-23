# Compare source and generated PMX structures inside the restricted environment

Any component that uses `source` is marked
`"restricted_not_releasable"`. A fitted private model does not make a
new source-derived comparison private; releasing such a diagnostic
requires a separate public justification or budgeted DP mechanism.

## Usage

``` r
compare_pmx(source, synthetic, roles, endpoints = NULL)
```

## Arguments

- source:

  Source PMX data.

- synthetic:

  Generated synthetic PMX data.

- roles:

  Explicit roles from
  [`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md).

- endpoints:

  Optional endpoint declarations.

## Value

A `pmx_comparison` containing component-level release metadata.
