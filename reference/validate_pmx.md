# Validate a pharmacometric event dataset

Checks schema usability, chronological event logic, explicit endpoint
semantics, derived timing fields, censoring conventions, baseline
constancy, subject properties, and occasion-assigned dose coherence. It
does not assess scientific or inferential validity.

## Usage

``` r
validate_pmx(data, roles, endpoints = NULL, strict = FALSE)
```

## Arguments

- data:

  A PMX event data frame or tibble.

- roles:

  Explicit roles from
  [`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md).

- endpoints:

  Optional named endpoint declarations from
  [`pmx_endpoint()`](https://iamstein.github.io/synpmx/reference/pmx_endpoint.md).

- strict:

  Stop when any error-level check fails.

## Value

A `pmx_validation` report with `valid`, `checks`, and `summary`.
