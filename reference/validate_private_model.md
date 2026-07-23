# Validate a fitted private PMX population model

Checks structural integrity, accounting, prohibited payload names,
backend status, and the absence of direct patient records. It cannot
independently prove that runtime configuration or the external backend
was honest.

## Usage

``` r
validate_private_model(private_model, strict = FALSE)
```

## Arguments

- private_model:

  An object returned by
  [`fit_private_pmx()`](https://iamstein.github.io/synpmx/reference/fit_private_pmx.md).

- strict:

  Stop on any failed check when `TRUE`.

## Value

A `pmx_private_validation` report.
