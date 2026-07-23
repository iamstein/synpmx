# Declare a public structural model

The structural model and its typical parameter values are public inputs:
they must be established without inspecting the confidential dataset.
For a first-in-human compound they normally come from preclinical
allometric scaling, which is also what selected the starting dose.

## Usage

``` r
pmx_structural_model(
  pk,
  typical,
  pd = "none",
  source,
  rx = NULL,
  iiv = c(cl = 0.3, v = 0.2),
  residual_cv = 0.15
)
```

## Arguments

- pk:

  One of `"1cmt_iv"`, `"1cmt_oral"`, `"1cmt_infusion"`, `"2cmt_iv"`,
  `"2cmt_oral"`.

- typical:

  Named numeric vector of typical parameter values, interpreted as the
  median of a lognormal population. Requires `cl` and `v` (the central
  volume), plus `ka` for oral models and `q` and `v2` for
  two-compartment models. Optional `f` defaults to 1. PD shapes
  additionally require `baseline`, plus `slope` for `"linear"` and
  `plateau` and `rate` for `"exponential"`.

- pd:

  The PD time course, with no exposure dependence. One of `"none"` (PK
  only), `"constant"`, `"linear"` (needs `slope`), or `"exponential"`
  (needs `plateau` and `rate`, covering both decay and growth). A simple
  shape is adequate for exercising longitudinal code and calibrates
  through a well-conditioned level correction.

- source:

  Required provenance string recording where the model and its typical
  values came from. Recorded in the release ledger.

- rx:

  Reserved for an `rxode2` model. **Not yet implemented**: the value is
  stored on the returned object but the generator always evaluates the
  built-in analytic solution, so supplying it warns.

- iiv:

  Named vector of between-subject variability, as CV on the log scale. A
  public assumption; it consumes no privacy budget.

- residual_cv:

  Proportional residual error, as a CV.

## Value

A `pmx_structural_model`.
