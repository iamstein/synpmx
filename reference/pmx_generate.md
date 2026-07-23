# Generate a synthetic PMX event table from a structural model

Works in two modes. Supplied a
[`pmx_structural_model()`](https://iamstein.github.io/synpmx/reference/pmx_structural_model.md)
it generates purely from public inputs, reads no confidential data, and
makes no privacy claim. Supplied a
[`fit_calibrated_pmx()`](https://iamstein.github.io/synpmx/reference/fit_calibrated_pmx.md)
result it uses the privately corrected parameters; that is
post-processing and consumes no further budget.

## Usage

``` r
pmx_generate(
  x,
  design = NULL,
  n_subjects = NULL,
  seed = NULL,
  dropout = 0,
  lloq = NULL,
  covariates = NULL
)
```

## Arguments

- x:

  A `pmx_structural_model` or a `pmx_calibrated_model`.

- design:

  A
  [`pmx_trial_design()`](https://iamstein.github.io/synpmx/reference/pmx_trial_design.md).
  Taken from `x` when it is a calibrated model.

- n_subjects:

  Number of subjects. Defaults to the planned cohort total, or to the
  released private count for a calibrated model.

- seed:

  Ordinary generation seed. Unrelated to privacy noise.

- dropout:

  Fraction of subjects who discontinue early. A public assumption from
  the protocol.

- lloq:

  Lower limit of quantification. Observations below it are flagged
  `CENS = 1` with `DV` at the limit, following the Monolix convention.

- covariates:

  Optional
  [`pmx_covariates()`](https://iamstein.github.io/synpmx/reference/pmx_covariates.md)
  for prior-mode generation. Ignored for a calibrated model, which
  carries its own released covariate summaries.

## Value

A data frame in PMX event-table form.
