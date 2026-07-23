# Calibrate a public structural model to confidential data

The only stage that reads source data. Each subject is reduced to
bounded multiplicative corrections against the structural model's own
prediction, clipped to public prior ranges, and released through a
validated differential-privacy backend.

## Usage

``` r
fit_calibrated_pmx(
  data,
  roles,
  model,
  design,
  priors,
  epsilon,
  covariates = NULL,
  backend = "opendp",
  public_source = FALSE
)
```

## Arguments

- data:

  Confidential PMX event data.

- roles:

  Column roles from
  [`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md).

- model:

  A public
  [`pmx_structural_model()`](https://iamstein.github.io/synpmx/reference/pmx_structural_model.md).

- design:

  A public
  [`pmx_trial_design()`](https://iamstein.github.io/synpmx/reference/pmx_trial_design.md).

- priors:

  Public
  [`pmx_priors()`](https://iamstein.github.io/synpmx/reference/pmx_priors.md)
  for each released correction.

- epsilon:

  Requested subject-level privacy budget.

- covariates:

  Optional public
  [`pmx_covariates()`](https://iamstein.github.io/synpmx/reference/pmx_covariates.md).
  Each declared covariate is released privately and adds one to the
  released dimension.

- backend:

  `"opendp"`, or `"public"` for an explicitly public fixture.

- public_source:

  Logical assertion that the input is already public.

## Value

A `pmx_calibrated_model`, carrying corrected typical parameters,
accounting, provenance, and a release ledger. It contains no raw
records.
