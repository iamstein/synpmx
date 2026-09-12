# Generate a dataset from a privately calibrated structural model

Keeps a public structural model's *shape* and spends a small privacy
budget correcting only its *magnitude*. Each subject is reduced to a
bounded multiplicative correction of the model's own prediction, clipped
to a public prior range, and released with calibrated noise. Only a
handful of numbers leave the data, which is why this mode stays usable
at 20 to 60 subjects.

## Usage

``` r
synpmx_calibrated(
  data,
  roles,
  model,
  design,
  priors,
  epsilon,
  n_subjects = NULL,
  seed = 123,
  n_datasets = 1L,
  covariates = NULL,
  backend = "opendp",
  public_source = FALSE
)
```

## Arguments

- data:

  The confidential dataset.

- roles:

  A
  [`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md)
  declaration for `data`. It names the columns of the output too, so the
  synthetic table wears the source study's schema, the way
  [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
  and
  [`synpmx_model()`](https://iamstein.github.io/synpmx/reference/synpmx_model.md)
  output does. A role left undeclared gets no column.

- model:

  A public
  [`pmx_structural_model()`](https://iamstein.github.io/synpmx/reference/pmx_structural_model.md).

- design:

  A public
  [`pmx_trial_design()`](https://iamstein.github.io/synpmx/reference/pmx_trial_design.md).

- priors:

  A
  [`pmx_priors()`](https://iamstein.github.io/synpmx/reference/pmx_priors.md)
  giving a public range per released correction.

- epsilon:

  The privacy budget. A governance decision, not a default;
  [`pmx_preflight()`](https://iamstein.github.io/synpmx/reference/pmx_preflight.md)
  reports the expected fold-error before any is spent.

- n_subjects:

  Number of subjects to generate. Defaults to the released noisy count.

- seed:

  Ordinary generation seed. Unrelated to privacy noise, which is
  controlled by the backend and is never user-seeded.

- n_datasets:

  Number of datasets to draw from the single release. One dataset is
  returned directly; several are returned as a list.

- covariates:

  Optional
  [`pmx_covariates()`](https://iamstein.github.io/synpmx/reference/pmx_covariates.md).

- backend:

  Privacy backend. Defaults to the validated OpenDP adapter and fails
  closed if it is unavailable.

- public_source:

  Assert that `data` is genuinely public. Required by, and only
  meaningful for, `backend = "public"`, which makes no DP claim.

## Value

A data frame under the names in `roles`, carrying its release so that
[`privacy_report()`](https://iamstein.github.io/synpmx/reference/privacy_report.md)
and
[`synpmx_generate()`](https://iamstein.github.io/synpmx/reference/synpmx_generate.md)
can read it. A list of such data frames when `n_datasets > 1`.

## Details

This is the recommended differentially private path for early-phase
cohorts. The tradeoff is that everything not calibrated is *asserted*:
curve shape, variability, residual error, and covariate relationships
come from the public model, so the output is only as realistic as that
model.

**This function spends privacy budget.** Calling it again spends the
budget again. To draw further datasets from the release you have already
paid for, use
[`synpmx_generate()`](https://iamstein.github.io/synpmx/reference/synpmx_generate.md)
or ask for several at once with `n_datasets`.

## Maintenance status

A secondary, provided-as-is path.
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
is the primary, maintained method. The differentially private modes are
complete and tested but not under active development, carry known open
findings, and have not been independently privacy-audited. Use them to
demonstrate the privacy/utility tradeoff, not as a production release
mechanism; a real regulated release needs specialist review and the
external OpenDP backend. Requires a session-level acknowledgment: call
[`synpmx_enable_dp_engines()`](https://iamstein.github.io/synpmx/reference/synpmx_enable_dp_engines.md)
once before this function will run (unless `backend = "public"`, which
makes no DP claim).

## See also

[`synpmx_generate()`](https://iamstein.github.io/synpmx/reference/synpmx_generate.md)
to draw more datasets for free,
[`privacy_report()`](https://iamstein.github.io/synpmx/reference/privacy_report.md)
for the realized accounting.
