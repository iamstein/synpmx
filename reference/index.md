# Package index

## Generating data

The four generation modes.
[`synthesize_pmx()`](https://iamstein.github.io/synpmx/reference/synthesize_pmx.md)
is the default AVATAR path; the others require a public structural
model.

- [`synthesize_pmx()`](https://iamstein.github.io/synpmx/reference/synthesize_pmx.md)
  : Synthesize a structurally faithful PMX dataset (AVATAR-style)
- [`pmx_generate()`](https://iamstein.github.io/synpmx/reference/pmx_generate.md)
  : Generate a synthetic PMX event table from a structural model
- [`fit_calibrated_pmx()`](https://iamstein.github.io/synpmx/reference/fit_calibrated_pmx.md)
  : Calibrate a public structural model to confidential data
- [`fit_private_pmx()`](https://iamstein.github.io/synpmx/reference/fit_private_pmx.md)
  : Fit a subject-level differentially private PMX population generator
- [`generate_pmx()`](https://iamstein.github.io/synpmx/reference/generate_pmx.md)
  : Generate a new PMX event dataset from a fitted private model

## Declaring the data

Roles, endpoints, schema, and bounds that describe an event table.

- [`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md)
  : Declare pharmacometric column roles

- [`pmx_endpoint()`](https://iamstein.github.io/synpmx/reference/pmx_endpoint.md)
  : Declare endpoint scientific-clock behavior

- [`pmx_schema()`](https://iamstein.github.io/synpmx/reference/pmx_schema.md)
  : Capture a schema asserted to be public

- [`pmx_bounds()`](https://iamstein.github.io/synpmx/reference/pmx_bounds.md)
  : Declare public numeric domains for private PMX fitting

- [`pmx_generated_roles()`](https://iamstein.github.io/synpmx/reference/pmx_generated_roles.md)
  :

  Roles for tables produced by
  [`pmx_generate()`](https://iamstein.github.io/synpmx/reference/pmx_generate.md)

## Public model and design inputs

Data-independent inputs for the model-based modes. See the model and
data elicitation articles for how to produce these without reading data.

- [`pmx_structural_model()`](https://iamstein.github.io/synpmx/reference/pmx_structural_model.md)
  : Declare a public structural model
- [`pmx_trial_design()`](https://iamstein.github.io/synpmx/reference/pmx_trial_design.md)
  : Declare a public trial design
- [`pmx_public_design()`](https://iamstein.github.io/synpmx/reference/pmx_public_design.md)
  : Declare public event-design information
- [`pmx_prior()`](https://iamstein.github.io/synpmx/reference/pmx_prior.md)
  : Declare one public prior range
- [`pmx_priors()`](https://iamstein.github.io/synpmx/reference/pmx_priors.md)
  : Collect public priors for the released corrections
- [`pmx_covariate()`](https://iamstein.github.io/synpmx/reference/pmx_covariate.md)
  : Declare one public baseline covariate
- [`pmx_covariates()`](https://iamstein.github.io/synpmx/reference/pmx_covariates.md)
  : Collect public covariate declarations
- [`pmx_covariates_auto()`](https://iamstein.github.io/synpmx/reference/pmx_covariates_auto.md)
  : Declare bootstrap-resampled covariates by column name

## Privacy accounting

Contribution limits, budget, preflight, and the release ledger.

- [`pmx_contribution_limits()`](https://iamstein.github.io/synpmx/reference/pmx_contribution_limits.md)
  : Declare subject contribution limits
- [`pmx_budget_allocation()`](https://iamstein.github.io/synpmx/reference/pmx_budget_allocation.md)
  : Allocate an epsilon budget across private summary groups
- [`pmx_preflight()`](https://iamstein.github.io/synpmx/reference/pmx_preflight.md)
  : Check whether a private release is worth its budget, before spending
  it
- [`privacy_report()`](https://iamstein.github.io/synpmx/reference/privacy_report.md)
  : Summarize a fitted model's privacy contract
- [`validate_private_model()`](https://iamstein.github.io/synpmx/reference/validate_private_model.md)
  : Validate a fitted private PMX population model
- [`dp_backend_status()`](https://iamstein.github.io/synpmx/reference/dp_backend_status.md)
  : Inspect the differential-privacy backend
- [`run_dp_backend_tests()`](https://iamstein.github.io/synpmx/reference/run_dp_backend_tests.md)
  : Run canonical checks against the configured DP backend

## Validation and diagnostics

Structural checks on generated data, and restricted comparisons.

- [`validate_pmx()`](https://iamstein.github.io/synpmx/reference/validate_pmx.md)
  : Validate a pharmacometric event dataset
- [`compare_pmx()`](https://iamstein.github.io/synpmx/reference/compare_pmx.md)
  : Compare source and generated PMX structures inside the restricted
  environment
- [`sampling_summary()`](https://iamstein.github.io/synpmx/reference/sampling_summary.md)
  : Summarize the fitted sampling design
- [`subject_property_summary()`](https://iamstein.github.io/synpmx/reference/subject_property_summary.md)
  : Summarize fitted subject properties and associated regimens

## Fixtures

Fully public example datasets for testing and demonstration.

- [`pmx_censoring_fixture()`](https://iamstein.github.io/synpmx/reference/pmx_censoring_fixture.md)
  : Public PMX censoring fixture
- [`pmx_simulated_fixture()`](https://iamstein.github.io/synpmx/reference/pmx_simulated_fixture.md)
  : Fully simulated public repeated-dose fixture
