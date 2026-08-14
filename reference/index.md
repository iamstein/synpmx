# Package index

## Generating data

The four generation modes.
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
is the default AVATAR path; the others require a public structural
model.

- [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
  : Synthesize a structurally faithful PMX dataset (AVATAR-style)
- [`synpmx_prior()`](https://iamstein.github.io/synpmx/reference/synpmx_prior.md)
  : Generate a dataset from public inputs only
- [`synpmx_calibrated()`](https://iamstein.github.io/synpmx/reference/synpmx_calibrated.md)
  : Generate a dataset from a privately calibrated structural model
- [`synpmx_empirical()`](https://iamstein.github.io/synpmx/reference/synpmx_empirical.md)
  : Generate a dataset from a dense differentially private release
- [`synpmx_generate()`](https://iamstein.github.io/synpmx/reference/synpmx_generate.md)
  : Draw another dataset from a release already paid for

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
  [`.generate_structural()`](https://iamstein.github.io/synpmx/reference/dot-generate_structural.md)

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

- [`synpmx_enable_dp_engines()`](https://iamstein.github.io/synpmx/reference/synpmx_enable_dp_engines.md)
  : Acknowledge the DP engines' unaudited status for this session

- [`synpmx_disable_dp_engines()`](https://iamstein.github.io/synpmx/reference/synpmx_disable_dp_engines.md)
  :

  Withdraw the acknowledgment from
  [`synpmx_enable_dp_engines()`](https://iamstein.github.io/synpmx/reference/synpmx_enable_dp_engines.md)

## Validation and diagnostics

Structural checks on generated data, and restricted comparisons.

- [`synpmx_scorecard()`](https://iamstein.github.io/synpmx/reference/synpmx_scorecard.md)
  : The scorecard for one synthetic dataset

- [`validate_pmx()`](https://iamstein.github.io/synpmx/reference/validate_pmx.md)
  : Validate a pharmacometric event dataset

- [`pmx_endpoint_types()`](https://iamstein.github.io/synpmx/reference/pmx_endpoint_types.md)
  : What kind of values each endpoint takes

- [`compare_pmx()`](https://iamstein.github.io/synpmx/reference/compare_pmx.md)
  : Compare source and generated PMX structures inside the restricted
  environment

- [`compare_pmx_distributions()`](https://iamstein.github.io/synpmx/reference/compare_pmx_distributions.md)
  : Compare per-covariate and per-endpoint distributions of source and
  synthetic

- [`compare_pmx_rare_levels()`](https://iamstein.github.io/synpmx/reference/compare_pmx_rare_levels.md)
  : Which rare source levels reached the synthetic output

- [`compare_pmx_strata_sizes()`](https://iamstein.github.io/synpmx/reference/compare_pmx_strata_sizes.md)
  : Stratum sizes, source against synthetic

- [`skeleton_uniqueness()`](https://iamstein.github.io/synpmx/reference/skeleton_uniqueness.md)
  : Score how many patients share each patient's event skeleton

- [`plot_pmx_schedule()`](https://iamstein.github.io/synpmx/reference/plot_pmx_schedule.md)
  : Draw a cohort's dosing and observation schedule

- [`unmaskable_strata()`](https://iamstein.github.io/synpmx/reference/unmaskable_strata.md)
  : Which strata can mask their own avatars

- [`pmx_masking_report()`](https://iamstein.github.io/synpmx/reference/pmx_masking_report.md)
  : Report what each masking mechanism did, and what it cost

- [`compare_pmx_proximity()`](https://iamstein.github.io/synpmx/reference/compare_pmx_proximity.md)
  : Are synthetic subjects sitting too close to real ones?

- [`flag_identifiable_subjects()`](https://iamstein.github.io/synpmx/reference/flag_identifiable_subjects.md)
  : Flag structurally unusual – and so easily identifiable – subjects

- [`remediate_identifiable_subjects()`](https://iamstein.github.io/synpmx/reference/remediate_identifiable_subjects.md)
  :

  Remove or shorten the subjects
  [`flag_identifiable_subjects()`](https://iamstein.github.io/synpmx/reference/flag_identifiable_subjects.md)
  flags

- [`sampling_summary()`](https://iamstein.github.io/synpmx/reference/sampling_summary.md)
  : Summarize the fitted sampling design

- [`strata_summary()`](https://iamstein.github.io/synpmx/reference/strata_summary.md)
  : Summarize fitted strata and associated regimens

## Fixtures

Fully public example datasets for testing and demonstration.

- [`pmx_censoring_fixture()`](https://iamstein.github.io/synpmx/reference/pmx_censoring_fixture.md)
  : Public PMX censoring fixture
- [`pmx_simulated_fixture()`](https://iamstein.github.io/synpmx/reference/pmx_simulated_fixture.md)
  : Fully simulated public repeated-dose fixture
