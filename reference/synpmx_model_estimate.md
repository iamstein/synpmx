# Estimate a population model from a trial

The only stage that reads patient data, and the only one that needs
`nlmixr2`. It works out which endpoint is the drug concentration and
what design produced it, fits the candidate models that design admits,
picks one on AIC, and returns that fit alongside the dosing and visit
models the generated subjects are built from. No patient row survives
it.

## Usage

``` r
synpmx_model_estimate(
  data,
  roles,
  pk = NULL,
  pd = NULL,
  endpoint_roles = NULL,
  covariate_effects = "auto",
  min_subjects = 20L,
  min_arm_patients = 3L,
  min_time_bins = 6L,
  estimation = "focei",
  seed = NULL,
  quiet = FALSE
)
```

## Arguments

- data:

  Source PMX event data.

- roles:

  Explicit column roles from
  [`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md),
  including `nominal_time`.

- pk:

  One of the five built-in structural models, forcing it and skipping
  the search. `NULL` searches the candidates the design admits.

- pd:

  Named character vector of PD shapes per endpoint, skipping that
  search. One of `"constant"`, `"linear"` or `"exponential"` each.

- endpoint_roles:

  Named character vector naming which endpoint is the drug
  concentration, as `c(pk = "cp")`, overriding the inference.

- covariate_effects:

  `"auto"` fits allometric scaling on clearance and volume where a
  weight-like covariate is declared and keeps it where it improves AIC.
  `"none"` fits nothing.

- min_subjects:

  Cohort floor. Below it the covariance matrix describes the subjects it
  was fitted to rather than a population.

- min_arm_patients:

  Minimum patients in every arm, as
  [`synpmx_pca_summarize()`](https://iamstein.github.io/synpmx/reference/synpmx_pca_summarize.md)
  uses.

- min_time_bins:

  Minimum distinct nominal times after a dose across the cohort. Below
  it no linear model is identifiable.

- estimation:

  Passed to `nlmixr2`. `"focei"` by default because the selection
  criterion is AIC and `"saem"` does not reliably produce one at these
  cohort sizes.

- seed:

  Seed for the one random step, which is imputing censored values before
  the fit.

- quiet:

  Suppress the per-candidate progress messages.

## Value

A `pmx_fitted_model`.

## Details

**The fitted parameters are not estimates to report.** They exist to
make simulated profiles look like the source study. The candidate set is
five linear models and the covariate model is allometric scaling or
nothing, which is too little to answer a scientific question, and the
object prints that warning with itself because its contents look exactly
like the output of a real population analysis.

`nominal_time` is required, for two reasons. The dosing and visit models
sit on the nominal grid, and a grid inferred from recorded times is a
statement about the protocol only the caller can make. Estimation,
separately, reads the recorded times and the recorded dosing history,
because a population fit is a statement about the dose that was actually
given.

No formal privacy guarantee is offered. No patient's measured value
reaches the output, which is the claim
[`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md)
makes and is stronger than
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)'s,
but the fixed effects and the covariance matrix are functions of the
individuals in the source and neither is noised. The cohort floor is the
whole defence and it is a threshold rather than an accounting.

## See also

[`synpmx_model_generate()`](https://iamstein.github.io/synpmx/reference/synpmx_model_generate.md),
[`synpmx_model()`](https://iamstein.github.io/synpmx/reference/synpmx_model.md),
[`model_report()`](https://iamstein.github.io/synpmx/reference/model_report.md),
[`model_candidates()`](https://iamstein.github.io/synpmx/reference/model_candidates.md),
[`model_parameters()`](https://iamstein.github.io/synpmx/reference/model_parameters.md).
