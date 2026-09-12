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
  pd_by_arm = FALSE,
  endpoint_roles = NULL,
  start_param = NULL,
  covariate_effects = "none",
  min_subjects = 20L,
  min_arm_patients = 3L,
  min_time_bins = 6L,
  max_fit_subjects = 60L,
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

- pd_by_arm:

  Fit each PD endpoint's shape per arm rather than once over the pooled
  cohort. `FALSE`, the default, gives every arm the same time course, so
  a synthetic patient's dose does not reach their response. `TRUE` gives
  each arm its own shape, selected on AIC within that arm and asserting
  no dose-response form; an arm holding too little of an endpoint to fit
  keeps the pooled shape. Costs nothing measurable, since these are
  least-squares fits.

- endpoint_roles:

  Which endpoint is the drug concentration, as `c(pk = "cp")`,
  overriding the inference.

  **More than one may be named**, for a study that measures two
  concentrations — two drugs, or a parent and its metabolite. Each gets
  its own structural model, its own parameters and its own residual
  error, evaluated against the one dose schedule they share; no
  correlation between their random effects is estimated. Write it as
  `list(pk = c("parent", "metabolite"))`, or
  `c(pk = c("parent", "metabolite"))`, which R renames to `pk1`/`pk2`
  and which is read the same way.

  Naming several is a declaration, never an inference. Left to itself
  the classification picks one concentration and treats every other
  continuous endpoint as a pharmacodynamic time course, because a second
  endpoint that passes the concentration signals is at least as often a
  biomarker as a metabolite — `onc_sim`'s tumour size passes them. Where
  a demoted endpoint really is a concentration, that time course has no
  dose term in it and the generated values lose their dose ordering,
  which is why naming both is worth doing.

- start_param:

  Starting values for the population PK parameters, as
  `start_param = c(cl = 4, v = 40, ka = 0.5)`. Anything not named is
  read off the curve as usual, so a caller who knows the clearance and
  not the absorption says only the clearance.

  Where the study fits **more than one** concentration endpoint, key it
  by endpoint —
  `start_param = list(parent = c(cl = 4), metabolite = c(cl = 9))` —
  because a flat vector would not say which one it describes, and is
  refused with that message. With one concentration the flat form is
  what to write. The PD endpoints are least-squares time courses that
  this does not reach.

  The automatic starting values are a non-compartmental read of the
  cohort's median profile. On a study that read cannot describe —
  sampled only at troughs, dosed by two routes whose reads disagree, or
  recorded in units that are not what they appear — it can start the
  search somewhere the optimizer cannot leave, which
  [`model_report()`](https://iamstein.github.io/synpmx/reference/model_report.md)
  then reports as a fit that did not move. This is how a caller who
  knows the compound fixes that in advance, and
  [`model_report()`](https://iamstein.github.io/synpmx/reference/model_report.md)
  says which values were declared.

- covariate_effects:

  `"none"`, the default, puts no covariate in the structural model.
  `"auto"` fits allometric scaling on clearance and volume where a
  weight-like covariate is declared and keeps it where it improves AIC.

  The default is `"none"` because a synthetic study does not need the
  relationship: covariates are generated from the source's own
  study-wide distributions either way, and a clearance that moves with
  weight buys nothing the generator spends. `fit$correlations` still
  reports where a covariate moves with a random effect, so the
  relationship the model is not carrying is still visible.

- min_subjects:

  Cohort size the fit should have. Below it the covariance matrix
  describes the subjects it was fitted to rather than a population,
  which warns rather than refuses: the fit runs on whatever the study
  has.

- min_arm_patients:

  Minimum patients in every arm, as
  [`synpmx_pca_summarize()`](https://iamstein.github.io/synpmx/reference/synpmx_pca_summarize.md)
  uses. Patients in a shorter arm are dropped with a warning before
  anything is fitted, so that arm is absent from the fitted model and
  from the data generated from it.

- min_time_bins:

  Distinct nominal times after a dose the cohort should hold. Below it a
  one-compartment model is not identifiable, which warns rather than
  refuses: the fit runs and its parameters sit close to their starting
  values. No post-dose observation at all is an error.

- max_fit_subjects:

  Most subjects the population model is fitted to, 60 by default. A
  study above it has its PK parameters estimated on a subset drawn in
  proportion to the arms, using `seed`; the dosing, visit and covariate
  models still read every subject. Fit time is linear in subjects and
  the parameters a synthetic study needs are settled well before the
  sixtieth, so this is where a twenty-minute fit becomes a five-minute
  one.
  [`model_report()`](https://iamstein.github.io/synpmx/reference/model_report.md)
  states the count. Cannot be below `min_subjects`.

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
