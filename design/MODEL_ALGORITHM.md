# The PMX model generator

Design for a third generator developed in this package, alongside AVATAR
blending and the principal-component (PCA) generator. It infers what each
endpoint is and how it was sampled, estimates a small set of linear
pharmacokinetic (PK) models with `nlmixr2`, picks one, and generates new
subjects by simulating from the fitted population model.

**Out of scope.** The fitted parameters are not estimates to report. They exist
to make simulated profiles look like the source study, and the candidate set is
too small and the covariate model too thin for any of it to answer a scientific
question. This is the same restriction the rest of the package carries, and it
binds harder here because the output of the fit looks exactly like the output of
a real population analysis.

This document is the record of what the algorithm should do and why, written
before the code exists. Nothing here is implemented.

## Where this sits among the three generators

AVATAR blending, the PCA generator and this one are developed in parallel and
none of them is retired. The maintainer's stated expectation on 2026-08-29 is
that this generator is the one he ends up using, which is a direction rather
than a decision: it does not make AVATAR secondary, and nothing in the package
should be written as though it had.

The practical consequence is O17 in `WRITING_FOR_ANDY.md`. Three generators that
get chosen between are read by comparison, so their documents have to hold the
same shape. They do not today: AVATAR has an algorithm document, a demo, a
scorecard document and a public-data survey, PCA has an algorithm document, a
demo and a trial-summary document, and this one would arrive with two. Either
the missing slots get filled or the shape gets stated somewhere as deliberate.

## Status of the decisions

Three interface questions were settled on 2026-08-29 and the rest of the
document assumes the answers.

| Question | Answer |
|---|---|
| What happens to endpoints that are not fittable PK | PK is fitted, the simple pharmacodynamic (PD) time-course shapes are fitted, and everything else — dose reductions, interrupted and missed cycles, discontinuation, visit attendance, covariates, censoring — comes from the models the PCA generator already builds |
| Who picks among candidate models | Automatic on Akaike information criterion (AIC), with a `pk` argument that forces one and skips the search |
| Where the `nlmixr2` dependency sits | Estimation only. Generation runs on base R, from a stored fit |

## The API

Three exported functions, in the shape `synpmx_pca_summarize()`,
`synpmx_pca_generate()` and `synpmx_pca()` already hold.

| Function | Reads patient data | Needs `nlmixr2` | Returns |
|---|---|---|---|
| `synpmx_model_estimate(data, roles, ...)` | yes | yes | `pmx_fitted_model` |
| `synpmx_model_generate(fitted_model, n_subjects, seed)` | no | no | data frame |
| `synpmx_model(data, roles, n_subjects, seed, ...)` | yes | yes | data frame |

Accessors on the fitted model, matching `pca_report()` and its siblings:
`model_report()` inventories everything the object carries,
`model_candidates()` returns the comparison table the selection was made from,
and `model_parameters()` returns fixed effects, the between-subject covariance
matrix and the residual error.

Arguments beyond the roles and the seed:

| Argument | Default | Effect |
|---|---|---|
| `pk` | `NULL` | One of the five built-in models. Skips the search. |
| `pd` | `NULL` | Named vector of PD shapes per endpoint. Skips that search. |
| `endpoint_roles` | `NULL` | Names which endpoint is the drug concentration, overriding inference. |
| `covariate_effects` | `"auto"` | `"auto"`, `"none"`, or a named list of parameter/covariate pairs. |
| `min_subjects` | `20L` | Refuse below this cohort size. |
| `min_arm_patients` | `3L` | Refuse below this many patients in any arm, as `synpmx_pca_summarize()` does. |
| `estimation` | `"saem"` | Passed to `nlmixr2`, with `"focei"` as the documented fallback. |

## Step 1: Classify the Endpoints

`pmx_endpoint_types()` already infers continuous, binary and ordinal from the
values, and that inference is reused unchanged. What this generator adds is a
second classification on top of it: which continuous endpoint is the drug
concentration.

Four signals, computed per endpoint:

1. **Compartment.** The `cmt` role puts the endpoint in the compartment the
   doses go to, or in the compartment one above a dosing compartment that is
   never observed.
2. **Post-dose only.** Observations before the first dose are absent, or sit at
   or below the censoring limit, in most subjects.
3. **Shape.** The median profile within the richest dose interval rises to a
   maximum and then declines without rising again.
4. **Dose proportionality.** Across dose levels, the median subject maximum
   scales with dose. The test is that the ratio of medians between the highest
   and lowest arm falls within a factor of two of the dose ratio.

Signals 2 and 4 are required, 1 and 3 break ties when more than one endpoint
passes. A single-arm study cannot compute signal 4, and there signals 1 to 3
decide. When no endpoint passes, or two pass and the tie does not break,
`synpmx_model_estimate()` errors and names `endpoint_roles` as the way through.

This is the inference-versus-declaration fork, answered the same way
`dose_covariate` answers it: infer where the data settles the question, and
offer the declaration as an override rather than as a requirement.

## Step 2: Detect the Design

Two properties of the design narrow the candidate set. Neither one decides the
model on its own: detection prunes candidates, and the fit picks among what
survives.

**Route.** A `rate` role with nonzero values, or dose records carrying a
duration, means infusion, and `1cmt_infusion` is the only candidate. Otherwise
the time of the subject maximum after the first dose separates the two cases.
Where the median subject's maximum falls on the first post-dose sample and the
profile declines from there, the intravenous (IV) candidates are offered; where
it falls on the second or later, the oral candidates are. A study whose first
sample is drawn late looks oral whichever route it used, so where fewer than
half the subjects agree with the median verdict, both routes are offered and
AIC settles it.

**Sampling richness.** Two counts, one per subject and one across the cohort.

- Per subject: distinct times after dose within the richest dose interval. At
  four or more, with at least one before and two after the median peak, the
  two-compartment candidates are offered. Below that they are not, because a
  distribution phase cannot be identified from troughs, and a two-compartment
  model fitted to troughs will report a `q` and a `v2` that mean nothing.
- Across the cohort: distinct time-after-dose bins holding an observation. Below
  six, no model is offered and the function errors, pointing at
  `synpmx_avatar()` and `synpmx_pca()`, which need no identifiable structure.

`sampling_summary()` computes something adjacent for the differentially private
path, but it reads a fitted private model rather than the source, so the counts
here are new code.

## Step 3: Estimate the Candidates

The candidate set is at most the five models in `.pk_models`, and never more
than that, for the reason in the section on the `nlmixr2` boundary below.

Each candidate is fitted with log-normal between-subject variability on `cl` and
`v`, plus `ka` for the oral models and `q` and `v2` for the two-compartment
ones, and a proportional residual error. An endpoint carrying values at or below
zero after censoring is handled gets an additive error instead, and the
substitution is recorded in the fit. Where `cens` is declared, censored
observations go through the M3 likelihood; where the estimation method in use
does not offer it, they are imputed at half the limit and that is recorded too.

Selection is on AIC over the candidates that converged. A candidate that fails
to converge drops out carrying its reason, and `model_candidates()` shows it
alongside the ones that fitted, so a search that came down to one survivor does
not look like a search that had one candidate. Where nothing converges, the
function errors rather than returning the least bad fit.

**Covariates.** The default `"auto"` fits allometric scaling on clearance and
volume when a weight-like covariate is declared, keeps it when it improves AIC,
and fits nothing else. Everything beyond that is opt-in through
`covariate_effects`. The cost of this default is explicit: covariates that
influence the real profiles and are not in the model are generated independently
of the profiles, so the synthetic data carries no relationship between them.
AVATAR preserves those relationships without modelling them, because a blended
subject's covariates and profile come from the same donors. `model_report()`
reports the correlation between each declared covariate and the individual
random effects, which is where an unmodelled relationship shows up.

**PD.** Every remaining continuous endpoint is fitted against the `constant`,
`linear` and `exponential` shapes in `.pd_models`, with between-subject
variability on the baseline, and selected on AIC. These shapes have no exposure
dependence, so a PD endpoint driven by concentration is reproduced as a time
course that happens to resemble the average subject's response. A dataset whose
point is the exposure-response relationship is not served by this generator.

Binary and ordinal endpoints are not fitted. They are generated from per-visit
marginals by arm, as the schedule machinery already does.

## Step 4: The Dosing Model and the Visit Model

Everything that is not the concentration-time curve comes from the models
`synpmx_pca_summarize()` already builds, unchanged in what they represent. Dose
reductions, interrupted cycles, discontinuation and missed visits are modelled
here exactly as they are for PCA, and the same is true of arm sizes, the
covariate distributions, the censoring boundary and the schema.

`.pca_dose_model()` and `.pca_arm_models()` build all of it. The work is to lift
them out of `pca.R` into a shared file rather than to write them a second time,
and that lift is the first commit, on its own, with the PCA tests unchanged as
its check. What the lifted code produces, per arm:

| Model | What it holds | Drawn at generation as |
|---|---|---|
| Planned schedule | The nominal dose times enough of the arm reached, and the modal amount at each cycle among patients still on their starting dose | The cycle grid every subject starts from |
| Dose ladder | The levels patients dropped to, as ratios to their own starting dose, built from within-patient decreases | The amount multiplier in force at a cycle |
| Reduction rate | Discrete-time hazard of stepping down a level | Decided before the cycle is dosed |
| Interruption rate | Discrete-time hazard of skipping a cycle without ending treatment | Decided at the cycle |
| Discontinuation rate | Discrete-time hazard of stopping treatment | Decided after the cycle is dosed |
| Visit model | Per endpoint and per retained nominal time, the fraction of the arm holding an observation there | Attendance, drawn per visit |

The reduction ladder is read from within-patient decreases rather than from the
spread of amounts across patients, which is what lets it work on a study dosed
by body weight, where no two patients share an amount. That reasoning is already
written down in `.pca_dose_model()` and does not change here.

**Two time axes.** PCA replaces `roles$time` with the nominal grid for its whole
pipeline, because every feature it fits is a cell on that grid. This generator
cannot do that. A population PK fit is a statement about time after the dose
that was actually given, so estimation reads recorded times and recorded dosing
histories, including the reductions and the skipped cycles as they happened. The
dosing and visit models read the nominal grid, as they do for PCA. `nominal_time`
is therefore required here too, and for the same reason: the grid is a statement
about the protocol that only the caller can make.

**Estimation reads actuals, not the plan.** Fitting against the planned schedule
where a patient's dose was reduced would push the resulting drop in
concentration into clearance, and the model would report a population that
eliminates the drug faster than the real one. The subject-level dosing records
go to `nlmixr2` as recorded.

## Step 5: What the Fitted Model Carries

`pmx_fitted_model` holds the estimated part and the models from Step 4.

Estimated: the selected structural model, the candidate comparison table, fixed
effects, the between-subject covariance matrix, residual error, the covariate
effects that survived, and the PD shape and parameters per endpoint.

Not estimated by `nlmixr2`: the per-arm dosing model and its three rates, the
visit model, arm sizes, the covariate distributions, the censoring boundary, the
schema and the roles.

**No individual estimates.** Empirical Bayes estimates are per-subject
quantities, and a fitted model that carried them would be writing out a
description of each real patient. They are not stored on the object, and
generation draws random effects from the covariance matrix instead. The
correlation report named above computes from them and keeps only the
correlation.

## Step 6: Generate

Per synthetic subject: assign an arm keeping the source arm shares, draw
covariates from the arm's covariate model, draw random effects from the
covariance matrix, apply the covariate effects to the typical parameters, draw
the dose schedule from the arm's dosing model with `.pca_draw_schedule()`, draw
the visits attended from the visit model, evaluate `.pk_profile()` at the
attended times against the drawn schedule, add residual error, apply the
censoring boundary and emit.

`.pk_profile()` already superposes linear doses, accepts a per-subject parameter
vector, and handles infusion duration, so generation needs no solver and no
compiler.

**Dose changes reach the concentrations.** The drawn schedule is the input to
`.pk_profile()`, so a synthetic subject who steps down a dose level has a lower
exposure from that cycle on, and one who skips a cycle has the trough that
implies. In PCA the schedule and the profile are drawn from separate models and
a reduction appears in the dosing records without appearing in the values. This
is the largest fidelity gain the generator offers over PCA. The regression
check for it compares exposure before and after a drawn reduction, rather than
only checking that the dosing records show one.

The order of the two draws is fixed: the schedule is drawn first and the profile
is computed from it. Drawing a profile and then a schedule would
reproduce PCA's disconnect with extra steps.

## The `nlmixr2` boundary

`nlmixr2` appears in `Suggests`, is loaded only by `synpmx_model_estimate()`,
and the candidate set is exactly the models `.pk_single_dose()` can evaluate in
closed form. A candidate the fitter could estimate and the generator could not
simulate would be a model that fits and then generates nothing, so the two lists
are one list.

Vignettes ship a stored `pmx_fitted_model` built by a script under `scripts/`,
so `R CMD check` and the pkgdown site never compile a model. A chunk that calls
`synpmx_model_estimate()` is guarded on `requireNamespace("nlmixr2")` and is
there to show the call, not to produce the object the rest of the document uses.

## Gates

Conditions under which the function refuses, each of which needs a registry row
and a test before the code exists.

| Gate | Threshold | Reason |
|---|---|---|
| Cohort size | 20 subjects | A covariance matrix fitted to a handful of subjects describes those subjects. PCA's floor is 10 and this one is higher because a parameter estimate concentrates faster than a score. |
| Cohort time coverage | 6 distinct time-after-dose bins | Below it no linear model is identifiable and the fit reports parameters that came from the starting values. |
| Arm size | 3 patients | Inherited from the dosing and visit models, which are summaries of an arm. An arm of one or two has no rates to pool. |
| `nominal_time` undeclared | — | The dosing and visit models sit on the nominal grid, and inferring that grid is a statement about the protocol only the caller can make. |
| No PK endpoint identified | — | Errors and names `endpoint_roles`. |
| No candidate converged | — | Errors rather than returning the least bad fit. |
| Fitted model used as a public input | — | `pmx_fitted_model` is a distinct class, and `pmx_prior()` and `synpmx_calibrated()` must refuse it. A model estimated from the confidential data that entered the differentially private path as a public structural input would consume no budget for information it took from the data. This is the one gate whose absence is a privacy defect rather than a quality defect. |

## Privacy posture

No patient's measured value reaches the output, which is the same claim
`synpmx_pca()` makes and stronger than AVATAR's. No formal guarantee is offered,
which is also the same. The fixed effects and the covariance matrix are
functions of the individuals in the source, and neither is noised.

As the cohort shrinks, a population parameter approaches the one subject it was
fitted to, and a one-subject cohort would produce a model that reproduces that
subject. The
cohort gate above is the whole defence, and it is a threshold rather than an
accounting.

## Documentation obligations

Adding this generator makes three existing statements false, and each one is a
defect in the commit that adds the code rather than a follow-up.

- `README.md` describes three secondary modes beside AVATAR, and calls AVATAR
  the primary maintained code. The count and the maintenance-status section
  change. Its voice belongs to the maintainer, so the restructuring is his to
  confirm rather than a defect fix taken without asking.
- `DESCRIPTION` was corrected on 2026-08-29, before this generator exists: the
  clause claiming the package "does not fit a pharmacokinetic,
  pharmacodynamic, or nonlinear mixed-effects model" is gone, PCA is named, and
  the out-of-scope statement is now about what the data is for. It still
  describes two generators, and gains a clause naming this one when the code
  lands.
- `vignettes/synpmx-methods.Rmd` introduces four modes and adds `synpmx_pca()`
  as a fifth in a paragraph. A sixth arriving the same way makes the document a
  list of exceptions, so the mode inventory needs restructuring rather than
  another paragraph.
- `_pkgdown.yml` needs the new functions and documents.

New documents, following the shape the AVATAR and PCA pairs hold:
`vignettes/model-algorithm.Rmd` for the mechanism and `vignettes/model-demo.Rmd`
for the first run. The cross-document contract in `AGENTS.md` grows a fourth
entry: each step in `model-algorithm.Rmd` describes an operation the function
performs, and the functions carry the pointer comment naming their section.

`synpmx_scorecard()` scores a synthetic dataset against its source and does not
know which generator produced it, so it needs nothing new. Whether it should
gain a row for goodness of fit is the open question below.

## Open questions

1. **Goodness of fit as a scorecard row.** The generator is usable only where
   the selected model describes the source, and nothing currently scores
   that. A visual predictive check belongs in `model_report()`, but a pass or
   fail row on the scorecard would put it where a user looks before deciding to
   use a dataset. The scorecard's row set is a cross-document contract, so this
   is a change to `avatar-scorecard.Rmd` as well.
2. **Whether the arm structure should come from the fit or the apparatus.**
   Dose and arm are currently apparatus, which means a dose-dependent PK
   nonlinearity in the source is generated away rather than reported. Fitting a
   dose effect on clearance would catch it, at the cost of a candidate set that
   is no longer just the five linear models.
3. **Occasion-varying dosing.** `mavoglurant` resets time within occasion and
   carries an occasion-varying assigned dose. The apparatus handles it; whether
   the estimation step should pool occasions or fit them separately is
   undecided.
