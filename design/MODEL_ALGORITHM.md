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

This document was the record of what the algorithm should do and why, written
before the code existed. **It is now implemented**, in `R/model.R`,
`R/model-design.R`, `R/model-estimate.R` and `R/model-generate.R`, and
documented for readers in `vignettes/pmxmodel-algorithm.Rmd` and
`vignettes/pmxmodel-demo.Rmd`. Where the two disagree the vignette is right:
it is generated from the code and the code is under test. What this document
keeps is the reasoning, the rejected alternatives, and the record of what moved
between the design and the implementation, which is the section below.

## Where this sits among the three generators

AVATAR blending, the PCA generator and this one are developed in parallel and
none of them is retired. The maintainer's stated expectation on 2026-08-29 is
that this generator is the one he ends up using, which is a direction rather
than a decision: it does not make AVATAR secondary, and nothing in the package
should be written as though it had.

The practical consequence is O17 in `WRITING_FOR_ANDY.md`. Three generators that
get chosen between are read by comparison, so their documents have to hold the
same shape. PCA's survey lives at
`vignettes/articles/pca-public-data-examples.Rmd`, so the surveys now cover two
of the three generators.

A second survey for PCA was written on this branch on 2026-08-29 and discarded
on 2026-09-01: one already existed on `claude/synpmx-package-work-cgcrue`, and
it was the better document. The difference is the whole lesson. The discarded
one copied the AVATAR survey's role declarations verbatim and let the generator
refuse whatever it refused, which ran two of eight datasets. The one that
survives declares or constructs a nominal grid per dataset, saying each time
what that declaration asserts about the protocol, and runs seven of eight. The
gate refuses to let the *package* infer a grid. It does not stop a *document*
from declaring one, and a survey that reads it as though it did is a survey
about the gate rather than about the generator.

This generator arrives with two documents, `pmxmodel-algorithm.Rmd` and
`pmxmodel-demo.Rmd`. Its public-data survey is developed after the maintainer is
satisfied with the demo, not alongside it: a survey written against a demo that
is still moving would be rewritten once for every change to it. That is a
sequencing decision rather than an exemption, and the contract in `AGENTS.md`
already says the third survey holds the same shape as the first two.

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

The lift happened on 2026-08-29: `.pca_dose_model()`, `.pca_draw_schedule()`
and `.pca_arm_models()` moved out of `R/pca.R` into `R/dose-visit-models.R` as
`.dose_model()`, `.draw_schedule()` and `.arm_models()`, so this generator calls
the same code rather than a copy of it. Nothing there reads a fitted basis:
`.arm_models()` takes the grid attendance is measured on, as a data frame of
`index`, `name`, `endpoint` and `time`, and `.pca_cells()` is the adapter that
reads that grid off a PCA fit. This generator writes its own adapter over the
nominal grid and calls the same function. What the shared code produces, per
arm:

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
written down in `.dose_model()` and does not change here.

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
the dose schedule from the arm's dosing model with `.draw_schedule()`, draw
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

Conditions under which the function refuses. Each needs a deterministic test
before the code exists. Only one needs a registry row: `design/ISSUES.md`
records defects and findings, and a gate built to this specification is
neither, so seven near-identical rows would be seven rows nobody reads. The
exception is the last one, which is a privacy defect if it is ever absent and
is `REV-042`. The cohort floor's own limitation is `SIM-056`.

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

New documents, following the shape the AVATAR and PCA sets hold:
`vignettes/pmxmodel-algorithm.Rmd` for the mechanism and
`vignettes/pmxmodel-demo.Rmd` for the first run — the maintainer's names, chosen
on 2026-09-01 — with `vignettes/articles/pmxmodel-public-data-examples.Rmd`
following once the demo has settled. The survey is an article, as both existing
surveys became on 2026-09-01, so `R CMD check` never runs it and
`./build.sh articles` is what proves it. The cross-document contract in
`AGENTS.md` grows a fourth entry: each step in `model-algorithm.Rmd` describes
an operation the function performs, and the functions carry the pointer comment
naming their section.

`synpmx_scorecard()` scores a synthetic dataset against its source and does not
know which generator produced it, so it needs nothing new. Whether it should
gain a row for goodness of fit is the open question below.

## What moved between the design and the code

Five things, each found by running the algorithm rather than by reading it.

1. **The route detector was rewritten.** The design put it on a per-subject vote
   with a "fewer than half agree" guard. That counts the sampling schedule
   rather than the drug: on `warfarin`, 22 of 32 patients are first sampled at
   24 hours, every one of their profiles declines from its first point, and the
   vote reads 69% intravenous on a study that is oral — above the guard, so it
   commits to the wrong answer. What actually separates the routes is that an
   oral dose cannot be in the blood at the moment it is swallowed. So: a peak at
   time zero after a dose is intravenous, a median profile that rises before it
   falls is oral, and a profile declining from a first sample drawn later is
   both at once and both candidates are offered. The rising share is still
   reported and decides nothing.

2. **`estimation` defaults to `"focei"`, not `"saem"`.** Selection is on AIC and
   SAEM does not reliably produce one: its log-likelihood is a Gaussian
   quadrature run after the fit, and at phase 1 cohort sizes it returns a
   non-finite value. `theo_sd` fits perfectly well under SAEM — clearance 2.75,
   volume 32.3, absorption 1.51 — and reports `AIC = Inf`.

3. **Signal 4 needs a definition of a dose level, and the design had none.**
   Reading each distinct amount as its own level makes a study dosed by body
   weight compare a median of one patient against a median of one patient over a
   dose ratio near 1, which almost any endpoint passes. Where the distinct
   amounts outnumber half the cohort the signal is not computable, as on a
   single-arm study.

4. **The richness rule refused every intravenous study.** "At least one before
   and two after the median peak" cannot be met by a profile that peaks at the
   dose. The before-peak requirement now applies only where there is an
   ascending limb to sample.

5. **`integer` endpoints are fittable.** The design said continuous. `warfarin`'s
   prothrombin activity is typed `integer` and is a continuous quantity that was
   rounded, which is what `.snap_endpoint_values()` puts back at emit.

Two more things the design did not anticipate. The PD shapes are fitted by
least squares rather than through `nlmixr2`, because they are three-parameter
curves on one endpoint and routing them through a population fitter would put a
compiler in the path of every PD endpoint. And a concentration is floored at
zero, because a proportional error near the assay limit produces negatives often
enough to matter.

`SIM-055` moved from a PCA row to a shared one: `synpmx_model()` uses the same
visit model and inherits the same missing attendance floor. On `warfarin` it
reads B4a = 2 FAIL, and the demo shows the failing row rather than hiding it.

## Build Order

Eight commits, ordered so that everything not needing `nlmixr2` is written and
under test before any estimation exists. Commit 4 produces a working generator
with no fitter in it at all, which is what makes the base-R half testable on its
own and keeps a slow, compiling dependency out of the loop until it has to be
there.

| # | Commit | Delivers | Checked by |
|---|---|---|---|
| 1 | The object and the gates | `pmx_fitted_model`, and a `synpmx_model_estimate()` that validates its inputs and refuses everything the Gates table says it should. No fitting. | Deterministic tests, one per gate, on fixtures. Needs no suggested package. |
| 2 | Endpoint classification and design detection | Steps 1 and 2: which endpoint is the drug concentration, the route, the richness counts, and the candidate set each implies. Returns the candidate set rather than fitting it. | Fixtures built to land on each branch, including the weak-route case that offers both. |
| 3 | The apparatus adapter | The nominal-grid `cells` adapter over `.arm_models()`, so a fitted model carries the dosing and visit models. | The PCA survey's own outputs are the reference: the same source through both adapters gives the same dosing model. |
| 4 | Generation | `synpmx_model_generate()` against a hand-constructed `pmx_fitted_model`, through `.pk_profile()`. A complete generator with no fitter in it. | End to end on fixtures, base R only. Includes the exposure-before-and-after-a-reduction check from Step 6. |
| 5 | Estimation | `nlmixr2` behind `requireNamespace()`: the candidate fits, the AIC table, the `pk` override, the failure-to-converge paths. | Guarded tests, skipped where `nlmixr2` is absent. |
| 6 | Covariates and PD | Allometric scaling under `covariate_effects = "auto"`, the three PD shapes, and the random-effect correlation report. | Guarded, as commit 5. |
| 7 | The documents | `vignettes/pmxmodel-algorithm.Rmd` and `vignettes/pmxmodel-demo.Rmd`, plus the stored fit `scripts/build-model-fits.R` writes and both knit against. | `./build.sh`. |
| 8 | The documentation sweep | The three false statements in the Documentation obligations section, and `_pkgdown.yml`. | Search for the affected names, as `AGENTS.md` requires. |

Commits 1 to 4 need nothing beyond base R, so they run in this repository's
ordinary loop. Commits 5 and 6 need a machine with `nlmixr2` installed and are
the first point at which a compiler is in the path.

**The order was changed on 2026-09-01, by the maintainer:** the two documents
were pulled forward and everything they need was built with them, so commits 2
to 7 landed together. The base-R half is still tested on its own — the
generation tests construct a `pmx_fitted_model` by hand and never load
`nlmixr2` — which is the property the order existed to protect. The estimation
tests skip where the fitter cannot build a model, which is a stronger guard than
`requireNamespace()`: on macOS, R's `FLIBS` points at `/opt/gfortran`, and where
that directory is absent every model compiles, none of them links, and the
failure reports itself as a missing C compiler.

The public-data survey is not in this list. It follows once the demo from
commit 7 has settled, for the reason in "Where this sits among the three
generators".

## Questions settled on 2026-09-01

The three questions this document left open were put to the maintainer and
answered. Two of them are answered "no", and the reason in both cases is the
same one: this generator is meant to be simple, and a candidate set that grows
to catch a case is how it stops being.

1. **Goodness of fit as a scorecard row.** Still open as a row, but the shape it
   would take is now settled, and it is not a shape specific to this generator.
   There is one scorecard function, `synpmx_scorecard()`, which already does not
   know what produced the dataset it scores. A row that does not apply to a
   given generator is marked as not applicable and says so in its result, which
   is what the card already does for the three rows that need a `pmx_settings`
   attribute and read `"unavailable"` without one. So a goodness-of-fit row
   would be a row on the one card, reading not applicable for AVATAR and PCA,
   rather than a per-generator card. The document follows: `avatar-scorecard.Rmd`
   becomes a generic scorecard document rather than gaining a sibling. Recorded
   as `REV-043`, with `REV-040` as the function half; it is not work this
   generator blocks on.

2. **A dose effect on clearance: no.** The candidate set stays the five linear
   models in `.pk_models`, and dose and arm stay apparatus. A source with
   dose-dependent PK is generated as though it were linear, and that is a
   limitation the documents state rather than a case the candidate set grows to
   cover. Fitting a dose effect would mean a candidate set that is no longer
   five closed-form models, which is the boundary the `nlmixr2` section draws
   and the reason generation needs no solver.

3. **Occasion-varying dosing: no.** Occasions are pooled. Estimation reads the
   recorded dosing history and the recorded times, and an `occasion` role is not
   read at all: no occasion-varying parameter is fitted and no per-occasion fit
   is offered. A dataset like `mavoglurant`, whose clock resets within occasion,
   is a dataset whose recorded times are already dose-relative, and putting it
   on a single cumulative axis before calling is the caller's work rather than
   the function's. This is the same fork as `nominal_time`: a statement about
   the protocol only the caller can make.
