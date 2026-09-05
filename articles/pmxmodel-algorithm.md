# The synpmx PMX Model Algorithm

This article is the full specification of
[`synpmx_model()`](https://iamstein.github.io/synpmx/reference/synpmx_model.md),
which fits a small set of linear pharmacokinetic (PK) models to a study
and generates new subjects by simulating from the fit. Use of the
algorithm is demonstrated in
[`vignette("pmxmodel-demo")`](https://iamstein.github.io/synpmx/articles/pmxmodel-demo.md).

**The aim is a simple model that roughly reproduces the properties of
the study, not one that characterises it.** The synthetic profiles it
generates are for prototyping analysis code. The goodness of fit of the
selected model is not assessed:
[`model_report()`](https://iamstein.github.io/synpmx/reference/model_report.md)
says what the search had to go on, not whether the model describes the
source. The candidate set is five linear models and the covariate model
is allometric scaling or nothing, which is too little to answer a
scientific question, so the fitted parameters are not estimates to
report. The fitted object prints that warning with itself, because its
output looks exactly like the output of a real population analysis.

No formal privacy guarantee is offered, although no patient’s measured
value reaches the output.

## Three Exported Functions

| Function | Reads patient data | Needs `nlmixr2` | Returns |
|----|----|----|----|
| `synpmx_model_estimate(data, roles, ...)` | yes | yes | `pmx_fitted_model` |
| `synpmx_model_generate(fitted_model, n_subjects, seed)` | no | no | data frame |
| `synpmx_model(data, roles, n_subjects, seed, ...)` | yes | yes | data frame |

[`model_report()`](https://iamstein.github.io/synpmx/reference/model_report.md)
inventories everything the fitted object carries,
[`model_candidates()`](https://iamstein.github.io/synpmx/reference/model_candidates.md)
returns the comparison table the selection was made from, and
[`model_parameters()`](https://iamstein.github.io/synpmx/reference/model_parameters.md)
returns the fixed effects, the between-subject covariance matrix and the
residual error.

## Two Time Axes, and `nominal_time` Is Required

**Estimation reads recorded times.** A population PK fit is a statement
about time after the dose that was actually given, so the subject-level
dosing records go to the fitter as recorded, reductions and skipped
cycles included. Fitting against the planned schedule where a patient’s
dose was reduced would push the drop in concentration that followed into
clearance, and the model would report a population that eliminates the
drug faster than the real one.

**Everything else reads the nominal grid.** The dosing model, the visit
model, and every signal in Steps 1 and 2 sit on the protocol’s planned
times. On recorded times a “median profile” is one point per time,
because no two patients share a clock reading, and a study sampled at
six times looks like a study sampled at two hundred. A shape, a peak
position and a count of sampling times are all statements about the
protocol.

Inferring that grid is a statement about the protocol only the caller
can make, so `nominal_time` is required rather than derived.

## Overview of Algorithm

1.  **Classify the endpoints.** Four signals decide which endpoint is
    the drug concentration. Everything else continuous becomes a
    pharmacodynamic (PD) endpoint.
2.  **Detect the design.** The route of administration and the sampling
    richness prune the set of PK models that will be offered to the
    fitter.
3.  **Estimate the parameters.** Fit each surviving PK model as a
    population nonlinear mixed-effects model through `nlmixr2` and
    select on AIC, with allometric scaling folded in where a weight-like
    covariate is declared. Fit each PD endpoint’s time course by least
    squares.
4.  **Fit a dosing model and a visit model per arm**: a planned dose
    schedule with rates for reduction, interruption and discontinuation,
    and the probability of a visit at each nominal time.
5.  **Store the fit**: the structural model, the parameters, the arm
    models and the schema, and no per-subject quantity.
6.  **Generate new subjects.** Draw covariates and random effects, draw
    the dose schedule, evaluate the profile against that schedule, add
    residual error and put the assay limit back.

Steps 1 to 4 read the source. Step 6 reads only what Step 5 stored,
which is why generation runs on base R and needs neither `nlmixr2` nor a
compiler.

## Step 1: Identify the Drug Concentration Endpoint

[`pmx_endpoint_types()`](https://iamstein.github.io/synpmx/reference/pmx_endpoint_types.md)
already infers continuous, integer, binary and ordinal from the values,
and that inference is reused unchanged. What this generator adds is a
second classification on top of it. Four signals, computed for every
endpoint that is a time course:

1.  **Compartment.** The `cmt` role puts the endpoint in the compartment
    the doses go to, or one above a dosing compartment nobody observes,
    which is the depot-and-central convention.
2.  **Post-dose only.** Observations before the first dose are absent,
    or sit at or below the censoring limit, in most subjects.
3.  **Shape.** The median profile within the richest dose interval rises
    to a maximum and declines without rising again.
4.  **Dose proportionality.** Between the highest and lowest dose level,
    the ratio of median subject maxima is within a factor of two of the
    dose ratio.

**Signals 2 and 4 are required; 1 and 3 break ties.** Being absent
before the first dose and scaling with the dose are properties only a
drug concentration has. A compartment number is a convention the
dataset’s author chose, and a rise-and-fall shape is one many biomarkers
also have.

Signal 4 is not computable on a single-level study, nor on one dosed by
body weight. A study that gives every patient a slightly different
amount has no dose levels: reading each amount as its own level compares
a median of one patient against a median of one patient over a dose
ratio near 1, which almost any endpoint passes. Where the distinct
amounts outnumber half the cohort the signal reports “not computable”
and signals 1 to 3 decide.

Where no endpoint passes, or two pass and neither tie-break separates
them, the function errors and names `endpoint_roles` as the way through.
That is the inference-versus-declaration fork, answered the way
`dose_covariate` answers it: infer where the data settles the question,
offer the declaration as an override.

``` r

model_report(fit)$endpoints$signals
#>   endpoint compartment post_dose shape proportional
#> 1       cp          NA      TRUE  TRUE           NA
#> 2      pca          NA      TRUE FALSE           NA
```

**What this cannot do.** On a study with one continuous endpoint and one
dose level the classification rests on signal 2 alone, and signal 2
passes for any endpoint with no pre-dose observation, including a
white-cell count. The remedy is `endpoint_roles`, and the evidence is in
the table above: a row reading “not computable” under `proportional` is
the classification saying how little it had to go on.

## Step 2: Detect the Design

Detection prunes the candidate set; the fit picks among what survives.
Neither property below decides the model on its own.

### Route

A `rate` role carrying a nonzero value is an infusion, and
`1cmt_infusion` is the only candidate. Otherwise one property separates
the two remaining routes, and it is not how many patients peak early.

A drug given by mouth cannot be in the blood at the moment it is
swallowed. So a concentration observed at time zero after a dose is
intravenous, and a median profile that rises before it falls is oral. A
profile that declines from a first sample drawn some time after the dose
is both of those at once, an intravenous bolus or an oral dose whose
absorption finished before anybody looked, and no amount of reading the
source separates them. There both routes are offered and AIC settles it.

Counting subjects instead is the mistake this replaces. In a study whose
sampling starts well past the peak, every subject’s profile declines
from its first point, and a per-subject vote reads an oral study as
intravenous: it counts the sampling schedule rather than the drug. The
share of subjects whose own profile rises is still reported, because it
says how much of the cohort the median profile speaks for, but it
decides nothing.

``` r

model_report(fit)$design$reason
#> [1] "the median profile rises to a peak at 9 before declining, and 31% of subjects do too"
```

### Sampling Richness

**Across the cohort**, distinct nominal times after a dose holding an
observation. Below six the fit still runs and warns: a one-compartment
model is not identifiable from that sampling, so its parameters stay
near their starting values and the simulated profiles are a plausible
shape rather than this study’s. There is no simpler structural model to
drop to — the default candidate set is already one compartment.
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
and
[`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md)
carry sparse sampling without fitting a structure to it. No post-dose
observation at all is an error, since there is then no curve to fit.

**Per subject**, distinct nominal times after a dose within the richest
dose interval. Four or more, with two after the median peak, is enough
to identify a distribution phase. A sample *before* the peak is required
only where there is an ascending limb to sample, since an intravenous
bolus peaks at the dose.

This count decides nothing. It is reported, and
[`model_report()`](https://iamstein.github.io/synpmx/reference/model_report.md)
says so when the sampling would support a two-compartment model, because
asking for one is your move rather than the search’s.

### One Compartment by Default

The default is one-compartment: `1cmt_iv`, `1cmt_oral` or
`1cmt_infusion`. Where the route is ambiguous both are fitted and AIC
chooses, which is the only case the default runs more than one fit.

A distribution phase is a refinement of a shape the one-compartment
model already has, and this generator exists to make simulated profiles
resemble the source study rather than to characterise it. The cost is
not small: a two-compartment fit takes several times as long, and on
data a one-compartment model describes it spends that time against a
flat likelihood and reports the worse AIC anyway.

**It remains available.** `pk = "2cmt_oral"` or `pk = "2cmt_iv"` forces
the model and skips the search, and the generator simulates from it
exactly as it does from any other, since the closed-form solution is
already there.

## Step 3: Estimate the PK and PD Parameters

**The candidates are the structural PK models Step 2 left standing.**
Each one is fitted as a population nonlinear mixed-effects model through
`nlmixr2`, which estimates the fixed effects, the between-subject
covariance matrix and the residual error, and where more than one
candidate is offered they are compared on AIC. The PD endpoints are
estimated separately, and not as mixed-effects models.

The candidate set is at most the five models in the package’s
closed-form set, and never more. A candidate the fitter could estimate
and the generator could not simulate would be a model that fits and then
generates nothing, so the two lists are one list. That is also why no
dose effect on clearance is fitted: it would mean a candidate set that
is no longer five closed-form models, and generation would need a
solver.

**A subject the endpoint was never measured in is not fitted.** It adds
no term to the likelihood, and a record carrying doses and no samples is
at best ignored by the solver. Its dosing and its visits are still real
and still reach the models of Step 4, which read the study rather than
the table handed to the fitter. How many subjects that left out is said
before the fit, and it warns where it takes the fitted cohort under
`min_subjects`.

Each candidate carries log-normal between-subject variability on every
structural parameter and a proportional residual error. An endpoint
holding values at or below zero after censoring is handled gets an
additive error instead, and the substitution is recorded on the fit.

**Censored observations are imputed, and the fit is told how many.** A
value below the assay limit is replaced by a uniform draw inside the
censoring region, a draw rather than a fixed LLOQ/2, which would swap
one artificial spike for another, and the boundary is put back when the
synthetic data is emitted. The same imputation is what lets the visit
model, the PD shapes and the covariate model read a latent value instead
of a stack of identical boundary substitutions. It is an assumption all
the same, and its weight is the share of the endpoint that carries it,
so that share is reported with the fit. Where much of an endpoint sits
below the limit, the fitted parameters are substantially a statement
about the draw rather than about measurements.

Starting values come from a non-compartmental reading of the source
rather than from a guess, because a population fit started far from the
answer either converges slowly or reports the starting values back. The
terminal slope is fitted by log-linear regression over the points after
the peak, and everything leans on it. Clearance is the dose over the
area **extrapolated to infinity**, not over the trapezoid alone, which
on a study sampled to four half-lives understates the area by about a
tenth and overstates clearance by the same. Volume is clearance over the
terminal slope, which is the quantity the terminal phase identifies,
since dose over the peak is not a volume of any kind for an oral dose.
Absorption is solved from where the peak falls: for a one-compartment
oral model the peak sits at `log(ka/ke)/(ka - ke)`, one equation in one
unknown once the slope has given `ke`. For a two-compartment model the
steady-state volume comes from the mean residence time, with the
absorption mean taken back off for an oral dose, and is split between
the central and peripheral compartments.

**The default fits one model and stops.** Not a search that happens to
have one candidate: one model, chosen by the route detection above, with
allometric scaling folded into it where a weight-like covariate is
declared. `pk` is how you buy accuracy with time.

``` r

# One fit, the default.
synpmx_model_estimate(data, roles)
# One fit, your choice of model.
synpmx_model_estimate(data, roles, pk = "2cmt_oral")
# Two fits, compared on AIC.
synpmx_model_estimate(data, roles, pk = c("1cmt_oral", "2cmt_oral"))
```

**How long it takes is set by the dose records, not the patients.**
Every likelihood evaluation superposes one contribution per dose per
subject, so the cost scales with the number of dose events rather than
with cohort size. A single-dose study of a few dozen patients fits in
seconds. A daily regimen carrying scores of dose records per subject
over a cohort several times larger is two orders of magnitude more work,
and a call that has not returned after half an hour is that cost rather
than a fault.

**Where there is a search, selection is on AIC**, and the estimation
method is `focei` for that reason. SAEM’s log-likelihood is a
Gaussian-quadrature step run after the fit, and at phase 1 cohort sizes
it returns a non-finite value, so a search over two candidates has
nothing to compare. `estimation = "saem"` remains available for a study
large enough to give it a likelihood. A candidate whose AIC is not
finite is recorded as not converged whichever method produced it.

A candidate that fails drops out carrying its reason and stays in the
table, so that a search which came down to one survivor does not look
like a search that had one candidate. Where nothing converges the
function errors rather than returning the least bad fit.

``` r

model_candidates(fit)
#>       model converged      aic note
#> 1 1cmt_oral      TRUE 895.9053
```

### Covariates

The default `"auto"` applies allometric scaling on clearance and volume
where a weight-like covariate is declared, and fits nothing else. The
exponents are the standard 0.75 and 1 rather than estimated ones, and
the effect is **asserted rather than tested**: it is folded into the one
fit, not compared against a model without it. Testing it would double
the cost of the whole call. `covariate_effects = "none"` switches it
off.

The covariate is recognised by name, `wt`, `weight`, `bw` and the like,
and must be numeric and positive. There is no way to recognise a body
weight from its values alone.

**The cost of this default is explicit.** A covariate that influences
the real profiles and is not in the model is generated independently of
them, so the synthetic data carries no relationship between the two.
[`model_report()`](https://iamstein.github.io/synpmx/reference/model_report.md)
reports the correlation between each declared covariate and the
individual random effects, which is where an unmodelled relationship
shows up.

### Pharmacodynamics

Every remaining continuous endpoint is fitted against a constant, a
linear and an exponential time course and selected on AIC, with
between-subject variability on the baseline. These shapes have no
concentration term, so a PD endpoint driven by exposure is reproduced as
a time course that happens to resemble the average subject’s response.

**The between-subject spread and the residual are read around each
subject’s own curve, not around the population one.** Taking the
residual as the spread of every point about the typical curve puts all
of the between-subject variation into it, and generation then emits that
as independent noise per observation. On a study where subjects share a
shape and sit at different levels, almost the whole of the structure
comes back as scatter and no synthetic subject has a profile at all.
`.pd_profile()` is linear in `baseline` for every shape, so a subject’s
own baseline is a least-squares projection needing no second optimizer,
and it is the same quantity generation draws: the between-subject term
is the spread of those baselines on the log scale, and the residual is
what is left around each subject’s own curve.

**A shape that fails to fit stays in the table.** The exponential is
tried from four start sets rather than one, because the single
median-based set fails on any response that falls and then recovers,
which is the shape a turnover endpoint has. Where none of the four
converges the candidate is reported with `converged` false and the
solver’s message rather than dropped, so a flat line winning on AIC
against two other flat lines is visible as that.

All three shapes are fitted by least squares rather than through
`nlmixr2`, so all three together cost no measurable time, and **the one
compiled fit in the whole call is the PK one.**

Binary and ordinal endpoints are not fitted at all. They are drawn from
the level frequencies their arm holds at each nominal time, which is
what the visit model already does for attendance.

## Step 4: Fit the Dosing Model and the Visit Model

Everything that is not the concentration-time curve comes from the
models
[`synpmx_pca_summarize()`](https://iamstein.github.io/synpmx/reference/synpmx_pca_summarize.md)
already builds, unchanged in what they represent. Dose reductions,
interrupted cycles, discontinuation and missed visits are modelled here
exactly as they are there, and so are arm sizes, the covariate
distributions, the censoring boundary and the schema. The two generators
call the same code rather than a copy of it; what differs is only the
adapter that hands it the grid, and this generator writes its own over
the nominal times the source holds.

Per arm:

| Model | What it holds | Drawn at generation as |
|----|----|----|
| Planned schedule | The nominal dose times enough of the arm reached, and the modal amount at each cycle among patients still on their starting dose | The cycle grid every subject starts from |
| Dose ladder | The levels patients dropped to, as ratios to their own starting dose, built from within-patient decreases | The amount multiplier in force at a cycle |
| Reduction rate | Discrete-time hazard of stepping down a level | Decided before the cycle is dosed |
| Interruption rate | Discrete-time hazard of skipping a cycle without ending treatment | Decided at the cycle |
| Discontinuation rate | Discrete-time hazard of stopping treatment | Decided after the cycle is dosed |
| Visit model | Per endpoint and per retained nominal time, the fraction of the arm holding an observation there | Attendance, drawn per visit |

A study where nobody reduces, skips or stops early has all three rates
at zero and one level, and the model then reproduces the planned
schedule exactly. No detector decides which kind of study this is,
because the rates already say it.

A grid cell is kept only where at least `min_arm_patients` distinct
patients hold an observation there. A nominal time one patient attended
is that patient, and generating from it would put them back.

## Step 5: Store the Fit

Estimated by `nlmixr2`: the selected structural model, the candidate
comparison table, fixed effects, the between-subject covariance matrix,
the residual error, the covariate effects that survived, and the PD
shape and parameters per endpoint.

Not estimated: the per-arm dosing model and its three rates, the visit
model, arm sizes, the covariate distributions, the censoring boundary,
the schema and the roles.

**No individual estimates.** Empirical Bayes estimates are per-subject
quantities, and a fitted model that carried them would be writing out a
description of each real patient. They are not stored on the object.
Generation draws random effects from the covariance matrix instead, and
the correlation report above computes from them and keeps only the
correlation.

``` r

model_report(fit)
#> What this fitted model carries
#> 
#> Estimated by nlmixr2
#>   structural model   1cmt_oral 
#>   fixed effects      cl 0.1366, v 8.176, ka 0.6119 
#>   between-subject    cl 0.243, v 0.0868, ka 0.685 (as SD on the log scale)
#>   residual error     proportional 0.21 
#>   covariate effects  cl ~ (wt/70)^0.75, v ~ (wt/70)^1.00 
#>   pd shapes          pca: exponential 
#>   not emitted below  cp 0.3; pca 4.5 (half the smallest value reported)
#> 
#> Summarized from the source, not estimated
#>   cohort             32 patients in 1 arm(s)
#>   visit model        22 grid cells over 2 endpoint(s)
#>   dosing model       1 planned cycle(s) per arm | no reductions, skips or early stops 
#> 
#> How the concentration endpoint was decided
#>   endpoint           cp (inferred) 
#>  endpoint compartment post_dose shape proportional
#>        cp          NA      TRUE  TRUE           NA
#>       pca          NA      TRUE FALSE           NA
#>   design             the median profile rises to a peak at 9 before declining, and 31% of subjects do too 
#> 
#> Covariate against the individual random effects
#>  covariate parameter correlation
#>        age        ka       -0.29
#>        sex         v       -0.21
#>        age        cl        0.19
#>         wt        ka        0.12
#>         wt        cl       -0.10
#> 
#>   A covariate that moves with a random effect and is not in the model above
#>   is generated independently of the profiles, so the synthetic data carries
#>   no relationship between them. `synpmx_avatar()` keeps those relationships
#>   without modelling them.
```

## Step 6: Generate New Subjects

Per synthetic subject: assign an arm keeping the source arm shares, draw
covariates from the arm’s covariate model, draw random effects from the
covariance matrix, apply the covariate effects to the typical
parameters, draw the dose schedule from the arm’s dosing model, draw the
visits attended from the visit model, evaluate the profile at the
attended times against the drawn schedule, add residual error, apply the
censoring boundary and emit.

**The order of the two draws is fixed.** The schedule is drawn first and
the profile is computed from it, so a synthetic subject who steps down a
dose level has a lower exposure from that cycle on, and one who skips a
cycle has the trough that implies. This is the one thing the generator
does that neither
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
nor
[`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md)
does. Drawing a profile and then a schedule would put a subject’s dose
reduction in the dosing records and not in the concentrations.

**The proportional multiplier is lognormal.** Written as `1 + N(0, cv)`
it is non-positive with probability `pnorm(-1 / cv)`, which is 2.3% of
draws at a coefficient of variation of 0.5 and 7.7% at 0.7, and each of
those is a concentration clamped to zero. A coefficient of variation
that high is a statement that the structural model does not describe the
data, and it belongs in the output as a wide band rather than as a
scatter of zeros on the floor of a log axis.
`exp(N(0, sqrt(log(1 + cv^2))))` has the same coefficient of variation
and a median of one, so the spread the fit estimated is preserved and no
draw reaches zero.

**A value is floored at the smallest one the study reported, halved.**
The residual is not the only thing that can put a synthetic value below
anything the assay could return: a one-compartment profile evaluated
late in a long dose interval underflows on its own. A study that
declares a censoring column says where its assay stopped, and that
boundary is put back. A study that declares none still had an assay, and
its smallest reported value is the only evidence of where the limit sat,
so half that value becomes a floor. An endpoint that reports a zero is
given no floor, since a zero is a value a floor would contradict.
[`model_report()`](https://iamstein.github.io/synpmx/reference/model_report.md)
carries the floor per endpoint.

**A floor that catches a large share of the output is a warning, not a
repair.** Raising a value is still changing it, and a floor doing that
to much of a dataset is hiding a fitted model that does not describe the
low end of the data. What the floor caught is recorded on the generated
dataset as the `pmx_floored` attribute, and above one observation in
twenty the generator says so.

## Gates

Only a study the generator cannot fit at all is refused. A study it can
fit badly is fitted, and told about.

| Gate | Threshold | Result |
|----|----|----|
| Cohort size | 20 subjects | Errors. A covariance matrix fitted to a handful of subjects describes those subjects. [`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md)’s floor is 10 and this one is higher because a parameter estimate concentrates on its cohort faster than a score does. |
| Cohort time coverage | 6 distinct nominal times after a dose | Warns and fits. Below it a one-compartment model is not identifiable and the parameters sit close to their starting values, which is the caller’s call to accept. |
| No observation after a dose | — | Errors, naming which of the role columns the empty count came from: nothing selected as an observation, nothing selected as a dose, or the two never meeting in one subject. |
| Arm size | 3 patients | Warns and drops those patients before anything is fitted, so the arm is absent from the model and from the data generated from it. Inherited from the dosing and visit models, which are summaries of an arm: an arm of one or two has no rates to pool. |
| `nominal_time` undeclared | — | Errors. The grid is a statement about the protocol only the caller can make. |
| No grid cell shared | `min_arm_patients` | The cell is dropped. A nominal time one patient attended is that patient. |
| No PK endpoint identified | — | Errors and names `endpoint_roles`. |
| No candidate converged | — | Errors rather than returning the least bad fit. |
| A fitted model as a public input | — | [`pmx_prior()`](https://iamstein.github.io/synpmx/reference/pmx_prior.md), [`synpmx_prior()`](https://iamstein.github.io/synpmx/reference/synpmx_prior.md) and [`synpmx_calibrated()`](https://iamstein.github.io/synpmx/reference/synpmx_calibrated.md) refuse a `pmx_fitted_model`, because its parameters were estimated from the confidential study. |
