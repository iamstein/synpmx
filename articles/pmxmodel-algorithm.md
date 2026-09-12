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
report, even though the printed object looks exactly like the output of
a real population analysis.

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

`nlmixr2` is in `Suggests`, is loaded only by
[`synpmx_model_estimate()`](https://iamstein.github.io/synpmx/reference/synpmx_model_estimate.md),
and the candidate set is exactly the models the generator can evaluate
in closed form. A candidate the fitter could estimate and the generator
could not simulate would be a model that fits and then generates
nothing, so the two lists are one list. The vignettes read a stored fit
built by `scripts/build-model-fits.R`, so `R CMD check` and the pkgdown
site never compile a model.

## The Arguments

Beyond the roles and the seed:

| Argument | Default | Effect |
|----|----|----|
| `pk` | `NULL` | One of the five built-in models, forcing it; or several, which is how a search is asked for. |
| `pd` | `NULL` | Named vector of PD shapes per endpoint. Skips that search. |
| `pd_by_arm` | `FALSE` | Fit each PD endpoint’s shape per arm rather than once over the pooled cohort. |
| `endpoint_roles` | `NULL` | Names which endpoint is the drug concentration, overriding inference. More than one may be named. |
| `start_param` | `NULL` | Starting values for the population fit, keyed by endpoint where more than one concentration is fitted. The escape hatch where the non-compartmental read of the median profile starts the optimizer somewhere it cannot move from. |
| `covariate_effects` | `"none"` | `"none"` puts no covariate in the structural model; `"auto"` applies allometric scaling on clearance and volume where a weight-like covariate is declared. |
| `min_subjects` | `20L` | Warn and fit anyway below this cohort size. |
| `min_arm_patients` | `3L` | Warn and drop the patients in any arm below this, as [`synpmx_pca_summarize()`](https://iamstein.github.io/synpmx/reference/synpmx_pca_summarize.md) does. |
| `min_time_bins` | `6L` | Warn and fit anyway below this many distinct nominal times after a dose; no post-dose observation at all still refuses. |
| `max_fit_subjects` | `60L` | Fit the population model to this many subjects, drawn in proportion to the arms. The dosing, visit and covariate models read every subject. |
| `estimation` | `"focei"` | Passed to `nlmixr2`. |

**The default path performs exactly one population fit.** The arguments
above are the ways to spend more time for more accuracy: a search costs
one fit per candidate, and a second concentration endpoint costs one
more. Neither happens unless it is asked for.

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

**Repeated doses are written out before anything reads them.** A regimen
recorded as one record plus `ADDL` and `II` — “and 364 more like it,
every day” — is expanded by
[`pmx_expand_doses()`](https://iamstein.github.io/synpmx/reference/pmx_expand_doses.md)
on the way in, so that every dose the patient received is a row to the
dosing model, the derived time after dose and the estimation table
alike. Before this, every reader counted rows: `onc_sim` showed the
fitter 237 dose records standing for 71,180 doses and got back a
clearance six orders of magnitude off, from a fit that had moved and so
tripped no gate (`REV-051`). The fitter is handed the compact form
again, because `nlmixr2` reads it natively and is fastest on it, and the
generated study is folded back by
[`pmx_compress_doses()`](https://iamstein.github.io/synpmx/reference/pmx_expand_doses.md)
on the way out, so it comes back in the encoding its source used.

## Overview of Algorithm

1.  **Classify the endpoints.** Four signals decide which endpoint or
    endpoints are drug concentrations. Everything else continuous
    becomes a pharmacodynamic (PD) endpoint.
2.  **Detect the design.** The route of administration and the sampling
    richness prune the set of PK models that will be offered to the
    fitter.
3.  **Estimate the parameters.** Fit each surviving PK model as a
    population nonlinear mixed-effects model through `nlmixr2` and
    select on AIC, once per concentration endpoint. No covariate enters
    the structural model unless `covariate_effects = "auto"` asks for
    it. Fit each PD endpoint’s time course by least squares.
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

### More than one endpoint can be a concentration

A study measuring two drugs, or a parent and its metabolite, has two
drug concentrations, and each is fitted its own population model — its
own structural model, its own parameters, its own residual error. No
correlation between their random effects is estimated, which is a
limitation rather than a claim: a patient with a high parent clearance
is given an independently drawn metabolite clearance.

**Which doses drive which endpoint is a declaration.** By default every
dose record drives every concentration: one administration, two
analytes, which is what a parent and its metabolite are. Two *different*
drugs given together are the opposite case — each model has to see its
own drug’s doses and not the other’s, or a clearance and a volume come
back having absorbed an input the patient never received of that drug.

The data cannot settle it. A metabolite has no dose records of its own,
so “this endpoint’s doses are the ones marked `2`” and “this endpoint
has no doses” are the same table, and guessing wrong leaves a model with
no input at all. So it is declared, on the argument `routes` is declared
on:

``` r

pmx_roles(
  adm            = "ADM",
  routes         = c("1" = "extravascular", "2" = "extravascular"),
  dose_endpoints = c("1" = "drug A PK", "2" = "drug B PK")
)
```

In a NONMEM dataset the compartment already does this job, so declare
the column twice — `cmt = "CMT", adm = "CMT"` — and key the mapping on
compartment numbers.

Declared, each drug gets its own dose schedule per arm — its own
amounts, its own interval, its own ladder of reductions — and each
endpoint’s fit and generated profile read only that drug’s doses.
**Including the dosing history each observation is read against**: time
after dose, the dose interval and the first dose amount follow the
endpoint’s own drug, which is what the starting values, the route
reading, the sampling richness and the dose-proportionality signal are
all computed from. Without that half, a fit table split correctly still
starts from a volume scaled by the wrong drug’s milligrams. The
generated study then writes both drugs’ dose records back with their own
`adm` value and compartment.

Measured on a simulated combination — 100 mg weekly at `cl` 5, `v` 50,
`ka` 1.2 beside 900 mg fortnightly at `cl` 20, `v` 200, `ka` 0.5,
sixteen subjects — the second drug comes back `cl` 14.4, `v` 3.26, `ka`
9.34 with the doses pooled, and `cl` 19.3, `v` 166, `ka` 0.57 with them
separated, which is what that data fitted on its own gives.

Undeclared,
[`model_report()`](https://iamstein.github.io/synpmx/reference/model_report.md)
says so on a study with more than one concentration, which is the case
where it matters.

**Inference promotes a second endpoint only on signal 4.** Dose
proportionality is the positive evidence that something is a drug: a
concentration scales with the dose and a biomarker does not. Signal 2 is
necessary and not sufficient, and it is satisfied by any endpoint with
no pre-dose observation. The distinction matters because signal 4 is
often *not computable*, and the required test treats that as passing —
which is right for choosing one endpoint and wrong for promoting a
second. `warfarin`’s prothrombin activity and `onc_sim`’s tumour size
both pass signals 2 and 4 that way, and neither is a drug concentration;
requiring signal 4 to be `TRUE` rather than merely not `FALSE` leaves
both where they belong.

So a study whose two endpoints both scale with dose is classified with
two concentrations and no declaration. A study where the evidence is
absent on both still refuses and asks, which is what `onc_sim` does.

`endpoint_roles = list(pk = c("parent", "metabolite"))` declares it
where inference cannot. Declaring several is worth doing wherever it is
true: a concentration demoted to a pharmacodynamic endpoint is fitted a
time course with no dose term in it, so the generated values lose their
dose ordering entirely — measured at 1.05 times across a three-fold dose
range, against 2.7 in the source (`SIM-080`).

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

**A declared route is not detected at all.**
`pmx_roles(adm = "ADM", routes = c("1" = "iv", "2" = "extravascular"))`
names the administration column and says what its values mean, and the
two are declared together or not at all: which id is which route is a
convention of the dataset and cannot be read off the numbers, so a guess
would put every dose in the wrong compartment without failing.

Where a study declares both routes, the route is a property of each dose
record rather than of the study, and the candidate is `1cmt_mixed` — the
intravenous form for the doses given intravenously and the extravascular
form for the rest, summed. That sum is exact rather than an
approximation, because every model in the closed-form set is linear in
dose, and it is what lets one patient receive both. Bioavailability
scales the extravascular doses and not the intravenous ones, which is
what bioavailability means and why it is identifiable in a study dosed
both ways and not in one dosed a single way. It is fitted as a fixed
effect with no between-subject term: a design with one dose each way
identifies the contrast between the routes, not a per-subject
distribution over it.

**A mixed model is fitted with explicit compartments rather than the
closed form**, and bioavailability is the whole reason. Under `linCmt()`
the `f(depot)` directive does not reach the doses the data sends to the
depot: evaluated at known parameters it scaled the intravenous doses
instead, so `f` was estimated as its own reciprocal and clearance and
volume as `cl/f` and `v/f` — a fit whose own predictions tracked the
data while reporting parameters that meant something else (`SIM-077`).
Written as `d/dt(depot)` and `d/dt(central)`, the compartments are real,
`f()` binds to the one it names, and the predictions match the analytic
profile on both routes. Only a mixed study pays the solver cost; a study
dosed one way has no `f` to place and keeps the closed form.

Starting values for a mixed study are read one route at a time, for the
same reason the fit has to be: pooled, the cohort’s median profile is an
intravenous decline and an extravascular rise averaged together, and no
non-compartmental quantity read off it describes either route.
Disposition comes from the intravenous records, where `f` is not in the
way; absorption from the extravascular ones, where there is a peak to
find; and their two clearances give `f` itself a starting value, since
an extravascular read returns `cl/f`.

Detection below runs only where nothing was declared. It answers with
one route for the whole study, which is the right answer for a study
dosed one way and no answer at all for a study dosed two ways — the
`rate` test reads a single nonzero rate anywhere as “this is an infusion
study”, so a mixed study came back as an infusion model with no
absorption in it.

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

**The population model is fitted to at most `max_fit_subjects`
patients**, 60 by default, drawn in proportion to the arms with the
run’s seed. Fit time is linear in subjects and goes into the per-subject
inner loop — on `onc_sim`, 40 patients fit in four minutes and 200 in
twenty-five — while the parameters a synthetic study needs are settled
long before the sixtieth patient. Every other model in Step 4 still
reads the whole study, and the report states the count.

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

**The concentration is fitted with its censored observations censored.**
Where the study declares a censoring column, the rows below the assay
limit reach `nlmixr2` as censored – the limit in `DV`, the flag in
`CENS`, and the other bound in `LIMIT` where the study reports one – so
each contributes the probability of falling below the limit rather than
a value nobody measured. On 30 subjects of `case1_pkpd`, 45% of them
below the limit, that moves the proportional residual from 0.517 to
0.290 against imputing the same rows, and the residual is what sets the
scatter of every generated observation. It costs about a factor of 1.7
in fitting time there, and less where less of the endpoint is censored.

**Everything else reads an imputed value**, because nothing else has a
likelihood to put censoring in: the visit model, the covariate model and
the PD shapes are least-squares fits, and a stack of identical boundary
substitutions would bend each of them toward the limit. Those rows are
replaced by a uniform draw inside the censoring region – a draw rather
than a fixed LLOQ/2, which would swap one artificial spike for another –
and the boundary is put back when the synthetic data is emitted. The
share of each endpoint sitting below the limit is reported with the fit
either way.

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

**What it did take is reported.** Each candidate’s fitting time is a
column of the candidate table, the fit prints the total with itself, and
the run says so as it finishes, because the number a caller weighs a
rerun against is how long the last one took.

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
#>       model converged      aic seconds note
#> 1 1cmt_oral      TRUE 926.9152  10.507
```

### Covariates

The default `"none"` puts no covariate in the structural model. A
synthetic study does not need the relationship: covariates are generated
from the source’s own study-wide distributions either way, so a
clearance that moves with weight buys the generator nothing.

`covariate_effects = "auto"` applies allometric scaling on clearance and
volume where a weight-like covariate is declared, and fits nothing else.
The exponents are the standard 0.75 and 1 rather than estimated ones,
and the effect is **asserted rather than tested**: it is folded into the
one fit, not compared against a model without it. Testing it would
double the cost of the whole call.

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
between-subject variability on the baseline.
`pd = c(endpoint = "linear")` names the shape and skips the search for
that endpoint; a named shape that will not fit is an error rather than a
silent fallback. There is no minimum number of observations: a shape is
a candidate where it has **one residual degree of freedom**, which is
what stops a line being fitted exactly through two points and then
winning the comparison it was never tested by. Below any candidate the
endpoint is a constant at its mean, which is what an endpoint measured
once is, and the report says that rather than naming a search that did
not run.

These shapes have no concentration term, so a PD endpoint driven by
exposure is reproduced as a time course that happens to resemble the
average subject’s response.

**One shape for the whole cohort by default, and one per arm on
request.** The pooled fit predicts a single number for every arm, so a
synthetic patient’s dose does not reach their response: on `case1_pkpd`
the source’s mean PD after 1500 h runs 83, 68, 70, 126, 237, 341 across
placebo and five ascending doses, and the pooled shape answers 149 to
all six. `pd_by_arm = TRUE` fits the shape within each arm instead,
which is the same per-arm summary the dosing and visit models already
are and asserts no dose-response form – an arm is fitted on its own
observations, on the same three-candidate ladder, or not at all. On that
study it answers 47, 76, 77, 126, 217, 363. It is not the default
because it is a different claim about the study: the pooled shape says
these endpoints are a time course, the per-arm shape says each arm has
its own. An arm holding too little of an endpoint to fit keeps the
pooled shape, so no arm is left with nothing to generate.

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
`nlmixr2`, so all three together cost no measurable time, and **the only
compiled fits in the whole call are the PK candidates.**
[`model_report()`](https://iamstein.github.io/synpmx/reference/model_report.md)
breaks the wait down by fit: each concentration’s block carries what
that fit took, one line per candidate where the search compiled more
than one. A search that spent nine of its ten minutes in a second
candidate is a run whose answer is to name `pk`.

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
the residual error, the covariate effects that survived, the PD shape
and parameters per endpoint, and how long the fitting and the whole call
took.

Not estimated: the per-arm dosing model and its three rates, the visit
model, arm sizes, the covariate distributions, the censoring boundary,
the schema and the roles.

**Everything on the object is an input to Step 6, so printing it reports
all of them** — the dose reduction, skipped cycle and early stop rates
per arm, the attendance frequency behind a missed observation, how each
covariate is drawn, which cells are drawn from recorded values rather
than simulated, the assay limit and how many observations sit below it,
the floor nothing is emitted below, and the columns the generated table
will carry.
[`model_report()`](https://iamstein.github.io/synpmx/reference/model_report.md)
returns the same account as a list.

**No individual estimates.** Empirical Bayes estimates are per-subject
quantities, and a fitted model that carried them would be writing out a
description of each real patient. They are not stored on the object.
Generation draws random effects from the covariance matrix instead, and
the correlation report above computes from them and keeps only the
correlation.

``` r

model_report(fit)
#> Summarized from the source, not estimated
#>   cohort             32 patients in 1 arm(s)
#>                      all (32)
#>   dose changes       none
#>   visit grid         2 endpoint(s) at 16 nominal time(s), 22 slot(s) in all
#>   visit attendance   median 95% of an arm attends a slot (9% to 100%)
#>   covariates         wt lognormal, age lognormal, sex categorical, drawn
#>                      once for the whole study, independently of the
#>                      profiles
#>   columns emitted    id, time, ntime, dv, amt, evid, dvid, wt, age, sex
#> 
#> Values at the lower limit of what was observed
#>     cp                 0.3, half the smallest value seen, no assay limit
#>     pca                4.5, half the smallest value seen, no assay limit
#> 
#> Each non-PK continuous endpoint, fitted as constant, linear, or exponential
#>   pca                exponential
#>                        plateau          27.34
#>                        baseline         96.3
#>                        rate             0.09877
#>                        between-subject  0.146 (SD on the log baseline)
#>                        residual         additive 12.4
#>                        chosen on AIC from constant, linear, exponential
#> 
#> PK endpoint for the PopPK model
#>   cp                 inferred from the following data characteristics:
#>                        absent before the first dose
#>                        rises to one peak and comes back down within one
#>                          dose interval
#>   route              oral: the median profile rises to a peak at 9 before
#>                      declining, and 31% of subjects do too
#> 
#> The PopPK model
#> 
#> Estimated by nlmixr2
#>   structural model   1cmt_oral 
#>   fitted on          all 32 patients with a concentration
#>   fixed effects      cl 0.1353, v 8.115, ka 0.5796 
#>   between-subject    cl 0.267, v 0.204, ka 0.68 (as SD on the log scale)
#>   residual error     proportional 0.211 
#>   time to fit        10.5 s
#>   whole call         10.6 s, against 10.5 s in the fitter
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
| Cohort size | 20 subjects | Warns and fits. A covariance matrix fitted to a handful of subjects describes those subjects, and nothing downstream will say so: the scorecard asks whether the output copies anybody or changed the study’s shape, and a small-cohort fit does neither. [`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md)’s floor is 10 and this one is higher because a parameter estimate concentrates on its cohort faster than a score does. |
| Fit cap | `max_fit_subjects` = 60 | The PK model is fitted to a subset drawn in proportion to the arms; the dosing, visit and covariate models read every subject. Reported in [`model_report()`](https://iamstein.github.io/synpmx/reference/model_report.md). Cannot be set below `min_subjects`. |
| Cohort time coverage | 6 distinct nominal times after a dose | Warns and fits. Below it a one-compartment model is not identifiable and the parameters sit close to their starting values, which is the caller’s call to accept. |
| No observation after a dose | — | Errors, naming which of the role columns the empty count came from: nothing selected as an observation, nothing selected as a dose, or the two never meeting in one subject. |
| Arm size | 3 patients | Warns and drops those patients before anything is fitted, so the arm is absent from the model and from the data generated from it. Inherited from the dosing and visit models, which are summaries of an arm: an arm of one or two has no rates to pool. |
| `nominal_time` undeclared | — | Errors. The grid is a statement about the protocol only the caller can make. |
| No grid cell shared | `min_arm_patients` | The cell is dropped. A nominal time one patient attended is that patient. |
| Administration column | `adm` and `routes` together | Errors on either alone. What an administration id means is a convention of the dataset, and reading it wrong routes every dose to the wrong compartment silently. |
| PD shape candidacy | 1 residual degree of freedom | The shape is dropped from the comparison, not the endpoint from the study. Below every candidate the endpoint is generated as a constant at its mean. A shape fitted exactly through its own points has `AIC` `-Inf` and would win any comparison it entered. |
| No PK endpoint identified | — | Errors and names `endpoint_roles`. |
| No candidate converged | — | Errors rather than returning the least bad fit. |
| The between-subject terms did not move | 1% of the eta init, fixed effects moved | A line under the parameter block in [`model_report()`](https://iamstein.github.io/synpmx/reference/model_report.md), not the banner: the fixed effects are estimates, but the spread synthetic subjects get is the starting value rather than this study’s. |
| Two endpoints tie as concentrations | neither tie-break separates them | Errors, naming both and `endpoint_roles`. Naming several is accepted, and each is fitted its own model. |
| The fit did not move | 1% of every starting value | Warns, and prints `!! THE FIT DID NOT MOVE !!` above the parameter block in [`model_report()`](https://iamstein.github.io/synpmx/reference/model_report.md). A converged-looking fit is not necessarily an estimated one: where the optimizer takes no effective step, `nlmixr2` returns an objective, an AIC and a full table of starting values, and nothing else downstream contradicts them — the generator simulates from them and the scorecard passes, because it asks whether the output copies anybody or changed the study’s shape and a fit that never moved does neither. Each parameter’s start, estimate and percent change are listed, because “did not move” is a claim the reader has to be able to check. |
| A fitted model as a public input | — | [`pmx_prior()`](https://iamstein.github.io/synpmx/reference/pmx_prior.md), [`synpmx_prior()`](https://iamstein.github.io/synpmx/reference/synpmx_prior.md) and [`synpmx_calibrated()`](https://iamstein.github.io/synpmx/reference/synpmx_calibrated.md) refuse a `pmx_fitted_model`, because its parameters were estimated from the confidential study. |
