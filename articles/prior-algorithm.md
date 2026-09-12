# The synpmx_prior Algorithm

[`synpmx_prior()`](https://iamstein.github.io/synpmx/reference/synpmx_prior.md)
simulates a dataset from a public structural model and a public trial
design. It reads no confidential data, so there is nothing to protect
and no budget to spend: this is `epsilon = 0`, the strongest guarantee
the package offers, and the weakest fidelity.

**Out of scope, stated once and bluntly.** The output is exactly as good
as the model asserted. It is not an estimate of anything, it contains no
information about any real study, and the catalogue of models it can
express is small and fixed — the section on what it cannot say is the
one to read before relying on the output. If your question is “does my
code handle this dataset shape?”, this is the right tool. If it is “what
would this study look like under my model?”, write the model in `rxode2`
and simulate it there.

Each numbered step below is an operation the function performs, in
order.

## Overview of Algorithm

| Step | Operation                                    | Reads data? |
|------|----------------------------------------------|-------------|
| 1    | Validate the public structural model         | No          |
| 2    | Validate the public trial design             | No          |
| 3    | Check the output declaration can be filled   | No          |
| 4    | Assign subjects to cohorts                   | No          |
| 5    | Draw one parameter vector per subject        | No          |
| 6    | Place dose and observation events on a clock | No          |
| 7    | Evaluate the profile by superposition        | No          |
| 8    | Apply residual error                         | No          |
| 9    | Apply dropout and the assay limit            | No          |
| 10   | Draw and merge baseline covariates           | No          |
| 11   | Name the columns from the declaration        | No          |

Every row reads no data, which is the whole point: there is no step at
which a confidential value could enter.

## Step 1: Validate the Public Structural Model

[`pmx_structural_model()`](https://iamstein.github.io/synpmx/reference/pmx_structural_model.md)
takes a PK form, a named vector of typical parameter values, an optional
PD form, and a required `source`. The `source` is not decoration: a
model without data-independent provenance cannot be treated as a public
input, and the constructor refuses without one.

The catalogue is small and fixed. Rather than restate it here, where it
could drift, these two chunks ask the package:

``` r

pmx_structural_model(pk = "list the options", typical = c(cl = 1, v = 1),
                     source = "illustrative")
#> Error in `match.arg()`:
#> ! 'arg' should be one of "1cmt_iv", "1cmt_oral", "1cmt_infusion", "2cmt_iv", "2cmt_oral", "1cmt_mixed", "2cmt_mixed"
```

``` r

pmx_structural_model(pk = "1cmt_iv", typical = c(cl = 1, v = 1),
                     pd = "list the options", source = "illustrative")
#> Error in `match.arg()`:
#> ! 'arg' should be one of "none", "constant", "linear", "exponential"
```

Every form is a closed-form analytic solution. There are no ODEs: `rx`
accepts an `rxode2` model and warns that it is never used, because
silently returning an analytic curve for a user-supplied ODE model would
be a fidelity claim the package cannot keep.

Each form requires its own parameters, and the error says which:

``` r

pmx_structural_model(pk = "2cmt_oral", typical = c(cl = 1), source = "x")
#> Error:
#> ! `typical` is missing required parameters: v, q, v2, ka.
```

The two `_mixed` forms exist for a study dosed by two routes. They
resolve to a single-route form per dose — intravenous or extravascular —
and `f` is the bioavailability applied to the extravascular one. An
intravenous dose carries no bioavailability term, because it is all of
the dose.

Two further inputs are assertions about variability:

- `iiv`, a CV per parameter, applied as lognormal between-subject
  variability. The default is `c(cl = 0.3, v = 0.2)`.
- `residual_cv`, proportional residual error. The default is `0.15`.

## Step 2: Validate the Public Trial Design

[`pmx_trial_design()`](https://iamstein.github.io/synpmx/reference/pmx_trial_design.md)
describes the protocol, and also requires a `source`. Dose can be given
in one of two ways, and never both:

- `dose_levels` with `cohort_sizes`, for parallel groups.
- `dose_escalation`, for within-subject escalation. One sequence applies
  to every subject; a list of sequences gives one per cohort.

`sampling` is where most of the expressive power sits. A bare vector
applies the same profile after every dose. A **list** with one entry per
dose says what each dose carries, and `NULL` means no samples after that
one — which is how a repeated-dose study that samples richly on the
first and last days and takes a single trough between them is written
down. `n_doses` and `dose_interval` place the doses, `duration` makes
them infusions, and `visit_window` is the fraction by which actual
sample times deviate from nominal.

``` r

rich <- c(0, 0.5, 1, 2, 4, 8, 24)
design <- pmx_trial_design(
  dose_levels = c(100, 300), cohort_sizes = c(12, 12),
  sampling = list(rich, 0, NULL, rich),
  n_doses = 4, dose_interval = 24,
  source = "illustrative protocol"
)
design
#> Public trial design
#>   doses: 100, 300  (n = 12, 12)
#>   dose times: 0, 24, 48, 72
#>   sampling: 
#>     dose 1: 0, 0.5, 1, 2, 4, 8, 24
#>     dose 2: 0
#>     dose 3: none
#>     dose 4: 0, 0.5, 1, 2, 4, 8, 24
#>   source: illustrative protocol
```

## Step 3: Check the Output Declaration Can Be Filled

`roles` is a
[`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md)
declaration naming the output columns, defaulting to
[`pmx_generated_roles()`](https://iamstein.github.io/synpmx/reference/pmx_generated_roles.md).
Every generator in the package takes one and writes its output under
those names, so a single declaration serves the generator,
[`compare_pmx_distributions()`](https://iamstein.github.io/synpmx/reference/compare_pmx_distributions.md)
and
[`synpmx_scorecard()`](https://iamstein.github.io/synpmx/reference/synpmx_scorecard.md).
The three that read a study take it as a description of the input; this
one reads nothing, so the same object prescribes the schema to generate
into.

Taking a declaration means honouring it, so a declaration naming a
column this generator cannot produce is refused rather than answered
with a table that fails
[`validate_pmx()`](https://iamstein.github.io/synpmx/reference/validate_pmx.md):

``` r

model <- pmx_structural_model(pk = "1cmt_oral",
                              typical = c(cl = 6, v = 35, ka = 1.5),
                              source = "illustrative")
synpmx_prior(model, design, pmx_roles(
  id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
  covariates = "WT"
))
#> Error:
#> ! `roles` declares columns this generator cannot fill.
#>   The covariate(s) WT are declared but not generated. These modes draw a
#>   covariate only from a public distribution, so pass `covariates =
#>   pmx_covariates(WT = pmx_covariate(...))` as well.
```

A role left undeclared simply gets no column.

## Step 4: Assign Subjects to Cohorts

Subjects are split across the declared cohorts in proportion to
`cohort_sizes`, with the rounding remainder distributed so the requested
total is exact. A cohort can end up empty if `n_subjects` is smaller
than the number of cohorts.

## Step 5: Draw One Parameter Vector per Subject

Each subject’s parameters are the typical values times a lognormal draw:

    p_i = typical * exp(N(0, sqrt(log(1 + cv^2))))

`typical` is the **median**, equivalently the geometric mean, following
the usual population-PK convention for a lognormal parameter. This also
has to match what the calibrated mode estimates: the released correction
is a mean on the log scale, so it targets the geometric mean. Centring
the arithmetic mean instead would leave a systematic `exp(sigma^2 / 2)`
gap between what is fitted and what is generated, which does not shrink
with the number of subjects or with epsilon.

One vector per subject, for the whole study. There is no inter-occasion
variability.

## Step 6: Place Dose and Observation Events on a Clock

Doses go at the times `n_doses` and `dose_interval` imply. Each
observation is assigned to the dose that qualifies it, and its time
after that dose is jittered by `visit_window`.

The jitter is applied to time **since the qualifying dose**, never to
absolute time. That keeps a predose sample exactly predose, keeps every
sample inside its own occasion, and guarantees a nonnegative time after
dose. Jittering absolute time does none of those things, and doing so
was a defect (`SIM-005`). A jittered time is also clipped so it cannot
reach the next dose.

## Step 7: Evaluate the Profile by Superposition

Every built-in model is linear in dose, so a multiple-dose profile is
the sum of single-dose profiles from each dose already given. That is
exact here rather than an approximation, and it is also what lets one
subject receive two routes.

The consequence is worth seeing, because it is the one fidelity property
this mode does hold exactly:

``` r

low  <- synpmx_prior(model, pmx_trial_design(
  dose_levels = 10, cohort_sizes = 1, sampling = c(1, 4, 12),
  source = "illustrative"), n_subjects = 1, seed = 1)
high <- synpmx_prior(model, pmx_trial_design(
  dose_levels = 1000, cohort_sizes = 1, sampling = c(1, 4, 12),
  source = "illustrative"), n_subjects = 1, seed = 1)
obs <- low$EVID == 0
round(high$DV[obs] / low$DV[obs], 6)
#> [1] 100 100 100
```

A hundred times the dose is exactly a hundred times the concentration.

## Step 8: Apply Residual Error

Each observation is multiplied by a lognormal draw with CV
`residual_cv`. Proportional only: there is no additive component and no
combined error model, so low concentrations are less noisy than a real
assay would make them.

## Step 9: Apply Dropout and the Assay Limit

`dropout` is the fraction of subjects who discontinue early — a protocol
assumption, not something learned. A subject who drops out has their
follow-up truncated at a uniform point, keeping at least two
observations.

`lloq` flags observations below the limit with `CENS = 1` and sets the
value to the limit, following the Monolix convention. It applies to the
PK endpoint only: a signed PD endpoint has no lower limit of
quantification, and treating one as censored would report a value above
the uncensored ones.

## Step 10: Draw and Merge Baseline Covariates

Declared covariates are drawn once per subject and merged onto the
finished table.
[`pmx_covariate()`](https://iamstein.github.io/synpmx/reference/pmx_covariate.md)
takes a public `range` for clipping plus a `median` and `cv` for the
distribution inside it;
[`pmx_covariates_reference()`](https://iamstein.github.io/synpmx/reference/pmx_covariates_reference.md)
supplies those for the usual covariates, and
[`pmx_covariate_reference_table()`](https://iamstein.github.io/synpmx/reference/pmx_covariate_reference_table.md)
says on what basis.

**The merge happens after the profile is computed, and that is the whole
story.** By the time a weight exists, the concentration beside it has
already been evaluated. Nothing links them.

## Step 11: Name the Columns from the Declaration

The table is built under the names in
[`pmx_generated_roles()`](https://iamstein.github.io/synpmx/reference/pmx_generated_roles.md)
and renamed to whatever `roles` declared, keeping the generator’s column
order. One collision is permitted, because NONMEM’s `CMT` routinely does
two jobs: a declaration naming the same column as both `cmt` and `dvid`
takes the compartment column, whose values already tell the endpoints
apart, and drops the separate endpoint key rather than writing it twice.

## What This Algorithm Cannot Express

The catalogue is a scope decision rather than a gap waiting to be
filled. The space of pharmacometric models is effectively unbounded and
every study wants something idiosyncratic; chasing that would slowly
produce a worse `rxode2`, which already does the job properly.

- **No covariate-parameter relationship.** No allometric exponent on
  clearance, no sex or biomarker effect. Covariates are columns beside
  the data, not inputs to it.
- **No exposure-driven PD.** This one is worth demonstrating, because
  the output looks plausible and is not. The PD forms are functions of
  time alone, so the PD endpoint does not respond to dose:

``` r

pkpd <- pmx_structural_model(
  pk = "1cmt_oral", pd = "linear",
  typical = c(cl = 6, v = 35, ka = 1.5, baseline = 100, slope = 2),
  source = "illustrative"
)
two_doses <- pmx_trial_design(dose_levels = c(10, 1000),
                              cohort_sizes = c(60, 60),
                              sampling = c(0, 4, 12, 24),
                              source = "illustrative")
generated <- synpmx_prior(pkpd, two_doses, n_subjects = 120, seed = 7)
rows <- generated$EVID == 0
medians <- function(endpoint) {
  keep <- rows & generated$DVID == endpoint
  round(tapply(generated$DV[keep],
               list(dose = generated$DOSE[keep],
                    time = generated$NTIME[keep]),
               stats::median), 2)
}
medians("cp")
#>       time
#> dose   0     4   12   24
#>   10   0  0.15 0.04 0.00
#>   1000 0 14.89 3.81 0.41
medians("pd")
#>       time
#> dose       0     4     12     24
#>   10   92.78 107.3 122.77 145.40
#>   1000 97.63 106.6 123.65 150.94
```

The concentration column separates the two dose groups by a factor of a
hundred. The PD column does not separate them at all: both arms rise
from the baseline on the same time course, and the differences are
residual error. A dose-response plot of this data would show nothing,
and an exposure-response analysis of it would be meaningless.

- **No inter-occasion variability**, no absorption lag or transit
  compartments, no time-varying parameters or enzyme induction, no
  additive or combined residual error, no ODE model.

## What Leaves the Source Data

Nothing. There is no source data. Every input is a declaration with a
`source` field recording where it came from, and the only way a real
study can reach the output is if someone fits a model to it and pastes
the estimates into
[`pmx_structural_model()`](https://iamstein.github.io/synpmx/reference/pmx_structural_model.md).
The function cannot detect that, and neither can anything downstream —
the `source` field is where that discipline is recorded, and keeping it
honest is the caller’s job.

## Where to go next

- [`vignette("prior-demo")`](https://iamstein.github.io/synpmx/articles/prior-demo.md)
  — one specification worked end to end, including the parts of it that
  could not be expressed.
- [Calibration](https://iamstein.github.io/synpmx/articles/synpmx-methods.html)
  — the same public model with a small privacy budget spent correcting
  its magnitude against a real study.
- [Building the structural
  model](https://iamstein.github.io/synpmx/articles/model-elicitation.html)
  and [Describing the
  trial](https://iamstein.github.io/synpmx/articles/data-elicitation.html)
  — producing the two public inputs without reading data.
- [The synpmx data generation
  algorithms](https://iamstein.github.io/synpmx/articles/synpmx-methods.html)
  — the other generation modes, and which one to use when.
