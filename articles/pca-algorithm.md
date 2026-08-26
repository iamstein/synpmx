# The synpmx PCA Algorithm

This article is the full specification of
[`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md),
which generates a synthetic pharmacometric (PMX) dataset from a fitted
model rather than from real patients’ values. Use of the algorithm is
demonstrated in
[`vignette("pca-demo")`](https://iamstein.github.io/synpmx/articles/pca-demo.md).

[`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md)
makes no formal privacy claim. It does not estimate a pharmacokinetic
(PK), pharmacodynamic (PD) or nonlinear mixed-effects model, so it
reports no clearance, no volume and no between-subject variability on
either. Use the synthetic datasets to develop and debug analysis code.
Do not use them to estimate parameters.

## What the principal components are used for

[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
also runs a principal component analysis (PCA), and uses it for
something else. There the components are a **distance metric**: the
coordinates decide which real patients are near which, those neighbours
donate, and the numbers written out are weighted averages of real
patients’ measurements. The output is assembled from real values.

Here the components are a **basis**. Each subject’s trajectories become
a short vector of scores, the scores become the response of a fitted
linear model, and new subjects come from drawing that model’s residual.
No donor index exists in the algorithm, and no number a patient measured
reaches the output. What leaves the source data is a mean, a scale, a
set of loadings, one mean score vector per arm, and a residual
covariance.
[`pmx_pca_report()`](https://iamstein.github.io/synpmx/reference/pmx_pca_report.md)
inventories all of it.

## Overview of Algorithm

1.  **Declare the columns.** The `roles` variable says which column is
    the id, time, dose, measurement, arm and covariate. Only named
    columns survive.
2.  **Build one feature vector per subject.** Place every row on the
    declared nominal grid, align each subject to their own first dose,
    and read each endpoint’s trajectory onto that grid, beside the
    baseline covariates.
3.  **Fill, standardize and decompose.** Replace censored values with a
    draw inside the censoring region, median-fill the cells a subject
    has no observation in, standardize every column, and take the
    principal components.
4.  **Model the scores against the arm.** Each arm gets its own mean
    score vector and its own residual covariance.
5.  **Draw new subjects** by adding a fresh residual to an arm’s mean,
    and invert the standardization to get a feature vector back.
6.  **Fit a dosing model and a visit model per arm**: the dose schedule
    the arm holds in common, and the probability of a visit at each
    nominal time.
7.  **Place the drawn values** on that skeleton, snapping a discrete
    endpoint back onto its levels and putting the assay limit back where
    the source had one.
8.  **Construct the full synthetic dataset** with new subject IDs and
    every column back in the shape the source had, and validate it.

Steps 1 to 4 and Step 6 read the source. Steps 5, 7 and 8 read only the
summaries those produced. The section below says why that split is
enforced rather than intended.

## Step 1: Declare the Meaning of the Columns

The package never guesses PMX roles from column names. Declare them with
[`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md),
as in
[`vignette("avatar-demo")`](https://iamstein.github.io/synpmx/articles/avatar-demo.md).

`strata` carries more weight here than it does in
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md).
It names the arms, and the arms are the groups Step 4 fits against, so a
column added to `strata` adds a group whose mean and covariance are
estimated separately.

**An arm of fewer than three patients is refused.** Its mean score
vector would be those patients and its covariance would be noise around
them.
[`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md)
stops and names the short arms rather than generating from one or two
people. `min_arm_patients` sets the floor. Pool the arm, drop the column
from `strata`, or exclude those patients before calling.

## Step 2: One Feature Vector per Subject

### `nominal_time` is required

[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
derives a visit grid when none is declared, and that derivation is a
masking mechanism: it decides where the visits are so that no patient’s
own time vector identifies them. Here the grid is the model’s axis.
Every feature is a cell on it, and dose rows and observation rows are
placed on it together. A grid inferred from elapsed time puts a sample
on the wrong side of a dose as soon as the doses were recorded as
actuals rather than as planned times, and choosing where a trough sample
belongs is a statement about the protocol rather than something an
algorithm can derive.

So
[`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md)
requires `nominal_time` and refuses without it. A study that does not
carry one needs the protocol’s planned times added as a column first.
Rows where it is missing are named in the error rather than filled.

Each subject is then aligned to their own first dose, so that time zero
means the same thing for everybody.

For each endpoint, the grid is **every nominal time in the study**.
There is no cap on its width.
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
reduces each endpoint to fifteen points because a distance metric gains
nothing from more, and a basis has no such ceiling: the number of
components is bounded by the subject count however many columns there
are.

Endpoint values are log-transformed where the endpoint is positive, by
the same rule
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
uses. Covariates join the same matrix, continuous ones as themselves and
categorical ones as one indicator column per level.

**Covariates and trajectories share one matrix on purpose.** A single
draw in Step 5 therefore carries the relationship between a covariate
and the profile, which a generator drawing covariates separately and
merging them on afterwards cannot do.

### The minimum patients per cell

A grid cell only a few patients ever reached is filled from the median
for everybody else. Such a cell describes those few patients and nobody
else, so it is dropped rather than modelled. The floor is
`min_column_patients`, defaulting to the larger of 3 and a tenth of the
cohort.

This is the only place the algorithm defends the late timepoints.
Follow-up thins out towards the end of a study, and the last visits are
where a cell held by one patient is most likely.

## Step 3: Fill, Standardize and Decompose

Missing cells take the column median, every column is centered and
divided by its standard deviation, and
[`prcomp()`](https://rdrr.io/r/stats/prcomp.html) gives the loadings and
the scores.

Median filling shrinks a subject towards the middle wherever they have
no observation. That is accepted rather than corrected: the cells where
it happens most are the ones Step 2 has already dropped, and an
iterative fill would add machinery for an effect the column guard has
bounded.

### How many components

The retained count is the smallest reaching `pca_variance` of the total,
capped at a fifth of the cohort. The cap is the part that matters. With
many components against few subjects the basis interpolates the source,
and a drawn subject can land on a real one; the variance rule alone does
not prevent that, because a small study reaches any variance target
quickly.

## Step 4: Model the Scores Against the Arm

Scores become the response of a linear model with the arm as the
predictor. Every arm gets its own mean score vector, and its own
residual covariance shrunk towards the pooled one in proportion to its
size.

`dose_term = "log"` is the alternative: one regression of the scores on
[`log1p()`](https://rdrr.io/r/base/Log.html) of each subject’s total
dose, which spends a single coefficient rather than a mean per arm and
generates at doses the study never ran.

**`"factor"` is the default, and the reason is measurable.** `"log"`
assumes the dose-response is log-linear and that every arm has the same
between-subject spread. A lower limit of quantification (LLOQ) breaks
both: the low arms sit on the floor, so the relationship is flat and
then rising rather than straight, and the arms on the floor have far
less spread than the top arm, so one pooled covariance smears them
upward and pulls the top arm in. On a six-arm study whose source PK
medians ran from the LLOQ to 0.53, `"log"` returned a range of 0.10 to
0.18 and `"factor"` recovered the arms. Use `"log"` where the arms are
genuinely log-linear and there are too many dose levels to spend a mean
on each.

## Step 5: Draw and Invert

A new subject’s scores are their arm’s mean plus a fresh multivariate
normal residual drawn through the Cholesky factor of that arm’s
covariance. Multiplying by the loadings, then by the column scales, then
adding the column centers, gives a feature vector in the original units.

The residual draw is what makes this generation rather than blending.
Nothing selects a real subject, and nothing weights one.

## Step 6: The Dosing Model and the Visit Model

Each generated subject is assigned an arm, keeping each arm’s share of
the cohort. The arm carries two models, and both are summaries of the
arm rather than facts about a patient.

### The dosing model

[`pmx_pca_dosing()`](https://iamstein.github.io/synpmx/reference/pmx_pca_dosing.md)
reports it: one row per arm and dose, with the time, the amount, and how
much of the arm stands behind it.

The schedule is the **modal** set of times and amounts among the arm’s
patients — the one the arm holds in common — never an individual’s. Two
columns say how safe that is. `share` is the fraction of the arm holding
it, so a schedule chosen by a bare plurality is visible as one rather
than passing as consensus. `distinct` is how many different schedules
the arm actually contained.

Every dosing event counts, including those carrying no drug. A placebo
arm records its administrations with an amount of zero, and reading only
the rows that carry drug would leave that arm with no dosing events at
all.

**A study recording actual dose times rather than planned ones will show
`distinct` close to its arm size and `share` close to `1/n`.** The
output then holds one schedule per arm where the source held one per
patient. The timing is right and the variety is gone. That is the
mechanism working, and it is the number to look at before trusting a
generated dataset for anything schedule dependent.

### The visit model

[`pmx_pca_visits()`](https://iamstein.github.io/synpmx/reference/pmx_pca_visits.md)
reports it: one row per arm, endpoint and modelled nominal time, giving
the probability that a generated subject in that arm has an observation
there. It is the fraction of the arm’s patients who did.

Attendance is drawn independently per visit from that probability. No
real patient’s visit set is reused, which is what the B1a row of
[`synpmx_scorecard()`](https://iamstein.github.io/synpmx/reference/synpmx_scorecard.md)
asks about for
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md).

### Why the timing survives

Dose rows and observation rows are placed on the **same declared nominal
grid**, which is what Step 2 requires `nominal_time` for. The
consequence is that the event structure making a dataset fittable comes
through at its source values: doses per patient, the time of the last
dose, the end of follow-up, the baseline sample before the first dose,
and the distribution of time after dose.

What does not come through is the variety of schedules, as above.

Neither B1a nor B1b nor C2 can be scored on this output. Those rows read
the AVATAR run record, and report `unavailable` here. The guarantee is
structural instead: no individual schedule and no individual visit set
exists in the model to be copied. B4a and B4b compare the finished
tables directly and are the rows that carry the claim.

## The Model Is the Only Thing Generation Reads

Fitting produces one object, and generation is a function of that object
and a subject count. No patient row is in scope while a synthetic
dataset is being written, and that is enforced rather than intended: the
generator’s only arguments are the model and the count, and a regression
test regenerates a dataset from the model alone and requires it to
match.

The object holds, and this is all of it:

| Group | What it is |
|----|----|
| Basis | Column means, column standard deviations, and the component loadings |
| Scores | One mean score vector per arm, and one residual covariance per arm |
| Dosing model | Per arm, the shared schedule of times and amounts |
| Visit model | Per arm, endpoint and time, the probability of an observation |
| Schema | Column order, column types as empty prototypes, the discrete endpoints’ level sets, the log-or-identity choice and assay limit per endpoint, and one value per arm for the strata and kept columns |
| Arms | The number of patients in each |

[`pmx_pca_report()`](https://iamstein.github.io/synpmx/reference/pmx_pca_report.md)
prints the inventory with the count of numbers in each group and, for
each, the smallest number of patients standing behind any one of them.

Two entries deserve naming because they are not averages. **The level
set of a discrete endpoint** is read from real values: a 0/1 endpoint
has no third value to emit, so a generated value is a source value. **A
categorical covariate’s levels** are likewise the ones the source held,
which is what the B5 row of
[`synpmx_scorecard()`](https://iamstein.github.io/synpmx/reference/synpmx_scorecard.md)
measures — a level held by one patient is a level that identifies them.

Subject identifiers are minted rather than derived. The only thing the
model keeps about the source’s identifiers is the largest one, so that a
synthetic ID cannot collide with a real one.

## Step 7: Place the Values

Each attended visit takes its drawn value for that endpoint, inverted
out of the log transform where one was applied, and snapped back onto
the source’s levels where the endpoint is discrete — the same correction
described under Step 11 of
[`vignette("avatar-algorithm")`](https://iamstein.github.io/synpmx/articles/avatar-algorithm.md).

Categorical covariates come out of the draw as one number per level, and
the level with the largest is the one written.

### The assay limit

The drawn value is the **latent** one: what the subject would have
measured with no assay limit. Where the source declared censoring, the
reported value, `CENS` and `LIMIT` are then reconstructed from it
together, so a value below the lower limit of quantification (LLOQ) is
reported *on* the boundary and flagged, as the source reports it.

This runs in both directions and both are needed.

Censored source values are replaced by a draw inside the censoring
region **before** the basis is fitted, back in Step 3. A grid cell where
most subjects sit exactly at the limit has almost no variance, so a
basis fitted on the reported values would learn the assay rather than
the patients, and the generated data would inherit a floor slightly
above the real one. This is the only random step in summarizing, which
is why
[`synpmx_pca_summarize()`](https://iamstein.github.io/synpmx/reference/synpmx_pca_summarize.md)
takes a seed of its own.

Reapplying the boundary at emit is what makes a low-dose arm come back
flat. Without it every value below the limit is emitted as itself, and
an arm that the study recorded as a straight line at the LLOQ comes back
as a spread of small numbers — visibly wrong on a log-scale
concentration plot, and wrong in a way that would carry into any BLQ
handling the dataset is used to develop.

## Step 8: Construct and Validate

New subject IDs, the strata and `keep` columns carried from the assigned
arm, every column back in the source’s order and type, and
[`validate_pmx()`](https://iamstein.github.io/synpmx/reference/validate_pmx.md)
over the result.

## Reading the Scorecard on This Output

[`synpmx_scorecard()`](https://iamstein.github.io/synpmx/reference/synpmx_scorecard.md)
was built to judge a synthetic dataset rather than a generator, and most
of it holds here. Fifteen of its eighteen rows read the source and
synthetic tables, compute the same way for any generator, and mean the
same thing.

Three do not, and they are the three that read the run’s own record
rather than the tables: **B1a**, **B1b** and **C2**. They report
`unavailable` on this output.

For B1a and B1b the guarantee is structural rather than unrecorded. No
individual’s visit set and no individual’s dose schedule exists in the
model to be copied, because the model holds a per-visit probability and
one shared schedule per arm. **B4a** and **B4b** ask the copy question
directly, on the finished tables, and are the rows that carry the claim.

C2 is the one to be careful about, because a blank row is at its most
misleading exactly where the loss is largest. C2 asks how many of the
study’s distinct dose-time schedules survived, and this algorithm
generates one schedule per arm. On a study whose recorded dose times are
actuals, that is a large loss and the card says nothing about it.
[`pmx_pca_dosing()`](https://iamstein.github.io/synpmx/reference/pmx_pca_dosing.md)
reports the same quantity as `distinct` and `share`, per arm, and is
what to read in its place.

## What This Algorithm Does Not Preserve

- **Within-visit detail.** The model holds a handful of components on a
  visit grid. Measurement noise from one visit to the next is not in
  them, so the generated profiles are smoother than the source.
- **Anything the components discarded.** Reconstructing the source
  subjects from their own retained scores shows how much that is, per
  endpoint. It is the diagnostic that says whether the retained count is
  large enough.
- **Interpretable pharmacology.** A score standard deviation is a
  variance decomposition, not a parameter. There is no clearance here to
  report a coefficient of variation on. A loading plotted against time
  can be described — flat and positive is overall magnitude, crossing
  zero separates early from late — but naming a mechanism for it is a
  claim the fit does not support.
- **Doses and times the study never ran**, under the default
  `dose_term`. `"log"` extrapolates; `"factor"` does not.
- **The variety of dose schedules.** Every patient in an arm receives
  the schedule that arm holds in common. The timing survives; the spread
  of schedules around it does not.
- **Within-subject dose escalation.** A subject’s total dose is one
  number and an arm’s schedule is one schedule, so a study escalating
  within a patient is represented by its arms rather than by the
  escalation.

## What Leaves the Source Data

The inventory is above, under “The Model Is the Only Thing Generation
Reads”, and
[`pmx_pca_report()`](https://iamstein.github.io/synpmx/reference/pmx_pca_report.md)
prints it for a given run. The column to read is the smallest number of
patients standing behind a quantity. A grid cell mean is backed by most
of the cohort and a per-arm quantity by one arm, while a rare covariate
level can be backed by a single patient.
