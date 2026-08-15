# The synpmx AVATAR Algorithm

This article is the full specification of the AVATAR synthetic data
generator applied to pharmacometrics modeling datasets. The original
method \[1\] was extended to population pharmacokinetic (PK) datasets in
\[2\], and here further extensions are provided so that it works well
with more realistic datasets. Use of the algorithm is demonstrated in:

- [`vignette("demo")`](https://iamstein.github.io/synpmx/articles/demo.md) -
  AVATAR demonstration on a single dataset
- [`vignette("public-data-examples")`](https://iamstein.github.io/synpmx/articles/public-data-examples.md) -
  AVATAR demonstration on a series of public datasets

AVATAR makes no formal privacy claim. It also does not necessarily fully
preserve the dose-exposure-response relationship. Use the synthetic
datasets to develop and debug analysis code. Do not use it to estimate
parameters.

Because modeling datasets are generally longitudinal, one set of
identifying features are the unique dose or observation times each
patient will have. The default setting of `synpmx_avatar` ensure that
the set of observation and dose times of the synthetic patient does not
match any one real patient data.

## Overview of Algorithm

Each synthetic patient is called an avatar. Each avatar is built from
one real **anchor** subject’s dosing event structure, then filled with
measurements blended from several similar real subjects, called
**donors**. The generator works on any number of endpoints, so one study
can carry multiple PK profiles alongside multiple pharmacodynamic (PD)
measures.

Every patient’s data can be thought of as containing two pieces of
information which are treated differently

- **event skeleton** is *when* the dose events and the DV observations
  happen.
- **values** are the DV numbers filling those observation slots, plus
  the baseline covariates.

The event skeleton is copied from one anchor patient with some masking
of the times, and the values are blended across several donors. The
`synpmx_avatar` function runs the following steps.

1.  **Declare the columns.** The `roles` variable say which column is
    the id, time, dose, measurement, etc. Only named columns survive.
2.  **Coarsen and align the times.** Snap every recorded time onto a
    shared visit grid (often based on nominal time), and align each
    subject to its own first dose.
3.  **Characterize each endpoint**: whether it is discrete, and whether
    it is log transformed. Both are applied later.
4.  **Build one numeric profile per source subject**, and reduce it to a
    handful of principal-component coordinates that measure similarity
    between patients.
5.  **Read each subject’s route of administration and dosing history**,
    the two facts that decide who may donate to whom.
6.  **Restrict the anchor pool and sample anchors.** A patient whose
    structure would give them away, an unusually long follow-up or a
    schedule nobody else shares, is not built upon, though their
    measurements still contribute as a donor. Each cohort (arm) keeps
    the same size (except for very small cohorts).
7.  **Copy the anchor’s event skeleton** and redraw the parts of it that
    identify the anchor: which visits carry an observation, and where
    dosing stopped.
8.  **Choose the donors**, nearest first within the anchor’s route.
9.  **Randomize and cap the donor weights**, so no one real patient is
    more than half of any avatar.
10. **Synthesis the baseline covariates by blending** from those donors.
11. **Synthesis the endpoint trajectories by blending**, then correct
    the generated value: snap a discrete endpoint back onto its levels,
    and reapply the censoring boundary.
12. **Recompute the dose** from the avatar’s own blended covariate,
    where dosing is proportional to a covariate (e.g. weight based
    dosing).
13. **Construct the full synthetic dataset**, with new synthetic subject
    IDs and every column back in the shape the source had, and validate
    it.

## The Masking Mechanisms

Across the steps above are six **masking mechanisms** of the individual
patient data.

- **M1. Coarsen the times.** Snap every time onto a shared visit grid
  (nominal time when available), then add deviations. This prevents
  identifying a patient from their TIME vector. *Step 2.*
- **M2. Restrict who may be an anchor.** Drop patients whose structure
  no amount of blending would hide, such as a patient followed twice as
  long as everyone else, and patients in a route arm too small to blend
  within. *Step 6.*
- **M3. Redraw which visits the avatar attended.** Reselect which DV
  values were collected at which time points, in case a patient has a
  unique set of missed visits. *Step 7.*
- **M4. Redraw where dosing stopped.** Where a patient’s dose schedule
  is one nobody else has, the avatar stops earlier, at the latest dosing
  time that several patients share. Dose *times* are never moved or
  invented, so this only ever truncates a real schedule earlier. *Step
  7.*
- **M5. Blend the values across several donors.** Covariates and DV,
  capped so no one donor dominates, plus noise. *Step 9.*
- **M6. Recompute the dose from the avatar’s own covariate.** Done when
  dosing is weight-based or body-surface-area based. Declared with
  `dose_covariate`, or inferred. *Step 12.*

## Step 1: Declare the Meaning of the Columns

The package never guesses critical pharmacometric (PMX) roles from
column names. The user must declare them with
[`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md):

``` r

roles <- pmx_roles(
  id = "ID",
  time = "TIME",
  dv = "DV",
  amt = "AMT",
  evid = "EVID",
  cmt = "CMT",
  dvid = "DVID",
  mdv = "MDV",
  rate = "RATE",
  nominal_time = "NTIME",  # protocol-scheduled time, if the study records it
  cens = "CENS",           # below-limit indicator, if the study has BLQ data
  limit = "LIMIT",         # the other end of an interval, when used
  covariates = c("WT", "AGE", "SEX"),
  dose_covariate = "WT",   # dose is mg/kg; say so rather than let it be inferred
  strata = "ARM",          # assigned arm; groups three mechanisms, see below
  keep = "STUDYID"         # carried verbatim, see the allowlist below
)
```

`id`, `time`, `dv`, and `evid` are required and
[`?pmx_roles`](https://iamstein.github.io/synpmx/reference/pmx_roles.md)
lists the rest. A column generally holds one role, with two exceptions.
NONMEM’s `CMT` routinely does both jobs, the dosing compartment on event
rows and the endpoint key on observation rows, so `cmt` and `dvid` may
name the same column. And `dose_covariate` must *also* be named in
`covariates`, because the avatar’s amount is rebuilt from its own
blended value of that column, so the column has to be one that gets
blended.

### The Role Declaration Specifies Which Columns Appear in the Synthetic Dataset

[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
keeps only columns a role names, and **drops every other column**,
reporting which ones are dropped. This is because AVATAR builds each
synthetic subject from one real anchor subject, and any column it
retains that it does not actively synthesize is copied from that anchor
verbatim. A secondary identifier such as `USUBJID`, a site number, or a
randomization date, left undeclared, would otherwise follow a real
subject’s real value straight into the synthetic data. Dropping by
default means a column you forget fails safe, toward removal, rather
than leaking.

That leaves five ways a column is treated:

| Declared as | What AVATAR does to it | Use for |
|----|----|----|
| `dv` (with `cens`/`limit`) | Blended across donors into a new trajectory | the measurement |
| `covariates` | Blended or resampled across donors into a new value | baselines you want *synthesized*: weight, age |
| `strata` | Copied verbatim from the anchor | assigned strata: treatment arm, dose group, cohort |
| `keep` | Copied verbatim from the anchor, and otherwise inert | assigned values you want kept faithful to that subject’s dosing, such as a study identifier |
| *(undeclared)* | Dropped | anything you do not need |

**`strata` carries a column the same way `keep` does, and then uses
it.** Both copy the anchor’s value verbatim, so the choice between them
is not about what appears in the output. It is about whether that column
should also *partition the cohort*, as will de desrcibed below. But,
`strata` is deliberately **not** a blending barrier. Donors are still
borrowed across strata. Declare a treatment arm as `strata` rather than
`keep` whenever the arms differ in dose rule or visit schedule, and use
`keep` for a label that carries no such structure.

Redundant endpoint labels have their own handling. A dataset that
carries both a numeric `YTYPE` and a character `NAME` for the same
endpoint declares **both** as `dvid`: `dvid = c("YTYPE", "NAME")`. The
first is the grouping key. Validation checks that the rest are a
consistent 1:1 mapping with it, catching a source where the two labels
disagree, which `keep` would carry through silently, and AVATAR carries
all of them into the output, aligned.

Immediately before generation,
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
calls `validate_pmx(..., strict = TRUE)`. Among other checks, validation
requires:

- nonmissing IDs and EVID values;
- finite numeric TIME values and nondecreasing TIME within subject;
- numeric DV and finite nonmissing observations;
- an endpoint value on every eligible observation when DVID is declared;
- multiple `dvid` columns to be a consistent 1:1 mapping;
- every declared baseline covariate to be constant within subject; and
- when `cens` is declared, each CENS flag to agree with its DV. A
  left-censored value at the limit cannot exceed an ordinary measurement
  of the same endpoint.

## Step 2: Coarsen the Times and Align to the First Event

### M1. Coarsening the Times, and Putting the Scatter Back

At the default `coarsen_time = TRUE`, source times are collapsed onto a
shared visit grid (`nominal_time` when provided, otherwise a grid is
computed automatically). This **snapping**, replaces each recorded time
with its visit’s grid time. The deviation removed is not discarded. It
goes into one pool shared by the whole cohort. **Restoring** occurs
after each avatar is built, when a deviation from that pool is drawn and
added back to the time.

### First-Event-Relative Time

Neighbor profiles and donor trajectories are constructed on an internal
aligned time. For subject $`i`$, the preferred origin is the first
nonzero-EVID row having a positive AMT:

``` math
t_{0i} = \min\{t_{ij}: \mathrm{EVID}_{ij}\ne0,\ \mathrm{AMT}_{ij}>0\}.
```

If no positive-AMT event exists, the first nonzero-EVID row is used. If
the subject has no event row, TIME is left unchanged. The internal
aligned time is

``` math
t^*_{ij}=t_{ij}-t_{0i}.
```

This removes a subject-specific calendar or study-time offset. It does
**not** compute time after each dose.

A declared `tad` column is recomputed in the output, as time since the
avatar’s own most recent dose.

## Step 3: Characterize Each Endpoint

Blending is a weighted mean of several donors, so before any blending
happens the endpoint must be categorized in two ways.

- What type is it: continuous, binary, ordinal, or integer.  
- Should the endpoint be log transformed (and offset)

Both categories are recorded on the output and used later.

### Endpoint types

- **`continuous`** when any recorded value is not a whole number.
  Nothing is constrained.
- **`binary`** when every value is a whole number and there are two
  distinct levels.
- **`ordinal`** when every value is a whole number and there are at most
  twelve distinct levels.
- **`integer`** when every value is a whole number and there are more
  than twelve.

The endpoint type is inferred from the data, but
`pmx_roles(endpoint_types = )` overrides the inference.

### Log-transformation

Transformations are selected per DVID endpoint, or once when the DVID
column is absent. An endpoint is **positive-like** when it has at least
one positive value, its median is positive, and at most 1% of values are
negative. The 1% tolerance lets a concentration endpoint with a couple
of small negative assay readings still count.

PK and PD data can have exact zeros and $`\log 0=-\infty`$ would make
blending impossible. So a small constant $`c_e`$ set to be half the
smallest positive value of that endpoint, giving the transformation of:

``` math
g_e(y)=\log\{\max(y,0)+c_e\},
\qquad
g_e^{-1}(z)=\max\{\exp(z)-c_e,0\}.
```
Transformations and offsets are stored in
`attr(synthetic, "pmx_settings")$endpoint_transforms`.

## Step 4: Build One Numeric Profile per Source Subject for Distance Calculation

Source subject $`i`$’s profile concatenates their baseline covariates
$`\mathbf{c}_i`$ and their longitudinal observations
$`\mathbf{y}_{ie}(G_e)`$ on each of the $`E`$ endpoints, interpolated
onto a common grid $`G_e`$ of at most fifteen times, built per endpoint
from the observation times pooled across subjects. This common grid is
not the nominal times, but rather a grid that is later used for
calculating the principal component analysis and distance between
subjects.

``` math
\mathbf{x}_i=
\left[
  \mathbf{c}_i,\
  \mathbf{y}_{i1}(G_1),\ldots,\
  \mathbf{y}_{iE}(G_E)
\right].
```

This $`\mathbf{x_i}`$ profile is used for distance calculations. It does
not replace the source data and is not returned as a synthetic patient.

### The Common Time Grid

For endpoint $`e`$, the package pools finite aligned observation times
and makes a common grid $`G_e`$. If there are at most 15 unique times,
all are retained. If there are more, 15 type-8 empirical quantiles are
used and duplicate quantiles are removed.

Each subject’s transformed trajectory is linearly interpolated onto
$`G_e`$. There is no extrapolation during profile construction, so grid
locations outside a subject’s observed window become missing profile
features. Tied times are averaged by `stats::approx(..., ties = mean)`.

For PK data a 15-point summary may compresses it. The grid points are
quantiles of the observed times, so they land where sampling is dense,
early and around absorption and peak, rather than spreading evenly
across a long tail. The grid is also built per endpoint, so PK and PD
each get their own grid times rather than sharing them. For rich,
multiple-dose data with several peaks and troughs, 15 points may not be
enough to well characterize the shape and when that happens the cost is
a worse choice of donors. Again, this 15 point grid is only for distance
calculating between patients and not for the final synthetic data
profiles.

### Standardization across Multiple Endpoints

Every endpoint’s columns go into one profile and one PCA.
Standardization of $`\mathbf{x_i}`$ equalizes…

The standardization below is what makes that possible, but it equalizes
**units**, not **endpoints**. Every grid point becomes one unit-variance
column, so an endpoint sampled at 15 times contributes 15 columns while
one sampled at 3 times contributes 3, and a Euclidean distance summing
over columns lets the first outvote the second. Donor selection is
therefore driven mostly by whichever endpoint was sampled most densely,
and baseline covariates are a small minority of the distance whenever
trajectories are rich. There is currently no way to reweight endpoints.

The profiles are stacked into a matrix with one row per source subject
$`i`$ and one column per **feature** $`l`$: one column for each
covariate, and one for each grid time of each endpoint, so a study with
two covariates and two endpoints on fifteen grid points each has 32
columns.

Missing entries are replaced by that column’s across-subject median
$`m_l`$:

``` math
x^{\mathrm{imp}}_{il}=
\begin{cases}
x_{il}, & x_{il}\text{ is finite},\\
m_l, & \text{otherwise}.
\end{cases}
```

Each column is then centred on its own across-subject mean and divided
by its own across-subject standard deviation $`s_l`$:

``` math
z_{il}=\frac{x^{\mathrm{imp}}_{il}-\bar{x}^{\mathrm{imp}}_l}{s_l}.
```

Every column therefore arrives at the PCA with mean 0 and variance 1.  
A column whose standard deviation is not greater than
$`\sqrt{\epsilon_{\mathrm{mach}}}`$ is treated as constant and is
dropped.

### Reducing the Profile to a Few Coordinates by PCA

The standardized matrix $`z_{il}`$ has $`n`$ rows, one per source
subject, and $`p`$ columns, one per surviving feature. Its columns are
heavily redundant: adjacent grid times on one trajectory carry nearly
the same information, so $`p`$ columns hold far fewer than $`p`$
independent facts about a patient, and a distance computed straight on
those columns counts the repeated ones several times over.

**PCA refresher.** Principal component analysis (PCA) rotates the $`p`$
columns into a new set of axes, the **principal components**. Each
component is one fixed weighted combination of the original columns, and
the weights are chosen so that the first component runs along the
direction in which subjects are most spread out, the second along the
most spread remaining while being uncorrelated with the first, and so
on. Those directions are the **eigenvectors** of the correlation matrix
of the columns, and the variance of the subjects along each one, written
$`\lambda_h`$, is the matching **eigenvalue**.

Because every column was standardized to variance 1 in the previous
section, the total variance is exactly $`p`$, and $`\lambda_h/p`$ reads
directly as the share of the whole profile’s spread that component $`h`$
accounts for. The retained dimension $`H`$ is the smallest number of
components whose shares add up to at least $`v`$, which is
`pca_variance`, default 0.90:

``` math
\frac{\sum_{h=1}^{H}\lambda_h}
     {\sum_{h=1}^{r_{\max}}\lambda_h}
\ge v .
```

Each subject is then projected into the PCA space of H dimensions, and
it is this space that is used to compute the distance between subjects
in Step 8:

``` math
\boldsymbol{\xi}_i=(\xi_{i1},\ldots,\xi_{iH}),
\qquad
\xi_{ih}=\text{subject } i\text{'s score on principal component } h .
```

**Rank-safe** is the bound on how many components can exist at all:

``` math
r_{\max}=\min(n-1,p).
```

There can be no more directions than columns, hence $`p`$; and centering
each column on its own mean spends one degree of freedom, hence $`n-1`$.
A twelve-subject study therefore has at most eleven components however
many grid points it carries, and a component beyond $`r_{\max}`$ would
be fitting numerical noise.

## Step 5: Read Each Subject’s Route and Dosing History

A subject given a single oral dose should not donate a trajectory to an
infusion template. So two facts are read off each source subject here,
once for the whole cohort, and used twice later: Step 6 needs the route
to size the arms it may anchor in, and Step 8 needs both to pick each
avatar’s donors.

**Route of administration is an absolute barrier.** A bolus, an infusion
and an oral dose produce categorically different concentration shapes,
so blending across them yields a trajectory no protocol could have
produced. Donors are never drawn from outside the anchor’s route. Each
subject’s route is the set of

``` math
\bigl(\mathrm{EVID},\ \mathrm{CMT},\ [\mathrm{RATE}\neq 0]\bigr)
```

triples on its dosing rows, so the compartment a dose enters and whether
it is delivered over time both bind, while dose *count* does not.
Because it is the whole set, **a patient given more than one route is
their own route class**: an infusion followed later by an oral dose
matches neither the oral-only nor the infusion-only patients, only
others with that same combination.

## Step 6: Restrict the Anchor Pool and Sample Anchors

### M2. Which Patients May Be Used as an Anchor

There are two reasons subjects may be excluded from being an anchor.

1.  A subject whose follow-up length or dose count exceeds twice the
    cohort’s 90th percentile is never used as an anchor (`screen`), so
    no avatar inherits a conspicuous skeleton.
2.  A subject in a route arm holding fewer than `k + 1` subjects has no
    legal donor set, so by default it is dropped from the pool too
    (`on_donor_shortfall`), where `k` is the number of subjects to be
    chosen for blending.

When an arm holds nobody who can be masked, individualised dosing in one
cohort of a trial for instance, its avatars keep the anchor they have,
the run raises an alert naming that arm, and rows B1a and B1b of
[`vignette("scorecard-synthetic-data-checks")`](https://iamstein.github.io/synpmx/articles/scorecard-synthetic-data-checks.md)
fail against it.
[`unmaskable_strata()`](https://iamstein.github.io/synpmx/reference/unmaskable_strata.md)
can assess whether any such strata exist.

### Sampling Anchors, and Keeping Each Stratum’s Share

Anchors are sampled with replacement from the surviving pool. As the
default `preserve_strata_balance = TRUE`, and so the target proportions
of the strata in the synthetic dataset match the source dataset.
However, **Strata holding fewer than three source patients are
deliberately not balanced**, because reproducing such a stratum’s size
exactly would say more about those few patients.

Setting `preserve_strata_balance = FALSE` allows strata sizes to vary.

## Step 7: Copy the Anchor’s Event Skeleton, and Redraw the Visit Schedule

The anchor’s **event skeleton** is a row-wise combination of TIME, EVID,
AMT, RATE, CMT, DVID.

### M3. Redrawing Which Visits the Avatar Attended

After coarsening, every subject sits on the same grid of visit times.
What still differs between them is which of those visits carry an
observation, and that combination can identify a patient on its own, as
shown in the example below.

    grid            wk1  wk2  wk3  wk4  wk5  wk6
    patient A        x    x    x    x    x    x     attended everything
    patient B        x    x    x    x                stopped after week 4
    patient C        x    x         x    x    x      missed week 3

Each row of x’s is a **visit set**. If patient C is the only one who
ever missed week 3, then “attended everything except week 3” points at
patient C as surely as an unusual sampling time would.

So the avatar is given a different visit set from its anchor’s. The
doses and the grid times are untouched; only which cells carry an
observation changes.

Two descriptions of a visit set are used to decide what the avatar may
be given. Its **arrangement** is where the misses sit, one of four
labels:

| Arrangement | Meaning |
|----|----|
| `complete` | nothing missed |
| `trailing` | every miss falls in the final visits, so the patient stopped: a discontinuation, or follow-up that has not reached those visits yet |
| `block` | the misses are consecutive but do not run to the end, so the patient was interrupted and came back |
| `scattered` | anything else |

With $`n`$ visit cells of which $`m`$ are missed, the arrangement is
`trailing` when every missed index exceeds $`n-m`$, `block` when the
missed indices are consecutive but not trailing, and `scattered`
otherwise. One missed visit that is not the last is a block of one. In
the picture above, patient B is `trailing` and patient C is `block`.

Its **shape** is that arrangement paired with the number of misses,
written `trailing|2` or `block|1`. The count is the only difference
between the two descriptions: patient B above is `trailing|2`, and a
patient who stopped one week later, missing week 6 alone, is
`trailing|1`. Same arrangement, different shape. Neither says *which*
visits were missed.

The two are separate because they are used at different points below. A
shape is specific enough that a set built to it resembles the patients
it came from, so it is tried first. An arrangement is coarser, and
groups patients a shape cannot, which is what makes it a usable last
resort.

With those two, the avatar’s visit set is chosen by trying three things
in order and using the first that works:

1.  **Reuse a set.** Take a visit set that at least `min_pattern_share`
    patients hold. The avatar then has a schedule several real patients
    had.
2.  **Build one with a shared shape.** Where no set is held that widely,
    take a shape that at least `min_pattern_share` patients share and
    build a fresh set with it.
3.  **Build one with a shared arrangement.** Where no shape is shared
    either, only the arrangement is used, and the number of misses is
    set to the median of the counts among the patients holding that
    arrangement. A study with staggered follow-up needs this: patients
    who each stopped at a different visit share no set and no shape,
    because every count has a single holder, while `trailing` is a label
    all of them carry.

Pairing the arrangement with the count is what makes step 2 work. Two
patients who each missed a single early visit hold different sets, so
each is a group of one and neither set can be reused. As `block|1` they
are one group of two, and a fresh block of one can be built.

When a fresh set is built, the arrangements that fit the wanted count
are enumerated and one that no single patient holds is taken. Trailing
misses have only one arrangement per count, so where that one is already
somebody’s, the count moves outward by a visit or two until a free
arrangement is found. Misses are only ever placed on endpoints that some
patient actually misses, since putting one on an endpoint everybody has
at every visit would invent missingness the study does not have.

What survives is how much missingness the study had and of what kind.
What is lost is which weeks any one patient missed.

### M4. Redrawing Where Dosing Stopped

A patient’s dose schedule is their list of dose times, and it identifies
them whenever nobody else has the same list.

Dose times are never moved or invented, because a drawn regimen could be
one no protocol permits. The one edit allowed is to shorten: dropping
the doses at the end of a schedule leaves a course some patient could
have been given, namely one who stopped early.

So when a patient’s schedule is one nobody else has, the avatar’s is
shortened. Shortening to $`d`$ doses means the avatar keeps that
patient’s own first $`d`$ dose times and stops there. A patient dosed at
0, 12, 24 and 36 hours shortened to 3 doses is dosed at 0, 12 and 24,
and their 36-hour dose row is dropped.

The search tries one dose fewer than the patient received, then one
fewer again, and stops at the first length that is safe. A length is
safe on either of two grounds:

- **Matched.** Those first $`d`$ times are exactly the whole schedule of
  at least `min_pattern_share` patients, so the avatar finishes where
  they finished.
- **Free.** No patient’s whole schedule is those $`d`$ times, and at
  least `min_pattern_share` patients were dosed at all of them and then
  carried on further. Nobody stopped there, so the stopping point
  identifies nobody, and every dose leading up to it is one that many
  patients received.

Worked through, with a floor of two:

    dose times (h)       0   12   24   36   48
    patient A            x    x    x                3 doses
    patient B            x    x    x                3 doses
    patient C            x    x    x                3 doses
    patient D            x    x    x    x    x      5 doses, shared by nobody

An avatar anchored on D would carry D’s schedule alone. Four doses is
tried first: {0, 12, 24, 36} is nobody’s whole schedule, but D is also
the only patient dosed at all four of those times, so it is neither
matched nor free. Three doses is {0, 12, 24}, which is exactly the whole
schedule of A, B and C, so it is matched and the avatar stops at 24
hours.

If it happens that a set of non-informative doses cannot be built for a
particluar cohort, then that cohort cannot be represented and this will
be shown in the reporting scorecard.

Two study shapes are served. In a protocol cohort with staggered
discontinuation, everyone is dosed on the same schedule, so any length
nobody stopped at is free and a target is available almost anywhere. The
other is routine clinical care, where every patient begins on the same
interval and diverges as their own clinician decides, so whole schedules
rarely repeat while the first several doses are common to many patients.

## Step 8: Choose the Donors

### M5. Blending

Blending subject data together across donors spans Steps 8 to 11.

### The Distance

For anchor $`a`$ and candidate donor $`r`$, the distance is the ordinary
straight-line (Euclidean) distance between their subject coordinates
from Step 4:

``` math
d_{ar}=\left\|\boldsymbol{\xi}_a-
                   \boldsymbol{\xi}_r\right\|_2
=\left(\sum_{h=1}^{H}(\xi_{ah}-\xi_{rh})^2\right)^{1/2},
```

where $`h`$ runs over the $`H`$ retained principal components.

### The Selection Algorithm

Every other subject in the study that shares the anchor’s route is a
candidate, whatever arm they were in. Strata never restrict donors, so
an avatar allocated to one arm can be blended from patients in another,
and only route excludes anybody. To choose $`k`$ donors for anchor
$`a`$, using the route and the dosing history of Step 5:

1.  **Same-history stage.** Among the other subjects sharing the
    anchor’s route that also share its dosing history, take the $`k`$
    nearest by $`d_{ar}`$, breaking ties by subject index. If there are
    at least $`k`$ of them, stop.
2.  **Fallback stage.** Otherwise fill the remainder from the nearest
    subjects on that route whatever their dosing history, and record a
    warning that donors were borrowed across dose groups.
3.  **Shortfall.** If the anchor’s route holds fewer than $`k`$ other
    subjects at all, the floor is unreachable, because the only
    remaining candidates are on another route and route is never
    crossed.

Selection is deterministic under a fixed seed. Note that `AMT` is not a
profile feature, so the choice of donors never directyl uses dose. It
compares baseline covariates and the shape and level of the DV
trajectory. Dose does influence the DV trajectory, but only
*indirectly*: a subject given a much larger dose than the anchor has
higher concentrations, and those concentrations are profile features, so
that potential donor lands further away. The ranking therefore tends to
prefer similar-dose donors without being told to.

### When the Donor Floor Cannot Be Reached

A route arm holding fewer than $`k+1`$ subjects can never supply a legal
donor set. There is no good answer in that situation, only a choice
between two bad ones, and which is worse depends on whether the sparse
arm matters more than the disclosure risk of reproducing it. That is the
caller’s judgement, so `on_donor_shortfall` makes it explicit:

| `on_donor_shortfall` | Behavior |
|----|----|
| `"drop"` *(default)* | Omit those subjects from the anchor pool. No avatar is built on them, and the synthetic cohort does not represent that arm. |
| `"noise"` | Keep them, blending however many same-route donors exist, possibly none, and relying on `subject_noise_sd` and `residual_noise_sd` for the rest. |
| `"error"` | Refuse to generate, naming the arm and both alternatives. |

`"noise"` is available but not recommended because a subject blended
from one donor, or from none, is a noised near-copy of a real patient.
If you use it because the sparse arm genuinely matters for the synthetic
data, screen the result and treat those subjects as individually
identifying. If *every* arm is below the floor, `"drop"` would leave
nothing to generate, so generation proceeds as if `"noise"` and says so.

## Step 9: Randomize and Cap the Donor Weights

Suppose $`K`$ donors have been selected. Each one draws an independent
$`E_r\sim\mathrm{Exponential}(1)`$, a random positive number from the
exponential distribution with rate 1. Such a draw has mean 1 and median
0.69, falls below 2 about 86% of the time, and occasionally comes out
several times larger, allowing a donor take a large share at times
regardless of distance.

Let $`R_r`$ be a random permutation of $`1,\ldots,K`$. It is a
randomized rank that stops the nearest neighbour from being
systematically the largest contributor. The raw donor weight is

``` math
q_r=\frac{E_r}{\max(d_{ar},\varepsilon)}2^{-R_r},
```

where $`d_{ar}`$ is the profile distance between the anchor and donor
$`r`$ from Step 8. Dividing by it is what makes a nearer donor count for
more.

The $`2^{-R_r}`$ term is **rank attenuation**. Since $`R_r`$ is a
permutation of $`1,\ldots,K`$, the $`K`$ donors are handed the
multipliers $`\tfrac12`$, $`\tfrac14`$, $`\tfrac18`$ and so on in random
order, so one donor counts for twice the next, four times the one after,
and the last counts for almost nothing.

It is there to stop an avatar being a cohort average. The $`k`$ nearest
neighbours sit at comparable distances from the anchor, so inverse
distance alone would give them near-equal weights, and an even average
of five patients is smoother and flatter than any one of them: every
avatar would land near the middle of the cohort and the synthetic
between-subject variability would collapse. Attenuating by rank leaves
each avatar mostly resembling one or two real patients, which is what
keeps that variability. It is inherited from the published AVATAR
method, and it is also the term that pushes hardest against privacy,
which is what the cap below exists to limit.

A donor whose profile is identical to the anchor’s has $`d_{ar}=0`$,
which would divide by zero, so the distance is floored at
$`\varepsilon`$:

``` math
\varepsilon=
\begin{cases}
\max\{10^{-8},10^{-6}\operatorname{median}(d_{ar}:d_{ar}>0)\},
  & \text{if a positive distance exists},\\
10^{-8}, & \text{otherwise}.
\end{cases}
```

The floor is scaled to the cohort’s own distances rather than fixed, so
it stays negligible against a real distance whatever units the profile
happens to be in.

The raw weights are then normalized to sum to 1, $`w_r=q_r/\sum_s q_s`$,
which is what makes each one a share of the avatar. If they are
nonfinite or have a nonpositive total, the implementation replaces them
with equal raw values before normalizing.

### Capping Every Donor

Each donor’s weight is its **share of the avatar**: five weights that
sum to 1. The floor $`k`$ says how many real patients go into the blend,
but it says nothing about how much of any one of them comes out. Five
donors weighted $`(0.95, 0.02,
0.02, 0.005, 0.005)`$ satisfy a floor of five while being, to any
practical purpose, a copy of one person. So the cap, rather than $`k`$,
is what controls how closely an avatar can resemble a single real
patient.

The raw formula above produces very uneven weights on purpose, which is
the rank attenuation doing its job. `max_donor_weight` (written $`c`$
below, default $`0.50`$) is the ceiling on any one share.

Enforcing it is less obvious than it sounds, because **weights have to
keep summing to 1**. Cutting the leader down does not remove that weight
from the blend; it has to be given to the others, in proportion to what
they already hold:

``` math
w_r \leftarrow \Bigl(1-\textstyle\sum_{s\in P}w_s\Bigr)
              \frac{w_r}{\sum_{s\notin P}w_s},
\qquad r\notin P,
```

where $`P`$ is the set of donors already pinned at $`c`$. And that
redistribution can itself push a donor over the ceiling. Starting from
$`(0.50, 0.25, 0.15, 0.06,
0.04)`$ with a tight $`c = 0.30`$: cutting the leader to $`0.30`$ frees
$`0.20`$, and the second donor’s share of it takes them from $`0.25`$ to
$`0.35`$, over the ceiling.

So the rule repeats, pinning everything over $`c`$, redistributing, and
looking again, which is the standard construction called
**water-filling**. Pinned donors are never topped up, so each pass pins
at least one more and it terminates. The tight example settles at
$`(0.30, 0.30, 0.24, 0.096, 0.064)`$.

## Step 10: Synthesis the Baseline Covariates by Blending

The covariates are blended, using the donors chosen in Step 8 and the
weights calculated in Step 9. A synthetic patient’s weight is a weighted
mean of its donors’ weights. Computing the mean requires:

1.  Choosing the scale to average on
2.  Choosing how much to perturb the estimate
3.  Choosing what to do with covariates that cannot be averaged (binary
    or categorical)

For a numeric covariate, let $`c_r`$ be donor $`r`$’s first nonmissing
value. When all available values are positive and

``` math
\frac{\max_r (c_r)}{\operatorname{median}_r(c_r)}>3,
```

the blend is formed on the log scale; otherwise it is formed on the
original scale. The working scale $`h`$ is $`h(c)=\log c`$ on the log
scale and $`h(c)=c`$ on the original one. The donor mean ($`\mu_c`$) is
then:

``` math
\mu_c=\sum_r w_r^* h(c_r).
```

The star marks a second normalization. In Step 9 the weights were
normalized across all the chosen donors, so $`w_r`$ already sums to 1.
But if one donor has no value for this covariate, the donors that do
have one carry weights summing to less than 1, and their mean would be
pulled toward zero. So the missing donor is dropped and the rest are
rescaled to sum to 1 among themselves. Where every donor has a value,
$`w_r^*=w_r`$.

Averaging several patients gives a value nearer the middle of the cohort
than any of them was, so scatter is added back in proportion to how
spread out the donors themselves were. Let $`s_c`$ be the ordinary
sample standard deviation of the available working-scale donor values,
with fallback $`\max(0.05|\mu_c|,0.01)`$ when it is unavailable or zero.
The generated working-scale covariate is

``` math
c^*=\mu_c+\eta_c,
\qquad
\eta_c\sim\mathrm{N}(0,(\sigma_{\mathrm{subj}}s_c)^2),
```

where `subject_noise_sd` is $`\sigma_{\mathrm{subj}}`$ (default 0.15).
The result is then inverse transformed by $`h^{-1}`$, so a log-scale
blend is exponentiated, and any all-positive covariates are floored
above zero. This final value is repeated on every row of that synthetic
subject.

For factor, character, and logical covariates, one available donor
category is sampled using the locally normalized weights. Categories are
not averaged, so a synthetic patient’s category is some real patient’s
actual category, copied.

### Covariates and Endpoint Are Blended Independently

Covariates and endpoint trajectories are generated in separate passes.
In most cases that is fine, but some datasets carry a baseline covariate
that is the same quantity as a longitudinal endpoint. For example: `B0`
baseline B-cell count beside a B-cell kinetic endpoint. AVATAR does not
know they are linked, so a synthetic subject’s baseline covariate need
not equal the baseline of its own generated trajectory. The two are
individually plausible but not mutually consistent.

This is a known limitation. It is usually harmless for
workflow-development use, where each column need only be realistic on
its own. If your analysis relies on a covariate agreeing with an
endpoint’s baseline, reconcile them after generation, or do not declare
the redundant covariate and derive the baseline covariate value from the
synthetic data trajectory instead.

## Step 11: Synthesize the Endpoint Trajectories by Blending

The anchor supplies the time of the observations, including its
missing-DV pattern, as redrawn by M3. For endpoint $`e`$, donor
trajectories are transformed with $`g_e`$ (from Step 3) and interpolated
to the anchor observation times. Let $`z_{rj}`$ be donor $`r`$’s
transformed value at target position $`j`$, where it has one. The blend
is

``` math
\bar z_j=\sum_r w_{rj}^{*}z_{rj}.
```

The weights carry a $`j`$ because the renormalization of Step 10 is done
separately at every target position: which donors have an observation
covering that time changes along the trajectory, so the blend runs over
whichever of them do and their weights are rescaled to sum to 1 among
themselves.

Where no donor covers a time at all, the row falls back to the cohort’s
median trajectory. The endpoint noise scale $`s_e`$ is 1 on the
offset-log scale. On the identity scale it is the source endpoint
standard deviation, with fallback
$`\max(0.1|\operatorname{median}(Y_e)|,0.01)`$ when necessary. A single
subject-level shift is shared across all positions of the endpoint:

``` math
b_e\sim\mathrm{N}(0,(\sigma_{\mathrm{subj}}s_e)^2).
```

Within-endpoint residual perturbations follow a stationary-scale
first-order autoregressive, AR(1), process in **observation order**:

``` math
\begin{aligned}
\epsilon_{e1} &\sim
  \mathrm{N}(0,(\sigma_{\mathrm{res}}s_e)^2),\\
\epsilon_{ej} &= \phi\epsilon_{e,j-1}+\nu_{ej},\\
\nu_{ej} &\sim
  \mathrm{N}\left(0,
    (\sigma_{\mathrm{res}}s_e)^2(1-\phi^2)\right).
\end{aligned}
```

Defaults are `residual_noise_sd = 0.05` and `residual_phi = 0.6`. The
generated DV is

``` math
Y^*_{ej}=g_e^{-1}(\bar z_j+b_e+\epsilon_{ej}).
```

Endpoints are processed separately, so PK and PD values cannot be
blended with one another. They nevertheless use the same donor subjects
and initial subject weights, preserving a limited form of cross-endpoint
subject coherence.

### Snapping a Discrete Endpoint Back onto Its Levels

$`Y^*_{ej}`$ above is a new number, which is what should happen for a
continuous endpoint, but not for a categorical endpoint. So for an
endpoint the value type of Step 3 made `binary` or `ordinal`, the
generated value is replaced by the **nearest observed level**, deciding
ties at the midpoints between adjacent levels. For an `integer` endpoint
it is rounded, and floored at zero where the source held no negative
value. A `continuous` endpoint is left alone.

### Below-the-Limit-of-Quantification (BLOQ) Data and a Monolix-Style CENS Column

If `cens` is declared in
[`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md),
AVATAR reconstructs censoring. The mechanism has three parts, and the
last of them happens here, alongside the snap above.

**The boundary is read from the source.** Under the Monolix convention a
censored row reports the boundary in `DV`, with `LIMIT` carrying the
other end of an interval when there is one. Left, right, and interval
censoring are each recognised. If a study used more than one assay
limit, the most conservative boundary is taken rather than inventing a
per-row rule.

**Censored source values are replaced before blending.** This happens to
the source, before profiles are built, so these values are *donor
inputs* rather than output. Every censored row in the source reports the
same number, the limit itself, so a donor’s censored stretch is a flat
line at the limit. Averaging those gives an avatar a latent value
sitting exactly on the boundary, and whether it comes out censored is
then decided by noise: roughly half the avatars are emitted just above
the limit, a value the assay could not have reported, and the synthetic
censoring rate is far below the study’s.

Each censored row is instead replaced by a draw inside its censoring
region, before the trajectory is transformed and blended. For a
left-censored row that is uniform between zero and the limit, and for an
interval-censored row uniform between the two ends.

**`DV`, `CENS`, and `LIMIT` are reconstructed together.** The
back-transformed blend is treated as a *latent* value: what the subject
would have measured with no assay limit. The boundary is then applied to
that latent mearurement, and all three columns are written from the
result. A generated value below the limit produces `CENS = 1` with `DV`
exactly at the limit.

## Step 12: Recompute the Dose from the Blended Covariate

### M6. Recomputing the Dose

When dosing is proportional to a baseline covariate (e.g. weight), the
multiplier is a protocol property the stratum shares and the covariate
is individual. Keeping the multiplier and recomputing the amount from
the avatar’s own covariate stops the dose being a verbatim real value.

**Declared or inferred.** Naming the covariate in
`pmx_roles(dose_covariate = )` states the relationship outright; leaving
it `NULL` makes the run infer one. Inference is deliberately
conservative: it requires the dose-to-covariate ratio to collapse onto a
handful of levels within 2%, because rewriting amounts on a study that
is *not* proportional would be worse than leaving them alone.

Two ordinary things defeat that test. The first is **dispensing in whole
units**. A protocol may say 2 mg/kg, but if the drug comes in 25 mg
vials the amount actually given is rounded to a multiple of 25, and
since every patient weighs something different that rounding is a
different fraction of each patient’s dose. A cohort with weights ranging
from 54 to 80 kg dosed that way records ratios from 1.84 to 2.21 mg/kg,
a spread of about 19%, so the ratios never collapse to within 2% and
inference correctly declines. The second is **patients escalating by
factors that are not identical**, which spreads the ratios the same way.
Either leaves a genuinely weight-based study with its amounts untouched.
Declaring the covariate skips inference and holds each dose row’s own
ratio, rather than snapping to shared levels, so intra-patient
escalation survives exactly.

## Step 13: Construct the Full Synthetic Dataset

The synthetic subjects are bound into one dataset, the original column
order is restored, and each result column is converted toward its source
class:

- factors recover source levels and ordering
- integer columns are rounded and converted to integer
- double, logical, character, Date, and POSIXct columns are restored

New IDs use labels such as `syn_001`, with a prefix added repeatedly if
needed to avoid collision.

If MDV is declared and the entire source obeys the standard relationship

``` math
[\mathrm{MDV}=0]=[\mathrm{EVID}=0\ \text{and DV is present}],
```

MDV is re-derived from that rule in the synthetic data and restored to
its original class.

## Reproducibility

Every random draw happens inside a local seed scope, so identical inputs
and arguments return an identical dataset. The caller’s random-number
state is restored on exit, whether or not one existed before the call,
so generating a dataset does not consume the random stream of the
analysis around it. `seed` must be an integer between zero and
`.Machine$integer.max`.

## Evaluating a Synthetic Dataset with the Scorecard

[`synpmx_scorecard()`](https://iamstein.github.io/synpmx/reference/synpmx_scorecard.md)
is the function that is used to assess a synthetic dataset. It is
explained in
[`vignette("scorecard-synthetic-data-checks")`](https://iamstein.github.io/synpmx/articles/scorecard-synthetic-data-checks.md),
and demonstrated in:

- [`vignette("demo")`](https://iamstein.github.io/synpmx/articles/demo.md) -
  AVATAR demonstration on a single dataset
- [`vignette("public-data-examples")`](https://iamstein.github.io/synpmx/articles/public-data-examples.md) -
  AVATAR demonstration on a series of public datasets

## Conclusions

This document describes an extension of the AVATAR algorithm from \[1,
2\]. This algorithm was designed to work with realistic PKPD datasets,
which can contain multiple PK and PD observations variables, multiple PD
observation types (binary, ordinal, or continuous observations), and
BLOQ data. The algorithm also provides additional mechanisms for masking
the event skeleton (dose and observation times), and handles
weight-based dosing.

This algorithm offers no formal privacy guarantees, nor does it
guarantee preservation of the true relationship between dose, exposure,
observation variables and covariates. The algorithm is desgined to
produce synthetic datasets that carries many of the properties of the
source data, for use in developing and debugging analysis code.

## References

1.  Guillaudeux M, Rousseau O, Petot J, et al. Patient-centric synthetic
    data generation, no reason to risk re-identification in biomedical
    data analysis. *npj Digital Medicine.* 2023;6. doi:
    [10.1038/s41746-023-00771-5](https://doi.org/10.1038/s41746-023-00771-5).

2.  Destere A, Lombardi R, Labriffe M, et al. *Can synthetic data
    overcome the privacy and fidelity bottleneck in Pharmacometrics? A
    comparative benchmark using a daptomycin population pharmacokinetic
    model.* medRxiv preprint, posted June 2, 2026. doi:
    [10.64898/2026.05.30.26354512](https://doi.org/10.64898/2026.05.30.26354512).
