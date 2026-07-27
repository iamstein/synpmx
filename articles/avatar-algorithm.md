# The AVATAR Algorithm

This article is the full specification of the default AVATAR generator:
every step from role declaration to the restored output schema, with the
mathematics, the edge cases, and a worked example.

It is the companion to
[`vignette("synpmx-method")`](https://iamstein.github.io/synpmx/articles/synpmx-method.md),
which introduces all four generation modes at a high level and says
which one to reach for. Read that first. Come here when you need to
defend, debug, or review what
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
actually did.

## The process at a glance

Before the step-by-step detail, here is the whole pipeline in one place.
Each synthetic subject is built from one real *anchor* subject’s event
structure, filled with measurements blended from several similar real
subjects.

The generator works on **any number of endpoints** — a PK concentration
alongside one or more PD measures, declared through the `dvid` role.
Each endpoint gets its own transformation and its own interpolation
grid, and they are then combined into a single per-subject profile, so
one donor set serves them all. Step 5 gives the details, including which
endpoint ends up dominating the choice of donors and why.

The pipeline has six **stages**, lettered (a) to (f). The rest of the
article walks through ten numbered **Steps**, and each stage below says
which Steps it covers — two stages span several Steps, which is why the
two are lettered and numbered separately rather than sharing one
sequence.

- **(a) Declare the columns** — *Step 1*. Roles say which column is the
  id, the time, the dose, the measurement, and so on. Only named columns
  survive.
- **(b) Split structure from values** — *Step 2*. Each subject’s rows
  divide into a fixed *event skeleton* — when doses and observations
  happen — and the measured DV values that fill it.
- **(c) Put every subject on a comparable footing** — *Steps 3, 4, and
  5*. Align time to the first event, choose an endpoint transformation,
  and reduce each subject to one numeric *profile* — baseline covariates
  plus its trajectory on a common grid — then to a handful of
  principal-component coordinates for measuring who resembles whom.
- **(d) Build each synthetic subject** — *Steps 6, 7, 8, and 9*. Copy a
  real anchor’s event skeleton — but by default do not anchor on a
  source subject whose structure is a gross outlier (a lone very long
  follow-up), so no avatar looks structurally extreme (`screen = TRUE`).
  Then pick the `k` (default 5) nearest donors, borrowing across
  dose/schedule groups when the anchor’s own group is too small — but
  never across route of administration; blend the donors’ trajectories
  onto the skeleton’s observation times, capping any one donor’s share
  at `max_donor_weight` (default 0.50); add subject-level and
  within-subject noise; and reconstruct any below-limit (BLOQ)
  censoring.
- **(e) Restore the original shape** — *Step 10*. Put back the source
  schema, column types, and conventions, and attach a record of what was
  done.
- **(f) Screen the result** — *after generation*, so it has no Step
  number. Because the event skeleton is copied from one anchor, a
  structurally unusual source subject yields a structurally unusual —
  and identifiable — avatar. A per-subject screen flags those, and a
  remediation step truncates, drops, or replaces them.

The two ideas that carry the most weight are the **anchor** (each avatar
keeps one real subject’s event structure) and the **donor blend** (its
measurements are averaged from several real subjects, never copied from
one). Keep those in mind and the rest is detail.

## Step 1: declare the meaning of the columns

The package never guesses critical PMX roles from column names. The user
must declare them with
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
  cens = "CENS",           # below-limit indicator, if the study has BLQ data
  limit = "LIMIT",         # the other end of an interval, when used
  covariates = c("WT", "AGE", "SEX"),
  keep = "ARM"             # carried verbatim; see the allowlist below
)
roles
#> Pharmacometric column roles:
#>   id: ID
#>   time: TIME
#>   nominal_time: <absent>
#>   tad: <absent>
#>   occasion: <absent>
#>   dv: DV
#>   amt: AMT
#>   evid: EVID
#>   cmt: CMT
#>   dvid: DVID
#>   mdv: MDV
#>   rate: RATE
#>   cens: CENS
#>   limit: LIMIT
#>   addl: <absent>
#>   ii: <absent>
#>   assigned_dose: <absent>
#>   covariates: WT, AGE, SEX
#>   keep: ARM
```

`id`, `time`, `dv`, and `evid` are required; the rest are optional.
`amt`, `cmt`, `mdv`, `rate`, `cens`, and `limit` each name one column.
`dvid` names one column, or several that encode the same endpoint (see
below). `covariates` and `keep` each take any number of columns. One
source column cannot be assigned to two roles. The roles
`subject_properties`, `assigned_dose`, and `exclude` belong to the
differentially private engines only;
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
rejects them and points at the role that does the job here.

### The role declaration is the allowlist

[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
keeps only columns a role names, and **drops every other column**,
reporting which ones. This is deliberate. AVATAR builds each synthetic
subject from one real *anchor* subject, and any column it retains that
it does not actively synthesize is copied from that anchor verbatim. A
secondary identifier — `USUBJID`, a site number, a randomization date —
left undeclared would otherwise ride a real subject’s real value
straight into the synthetic data. Dropping by default means a column you
forget fails safe, toward removal, rather than leaking.

That leaves four ways a column is treated, and the difference matters:

| Declared as | What AVATAR does to it | Use for |
|----|----|----|
| `dv` (with `cens`/`limit`) | Blended across donors into a new trajectory | the measurement |
| `covariates` | Blended/resampled across donors into a new value | baselines you want *synthesized* — weight, age |
| `keep` | Copied verbatim from the one anchor subject | assigned values you want kept faithful to that subject’s dosing — arm, dose group, a redundant endpoint label |
| *(undeclared)* | Dropped | anything you do not need |

`covariates` and `keep` are opposites, and choosing wrongly is the
common mistake. A **covariate is blended**: the synthetic subject gets a
genuinely new value interpolated from several neighbours, decoupled from
any one person. A **kept column is copied**: the synthetic subject gets
one real subject’s real value, unchanged. For a treatment arm, a dose
group, or a randomization sequence — anything that must stay consistent
with the doses AVATAR copied from the same anchor — `keep` is correct,
precisely because it never leaves that anchor’s side. Because a kept
value is a real subject’s real value, output containing it must stay
within the source data’s own access controls and obligations.

Redundant endpoint labels have their own handling. A dataset that
carries both a numeric `YTYPE` and a character `NAME` for the same
endpoint declares **both** as `dvid`: `dvid = c("YTYPE", "NAME")`. The
first is the grouping key; validation checks that the rest are a
consistent 1:1 mapping with it — catching a source where the two labels
disagree, which `keep` would carry through silently — and AVATAR carries
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
- when `cens` is declared, each CENS flag to agree with its DV — a
  left-censored value at the limit cannot exceed an ordinary measurement
  of the same endpoint. A pooled dataset with two assays at different
  limits will trip this legitimately; synthesize one study at a time, or
  drop the `cens` role.

Every error names the role and the column it maps to, and the count and
an example of the offending rows, so a failure says exactly what to fix.

Rows are eligible observations when EVID is zero and, if MDV is
declared, MDV is also zero. A nonmissing DV is additionally required
when a value is used to build a profile or donor trajectory. In symbols,
with bracket notation for a logical condition,

``` math
O_{ij} = [\mathrm{EVID}_{ij}=0]
         [\mathrm{MDV}_{ij}=0\;\text{if MDV is declared}]
         [Y_{ij}\;\text{is present}].
```

Rows with EVID equal to zero but missing DV, or with nonzero MDV, remain
in the sampled template but do not contribute a measurement value.

## Step 2: split structure from values

Let subject $`i`$ have source rows $`D_i`$. For every requested
synthetic subject, the generator samples an **anchor** $`a`$ with
replacement and begins with

``` math
S^{(0)} = D_a.
```

Thus $`S^{(0)}`$ contains the anchor’s full row count and row-wise
combination of TIME, EVID, AMT, RATE, CMT, DVID, and any other source
columns. This is the event skeleton. By copying these rows together, the
method avoids impossible records that can result from generating event
fields independently.

Only a few parts of this skeleton are subsequently changed:

- TIME may be jittered if `time_jitter > 0`;
- declared baseline covariates are replaced for all rows of the
  synthetic subject;
- eligible, originally present observation DVs are synthesized;
- standard MDV may be re-derived;
- ID is replaced; and
- rows are sorted by TIME, retaining original order within ties.

Event-row DV values and undeclared columns otherwise remain
anchor-template values. Consequently, users should not place a
semantically important column in the input and assume the package
understands it merely because it is present.

## Step 3: align time without pretending every subject was sampled identically

### First-event-relative time

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

### TIME is not TAD

A conventional time-after-dose (TAD) variable resets after every dose.
AVATAR does not use a declared `tad` role for alignment and never
recomputes a TAD column. With multiple doses, $`t^*`$ continues to
increase from the first qualifying dose; there is no later reset. A
declared `tad` column is retained and copied from the anchor template
unchanged; an undeclared one is dropped, as any undeclared column is.

This has two practical consequences:

1.  Use the actual analysis-time column as `time`. Do not expect a TAD
    column to drive neighbor construction or interpolation.
2.  If `time_jitter > 0`, a declared TAD column is copied, not
    recalculated, and can become inconsistent with the jittered TIME.
    Recompute TAD downstream if the workflow requires it.

### Different observation times across subjects

Two subjects sampled at different times cannot be compared as vectors
until they are put on shared coordinates. That is what the common grid
is for — and it is worth saying up front what it is *not* for, because
the 15-point cap below invites a misreading.

> **The grid is a ruler, not a resampler.** It is used only to build the
> profile that measures which subjects resemble each other (Step 5). It
> never limits the resolution of the generated data. A synthetic subject
> carries every one of its anchor’s observation times, however many that
> is, and the donor values placed there are interpolated from each
> donor’s **full** observed trajectory — not from its 15-point summary.
> A richly sampled source gives a richly sampled avatar.

For endpoint $`e`$, the package pools finite aligned observation times
and makes a common grid $`G_e`$. If there are at most 15 unique times,
all are retained. If there are more, 15 type-8 empirical quantiles are
used and duplicate quantiles are removed.

Each subject’s transformed trajectory is linearly interpolated onto
$`G_e`$. There is no extrapolation during profile construction, so grid
locations outside a subject’s observed window become missing profile
features. Tied times are averaged by `stats::approx(..., ties = mean)`.

The cap does bind on ordinary PK data: `theo_md` pools 156 distinct
aligned times, so donor ranking there runs on a 15-point summary of
trajectories with about 22 observations each, compressed further by the
PCA of Step 5. Two things limit the damage. The grid points are
*quantiles of the observed times*, so they land where sampling is dense
— early, around absorption and peak — rather than spreading evenly
across a long tail. And the grid is built per endpoint, so PK and PD
each get their own 15 rather than sharing them. The residual risk is
worst for multiple-dose data with several peaks and troughs, where 15
points across a long window can blur the shape that distinguishes
subjects. When that happens the cost is a *worse choice of donors*,
never a coarser output trajectory. `max_points` is currently fixed at 15
and is not an argument.

During final synthesis, donors are first interpolated directly to the
anchor’s aligned observation times. Suppose target time $`t`$ lies
outside donor $`r`$’s window. Rather than extrapolating indefinitely,
the method maps the target’s relative position in the whole anchor
observation window to the donor window:

``` math
q(t)=\frac{t-a_{\mathrm{target}}}
           {b_{\mathrm{target}}-a_{\mathrm{target}}},
\qquad
u_r(t)=a_r+q(t)(b_r-a_r),
```

where $`[a_{\mathrm{target}},b_{\mathrm{target}}]`$ and $`[a_r,b_r]`$
are the target and donor ranges. The donor trajectory is evaluated at
$`u_r(t)`$. This normalized-window fallback preserves progress through a
trajectory, but it is not a mechanistic warping model. A donor with only
one unique time supplies its mean at every requested target time.

### Optional coherent time jitter

The default `time_jitter = 0` copies template times exactly. Otherwise,
one normal perturbation is drawn for each unique template time $`u_l`$,

``` math
u'_l=u_l+\delta_l, \qquad
\delta_l\sim\mathcal{N}(0,\sigma_t^2),
```

and constrained between the midpoints of adjacent nominal times.
Nonnegative time is enforced at the first time. All rows tied at $`u_l`$
receive the same $`u'_l`$, so a dose and observation tied at a nominal
time remain tied and time ordering cannot cross.

## Step 4: choose an endpoint transformation

### Why there is a transformation at all

Blending happens *on the transformed scale*. So this step is really the
choice of **what “average” means** when five patients are mixed
together.

Concentrations span orders of magnitude. Averaging 0.5, 0.8, and 12.0 on
the raw scale gives 4.4 — a number dominated by the single high patient
and representative of none of them. Averaging on the log scale gives
about 1.6, the geometric mean, which is the pharmacometric convention
for precisely this reason. Blend concentrations on the raw scale instead
and one high-exposure donor drags every avatar it touches upward.

Logs have one problem with PK data: exact zeros. Pre-dose rows and
imputed BLQ values are genuinely 0, and $`\log 0=-\infty`$ would poison
the entire blend. So a small constant $`c_e`$ is added before taking
logs. Choosing $`c_e`$ to be **half the smallest positive value in that
endpoint** puts zeros just below the smallest real measurement — near
the bottom of the observed range, where a zero belongs — and scales it
to the data rather than fixing an arbitrary $`10^{-6}`$.

Not every endpoint wants logs. A PD score, a change from baseline, or a
temperature difference can legitimately be negative, where a log is
undefined. The “positive-like” test below is what separates the two
cases; its 1% tolerance lets a concentration endpoint with a couple of
small negative assay readings still count as positive-like.

### The rule

Transformations are selected separately for each DVID endpoint, or once
for the implicit endpoint `"DV"` when DVID is absent. After discarding
nonfinite values, an endpoint is considered positive-like when it has at
least one positive value, its median is positive, and at most 1% of
values are negative.

For a positive-like endpoint, define

``` math
c_e=\max\left\{\frac{1}{2}\min(Y_e[Y_e>0]),\sqrt{\epsilon_{\mathrm{mach}}}\right\}
```

and use the offset-log pair

``` math
g_e(y)=\log\{\max(y,0)+c_e\},
\qquad
g_e^{-1}(z)=\max\{\exp(z)-c_e,0\}.
```

Other endpoints use $`g_e(y)=y`$. The truncations mean a positive-like
endpoint cannot generate a negative DV, whatever the noise of Step 9
does — a guarantee that comes free with the back-transformation rather
than needing a separate clamp. Transformation choices and offsets are
stored in `attr(synthetic, "pmx_settings")$endpoint_transforms`, so the
scale an endpoint was blended on is always recoverable from the output.

## Step 5: build one numeric profile per source subject

A profile concatenates baseline and longitudinal information:

``` math
\mathbf{x}_i=
\left[
  \mathbf{c}_i,\
  \mathbf{y}_{i1}(G_1),\ldots,\
  \mathbf{y}_{iE}(G_E)
\right].
```

Here $`\mathbf{c}_i`$ contains the first nonmissing value of each
declared baseline covariate. Numeric covariates contribute one feature.
Factors and character covariates contribute one indicator for every
factor level or observed character value. The vector
$`\mathbf{y}_{ie}(G_e)`$ is subject $`i`$’s transformed endpoint-$`e`$
trajectory interpolated on its common grid.

This profile is used only for distance calculations. It does not replace
the source data and is never returned as a synthetic patient.

### Multiple endpoints

[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
handles any number of endpoints — a PK concentration and one or more PD
measures, declared through `dvid`. Each is processed *separately* first
and only then merged, in this order:

1.  **Transform, per endpoint.** Step 4 runs once per endpoint, so PK
    may be blended on the log scale while a PD score stays on the
    identity scale.
2.  **Interpolate, per endpoint.** Each endpoint gets its *own* common
    grid, built from the pooled observation times of that endpoint
    alone. A sparsely sampled PD is not padded out to match a rich PK.
3.  **Concatenate.** Covariate features, then each endpoint’s grid
    columns, are laid side by side into the single row $`\mathbf{x}_i`$
    above.
4.  **Impute, drop, standardize.** Median imputation, removal of
    constant columns, then column-by-column standardization.
5.  **One PCA over everything**, giving a single $`\boldsymbol{\xi}_i`$
    that mixes all endpoints — not one PCA per endpoint.

Consequently there is **one distance, one donor set, and one set of
weights**, shared by every endpoint. A donor is chosen for overall
similarity; a donor that matches PK closely but PD poorly still donates
its PD.

Standardization is what makes endpoints commensurable at all: a log
concentration and a raw neutrophil count cannot be compared until each
column is divided by its across-subject standard deviation.

But note carefully what standardization does *not* do. It equalizes
**units**, not **endpoints**. Every grid point becomes one unit-variance
column, so an endpoint sampled at 15 times contributes 15 columns while
one sampled at 3 times contributes 3, and a Euclidean distance summing
over columns lets the first outvote the second by roughly five to one.
On a worked two-endpoint example with 7 usable PK columns, 3 PD columns,
and 1 covariate, the retained variance splits

| block     | features | share of retained variance |
|-----------|----------|----------------------------|
| covariate | 1        | 4.4%                       |
| PD        | 3        | 28.1%                      |
| PK        | 7        | 67.5%                      |

so donor selection is driven mostly by whichever endpoint was sampled
most densely, and baseline covariates are a small minority of the
distance whenever trajectories are rich. **This weighting is an accident
of the sampling design, not a modeling choice**, and there is currently
no way to reweight endpoints. If donors must match on a sparse endpoint,
that is a limitation to know about before trusting the selection.

Two related behaviors follow. A subject missing an endpoint entirely is
median-imputed on those columns, so it looks *average* there and sits
nearer the middle of the space than it deserves to. And at generation,
if a chosen donor has no value for the endpoint being filled, the
weights renormalize locally over the donors that do — so a missing
endpoint thins the blend rather than punching a hole in it.

### Missing features and scaling

An entirely missing feature is removed. Otherwise, a missing profile
entry is replaced with its across-subject feature median $`m_l`$:

``` math
x^{\mathrm{imp}}_{il}=
\begin{cases}
x_{il}, & x_{il}\text{ is finite},\\
m_l, & \text{otherwise}.
\end{cases}
```

This imputation exists only inside neighbor finding. It does not fill
source or output DV rows. Features whose standard deviation is not
greater than $`\sqrt{\epsilon_{\mathrm{mach}}}`$ are removed. The
remaining columns are standardized as

``` math
z_{il}=\frac{x^{\mathrm{imp}}_{il}-\bar{x}^{\mathrm{imp}}_l}{s_l}.
```

### Rank-safe PCA

When there are at least two subjects and two usable features, PCA is
applied to the standardized matrix with maximum rank

``` math
r_{\max}=\min(n-1,p).
```

Let $`\lambda_h`$ be the variance of principal component $`h`$. The
retained dimension $`H`$ is the smallest value satisfying

``` math
\frac{\sum_{h=1}^{H}\lambda_h}
     {\sum_{h=1}^{r_{\max}}\lambda_h}
\ge v,
```

where $`v`$ is `pca_variance` (default 0.90).

The **subject coordinate** is the output of this whole step and the
input to every distance in Step 6, so it gets its own symbol:

``` math
\boldsymbol{\xi}_i=(\xi_{i1},\ldots,\xi_{iH}),
\qquad
\xi_{ih}=\text{subject } i\text{'s score on principal component } h .
```

Read $`\boldsymbol{\xi}`$ as the Greek letter *xi*. It is one point per
subject in an $`H`$-dimensional space, and $`H`$ is typically a handful
— the number of components needed to reach $`v`$. Everything upstream
exists to produce it: the raw profile $`\mathbf{x}_i`$ is imputed to
$`x^{\mathrm{imp}}_{il}`$, standardized to $`z_{il}`$ so that a weight
in kilograms and a log-concentration compete on equal footing rather
than by the size of their units, and finally rotated by PCA into
$`\boldsymbol{\xi}_i`$.

If PCA is not available but standardized features exist, those features
are used directly as $`\boldsymbol{\xi}_i`$. If none exist, all subjects
receive the same one-dimensional zero coordinate.

The term **rank-safe** means the attempted PCA rank cannot exceed either
the feature rank or $`n-1`$. It does not imply that sparse or nearly
duplicated profiles contain useful separation.

## Below-the-limit-of-quantification (BLOQ) data and a Monolix-style CENS column

Declare `cens` and, if the source uses one, `limit` in
[`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md),
and AVATAR reconstructs censoring rather than copying it. The mechanism
has three parts.

**The boundary is read from the source.** AVATAR makes no formal privacy
claim, so it may recover the assay limit from the data itself: under the
Monolix convention a censored row reports the boundary in `DV`, with
`LIMIT` carrying the other end of an interval when there is one. Left,
right, and interval censoring are each recognised. A differentially
private engine may *not* do this — there the limit must be a declared
public input, because a boundary inferred from confidential records is
itself a release. That asymmetry is the reason the two engines derive
boundaries differently.

If a study used more than one assay limit, the most conservative
boundary is taken rather than inventing a per-row rule the source cannot
support.

**Censored donor values are imputed before blending.** Every censored
source row reports the same substituted number, so blending those values
directly would put a floor under the synthetic data that the real study
does not have. Each censored row is instead replaced by a draw inside
its censoring region — uniform below the limit for a left-censored row,
uniform within the interval for an interval-censored one — before the
trajectory is transformed and blended. A uniform draw rather than a
fixed LLOQ/2 avoids replacing one artificial spike with another.

**`DV`, `CENS`, and `LIMIT` are reconstructed together.** The
back-transformed blend is treated as a *latent* value: what the subject
would have measured with no assay limit. The boundary is then applied to
that latent, and all three columns are written from the result, so they
cannot disagree. A generated value below the limit produces `CENS = 1`
with `DV` exactly at the limit, reproducing the spike of identical
values that characterises real censored data.

Two honest limitations. The censored *fraction* is reproduced
approximately, not exactly: blending shrinks variance toward the middle
while the imputation draws uniformly across the censoring region, so the
synthetic rate typically runs somewhat above the source’s. And if a
`cens` role is declared but no boundary can be read from it — rows
flagged with `DV` missing, for instance — AVATAR warns and carries the
flag through untouched rather than guessing.

## Step 6: choose the donors

Similarity in PCA space is not sufficient on its own. A subject
receiving a single oral dose should not donate a trajectory to an
infusion or repeat-dose template. But the two ways a donor can be
unsuitable are not the same kind of problem, and the generator treats
them differently.

### Two structural axes, one of them absolute

**Route of administration is an absolute barrier.** A bolus, an
infusion, and an oral dose produce categorically different concentration
shapes. Blending across them does not yield a noisier version of the
truth; it yields a trajectory no protocol could have produced. Donors
are therefore *never* drawn from outside the anchor’s route, at any
stage, for any shortfall — including when that leaves the anchor with no
legal donor at all. Each subject’s route key is read from its dosing
rows only, as the set of

``` math
\bigl(\mathrm{EVID},\ \mathrm{CMT},\ \mathbb{1}[\mathrm{RATE}\neq 0]\bigr)
```

triples appearing on them, so that the compartment a dose enters and
whether it is delivered over time both bind, while dose *count* does
not: three oral doses and five oral doses share a route. NONMEM’s
`RATE < 0` (modeled rate or duration) counts as an infusion.

**Everything else is a preference.** Each subject also receives an event
signature containing:

- the ordered EVID values on event rows;
- optional CMT and DVID values on those rows;
- the sign and rounded numeric magnitude of optional AMT and RATE;
- the number of positive-dose starts, or all events if no positive AMT
  exists;
- successive start-time gaps rounded to two significant digits; and
- the set of endpoints on eligible observation rows.

Notice what is not in the signature: the observation-time schedule.
Event values still come from the anchor template, while unequal
observation times are handled by interpolation. Dose magnitude *is* in
the signature, so same-dose subjects are preferred — but unlike route it
is not a barrier, and the fallback below will cross it.

### The distance

For anchor $`a`$ and candidate donor $`r`$, the distance is the ordinary
straight-line (Euclidean) distance between their subject coordinates
from Step 5:

``` math
d_{ar}=\left\|\boldsymbol{\xi}_a-
                   \boldsymbol{\xi}_r\right\|_2
=\sqrt{\sum_{h=1}^{H}(\xi_{ah}-\xi_{rh})^2},
```

where

- $`h`$ indexes the **retained principal components**, running from 1 to
  $`H`$ — not subjects and not time points;
- $`\xi_{ah}`$ is the **anchor’s** score on component $`h`$, and
  $`\xi_{rh}`$ is **candidate $`r`$’s** score on that same component,
  both defined in Step 5;
- each squared term is one component’s disagreement between the two
  subjects, and the sum over all $`H`$ components, square-rooted, is how
  far apart they sit.

It is Pythagoras in $`H`$ dimensions; all the work went into choosing
the coordinates. Because the features were standardized before PCA, no
single original measurement dominates the distance merely by being
recorded in larger units.

This is the only distance in the selection rule. There is no separate
structural metric trading dose differences against schedule differences:
structure enters as the two-stage ordering below, not as a score.

### The selection algorithm

Write $`\mathcal{R}_a=\{r\neq a: \mathrm{route}_r=\mathrm{route}_a\}`$
for the route-compatible pool and $`\mathcal{S}_a=\{r\in\mathcal{R}_a:
\mathrm{sig}_r=\mathrm{sig}_a\}`$ for the exact-signature subset. To
choose $`k`$ donors for anchor $`a`$:

1.  **Exact stage.** Sort $`\mathcal{S}_a`$ by increasing $`d_{ar}`$,
    breaking ties by subject index, and take the first $`k`$. If
    $`|\mathcal{S}_a|\ge k`$, stop.
2.  **Fallback stage.** Sort $`\mathcal{R}_a\setminus\mathcal{S}_a`$ the
    same way and take enough of it to bring the total to $`k`$, whatever
    the dose or schedule. Record a warning that donors were borrowed
    across dose groups.
3.  **Shortfall.** If $`|\mathcal{R}_a|<k`$ the floor is unreachable,
    because the only remaining candidates are on another route and route
    is never crossed.

Ties break by subject index, so selection is deterministic under a fixed
seed. Stage-2 donors’ measurements are mapped onto the anchor’s own
observation times by interpolation, so the avatar keeps its anchor’s
regimen while its values are averaged across at least $`k`$ real
patients. This borrowing matters most for datasets with individualized
dosing, where an exact dose magnitude can make nearly every subject its
own signature group: in `theo_md`, weight-based dosing yields 11
signature groups across 12 subjects, so stage 2 does essentially all the
work.

### What stage 2 costs, stated plainly

Be clear about what the fallback distance does and does not know. `AMT`
is not a profile feature, so the stage-2 ranking never compares doses
directly — it compares baseline covariates and the shape and level of
the DV trajectory. Dose does influence it, but only *indirectly*: a
subject given a much larger dose has higher concentrations, and those
concentrations are profile features, so it lands further away. The
ranking therefore tends to prefer similar-dose donors without being told
to.

Tendency is not a guarantee, and nothing rescales. A borrowed donor’s
values are interpolated onto the anchor’s times and blended **raw**, by
deliberate design (structural realism over statistical fidelity). The
consequence is worth stating directly rather than leaving implicit: a
synthetic subject carries its anchor’s `AMT` while its concentrations
may be blended from donors given a different dose, so **the
dose–exposure relationship in AVATAR output is not guaranteed to be
preserved**. For developing and debugging model code — the use this
generator exists for — that is acceptable and understood. For estimating
parameters, it is not: fitting a structural model to AVATAR output can
recover a dose–exposure relationship that is an artifact of blending.
Use
[`synpmx_calibrated()`](https://iamstein.github.io/synpmx/reference/synpmx_calibrated.md)
when the parameters themselves have to mean something.

### When the floor cannot be reached

Because route is absolute, a route arm holding fewer than $`k+1`$
subjects can never supply a legal donor set. There is no good answer in
that situation, only a choice between two bad ones, and which is worse
depends on whether the sparse arm matters more than the disclosure risk
of reproducing it. That is the caller’s judgement, so
`on_donor_shortfall` makes it explicit:

| `on_donor_shortfall` | Behavior |
|----|----|
| `"drop"` *(default)* | Omit those subjects from the anchor pool. No avatar is built on them, and the synthetic cohort does not represent that arm. |
| `"noise"` | Keep them, blending however many same-route donors exist — possibly none — and relying on `subject_noise_sd` and `residual_noise_sd` for the rest. |
| `"error"` | Refuse to generate, naming the arm and both alternatives. |

Every branch alerts loudly, because each one silently changes something
a reader would want to know: what the cohort covers, or how identifying
it is.

`"noise"` is **available but not recommended**, and the reason is the
defect this whole mechanism exists to prevent. A subject blended from
one donor, or from none, is a noised near-copy of a real patient —
exactly the behavior `REV-025` was raised about. If you use it because
the sparse arm genuinely matters more, screen the result with
[`flag_identifiable_subjects()`](https://iamstein.github.io/synpmx/reference/flag_identifiable_subjects.md)
and treat those subjects as individually identifying.

One case forces the issue: when *every* arm is below the floor, `"drop"`
would leave nothing to generate, so generation proceeds as if `"noise"`
and says so.

## Step 7: randomize and cap donor weights

Suppose $`K`$ donors have been selected. Let
$`E_r\sim\mathrm{Exponential}(1)`$ and let $`R_r`$ be a random
permutation of $`1,\ldots,K`$. Importantly, $`R_r`$ is a randomized
rank, not the distance order. The raw donor weight is

``` math
q_r=\frac{E_r}{\max(d_{ar},\varepsilon)}2^{-R_r},
```

where

``` math
\varepsilon=
\begin{cases}
\max\{10^{-8},10^{-6}\operatorname{median}(d_{ar}:d_{ar}>0)\},
  & \text{if a positive distance exists},\\
10^{-8}, & \text{otherwise}.
\end{cases}
```

The first normalization is $`w_r=q_r/\sum_s q_s`$. If the raw weights
are nonfinite or have a nonpositive total, the implementation replaces
them with equal raw values before normalizing.

### Capping every donor, not just the largest

Each donor’s weight is its **share of the avatar**: five weights that
sum to 1. The floor $`k`$ says how many real patients go into the blend,
but it says nothing about how much of any one of them comes out. Five
donors weighted $`(0.95, 0.02,
0.02, 0.005, 0.005)`$ satisfy a floor of five while being, to any
practical purpose, a copy of one person. So the cap, not $`k`$, is what
controls how closely an avatar can resemble a single real patient.

The raw formula above produces very uneven weights on purpose — that is
the rank attenuation doing its job. Left alone it hands one donor a
median 58% of the avatar. `max_donor_weight` (written $`c`$ below,
default $`0.50`$) is the ceiling on any one share.

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
$`(0.50, 0.25, 0.15,
0.06, 0.04)`$ with a tight $`c = 0.30`$: cutting the leader to $`0.30`$
frees $`0.20`$, and the second donor’s share of it takes them from
$`0.25`$ to $`0.35`$ — over the ceiling. Capping only the largest weight
and stopping, which is what the implementation did throughout the
$`0.80`$ era, leaves exactly that violation in place: the stated
maximum, broken by the donor the redistribution created.

So the rule repeats — *pin everything over $`c`$, redistribute, look
again* — the standard construction called **water-filling**. Pinned
donors are never topped up, so each pass pins at least one more and it
terminates. The tight example settles at
$`(0.30, 0.30, 0.24, 0.096, 0.064)`$.

### At the default cap, one pass is provably enough

Worth knowing, because it means the loop above is not the complexity it
looks like. A second pass is needed only when redistribution lifts the
runner-up over the ceiling, that is when

``` math
(1-c)\,\frac{w_2}{1-w_1}>c .
```

At $`c=0.5`$ this requires $`w_2>1-w_1`$, i.e. $`w_1+w_2>1`$, which no
pair of weights can satisfy. **At any cap of $`0.5`$ or above the loop
runs at most once**, and simulation agrees exactly: at $`c=0.5`$ the cap
is untouched 33% of the time and pinned once 67%, never twice. Only
tighter caps exercise the iteration — at $`c=0.3`$ two passes are needed
53% of the time and three passes 11%.

The general routine is kept because `max_donor_weight` is an argument
and a user may set it below $`0.5`$, where a single pass is genuinely
wrong. But at the shipped default, the simple thing and the correct
thing coincide.

### The infeasible-cap boundary

With $`K`$ donors the weights average $`1/K`$, so **no weight vector can
satisfy a ceiling below $`1/K`$** — asking for $`c=0.3`$ from three
donors is asking three numbers below $`0.3`$ to sum to 1. Such a cap
relaxes to exactly $`1/K`$: uniform weights, the flattest blend
available, rather than an error. Two donors under $`c=0.5`$ give
$`(0.5,0.5)`$, which meets the cap exactly. With one donor the weight is
necessarily 1; a cap cannot invent a second donor.

### Why a cap is needed at all

It is reasonable to ask whether the cap earns its place, given that
$`k`$ already forces five donors into every avatar. The answer is that
**without a cap, the floor is largely decorative.** Simulating the raw
weight formula at $`k=5`$ with representative donor distances:

| uncapped largest donor share        | value       |
|-------------------------------------|-------------|
| median                              | 0.58        |
| interquartile range                 | 0.47 – 0.72 |
| $`P(\text{largest} > 0.5)`$         | 0.67        |
| $`P(\text{largest} > 0.8)`$         | 0.14        |
| effective donors $`1/\sum_r w_r^2`$ | 2.37        |

So the unconstrained formula typically puts **58% of an avatar into a
single donor**, and one avatar in seven takes more than 80% from one
person. Nominally five patients are blended; effectively about two and a
half are. The cap is the only thing that closes that gap.

### Choosing the cap

The useful diagnostic is **how often the cap fires**, because that says
what role it is playing:

| `max_donor_weight` | binds on | effective donors |
|--------------------|----------|------------------|
| 0.30               | 99%      | 3.92             |
| 0.40               | 89%      | 3.28             |
| 0.50 *(default)*   | 68%      | 2.90             |
| 0.60               | 47%      | 2.64             |
| 0.80               | 15%      | 2.40             |
| none               | 0%       | 2.37             |

Both ends of that table are the wrong kind of parameter. A cap that
fires on 99% of subjects is not a guardrail — it *is* the weighting
scheme, and the inverse-distance term underneath it barely matters. A
cap at $`0.80`$ fires on 15%, trimming only the worst tail while still
permitting one patient to be four-fifths of an avatar; that was the old
default, and it was never really protecting anything.

The default is **$`0.50`$**, which states as a single checkable
sentence: *no single real patient is more than half of any synthetic
patient.* It fires on about two thirds of subjects, so it genuinely
constrains without wholly replacing the distance weighting.

Its honest cost is that “five donors are blended” is really “about three
effective donors” ($`2.90`$). Tightening to $`0.30`$ buys about one more
effective donor at the price of a cap that fires essentially always.
Both are defensible; $`0.50`$ is the one that can be explained in a
sentence.

Because these numbers depend on the cohort, `pmx_settings` records
`cap_binding_fraction` and `mean_effective_donors` for every run, so the
same question can be asked of real data rather than of a simulation.

### What the cap costs

Treating donors as independent, a weighted blend retains
$`\sum_r w_r^2`$ of individual variance, so the synthetic cohort’s
between-subject standard deviation shrinks by $`\sqrt{\sum_r w_r^2}`$.
The reciprocal $`1/\sum_r w_r^2`$ — the effective number of donors,
reported as `mean_effective_donors` in `pmx_settings` — reads as a
privacy floor and as the variability cost in one number. It equals $`K`$
for uniform weights and $`1`$ for a sole donor.

That independence assumption is, however, pessimistic in the direction
that matters. Donors are *nearest neighbours*, so they are strongly
positively correlated and averaging them destroys much less variance
than the formula implies; the subject-level perturbation of Step 9 then
restores more. Measured on `theo_md` — between-subject SD of
$`\log\mathrm{AUC}`$, averaged over 20 seeds, against a source SD of
$`0.273`$ — the cap turns out to be nearly free:

| `max_donor_weight` | effective donors | BSV retained |
|--------------------|------------------|--------------|
| 0.80               | 2.50             | 72%          |
| 0.50               | 2.93             | 75%          |
| 0.30               | 3.91             | 78%          |
| 0.25               | 4.43             | 74%          |
| 0.20 (uniform)     | 5.00             | 68%          |

Tightening the cap raises the effective donor count substantially — 2.4
with no cap, 2.9 at the default $`0.50`$, 3.9 at $`0.30`$ — while
between-subject variability stays flat within noise at a 12-subject
source. The independence formula would have predicted a fall from 81% to
49% over that range; correlation among neighbours is why it does not
happen. This is one dataset and one summary, so the table is evidence
that capping is affordable here, not a general result. It is also the
reason the cap can be chosen on privacy grounds: the fidelity side of
the trade is nearly flat.

The same selected subjects and weights are used for every declared
covariate and endpoint. When donor $`r`$ lacks a usable value at target
location $`j`$, the available weights are renormalized locally:

``` math
w_{rj}^{*}=\frac{w_r I_{rj}}
                 {\sum_s w_s I_{sj}},
```

where $`I_{rj}=1`$ when the value is finite and zero otherwise.

## Step 8: synthesize baseline covariates

For a numeric covariate, let $`c_r`$ be donor $`r`$’s first nonmissing
value. When all available values are positive and

``` math
\frac{\max_r c_r}{\operatorname{median}_r(c_r)}>3,
```

the blend is formed on the log scale; otherwise it is formed on the
original scale. With working-scale values $`h(c_r)`$, the donor center
is

``` math
\mu_c=\sum_r w_r^* h(c_r).
```

Let $`s_c`$ be the ordinary sample standard deviation of the available
working-scale donor values. If it is unavailable or zero, the fallback
is $`\max(0.05|\mu_c|,0.01)`$. The generated working-scale covariate is

``` math
c^*=\mu_c+\eta_c,
\qquad
\eta_c\sim\mathcal{N}(0,(\sigma_{\mathrm{subj}}s_c)^2),
```

where `subject_noise_sd` is $`\sigma_{\mathrm{subj}}`$ (default 0.15). A
log-scale result is exponentiated, and any all-positive covariate is
floored above zero. The final value is repeated on every row of that
synthetic subject.

For factor, character, and logical covariates, one available donor
category is sampled using the locally normalized weights. Categories are
not averaged.

### A covariate and an endpoint are blended independently

Covariates (this step) and endpoint trajectories (the next) are
generated in separate passes, each with its own donor weighting and its
own perturbation. That is fine when they are unrelated, but some
datasets carry a **baseline covariate that is the same quantity as a
longitudinal endpoint** — a `B0` baseline B-cell count beside a B-cell
kinetic endpoint, a baseline biomarker beside its own time course.
AVATAR does not know they are linked, so a synthetic subject’s baseline
covariate need not equal the baseline of its own generated trajectory.
The two are individually plausible but not mutually consistent.

This is a known limitation (`REV-022`). It is usually harmless for
workflow-development use, where each column need only be realistic on
its own. If your analysis relies on a covariate agreeing with an
endpoint’s baseline, reconcile them after generation, or do not declare
the redundant covariate and read the baseline from the trajectory
instead.

## Step 9: synthesize each endpoint trajectory

Only positions that were eligible, nonmissing observations in the anchor
template receive synthesized values. Thus the anchor supplies the number
and placement of observations, including its missing-DV pattern.

For endpoint $`e`$, donor trajectories are transformed with $`g_e`$ and
interpolated to the anchor observation times. Let $`z_{rj}`$ be donor
$`r`$’s transformed value at target position $`j`$. The deterministic
blend is

``` math
\bar z_j=\sum_r w_{rj}^{*}z_{rj}.
```

If every donor is unavailable at a target, the transformed dataset
median for that endpoint is used and a warning is stored.

The endpoint noise scale $`s_e`$ is 1 on the offset-log scale. On the
identity scale it is the source endpoint standard deviation, with
fallback $`\max(0.1|\operatorname{median}(Y_e)|,0.01)`$ when necessary.
A single subject-level shift is shared across all positions of the
endpoint:

``` math
b_e\sim\mathcal{N}(0,(\sigma_{\mathrm{subj}}s_e)^2).
```

Within-endpoint residual perturbations follow a stationary-scale
first-order autoregressive, AR(1), process in **observation order**:

``` math
\begin{aligned}
\epsilon_{e1} &\sim
  \mathcal{N}(0,(\sigma_{\mathrm{res}}s_e)^2),\\
\epsilon_{ej} &= \phi\epsilon_{e,j-1}+\nu_{ej},\\
\nu_{ej} &\sim
  \mathcal{N}\left(0,
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
subject coherence. This is a synthesis device, not a fitted multivariate
residual model.

## A transparent worked example

The small dataset below has six subjects, one dose each, one positive
endpoint, two baseline covariates, and deliberately different
observation times. It is large enough to have non-anchor neighbors but
small enough to inspect.

``` r

make_subject <- function(id, observation_time, wt, sex, peak) {
  time <- c(0, observation_time)
  evid <- c(1L, rep(0L, length(observation_time)))
  data.frame(
    ID = as.integer(id),
    TIME = time,
    DV = c(0, peak * exp(-0.55 * observation_time) + 0.03 * id),
    AMT = c(100, rep(0, length(observation_time))),
    EVID = evid,
    CMT = c(1L, rep(2L, length(observation_time))),
    MDV = ifelse(evid == 0L, 0L, 1L),
    WT = rep(c(58, 64, 70, 76, 82, 88)[id], length(time)),
    SEX = factor(
      rep(sex, length(time)),
      levels = c("female", "male")
    )
  )
}

worked_source <- do.call(rbind, list(
  make_subject(1, c(0.5, 1.0, 2.0, 4.0), 58, "female", 9.0),
  make_subject(2, c(0.25, 1.5, 3.0, 6.0), 64, "male", 10.0),
  make_subject(3, c(0.75, 1.25, 2.5, 5.0), 70, "female", 11.0),
  make_subject(4, c(0.4, 1.8, 3.5, 5.5), 76, "male", 12.0),
  make_subject(5, c(0.6, 1.4, 2.8, 4.5), 82, "female", 13.0),
  make_subject(6, c(0.3, 1.1, 3.2, 6.5), 88, "male", 14.0)
))
rownames(worked_source) <- NULL

worked_roles <- pmx_roles(
  id = "ID", time = "TIME", dv = "DV", amt = "AMT",
  evid = "EVID", cmt = "CMT", mdv = "MDV",
  covariates = c("WT", "SEX")
)

knitr::kable(
  worked_source[worked_source$ID %in% 1:2, ],
  digits = 3,
  caption = "The first two source subjects"
)
```

|  ID | TIME |    DV | AMT | EVID | CMT | MDV |  WT | SEX    |
|----:|-----:|------:|----:|-----:|----:|----:|----:|:-------|
|   1 | 0.00 | 0.000 | 100 |    1 |   1 |   1 |  58 | female |
|   1 | 0.50 | 6.866 |   0 |    0 |   2 |   0 |  58 | female |
|   1 | 1.00 | 5.223 |   0 |    0 |   2 |   0 |  58 | female |
|   1 | 2.00 | 3.026 |   0 |    0 |   2 |   0 |  58 | female |
|   1 | 4.00 | 1.027 |   0 |    0 |   2 |   0 |  58 | female |
|   2 | 0.00 | 0.000 | 100 |    1 |   1 |   1 |  64 | male   |
|   2 | 0.25 | 8.775 |   0 |    0 |   2 |   0 |  64 | male   |
|   2 | 1.50 | 4.442 |   0 |    0 |   2 |   0 |  64 | male   |
|   2 | 3.00 | 1.980 |   0 |    0 |   2 |   0 |  64 | male   |
|   2 | 6.00 | 0.429 |   0 |    0 |   2 |   0 |  64 | male   |

The first two source subjects {.table}

The next chunk uses internal helpers solely to expose a reproducible
teaching trace. These helpers are not public API. Ordinary analysis code
should call
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
rather than depend on them.

With `seed = 2026`, the anchor and its selected compatible donors are:

|     | role   | source_ID | profile_distance | weight |
|:----|:-------|----------:|-----------------:|-------:|
|     | anchor |         5 |               NA |     NA |
| 4   | donor  |         4 |           1.7462 | 0.4364 |
| 3   | donor  |         3 |           2.7489 | 0.0032 |
| 6   | donor  |         6 |           3.2730 | 0.5000 |
| 1   | donor  |         1 |           6.0532 | 0.0070 |
| 2   | donor  |         2 |           6.1283 | 0.0535 |

This trace shows the cap of Step 7 doing its work, and shows why the cap
is not optional. One donor is pinned at exactly `max_donor_weight`
$`=0.50`$; without the cap that donor would have taken more. Notice too
that the pinned donor is *not* the nearest one — the randomized rank
$`R_r`$, not the distance order, decides who dominates a given avatar,
which is what stops the nearest neighbour from being systematically the
largest contributor.

The recorded `mean_effective_donors`, $`1/\sum_r w_r^2`$, is 2.26
against a floor of 5 — the gap between “five donors were used” and “the
blend is worth about 2.3 independent donors”. That gap is the honest
cost of a cap chosen for simplicity, and it is why
`cap_binding_fraction` and `mean_effective_donors` are both recorded
rather than left to be assumed.

The anchor contributes the event skeleton. Donors contribute transformed
DVs after interpolation to the anchor times. The following table shows
the exact pre-noise blend used by the implementation.

| anchor_TIME | donor_4_z | donor_3_z | donor_6_z | donor_1_z | donor_2_z | blended_z | deterministic_DV | final_synthetic_DV |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.6 | 2.194 | 2.026 | 2.348 | 1.905 | 2.012 | 2.259 | 9.360 | 8.523 |
| 1.4 | 1.776 | 1.688 | 1.934 | 1.486 | 1.591 | 1.843 | 6.100 | 5.428 |
| 2.8 | 1.078 | 0.992 | 1.231 | 0.792 | 0.886 | 1.142 | 2.919 | 2.607 |
| 4.5 | 0.324 | 0.231 | 0.530 | 0.216 | 0.172 | 0.418 | 1.304 | 1.159 |

Interpolation and blending for the anchor endpoint; z is the endpoint
working scale {.table style="width:100%;"}

`deterministic_DV` is the back-transformed donor blend.
`final_synthetic_DV` also contains the one subject shift and AR(1)
residual sequence. The event rows, new ID, and subject-constant
synthesized covariates can be inspected directly:

``` r

knitr::kable(worked_synthetic, digits = 3)
```

|  ID | TIME |    DV | AMT | EVID | CMT | MDV |     WT | SEX  |
|----:|-----:|------:|----:|-----:|----:|----:|-------:|:-----|
|   7 |  0.0 | 0.000 | 100 |    1 |   1 |   1 | 80.579 | male |
|   7 |  0.6 | 8.523 |   0 |    0 |   2 |   0 | 80.579 | male |
|   7 |  1.4 | 5.428 |   0 |    0 |   2 |   0 | 80.579 | male |
|   7 |  2.8 | 2.607 |   0 |    0 |   2 |   0 | 80.579 | male |
|   7 |  4.5 | 1.159 |   0 |    0 |   2 |   0 | 80.579 | male |

The result records enough settings to audit the public generator call:

| setting | value |
|:---|:---|
| seed | 2026 |
| n_subjects | 1 |
| k | 5 |
| pca_variance | 0.9 |
| subject_noise_sd | 0.15 |
| residual_noise_sd | 0.05 |
| residual_phi | 0.6 |
| time_jitter | 0 |
| alignment | time relative to first positive dose within compatible schedules; normalized observation-window fallback |
| compatible_event_groups | 1 |
| routes | 1 |
| on_donor_shortfall | drop |
| max_donor_weight | 0.5 |
| cap_binding_fraction | 1 |
| mean_effective_donors | 2.256 |
| warnings |  |

Recorded generator settings {.table}

| endpoint | method     |  offset | positive |
|:---------|:-----------|--------:|:---------|
| DV       | log_offset | 0.21442 | TRUE     |

Recorded endpoint transformation {.table}

## Step 10: restore PMX schema and deterministic conventions

After all synthetic subjects are row-bound, the original column order is
restored. Each result column is converted toward its source class:

- factors recover source levels and ordering; new factor IDs add new
  levels;
- integer columns are rounded and converted to integer;
- double, logical, character, Date, and POSIXct columns are restored;
  and
- a non-`data.frame` source class, such as a tibble class vector, is
  reapplied.

New numeric IDs start above the source maximum. New character and factor
IDs use labels such as `syn_001`, with a prefix added repeatedly if
needed to avoid collision.

If MDV is declared and the entire source obeys the standard relationship

``` math
[\mathrm{MDV}=0]=[\mathrm{EVID}=0\ \text{and DV is present}],
```

MDV is re-derived from that rule in the synthetic data and restored to
its original class. For any nonstandard source-wide MDV convention,
anchor values are left as copied. Missing observation positions are not
invented or removed.

Finally, `validate_pmx(..., strict = TRUE)` is run on the assembled
result. A generation call cannot silently return an output that fails
the package’s own structural requirements.

## Reproducibility and diagnostics

All random operations occur inside a local seed scope. `seed` must be an
integer between zero and `.Machine$integer.max`. The caller’s RNG kind
and `.Random.seed` are restored on exit, whether or not `.Random.seed`
existed before the call. Consequently, identical inputs and arguments
produce an identical complete result without consuming the surrounding
analysis RNG stream.

`attr(synthetic, "pmx_settings")` records arguments, explicit roles,
endpoint transformations, the alignment description, the number of
compatible event groups and distinct routes, the donor weight cap and
how often it bound, the mean and minimum effective donor count, and
unique fallback warnings. It does not currently record the anchor ID,
donor IDs, distances, or realized weights for each synthetic subject.

Use the public checks for different questions:

- [`validate_pmx()`](https://iamstein.github.io/synpmx/reference/validate_pmx.md)
  asks whether one dataset is structurally usable under the declared
  roles.
- [`compare_pmx()`](https://iamstein.github.io/synpmx/reference/compare_pmx.md)
  reports source/synthetic row counts, subjects, event counts,
  endpoints, column classes, validations, and exploratory trajectories
  when `ggplot2` is installed.
- [`compare_pmx_distributions()`](https://iamstein.github.io/synpmx/reference/compare_pmx_distributions.md)
  summarizes the dependent variable per endpoint and each baseline
  covariate, source against synthetic, so their ranges can be compared
  at a glance.

None of these establishes distributional equivalence, privacy, or model
fidelity.

## After generation: screening for identifiable subjects

The blending in Steps 6–9 protects the *measurements*: every avatar’s DV
values are averaged from at least `k` real subjects, so no one person’s
numbers survive verbatim. It does **not** protect the *structure*. Stage
(d) copies the event skeleton — the number of doses, their sizes, and
the observation times — from a single anchor. So a source subject with
an unusual structure (a lone long-followed patient, a one-off dose
level) produces an avatar with that same unusual structure, and a
structural oddity is enough to single someone out even when their values
are blended.

Two functions address this after the data is generated.

[`flag_identifiable_subjects()`](https://iamstein.github.io/synpmx/reference/flag_identifiable_subjects.md)
scores every subject, one axis at a time, with a robust median/MAD
outlier statistic (the Iglewicz–Hoaglin modified z-score, default cutoff
3.5) on four structural features:

- **follow-up time** — the last observation time, which catches the lone
  long-followed subject;
- **number of doses** — an unusual dosing-history length;
- **dose magnitude** — a rare dose level; and
- **DV value** — an extreme peak measurement.

A subject is flagged when it is an outlier on any axis, and the report
names the axes. Scoring is per subject on one summary per axis, so a
single long-followed patient is one flag on the time axis — not one flag
per late row — and one extreme patient is one flag on the DV axis. The
screen is a heuristic, not a privacy proof: a flag is a candidate to
remove, and an empty list is not a certificate of anonymity.

[`remediate_identifiable_subjects()`](https://iamstein.github.io/synpmx/reference/remediate_identifiable_subjects.md)
acts on those flags. Its default policy splits by what a value-level
edit can honestly fix:

- a subject flagged **only** for a long follow-up is **truncated** back
  to the cohort’s longest ordinary follow-up — shortening a long
  timeline leaves a shorter but ordinary subject;
- a subject flagged for **any other** reason is **dropped**, because an
  extreme-DV subject is elevated throughout (trimming points only
  mangles it) and a rare dose cannot be trimmed without breaking the
  regimen;
- a subject flagged for both is dropped, since truncation would not
  resolve the other reason.

When a `source` is supplied, each dropped subject is **replaced**: fresh
avatars are generated, screened by the same policy, and appended with
new ids until the cohort is back to its original size — so screening
does not shrink the dataset. The `time` and `other` arguments expose the
policy.

This detect-and-remediate path is a deliberate alternative to preventing
structural outliers at generation time. A more thorough approach would
sample each avatar’s event skeleton from the cohort, so no avatar
carries any one real subject’s unique structure; that is a possible
future direction, not current behavior.

## Edge cases and documented fallbacks

The implementation favors an explicit warning and structurally valid
output over a silent incompatible blend.

| Situation | Behavior |
|:---|:---|
| Fewer than k same-signature donors | Borrow the nearest donors from other dose/schedule groups on the same route |
| Route arm smaller than the k-donor floor | Follow on_donor_shortfall: drop (default), noise, or error; alert loudly |
| Every route arm below the floor | Generate anyway and raise a loud alert; treat output as identifying |
| Cap below 1/K for K donors | Relax the cap to 1/K, i.e. uniform weights, rather than erroring |
| Source smaller than the k-donor floor | Blend all available donors and raise a loud alert |
| All selected distances essentially zero | Use epsilon-stabilized randomized weights |
| One donor time only | Repeat that donor’s mean at every target time |
| No donor value at a target | Use transformed endpoint median and record a warning |
| No variable profile feature | Use a shared zero coordinate; distance cannot separate subjects |
| Identical donors with all noise disabled | A source-shaped trajectory may be mathematically unavoidable; warn |
| Nonstandard MDV convention | Copy anchor MDV values rather than re-derive them |

Other boundaries are deliberate.
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
supports only `event_method = "template"` and
`dv_method = "avatar_blend"`. It does not fit a structural model,
reconstruct covariate-parameter relationships, simulate occasions,
enforce a dosing grammar beyond the sampled template, or provide a
censoring model. Fitting an explicit structural model is the job of
[`synpmx_calibrated()`](https://iamstein.github.io/synpmx/reference/synpmx_calibrated.md),
described below.

## Relationship to AVATAR

The method is **AVATAR-inspired**, not an implementation of published
AVATAR software. “AVATAR” is a method name rather than an initialism,
from the patient-centric *avatarization* literature in which each
synthetic record is built from the local neighborhood of real records.
The shared family resemblance is the use of standardized subject
features, PCA, local neighbors, inverse-distance influence, exponential
randomization, and rank attenuation. Destere and colleagues describe a
modified AVATAR benchmark for longitudinal population pharmacokinetic
(PopPK) data that first widened the data by subject and then used PCA,
K-nearest neighbors, inverse Euclidean distance, exponential stochastic
weights, and rank attenuation before reverse transformation \[1\]. The
original patient-centric AVATAR method was reported by Guillaudeux and
colleagues \[2\].

`synpmx` adapts those ideas to a PMX event table in several material
ways:

- it preserves a sampled longitudinal event template rather than
  generating all fields from a wide vector;
- route of administration is an absolute barrier, while dose and
  schedule differences only make a donor less preferred;
- endpoints are transformed and interpolated separately;
- the same donors support covariates and endpoints;
- every multi-donor weight, not only the dominant one, is capped at
  `max_donor_weight` (default 0.50);
- subject and AR(1) perturbations are added on the endpoint working
  scale; and
- there is no inverse-PCA reconstruction of a full synthetic patient
  vector.

Destere et al. also demonstrate why method labels are not validation:
their benchmark found materially different PK parameter and
residual-error biases across synthesis algorithms \[1\]. That paper is a
2026 medRxiv preprint and was not peer reviewed at the version supplied
with this package repository. Its privacy analyses do not transfer to
`synpmx`; this package has undergone no corresponding attack-based
privacy validation.

## Algorithm summary

What
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
actually executes, in order. These are implementation operations,
finer-grained than the ten Steps and not in one-to-one correspondence
with them, so the last column says where each is explained rather than
leaving two numberings to be guessed at.

| \# | Operation | Explained in |
|----|----|----|
| 1 | Validate source data, roles, and generator arguments | Step 1 |
| 2 | Compute first-event-relative time and endpoint-specific transformations | Steps 3–4 |
| 3 | Assemble baseline covariates and interpolate transformed endpoint trajectories into one profile per source subject | Step 5 |
| 4 | Median-impute profile features, remove unusable features, standardize, and retain enough PCA components to meet `pca_variance` when possible | Step 5 |
| 5 | Construct event/schedule signatures and route keys | Step 6 |
| 6 | Resolve anchors whose route arm cannot reach the donor floor per `on_donor_shortfall`, apply the structural screen, and sample the remaining anchors with replacement | Step 6 |
| 7 | For each anchor, copy its whole event template | Step 2 |
| 8 | Find `k` donors on the anchor’s route — exact-signature nearest-first, then the nearest remaining route-compatible subjects — and randomize their inverse-distance weights, water-filling every weight to `max_donor_weight` | Steps 6–7 |
| 9 | Generate subject-constant covariates from those donors | Step 8 |
| 10 | For each endpoint, interpolate donor trajectories to anchor times, blend, add a subject shift and AR(1) perturbations, and back-transform | Step 9 |
| 11 | Re-derive standard MDV when applicable, assign a new ID, preserve tie order, restore the source schema, and validate the assembled result | Step 10 |
| 12 | Attach `pmx_settings` and emit any collected fallback warnings | Step 10 |

## References

1.  Destere A, Lombardi R, Labriffe M, et al. *Can synthetic data
    overcome the privacy and fidelity bottleneck in Pharmacometrics? A
    comparative benchmark using a daptomycin population pharmacokinetic
    model.* medRxiv preprint, posted June 2, 2026. doi:
    [10.64898/2026.05.30.26354512](https://doi.org/10.64898/2026.05.30.26354512).

2.  Guillaudeux M, Rousseau O, Petot J, et al. Patient-centric synthetic
    data generation, no reason to risk re-identification in biomedical
    data analysis. *npj Digital Medicine.* 2023;6. doi:
    [10.1038/s41746-023-00771-5](https://doi.org/10.1038/s41746-023-00771-5).
