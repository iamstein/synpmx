# Declare pharmacometric column roles

Column roles are explicit: `synpmx` does not infer critical PMX
semantics from column names. The declaration is also the complete
manifest of what survives into synthetic data.
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
drops every column not named by some role, so a column you forget is
dropped rather than silently copied out of a real subject. Name a column
in `keep` to carry it through.

## Usage

``` r
pmx_roles(
  id,
  time,
  dv,
  amt = NULL,
  evid,
  cmt = NULL,
  dvid = NULL,
  mdv = NULL,
  rate = NULL,
  nominal_time = NULL,
  tad = NULL,
  occasion = NULL,
  cens = NULL,
  limit = NULL,
  addl = NULL,
  ii = NULL,
  covariates = NULL,
  strata = NULL,
  dose_covariate = NULL,
  assigned_dose = NULL,
  keep = NULL,
  exclude = NULL
)
```

## Arguments

- id, time, dv, evid:

  Required single column names for subject ID, actual time, dependent
  variable, and event ID.

- amt, cmt, mdv, rate:

  Optional single column names for amount, compartment, missing-DV
  indicator, and infusion rate.

- dvid:

  Endpoint-key column(s). Usually one column. A dataset that labels the
  same endpoint two ways — a numeric `YTYPE` beside a character `NAME` —
  may declare both, `dvid = c("YTYPE", "NAME")`. The first is the
  grouping key; validation checks the rest are a consistent 1:1 mapping
  with it and errors if they disagree, and
  [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
  carries all of them through.

  **`dvid` and `cmt` answer different questions.** `dvid` is which
  endpoint a measurement is; `cmt` is which compartment a dose enters,
  and it is read only on event rows. Nothing infers one from the other,
  so a source with more than one endpoint needs `dvid` — without it
  every measurement is treated as one endpoint, sharing a single value
  transform and a single censoring boundary.
  [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
  refuses rather than let that happen silently.

  A NONMEM `CMT` usually does both jobs, so **one column may be named as
  both roles**: `pmx_roles(..., cmt = "CMT", dvid = "CMT")`. This is the
  only permitted overlap; every other collision is an error.

- nominal_time, occasion:

  Optional time metadata columns.

- tad:

  Time-after-dose column. **It is an output, not an input.**
  [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
  recomputes it from the generated times and the generated dose rows and
  overwrites whatever the source held, so declaring the role says which
  column to overwrite and to carry through — the source's values are
  never used to generate anything.

  [`validate_pmx()`](https://iamstein.github.io/synpmx/reference/validate_pmx.md)
  does read them, and reports (as a non-fatal warning) where the
  declared column disagrees with time since the most recent dose row.
  That disagreement is worth knowing: a study may measure TAD from the
  end of an infusion, from a nominal dose time, or from an assigned
  dosing occasion rather than the most recent dose.
  [`nlmixr2data::nimoData`](https://nlmixr2.github.io/nlmixr2data/reference/nimoData.html)
  disagrees on 45% of its observation rows. Where it does, the synthetic
  column follows the derivation rather than the source's convention.

  Two limits. Samples taken before any dose are reported as 0, because
  [`validate_pmx()`](https://iamstein.github.io/synpmx/reference/validate_pmx.md)
  refuses a negative TAD — not because a baseline sample is genuinely
  zero hours after a dose it precedes. And where `addl` or `ii` is
  declared, the derivation cannot see the doses those imply, so the
  agreement check is skipped and says so.

- cens, limit:

  Optional Monolix-style censoring indicator and other interval-boundary
  columns.

- addl, ii:

  Optional additional-dose and interdose-interval columns.

- covariates:

  Baseline covariate column names, or `NULL`.

- strata:

  Treatment arm, dose group, cohort — any **assigned, subject-level
  stratum**, as opposed to a measured characteristic, which is a
  `covariate`. Must be constant within subject; subjects with no
  recorded value are grouped as their own stratum, with a warning rather
  than an error. Several columns may be named, and their combination
  defines the stratum.

  [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
  carries these verbatim from the subject that supplied the event
  skeleton and uses the stratum to group two things that are protocol
  properties rather than patient properties: the dose-to-covariate
  relationship, and the pool of attendance patterns an avatar may draw
  from. It is **not** a blending barrier — only route of administration
  is (see
  [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)),
  so donors are still borrowed across strata to reach the donor floor.

  The differential-privacy engines model the same columns jointly with
  the regimen as a released category domain.

- dose_covariate:

  The covariate the dose is a fixed multiple of, for weight-based or
  body-surface-area dosing: name the **covariate**, not the amount.
  Declaring it tells
  [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
  to recompute each avatar's `amt` from the avatar's own blended
  covariate at the anchor's own milligrams-per-unit, instead of copying
  the anchor's milligrams.

  This matters twice over. An avatar's covariates are blended while its
  `amt` would otherwise be copied verbatim, so under proportional dosing
  the avatar's own mg/kg comes out wrong — every generated patient
  violates the protocol it claims to follow. And a copied amount is one
  real patient's real dose, which under proportional dosing discloses
  that patient's weight exactly.

  Left `NULL`,
  [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
  tries to *infer* the relationship, and that inference deliberately
  fails closed: it requires the dose-to-covariate ratio to collapse onto
  a handful of protocol levels, so a study that rounds doses to vial
  sizes, or escalates within a patient by ratios that are not quite
  equal, is refused and the amounts are left alone. Declaring the
  covariate skips inference entirely and holds each dose row's own
  ratio, so **intra-patient escalation is preserved exactly** — three
  doses at 1, 2 and 4 mg/kg stay at 1, 2 and 4 mg/kg. The run report
  says which path was taken and, when inference declined, why.

- assigned_dose:

  Differential-privacy engines only. A nominal assigned-dose column
  reconstructed from the generated regimen.
  [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
  does not use this — carry the column with `keep`.

- keep:

  Columns to carry into synthetic data verbatim, copied from the same
  source subject that supplied the event skeleton, with no blending or
  synthesis. This is for assigned, subject-defining values you want kept
  faithful to a subject's dosing — a treatment arm, a dose group, a
  randomization sequence, or a redundant endpoint label such as a
  character `NAME` beside a numeric `dvid`. Because the value comes from
  the same anchor as the doses, it stays coherent with them. Contrast
  `covariates`, which are *blended* into new values across neighbours. A
  kept value is one real subject's real value, so use it only where the
  source data's own access controls and confidentiality obligations
  still apply.

- exclude:

  Differential-privacy engines only. Columns removed before private
  fitting, such as direct identifiers.
  [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
  does not use this — it drops every undeclared column by default, so
  not naming a column is how you drop it.

## Value

A `pmx_roles` object used by the fitting, generation, validation, and
comparison functions.

## Examples

``` r
roles <- pmx_roles(
  id = "ID", time = "TIME", dv = "DV", amt = "AMT",
  evid = "EVID", cmt = "CMT", tad = "TAD", covariates = "WT"
)

# Two columns for one endpoint, and a treatment arm carried through verbatim.
roles <- pmx_roles(
  id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
  dvid = c("YTYPE", "NAME"), covariates = "WT", keep = "ARM"
)
```
