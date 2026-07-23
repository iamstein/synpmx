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
  subject_properties = NULL,
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

- nominal_time, tad, occasion:

  Optional time metadata columns.

- cens, limit:

  Optional Monolix-style censoring indicator and other interval-boundary
  columns.

- addl, ii:

  Optional additional-dose and interdose-interval columns.

- covariates:

  Baseline covariate column names, or `NULL`.

- subject_properties:

  Differential-privacy engines only. Subject-level assignment or
  grouping columns (`ACTARM`, `TRT`, a nominal dose group) modeled
  jointly with the regimen as a released category domain.
  [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
  does not use this — carry such a column with `keep`, which copies it
  verbatim from the subject that supplied the doses.

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
  kept value is one real subject's real value, so use it only inside a
  trusted environment.

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
