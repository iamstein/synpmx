# Declare pharmacometric column roles

Column roles are explicit: `synpmx` does not infer critical PMX
semantics from column names. Columns listed in `exclude` are removed
before fitting and do not appear in generated data.

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
  exclude = NULL
)
```

## Arguments

- id, time, dv, evid:

  Required single column names for subject ID, actual time, dependent
  variable, and event ID.

- amt, cmt, dvid, mdv, rate:

  Optional single column names for amount, compartment, endpoint,
  missing-DV indicator, and infusion rate.

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

  Subject-level assignment or grouping columns, such as `ACTARM`, `TRT`,
  or a nominal dose group. These are treated as categorical, must be
  constant and nonmissing within a subject, and are modeled jointly with
  that subject's regimen rather than as independent baseline covariates.

- assigned_dose:

  Optional nominal assigned-dose column. It may vary by occasion but
  must be constant within subject and occasion and agree with the
  positive event amount. Generated values are derived from the generated
  regimen rather than sampled independently.

- exclude:

  Columns explicitly excluded before private fitting, such as direct
  identifiers. An ID role is still required as the privacy unit.

## Value

A `pmx_roles` object used by the fitting, generation, validation, and
comparison functions.

## Examples

``` r
roles <- pmx_roles(
  id = "ID", time = "TIME", dv = "DV", amt = "AMT",
  evid = "EVID", cmt = "CMT", tad = "TAD", covariates = "WT"
)
```
