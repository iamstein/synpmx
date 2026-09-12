# Generate a dataset from public inputs only

Simulates from a public structural model and a public protocol. No
confidential data is read, so there is nothing to protect and no budget
to spend: this is `epsilon = 0`, the strongest possible guarantee.

## Usage

``` r
synpmx_prior(
  model,
  design,
  roles = pmx_generated_roles(),
  n_subjects = NULL,
  seed = NULL,
  dropout = 0,
  lloq = NULL,
  covariates = NULL
)
```

## Arguments

- model:

  A
  [`pmx_structural_model()`](https://iamstein.github.io/synpmx/reference/pmx_structural_model.md).

- design:

  A
  [`pmx_trial_design()`](https://iamstein.github.io/synpmx/reference/pmx_trial_design.md).

- roles:

  A
  [`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md)
  naming the columns of the output. Defaults to
  [`pmx_generated_roles()`](https://iamstein.github.io/synpmx/reference/pmx_generated_roles.md),
  the package's own generated schema.

- n_subjects:

  Number of subjects. Defaults to the planned cohort total.

- seed:

  Ordinary generation seed. Unrelated to privacy noise.

- dropout:

  Fraction of subjects who discontinue early. A public assumption from
  the protocol.

- lloq:

  Lower limit of quantification. Observations below it are flagged
  `CENS = 1` with `DV` at the limit, following the Monolix convention.

- covariates:

  Optional
  [`pmx_covariates()`](https://iamstein.github.io/synpmx/reference/pmx_covariates.md).

## Value

A data frame in event-table form, under the names in `roles`.

## Details

The typical parameter values must come from somewhere that is not the
data – allometric scaling from preclinical work, a published model for
the compound class, or the reasoning that set the starting dose. The
output is exactly as good as that prior.

## Roles here name the output, not an input

Every generator in the package takes a
[`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md)
declaration and writes its output under those column names, so that one
declaration serves the generator,
[`compare_pmx_distributions()`](https://iamstein.github.io/synpmx/reference/compare_pmx_distributions.md)
and
[`synpmx_scorecard()`](https://iamstein.github.io/synpmx/reference/synpmx_scorecard.md).
For
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md),
[`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md)
and
[`synpmx_model()`](https://iamstein.github.io/synpmx/reference/synpmx_model.md)
the declaration describes the study being read. This function reads
nothing, so the same object does the opposite job here: it prescribes
the schema to generate into. Pass the roles of the study the synthetic
data has to stand in for and the output drops straight into code written
against it. A role left undeclared gets no column. The default,
[`pmx_generated_roles()`](https://iamstein.github.io/synpmx/reference/pmx_generated_roles.md),
names every column the generator can produce.

## See also

[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md),
[`synpmx_calibrated()`](https://iamstein.github.io/synpmx/reference/synpmx_calibrated.md),
[`synpmx_empirical()`](https://iamstein.github.io/synpmx/reference/synpmx_empirical.md)

## Examples

``` r
model <- pmx_structural_model(
  pk = "1cmt_oral", typical = c(cl = 6, v = 35, ka = 1.5),
  source = "illustrative allometric scaling"
)
design <- pmx_trial_design(
  dose_levels = 320, cohort_sizes = 12, sampling = c(0, 1, 2, 4, 9, 24),
  source = "illustrative protocol"
)
syn <- synpmx_prior(model, design, n_subjects = 12, seed = 202)
head(syn, 3)
#>   ID      TIME NTIME       TAD OCC       DV AMT RATE EVID CMT DVID MDV CENS
#> 1  1 0.0000000     0 0.0000000   1       NA 320    0    1   1 <NA>   1    0
#> 2  1 0.0000000     0 0.0000000   1 0.000000   0    0    0   2   cp   0    0
#> 3  1 0.9791034     1 0.9791034   1 6.823022   0    0    0   2   cp   0    0
#>   DOSE
#> 1  320
#> 2  320
#> 3  320

# Generating into a study's own schema, so the output needs no renaming.
study_roles <- pmx_roles(
  id = "SUBJID", time = "TIME", nominal_time = "NTIME", dv = "DV",
  amt = "AMT", evid = "EVID", cmt = "CMT", mdv = "MDV"
)
names(synpmx_prior(model, design, study_roles, n_subjects = 3, seed = 202))
#> [1] "SUBJID" "TIME"   "NTIME"  "DV"     "AMT"    "EVID"   "CMT"    "MDV"   
```
