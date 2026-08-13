# Compare per-covariate and per-endpoint distributions of source and synthetic

A numeric sanity check to run right after generating data. For each
baseline covariate and each endpoint (`dvid`), it summarizes the
distribution in the source and in the synthetic dataset side by side.
The dependent variable and continuous covariates get n, mean, standard
deviation, minimum, quartiles, and maximum; categorical covariates get
per-level counts and proportions.

## Usage

``` r
compare_pmx_distributions(source, synthetic = NULL, roles)
```

## Arguments

- source:

  Source PMX data.

- synthetic:

  Generated synthetic PMX data, or `NULL` to summarize `source` on its
  own.

- roles:

  Explicit roles from
  [`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md).

## Value

A `pmx_distribution_summary`: a list of `endpoints`,
`covariates_numeric`, and `covariates_categorical` data frames. Each is
`NULL` when the dataset declares no columns of that kind.

## Details

This is the distributional companion to
[`compare_pmx()`](https://iamstein.github.io/synpmx/reference/compare_pmx.md).
That function answers whether the *structure* matches — schema, event
grammar, row and event counts; this one answers whether the *numbers*
land in the same range. It is a diagnostic, not a validation of
statistical fidelity: AVATAR and the differentially private engines
deliberately do not reproduce source distributions exactly, so expect
the summaries to be close in magnitude and shape, not identical.

Every table is source-derived, so each is marked
`"restricted_not_releasable"`: it reads real covariate and endpoint
values and stays under the source data's access controls like any other
source-versus-synthetic diagnostic.

## See also

[`compare_pmx()`](https://iamstein.github.io/synpmx/reference/compare_pmx.md)
for the structural comparison.

## Examples

``` r
data <- pmx_simulated_fixture(20)
roles <- pmx_roles(
  id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
  cmt = "CMT", dvid = "DVID", covariates = c("WT", "SEX")
)
synthetic <- suppressWarnings(synpmx_avatar(data, roles, seed = 1))
#> synpmx_avatar(): dropped 8 undeclared column(s): NTIME, TAD, OCC, RATE, MDV, CENS, LIMIT, AGE.
#>   Declare a column in `keep` to carry it through verbatim.
compare_pmx_distributions(data, synthetic, roles)
#> Restricted PMX source-versus-synthetic distribution summary
#> 
#> Endpoints (dependent variable on observation rows):
#>  variable   dataset   n n_subjects mean   sd   min  q25 median  q75  max
#>        cp    source 160         20 3.75 2.88  0.88 1.33   2.66 5.29 9.52
#>        cp synthetic 160         20 3.95 3.18 0.628 1.33   2.76 5.98   14
#>        pd    source 120         20   60 14.3  35.2 48.6   59.6 71.3 89.6
#>        pd synthetic 120         20 57.9 16.8  21.5 43.6   55.7 69.7  101
#> 
#> Continuous covariates (baseline, per subject):
#>  variable   dataset  n mean   sd  min  q25 median  q75 max
#>        WT    source 20 71.4 9.01   58 62.6   73.8 79.5  82
#>        WT synthetic 20 72.9 7.92 61.1 64.7   76.7 79.8  81
#> 
#> Categorical covariates (baseline, per subject):
#>  variable   dataset  level  n proportion
#>       SEX    source female 10        0.5
#>       SEX    source   male 10        0.5
#>       SEX synthetic female 13       0.65
#>       SEX synthetic   male  7       0.35
#> 
#> Source-derived; not releasable unless separately public or privately budgeted.
```
