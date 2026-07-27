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
#>  variable   dataset   n n_subjects   mean     sd     min    q25 median    q75
#>        cp    source 160         20  3.750  2.879  0.8800  1.335  2.656  5.288
#>        cp synthetic 160         20  3.985  3.088  0.7063  1.378  2.633  5.829
#>        pd    source 120         20 60.000 14.280 35.2000 48.640 59.610 71.310
#>        pd synthetic 120         20 62.980 20.490 26.5600 48.030 60.490 75.420
#>     max
#>    9.52
#>   10.74
#>   89.60
#>  121.70
#> 
#> Continuous covariates (baseline, per subject):
#>  variable   dataset  n  mean    sd   min   q25 median   q75   max
#>        WT    source 20 71.42 9.010 58.01 62.59  73.77 79.53 81.97
#>        WT synthetic 20 72.85 7.355 61.13 64.45  76.34 79.18 81.72
#> 
#> Categorical covariates (baseline, per subject):
#>  variable   dataset  level  n proportion
#>       SEX    source female 10       0.50
#>       SEX    source   male 10       0.50
#>       SEX synthetic female 13       0.65
#>       SEX synthetic   male  7       0.35
#> 
#> Source-derived; not releasable unless separately public or privately budgeted.
```
