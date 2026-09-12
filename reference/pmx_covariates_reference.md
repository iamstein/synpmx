# Reference distributions for the usual baseline covariates

Builds
[`pmx_covariate()`](https://iamstein.github.io/synpmx/reference/pmx_covariate.md)
declarations for body size, age, sex and the common renal and hepatic
laboratory values, so that a public-model generator has plausible
covariate columns without each distribution being elicited by hand. The
values are round numbers for a broad adult population.

## Usage

``` r
pmx_covariates_reference(
  names,
  medians = NULL,
  source = paste("synpmx reference distribution for a",
    "broad adult population; a package", "default, not measured from any study")
)
```

## Arguments

- names:

  Covariate column names. A plain vector takes each name as the
  quantity, so `c("WT", "AGE")` is body weight and age. A named vector
  maps a column to a quantity, so `c(BWT = "WT")` puts the weight
  distribution in a column called `BWT`. A handful of common spellings
  resolve on their own: `WEIGHTB`, `WEIGHT`, `WGT` and `BW` to `WT`,
  `HGT` and `HEIGHT` to `HT`, `GENDER` to `SEX`, `CLCR` to `CRCL`.

- medians:

  Optional named numeric replacing the reference median for those
  columns, in the units listed by
  [`pmx_covariate_reference_table()`](https://iamstein.github.io/synpmx/reference/pmx_covariate_reference_table.md).
  Named by output column, so `c(BWT = 82)` goes with
  `names = c(BWT = "WT")`.

- source:

  Provenance recorded on every covariate built here. The default says
  these are package reference values rather than a study's own.

## Value

A `pmx_covariates` object.

## What these numbers are, and are not

They are the same kind of input as the built-in structural models in
[`pmx_structural_model()`](https://iamstein.github.io/synpmx/reference/pmx_structural_model.md):
illustrative on purpose, so that covariate-handling code has something
to run against. They are not a rendering of any particular trial and
were not measured from one.

The spread is the half that travels. Body weight holds much the same
coefficient of variation across adult populations whatever its median
is, so the shape and the CV are reusable and the median is not: an
ordinary adult cohort and a cohort selected for obesity differ in the
median alone. `medians` is how to move it, and a regression test holds
the reference median against the adult studies in this package's own
surveys.

Age, sex and race are not population constants. The protocol's inclusion
criteria and where it enrolled decide them, and a phase 1
healthy-volunteer cohort, a renal-impairment study and an oncology trial
have nothing in common here. The entries put a column in the data. State
your own if anything downstream reads it. No dataset in this package
carries a race column at all.

Nothing here fits a neonatal or paediatric cohort, where body size
differs from an adult by more than an order of magnitude. Declare those
with
[`pmx_covariate()`](https://iamstein.github.io/synpmx/reference/pmx_covariate.md).

Before generated data crosses a trust boundary, replace these with the
protocol's own criteria or a published description of the population,
and record that in each covariate's `source`.

## See also

[`pmx_covariate_reference_table()`](https://iamstein.github.io/synpmx/reference/pmx_covariate_reference_table.md)
for what is available,
[`pmx_covariate()`](https://iamstein.github.io/synpmx/reference/pmx_covariate.md)
to state one by hand,
[`pmx_covariates_auto()`](https://iamstein.github.io/synpmx/reference/pmx_covariates_auto.md)
to resample from the data instead.

## Examples

``` r
pmx_covariates_reference(c("WT", "AGE", "SEX"))
#> Covariates
#>   WT: continuous [35, 160]; lognormal, median 75, CV 18% (public, DP)
#>   AGE: continuous [18, 90]; normal, median 45, CV 30% (public, DP)
#>   SEX: categorical {M 50%, F 50%} (public, DP)

# A heavier cohort, in a column named the way the study names it.
pmx_covariates_reference(c(WEIGHTB = "WT"), medians = c(WEIGHTB = 117))
#> Covariates
#>   WEIGHTB: continuous [35, 234]; lognormal, median 117, CV 18% (public, DP)

# `WEIGHTB` resolves on its own, so the mapping is only needed to rename.
pmx_covariates_reference("WEIGHTB")
#> Covariates
#>   WEIGHTB: continuous [35, 160]; lognormal, median 75, CV 18% (public, DP)
```
