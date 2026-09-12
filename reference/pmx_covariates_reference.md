# Reference distributions for the usual baseline covariates

Builds
[`pmx_covariate()`](https://iamstein.github.io/synpmx/reference/pmx_covariate.md)
declarations for body size, age, sex and the common renal and hepatic
laboratory values, so that a public-model generator has plausible
covariate columns without each distribution being elicited by hand. The
values are round numbers for a broad adult population.

## Usage

``` r
pmx_covariates_reference(names, medians = NULL, cvs = NULL)
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

- cvs:

  Optional named numeric replacing the reference coefficient of
  variation, as a proportion. Named by output column, the same way as
  `medians`. A trial's eligibility criteria usually make this smaller
  than the survey's.

## Value

A `pmx_covariates` object.

## Where the numbers come from

Body size comes from NHANES, the CDC's continuous National Health and
Nutrition Examination Survey, for adults aged 20 and over of both sexes.
The laboratory values come from conventional adult reference intervals,
which is what a stated normal range is. Every covariate records its own
basis in its `source`, and
[`pmx_covariate_reference_table()`](https://iamstein.github.io/synpmx/reference/pmx_covariate_reference_table.md)
prints the basis alongside the value. Adults only, and deliberately: a
paediatric population needs its own declaration through
[`pmx_covariate()`](https://iamstein.github.io/synpmx/reference/pmx_covariate.md).

## A survey population is not a trial cohort

NHANES describes the general adult population. Eligibility criteria cut
the tails off that, so a trial runs both lighter and tighter than the
survey – by weight, a median around 10 per cent below the survey's and a
coefficient of variation closer to 15 per cent than 25. A cohort
selected for obesity runs the other way. `medians` and `cvs` are how to
move each, and a regression test holds this gap against the adult
studies in this package's own public-data surveys so the starting point
cannot drift silently.

Age, sex and race have no population source worth quoting here. The
protocol's inclusion criteria and where it enrolled decide them, and a
phase 1 healthy-volunteer cohort, a renal-impairment study and an
oncology trial have nothing in common. Those three entries exist to put
a column in the table and say so in their `source`. State your own if
anything downstream reads them. No dataset in this package carries a
race column at all.

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
#>   WT: continuous [35, 180]; lognormal, median 84, CV 25% (public, DP)
#>   AGE: continuous [18, 90]; normal, median 45, CV 30% (public, DP)
#>   SEX: categorical {M 50%, F 50%} (public, DP)

# What the survey holds, and what it is.
pmx_covariate_reference_table()
#>    covariate                             quantity          unit distribution
#> 1         WT                          body weight            kg    lognormal
#> 2         HT                               height            cm    lognormal
#> 3        BMI                      body mass index         kg/m2    lognormal
#> 4        BSA                    body surface area            m2    lognormal
#> 5        AGE                                  age         years       normal
#> 6       CRCL                 creatinine clearance        mL/min    lognormal
#> 7       EGFR estimated glomerular filtration rate mL/min/1.73m2    lognormal
#> 8        ALB                        serum albumin           g/L       normal
#> 9        SEX                                  sex          <NA>  categorical
#> 10      RACE                                 race          <NA>  categorical
#>                                                          typical   cv
#> 1                                                             84  25%
#> 2                                                            168   6%
#> 3                                                             29  24%
#> 4                                                              2  14%
#> 5                                                             45  30%
#> 6                                                            100  25%
#> 7                                                             95  25%
#> 8                                                             43  10%
#> 9                                                   M 50%, F 50% <NA>
#> 10 White 65%, Black or African American 12%, Asian 15%, Other 8% <NA>
#>         range
#> 1   35 to 180
#> 2  130 to 210
#> 3    15 to 60
#> 4    1.1 to 3
#> 5    18 to 90
#> 6   15 to 200
#> 7   15 to 180
#> 8    25 to 55
#> 9        <NA>
#> 10       <NA>
#>                                                                                                basis
#> 1                                                   NHANES 2015-2018, adults 20 and over, both sexes
#> 2                                                   NHANES 2015-2018, adults 20 and over, both sexes
#> 3                                                   NHANES 2015-2018, adults 20 and over, both sexes
#> 4                          Mosteller formula applied to the NHANES 2015-2018 adult height and weight
#> 5                                         placeholder: the protocol's inclusion criteria decide this
#> 6                                       conventional adult reference interval, normal renal function
#> 7                                          KDIGO stage G1, normal or high glomerular filtration rate
#> 8                                                conventional adult reference interval, 35 to 50 g/L
#> 9                                                placeholder: the protocol's enrollment decides this
#> 10 placeholder: the OMB and FDA reporting categories; where a trial enrolled decides the proportions

# A trial that enrolled lighter and narrower than the general population,
# in a column named the way the study names it.
pmx_covariates_reference(c(WEIGHTB = "WT"), medians = c(WEIGHTB = 75),
                         cvs = c(WEIGHTB = 0.16))
#> Covariates
#>   WEIGHTB: continuous [35, 180]; lognormal, median 75, CV 16% (public, DP)

# `WEIGHTB` resolves on its own, so the mapping is only needed to rename.
pmx_covariates_reference("WEIGHTB")
#> Covariates
#>   WEIGHTB: continuous [35, 180]; lognormal, median 84, CV 25% (public, DP)
```
