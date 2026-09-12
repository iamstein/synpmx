# What [`pmx_covariates_reference()`](https://iamstein.github.io/synpmx/reference/pmx_covariates_reference.md) holds

One row per available covariate, with its units, the distribution it
builds, and the `basis` that value came from. Read it before relying on
a reference value. Two columns earn the reading: `unit`, because the
units are the ones this package writes and a study recording height in
metres rather than centimetres needs its own declaration; and `basis`,
because a value marked as a placeholder is one the protocol is supposed
to decide.

## Usage

``` r
pmx_covariate_reference_table()
```

## Value

A data frame.

## See also

[`pmx_covariates_reference()`](https://iamstein.github.io/synpmx/reference/pmx_covariates_reference.md)

## Examples

``` r
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
```
