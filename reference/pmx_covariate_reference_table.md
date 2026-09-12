# What [`pmx_covariates_reference()`](https://iamstein.github.io/synpmx/reference/pmx_covariates_reference.md) holds

One row per available covariate, with its units and the distribution it
builds. Read it before relying on a reference value: the units are the
ones this package writes, and a study recording height in metres rather
than centimetres needs its own declaration.

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
#> 1                                                             75  18%
#> 2                                                            170   6%
#> 3                                                             26  18%
#> 4                                                            1.9  12%
#> 5                                                             45  30%
#> 6                                                            100  25%
#> 7                                                             95  25%
#> 8                                                             42  12%
#> 9                                                   M 50%, F 50% <NA>
#> 10 White 65%, Black or African American 12%, Asian 15%, Other 8% <NA>
#>         range
#> 1   35 to 160
#> 2  130 to 210
#> 3    15 to 55
#> 4  1.1 to 2.9
#> 5    18 to 90
#> 6   15 to 200
#> 7   15 to 180
#> 8    25 to 55
#> 9        <NA>
#> 10       <NA>
```
