# Fully simulated public repeated-dose fixture

Creates a deterministic two-endpoint PMX event table for privacy-utility
and workflow tests. The data are generated from fixed formulas and
contain no patient or proprietary information, so no random seed is
required.

## Usage

``` r
pmx_simulated_fixture(n_subjects = 60L)
```

## Arguments

- n_subjects:

  Positive number of fully simulated subjects. The default is 60 to
  provide a larger companion to six- and twelve-subject examples.

## Value

A deterministic data frame with two doses per subject, dose-relative
`cp`, global study-time `pd`, explicit nominal/TAD/occasion roles, and
fixed public factor levels.

## Details

`TIME`, `NTIME`, and `TAD` are in hours; `AMT` is in arbitrary dose
units; `WT` is in kg; `AGE` is in years; `cp` uses arbitrary
concentration units; and `pd` uses arbitrary response units.

## Examples

``` r
public_data <- pmx_simulated_fixture(12)
length(unique(public_data$ID))
#> [1] 12
table(public_data$DVID, public_data$EVID == 0)
#>     
#>      FALSE TRUE
#>   cp    24   96
#>   pd     0   72
```
