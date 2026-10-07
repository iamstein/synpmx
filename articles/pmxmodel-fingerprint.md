# The PMX model study fingerprint

A study’s **fingerprint** is the set of numbers a synthetic study is
generated from. For
[`synpmx_model()`](https://iamstein.github.io/synpmx/reference/synpmx_model.md)
it is a population PK model, a time course for each PD endpoint, and a
summary of the study’s design: its arms, dose schedules, attendance and
covariate distributions. Synthetic data is drawn from the fingerprint
and random numbers alone, so the fingerprint is everything about the
source that can reach the output.

[`synpmx_model_estimate()`](https://iamstein.github.io/synpmx/reference/synpmx_model_estimate.md)
returns a fit. The fit holds the fingerprint plus diagnostics for
whoever ran it, such as the candidate models and the correlations
between covariates and random effects.
[`model_release()`](https://iamstein.github.io/synpmx/reference/model_release.md)
extracts the fingerprint, and that is the object to carry out of the
environment that holds the study.
[`synpmx_model_generate()`](https://iamstein.github.io/synpmx/reference/synpmx_model_generate.md)
accepts either and reads only the fingerprint.

This vignette walks through the fingerprint of one study. Read
[`vignette("pmxmodel-demo")`](https://iamstein.github.io/synpmx/articles/pmxmodel-demo.md)
first for a run end to end, and
[`vignette("pmxmodel-algorithm")`](https://iamstein.github.io/synpmx/articles/pmxmodel-algorithm.md)
for how each quantity is estimated.

## The study

[`nlmixr2data::warfarin`](https://nlmixr2.github.io/nlmixr2data/reference/warfarin.html):
32 patients, one oral dose scaled to body weight, a warfarin
concentration and a prothrombin complex activity (PCA).

``` r

data("warfarin", package = "nlmixr2data")
warfarin <- as.data.frame(warfarin)
warfarin$ntime <- warfarin$time
roles <- pmx_roles(
  id = "id", time = "time", nominal_time = "ntime", dv = "dv", amt = "amt",
  evid = "evid", dvid = "dvid", covariates = c("wt", "age", "sex")
)
```

``` r

fit <- synpmx_model_estimate(warfarin, roles, seed = 1)
```

This document reads a stored copy of that fit, so it needs no compiler.

``` r

release <- model_release(fit)
names(release)
#>  [1] "structural"           "parameters"           "pk_models"           
#>  [4] "pd"                   "arms"                 "dosing"              
#>  [7] "visits"               "cells"                "discrete"            
#> [10] "covariates"           "covariate_effects"    "schema"              
#> [13] "roles"                "endpoints"            "quantification_floor"
#> [16] "n_source"             "settings"             "privacy"
setdiff(names(fit), names(release))  # diagnostics that stay with the fit
#> [1] "candidates"   "design"       "correlations" "censoring"    "timing"      
#> [6] "movement"     "fit_subjects" "start_param"  "dose_records"
```

[`model_report()`](https://iamstein.github.io/synpmx/reference/model_report.md)
summarizes the whole fit in one printout:

``` r

model_report(fit)
#> Summarized from the source, not estimated
#>   cohort             32 patients in 1 arm(s)
#>                      all (32)
#>   dose changes       none
#>   visit grid         2 endpoint(s) at 16 nominal time(s), 22 slot(s) in all
#>   visit attendance   share of scheduled slots attended, one rate per
#>                      endpoint: cp 55%, pca 91%
#>   covariates         wt lognormal, age lognormal, sex categorical, drawn
#>                      once for the whole study, independently of the
#>                      profiles
#>   columns emitted    id, time, ntime, dv, amt, evid, dvid, wt, age, sex
#> 
#> Values at the lower limit of what was observed
#>     cp                 0.2, no assay limit: half the lowest value several patients reached, rounded down
#>     pca                5, no assay limit: half the lowest value several patients reached, rounded down
#> 
#> Each non-PK continuous endpoint, fitted as constant, linear, or exponential
#>   pca                exponential
#>                        plateau          30
#>                        baseline         100
#>                        rate             0.1
#>                        between-subject  0.1 (SD on the log baseline)
#>                        residual         additive 10
#>                        chosen on AIC from constant, linear, exponential
#> 
#> PK endpoint for the PopPK model
#>   cp                 inferred from the following data characteristics:
#>                        absent before the first dose
#>                        rises to one peak and comes back down within one
#>                          dose interval
#>   route              oral: the median profile rises to a peak at 9 before
#>                      declining, and 31% of subjects do too
#>   sampling           median 6 distinct times after a dose, 6 after the
#>                      peak: sparse within-subject sampling
#> 
#> The PopPK model
#> 
#> Estimated by nlmixr2
#>   structural model   2cmt_oral 
#>   selected by        two-compartment model passed acceptance checks
#>   fitted on          all 32 patients with a concentration
#>   fixed effects      cl 0.1, v 7, q 0.1, v2 2, ka 0.4 
#>   between-subject    cl 0.265, v 0.2, ka 0.548, q 0.0775, v2 0.632 (as SD on the log scale)
#>   residual error     proportional 0.2 
#>   time to fit        50.9 s
#>   whole call         51.1 s, against 50.9 s in the fitter
#> 
#> Privacy
#>   one patient's pull pass: largest: cp: ka between-subject SD, 9.84 points;
#>                      `model_privacy_checks()` has every check
```

## Every Number in the Fingerprint

`warfarin` has one PK endpoint, fitted a two-compartment oral model, and
one PD endpoint, which makes it a typical small study. Each row below is
one number in its fingerprint that was estimated from patients, and
`model` says which part of the fingerprint it belongs to:

- **PK**: the population PK model, a typical value and a between-subject
  variance per parameter and one residual error.
- **PD**: each PD endpoint’s time course, its shape parameters, the
  spread of the patients’ baselines and a residual error.
- **Covariates**: a trimmed mean and SD per continuous covariate, and
  one fewer share than levels per categorical one.
- **Dose changes**: the rates of reducing, skipping and stopping, pooled
  over the arms, and each reduced dose level as a fraction of the
  starting dose.
- **Missed visits**: one attendance rate per endpoint. An endpoint every
  patient attended at every visit has a rate of 1 and is not listed.
- **Discrete endpoints**: one fewer level share than levels for each
  binary or ordinal endpoint.
- **Arms and floors**: the size of each arm and of the cohort, the
  lowest value generated for each endpoint without a declared assay
  limit, and the reference weight where body-weight scaling was
  requested.

``` r

items <- fingerprint_items(release)
show(items, "Every estimated number in the warfarin fingerprint", paged = TRUE)
```

``` r

counts <- table(factor(items$model, levels = models))
c(counts, total = sum(counts), `planned cycles (protocol)` =
    planned_cycles(release))
#>                        PK                        PD                Covariates 
#>                        11                         5                         5 
#>              Dose changes             Missed visits        Discrete endpoints 
#>                         3                         2                         0 
#>           Arms and floors                     total planned cycles (protocol) 
#>                         4                        30                         1
```

The planned dose times and amounts are the protocol rather than
estimates, and are counted apart. So are the nominal visit grid and the
description of the table to generate, which are not counted.

### Across the Public Studies

The same count for the stored fit of each study in the [public-data
evaluation](https://iamstein.github.io/synpmx/articles/pmxmodel-public-data-examples.html).
A study with more endpoints, arms or dose levels has more rows, and a PD
shape fitted per arm under `pd_by_arm = TRUE` counts once per arm.

| study | patients | PK | PD | Covariates | Dose changes | Missed visits | Discrete endpoints | Arms and floors | total | planned cycles (protocol) |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| warfarin | 32 | 11 | 5 | 5 | 3 | 2 | 0 | 4 | 30 | 1 |
| theo_md | 12 | 11 | 0 | 2 | 3 | 1 | 0 | 3 | 20 | 7 |
| mad | 60 | 11 | 10 | 3 | 3 | 0 | 3 | 9 | 39 | 36 |
| case1_pkpd | 180 | 11 | 5 | 2 | 3 | 0 | 0 | 7 | 28 | 510 |
| wbcSim | 45 | 0 | 3 | 0 | 3 | 1 | 0 | 3 | 10 | 1 |
| mavoglurant | 120 | 9 | 0 | 8 | 3 | 1 | 0 | 3 | 24 | 1 |
| nimoData | 12 | 9 | 0 | 6 | 3 | 1 | 0 | 3 | 22 | 10 |
| pheno_sd | 59 | 9 | 0 | 4 | 6 | 1 | 0 | 3 | 23 | 14 |
| mixroute_sim | 90 | 12 | 0 | 2 | 3 | 0 | 0 | 4 | 21 | 9 |
| onc_sim | 200 | 11 | 4 | 5 | 4 | 2 | 0 | 4 | 30 | 842 |

## The PK model

### Structural model

``` r

release$structural
#> [1] "2cmt_oral"
model_candidates(fit)
#>       model converged accepted      aic seconds note
#> 1 2cmt_oral      TRUE     TRUE 922.2948  50.869
```

By default the fit tries two compartments and falls back to one only if
two fail the acceptance checks. `pk` names a model, or several to
compare by AIC.

### Parameters

Typical values on the natural scale, the between-subject covariance on
the log scale, and the residual error. Every PK and PD estimate in the
fingerprint is rounded to one significant figure.

``` r

model_parameters(fit)
#> $fixed
#>  cl   v   q  v2  ka 
#> 0.1 7.0 0.1 2.0 0.4 
#> 
#> $omega
#>      cl    v  ka     q  v2
#> cl 0.07 0.00 0.0 0.000 0.0
#> v  0.00 0.04 0.0 0.000 0.0
#> ka 0.00 0.00 0.3 0.000 0.0
#> q  0.00 0.00 0.0 0.006 0.0
#> v2 0.00 0.00 0.0 0.000 0.4
#> 
#> $residual
#> $residual$kind
#> [1] "proportional"
#> 
#> $residual$cv
#> [1] 0.2
```

No covariate is in the model: `covariate_effects = "none"` is the
default, so `release$covariate_effects` is empty and every synthetic
patient’s parameters are drawn from the covariance matrix alone.

``` r

release$covariate_effects
#> list()
```

The fingerprint holds no patient’s own random effects (empirical Bayes
estimates). They are read during estimation and dropped.

`warfarin` reports four concentrations at or below zero, which is what
an assay returns near its limit. They are fitted at half the smallest
positive value, and the residual error stays proportional.

### Censoring

``` r

fit$censoring
#> NULL
```

`warfarin` declares no assay limit, so this is empty. Where a study
declares one, its below-limit concentrations enter the population fit as
censored, and this table gives the share of each endpoint below the
limit. On
[`xgxr::case1_pkpd`](https://rdrr.io/pkg/xgxr/man/case1_pkpd.html) it is
46% of the concentrations, so much of what the PK model knows there is a
bound rather than a measurement.

## The PD time courses

Each continuous endpoint other than the concentration gets a time course
with no exposure term, chosen by AIC from a constant, a linear and an
exponential shape. It reproduces the average response over time, which
is enough for testing analysis code, and is not an exposure-response
model.

``` r

lapply(release$pd, function(shape) c(shape = shape$pd, shape$typical))
#> $pca
#>         shape       plateau      baseline          rate 
#> "exponential"          "30"         "100"         "0.1"
```

``` r

show(fit$pd[[1L]]$candidates, "The three shapes, compared on AIC")
```

## Covariates

Each baseline covariate has one distribution for the whole study, and
synthetic patients draw their covariates from it independently of
everything else. A continuous covariate is a mean and SD on the log
scale, taken after the highest and lowest 5% of patients are set aside.
A categorical one is a set of level frequencies, keeping only levels at
least `min_category_patients` patients (default 3) hold.

``` r

str(release$covariates, max.level = 2)
#> List of 3
#>  $ wt :List of 3
#>   ..$ kind   : chr "lognormal"
#>   ..$ meanlog: num 4.23
#>   ..$ sdlog  : num 0.18
#>  $ age:List of 3
#>   ..$ kind   : chr "lognormal"
#>   ..$ meanlog: num 3.37
#>   ..$ sdlog  : num 0.34
#>  $ sex:List of 3
#>   ..$ kind       : chr "categorical"
#>   ..$ levels     : chr [1:2] "female" "male"
#>   ..$ probability: num [1:2] 0.156 0.844
```

Because no covariate is in the PK model, a relationship the real study
had, such as heavier patients having a larger volume, is absent from the
synthetic data. The fit reports where that happens: the correlation
between each covariate and each patient’s random effects, as a
diagnostic that stays with the fit.

``` r

show(fit$correlations[order(-abs(fit$correlations$correlation)), ],
     "Each covariate against the random effects")
```

Here weight and sex both move with volume. If the synthetic data needs
that relationship, `covariate_effects = "auto"` adds allometric scaling
of clearance and volume on body weight.
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
keeps such relationships without modelling them, because each synthetic
patient’s covariates and profile come from the same real patients.

## The study design

The design half of the fingerprint is the same as the one
[`synpmx_pca_summarize()`](https://iamstein.github.io/synpmx/reference/synpmx_pca_summarize.md)
builds, from the same code.

### Arms

``` r

data.frame(arm = release$arms$arms, patients = as.integer(release$arms$sizes))
#>   arm patients
#> 1 all       32
```

### Dosing

Each arm has its planned schedule. The levels doses were reduced to, and
the rates of reducing, skipping and stopping, are pooled over the arms.

``` r

arm <- release$arms$arms[[1L]]
show(release$dosing[[arm]]$planned, "The planned schedule")
```

``` r

dosing <- release$dosing[[arm]]
c(levels = length(dosing$levels), reduction = dosing$reduction,
  interruption = dosing$interruption, discontinuation = dosing$discontinuation)
#>          levels       reduction    interruption discontinuation 
#>               1               0               0               0
```

`warfarin` is a single dose everyone received, so all three rates are
zero and the planned schedule is reproduced exactly. Where patients do
reduce or skip doses, the synthetic schedule is drawn before the profile
is simulated, so a synthetic patient who steps down gets the lower
exposure.

### Visits

The visits each arm has, and one attendance rate per endpoint.
Attendance is drawn visit by visit at generation. The rates, then the
number of visits, per endpoint:

``` r

rate <- do.call(pmax, unname(lapply(release$visits, function(v) v$probability)))
signif(tapply(rate, release$cells$endpoint, max), 3)
#>    cp   pca 
#> 0.551 0.906
table(release$cells$endpoint)
#> 
#>  cp pca 
#>  14   8
```

A visit is kept only where at least `min_arm_patients` patients were
observed there.

### Schema

The shape of the table to generate: its columns and their classes, the
compartment numbers, and any assay limit per endpoint.

``` r

release$schema$columns
#>  [1] "id"    "time"  "ntime" "dv"    "amt"   "evid"  "dvid"  "wt"    "age"  
#> [10] "sex"
unlist(release$schema$cmt_obs)
#> NULL
```

## Generating from it

``` r

synthetic <- synpmx_model_generate(release, n_subjects = 32, seed = 7)
c(rows = nrow(synthetic),
  subjects = length(unique(synthetic$id)),
  valid = validate_pmx(synthetic, roles)$valid)
#>     rows subjects    valid 
#>      525       32        1
```

## Where to go next

- [`vignette("pmxmodel-demo")`](https://iamstein.github.io/synpmx/articles/pmxmodel-demo.md):
  this generator run on a study end to end.
- [`vignette("pmxmodel-algorithm")`](https://iamstein.github.io/synpmx/articles/pmxmodel-algorithm.md):
  how each quantity above is estimated.
- [`vignette("pca-fingerprint")`](https://iamstein.github.io/synpmx/articles/pca-fingerprint.md):
  the same walk for the PCA generator, which shares the study design and
  describes the profiles differently.
- [`vignette("scorecard")`](https://iamstein.github.io/synpmx/articles/scorecard.md):
  the checks that compare generated data with its source.
- [Privacy protections in the PMX model
  generator](https://iamstein.github.io/synpmx/articles/pmxmodel-privacy.html):
  what protects the fingerprint, and what risk remains.
