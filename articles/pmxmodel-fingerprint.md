# The PMX model study fingerprint

Every generator in this package reduces a study to a **fingerprint** — a
set of summaries small enough to carry out of the environment that holds
the real data, and complete enough to build a synthetic study from. The
fingerprint is the whole of what a synthetic dataset descends from.
Nothing else about the source reaches it, which is why looking at the
fingerprint is how you decide whether the output can be trusted, and how
you satisfy yourself that no patient left the room.

[`synpmx_model_estimate()`](https://iamstein.github.io/synpmx/reference/synpmx_model_estimate.md)
returns this generator’s fingerprint, as a `pmx_fitted_model`.
[`synpmx_model_generate()`](https://iamstein.github.io/synpmx/reference/synpmx_model_generate.md)
reads that object and no patient row.

This vignette walks through it group by group. It is a reference: read
[`vignette("pmxmodel-demo")`](https://iamstein.github.io/synpmx/articles/pmxmodel-demo.md)
first for a run end to end, and
[`vignette("pmxmodel-algorithm")`](https://iamstein.github.io/synpmx/articles/pmxmodel-algorithm.md)
for how each quantity is arrived at and why.
[`vignette("pca-fingerprint")`](https://iamstein.github.io/synpmx/articles/pca-fingerprint.md)
is the same walk through the PCA generator’s fingerprint, and the two
are worth reading side by side.

## Two halves, and they are not alike

This fingerprint divides more sharply than the PCA one does.

**The estimated half** is a population model: a structural form, a
handful of fixed effects, a covariance matrix, a residual error. It is
perhaps twenty numbers describing the shape of every profile in the
study, which is a very small fingerprint for a very strong claim. It
also looks exactly like the output of a real population analysis, and it
is not one — the candidate set is small and the covariate model is
allometric scaling or nothing.

**The summarized half** is the study’s apparatus: who was in which arm,
when doses were planned and how patients departed from that plan, which
visits were attended, what the covariates looked like, where the assay
limit sat. None of this is estimated. It is exactly the apparatus
[`synpmx_pca_summarize()`](https://iamstein.github.io/synpmx/reference/synpmx_pca_summarize.md)
builds, from the same code.

The two halves carry different risks. Twenty fitted numbers are a
compact summary of 32 people; a per-arm dosing model with a cycle grid
can be much larger and is backed by one arm’s patients rather than the
whole cohort.

## The study

[`nlmixr2data::warfarin`](https://nlmixr2.github.io/nlmixr2data/reference/warfarin.html)
— thirty-two patients, a single oral dose, a concentration endpoint and
a pharmacodynamic one.

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

The chunk above is shown rather than run: fitting compiles a model, so
this document reads a stored fit and `R CMD check` never needs a
compiler.

## The inventory

[`model_report()`](https://iamstein.github.io/synpmx/reference/model_report.md)
is the top-level account, and it separates the two halves because they
answer to different things. Printing the fitted object shows the same
account, because everything on it is an input to the generation step.

``` r

model_report(fit)
#> The PopPK model
#> 
#> Estimated by nlmixr2
#>   structural model   1cmt_oral 
#>   fitted on          all 32 patients with a concentration
#>   fixed effects      cl 0.1362, v 8.175, ka 0.603 
#>   between-subject    cl 0.246, v 0.0854, ka 0.68 (as SD on the log scale)
#>   residual error     proportional 0.21 
#>   time to fit        9.4 s (9.6 s for the whole call)
#>                      nlmixr2: 1cmt_oral 9.4 s
#>                      least squares: pca 0.0 s
#>   covariate effects  cl ~ (wt/70)^0.75, v ~ (wt/70)^1.00 
#> 
#> Each other continuous endpoint, fitted as a shape in time
#>   pca                exponential
#>                        plateau          27.34
#>                        baseline         96.3
#>                        rate             0.09877
#>                        between-subject  0.146 (SD on the log baseline)
#>                        residual         additive 12.4
#>                        chosen on AIC from constant, linear, exponential
#> 
#> Values at the bottom of the scale
#>   No assay limit declared, so nothing is generated below:
#>     cp                 0.3
#>     pca                4.5
#> 
#> Summarized from the source, not estimated
#>   cohort             32 patients in 1 arm(s)
#>                      all (32)
#>   dose changes       none
#>   visit attendance   22 cell(s) of the visit grid (one endpoint at one
#>                      nominal time) over 2 endpoint(s), attended at median
#>                      95%, from 9% to 100%
#>   covariates         wt lognormal, age lognormal, sex categorical, drawn
#>                      per arm, independently of the profiles
#>   columns emitted    id, time, ntime, dv, amt, evid, dvid, wt, age, sex
#> 
#> PK endpoint for the PopPK model
#>   cp                 inferred: absent before the first dose; dose
#>                      proportionality not computable here; rises to one peak
#>                      and comes back down
#>   route              oral: the median profile rises to a peak at 9 before
#>                      declining, and 31% of subjects do too
#> 
#> Covariate against the individual random effects
#>  covariate parameter correlation
#>        age        cl        0.37
#>        sex        cl       -0.23
#>         wt        ka        0.21
#>         wt        cl       -0.13
#>         wt         v       -0.12
```

Every section below expands one part of it. The object is a plain list
and the names used are its own, so a quantity can be reached directly
when no accessor covers it.

``` r

names(fit)
#>  [1] "structural"           "candidates"           "parameters"          
#>  [4] "arms"                 "dosing"               "visits"              
#>  [7] "schema"               "roles"                "settings"            
#> [10] "n_source"             "cells"                "pd"                  
#> [13] "covariate_effects"    "covariates"           "discrete"            
#> [16] "design"               "correlations"         "censoring"           
#> [19] "quantification_floor" "timing"               "movement"            
#> [22] "fit_subjects"         "start_param"          "pk_models"           
#> [25] "endpoints"
```

## The settings that produced it

``` r

unlist(fit$settings)
#>      min_subjects  min_arm_patients     min_time_bins  max_fit_subjects 
#>              "20"               "3"               "6"              "60" 
#>        estimation covariate_effects             error 
#>           "focei"            "auto"            "prop"
c(patients = fit$n_source, arms = length(fit$arms$arms))
#> patients     arms 
#>       32        1
```

## How the concentration endpoint was decided

Not a released quantity, but the first thing to read: everything
downstream is a model *of* this endpoint, and if the wrong one was
chosen nothing else matters.

``` r

show(fit$endpoints$signals, "The four signals, per endpoint")
```

`proportional` reads `NA` here because `warfarin` gives every patient
the same dose, so there are no dose levels to compare — the
classification rests on the remaining signals. A column of `NA` under
`proportional` is the object telling you how little evidence it had, and
`endpoint_roles` is how you overrule it.

``` r

fit$design$reason
#> [1] "the median profile rises to a peak at 9 before declining, and 31% of subjects do too"
fit$endpoints$decided_by
#> [1] "inferred"
```

## The structural model

``` r

fit$structural
#> [1] "1cmt_oral"
model_candidates(fit)
#>       model converged     aic seconds note
#> 1 1cmt_oral      TRUE 895.892   9.425
```

One row, because the default fits one model. `pk` is what asks for more.

## The fixed effects

The typical parameters, on the natural scale. These are the numbers that
look most like a result and are least entitled to be read as one.

``` r

model_parameters(fit)$fixed
#>        cl         v        ka 
#> 0.1362061 8.1753822 0.6029980
```

## Between-subject variability

A covariance matrix on the log scale, one row per parameter carrying a
random effect. Generation draws each synthetic subject’s deviations from
this matrix.

``` r

model_parameters(fit)$omega
#>            cl           v       ka
#> cl 0.06071257 0.000000000 0.000000
#> v  0.00000000 0.007300333 0.000000
#> ka 0.00000000 0.000000000 0.462376
sqrt(diag(model_parameters(fit)$omega))  # as CV on the log scale
#>         cl          v         ka 
#> 0.24639922 0.08544199 0.67998235
```

**No individual estimates.** Empirical Bayes estimates are per-subject
quantities, and an object carrying them would be a description of each
real patient. They are not here. This matrix is what stands in for them,
and it is a statement about the population rather than about anybody in
it.

## The residual error

``` r

model_parameters(fit)$residual
#> $kind
#> [1] "proportional"
#> 
#> $cv
#> [1] 0.2102603
```

Additive here rather than proportional, because `warfarin` holds
concentrations recorded as zero and a proportional error on zero is
zero. The substitution is recorded on the object rather than being
silent.

## What the assay limit cost

The concentration is fitted with its below-limit rows censored, so each
contributes the probability of falling below the limit rather than a
value nobody measured. Every other part of the fingerprint is a
least-squares fit with no likelihood to put censoring in, and reads a
uniform draw inside the censoring region instead — an assumption, whose
weight is the share of each endpoint carrying it. The boundary is put
back when data is generated.

``` r

fit$censoring
#> NULL
```

`warfarin` declares no censoring, so nothing here is below a limit. On a
study where it is, this is the line to read before trusting the
concentrations:
[`xgxr::case1_pkpd`](https://rdrr.io/pkg/xgxr/man/case1_pkpd.html) has
46% of its concentrations below the limit and 95% of the lowest dose
arm, and there the concentration curve is drawn mostly from where the
assay stopped reporting.

## The covariate effects

``` r

fit$covariate_effects
#> $cl
#> $cl$covariate
#> [1] "wt"
#> 
#> $cl$reference
#> [1] 70
#> 
#> $cl$exponent
#> [1] 0.75
#> 
#> 
#> $v
#> $v$covariate
#> [1] "wt"
#> 
#> $v$reference
#> [1] 70
#> 
#> $v$exponent
#> [1] 1
```

Allometric scaling on clearance and volume, with the standard exponents.
It is asserted rather than tested —
[`vignette("pmxmodel-algorithm")`](https://iamstein.github.io/synpmx/articles/pmxmodel-algorithm.md)
says why — and `covariate_effects = "none"` removes it.

### What is *not* modelled shows up here

``` r

show(fit$correlations[order(-abs(fit$correlations$correlation)), ],
     "Each covariate against the individual random effects")
```

This is the most important table in the document, and it is a diagnostic
rather than a released quantity. A covariate that moves with a random
effect and is not in the model above is generated **independently** of
the profiles, so the synthetic data carries no relationship between the
two.
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
preserves those relationships without modelling them, because a blended
subject’s covariates and profile come from the same donors.

The correlations are computed from the individual random effects and
only the correlations are kept, so no per-subject quantity survives into
the object.

## The pharmacodynamic shapes

``` r

lapply(fit$pd, function(shape) c(shape = shape$pd, round(shape$typical, 3)))
#> $pca
#>         shape       plateau      baseline          rate 
#> "exponential"      "27.335"      "96.296"       "0.099"
```

``` r

show(fit$pd[[1L]]$candidates, "The three shapes, compared on AIC")
```

A time course with no exposure term. A PD endpoint driven by
concentration is reproduced as a curve that resembles the average
subject’s response, which is adequate for exercising longitudinal code
and is not an exposure-response model.

## The covariate distributions

One distribution per covariate per arm, drawn from independently at
generation.

``` r

str(fit$covariates, max.level = 3)
#> List of 1
#>  $ all:List of 3
#>   ..$ wt :List of 4
#>   .. ..$ kind   : chr "lognormal"
#>   .. ..$ meanlog: num 4.23
#>   .. ..$ sdlog  : num 0.19
#>   .. ..$ median : num 71.7
#>   ..$ age:List of 4
#>   .. ..$ kind   : chr "lognormal"
#>   .. ..$ meanlog: num 3.39
#>   .. ..$ sdlog  : num 0.304
#>   .. ..$ median : num 27.5
#>   ..$ sex:List of 3
#>   .. ..$ kind       : chr "categorical"
#>   .. ..$ levels     : chr [1:2] "female" "male"
#>   .. ..$ probability: num [1:2] 0.156 0.844
```

## The arms

``` r

data.frame(arm = fit$arms$arms, patients = as.integer(fit$arms$sizes))
#>   arm patients
#> 1 all       32
```

## The dosing model

Per arm: a planned schedule, a dose ladder, and three discrete-time
hazards. This is the apparatus
[`synpmx_pca_summarize()`](https://iamstein.github.io/synpmx/reference/synpmx_pca_summarize.md)
builds, unchanged.

``` r

arm <- fit$arms$arms[[1L]]
show(fit$dosing[[arm]]$planned, "The planned schedule")
```

``` r

dosing <- fit$dosing[[fit$arms$arms[[1L]]]]
c(levels = length(dosing$levels), reduction = dosing$reduction,
  interruption = dosing$interruption, discontinuation = dosing$discontinuation)
#>          levels       reduction    interruption discontinuation 
#>               1               0               0               0
```

`warfarin` is a single dose that everybody received, so the ladder has
one level and all three rates are zero. The model then reproduces the
planned schedule exactly, which is the right answer and is what those
zeros say. On a study where patients reduce, skip cycles or come off
treatment, these are the numbers that carry it — and because the
schedule is drawn *before* the profile is computed from it, a synthetic
subject who steps down has the lower exposure that implies.

## The visit model

Per endpoint and per retained nominal time, the fraction of the arm
holding an observation there. Attendance is drawn per visit at
generation.

``` r

cells <- fit$cells
cells$probability <- as.numeric(fit$visits[[fit$arms$arms[[1L]]]]$probability)
show(cells[, c("endpoint", "time", "probability")],
     "Attendance, per grid cell", paged = TRUE)
```

``` r

library(ggplot2)
library(xgxr)
xgx_theme_set()

ggplot(cells, aes(time, probability, colour = endpoint)) +
  geom_line() + geom_point(size = 1.5) +
  ylim(0, 1) +
  xgx_scale_x_time_units("hours", breaks = seq(0, 120, by = 24)) +
  labs(x = "Nominal time (hours)", y = "Fraction of the arm observed",
       colour = NULL) +
  theme(legend.position = "top")
```

![](pmxmodel-fingerprint_files/figure-html/visits-plot-1.png)

A cell is kept only where at least `min_arm_patients` distinct patients
hold an observation there. A nominal time one patient attended is that
patient, and generating from it would put them back.

## The schema

What the generated table has to look like to be the same study: the
columns and their classes, the compartment numbers, the assay limit per
endpoint, and how identifiers are written.

``` r

fit$schema$columns
#>  [1] "id"    "time"  "ntime" "dv"    "amt"   "evid"  "dvid"  "wt"    "age"  
#> [10] "sex"
fit$schema$cmt_dose
#> NULL
unlist(fit$schema$cmt_obs)
#> NULL
fit$schema$censoring
#> $cp
#> NULL
#> 
#> $pca
#> NULL
```

## Generating from it

Everything above, and nothing else, is what the next line reads.

``` r

synthetic <- synpmx_model_generate(fit, n_subjects = 32, seed = 7)
c(rows = nrow(synthetic),
  subjects = length(unique(synthetic$id)),
  valid = validate_pmx(synthetic, roles)$valid)
#>     rows subjects    valid 
#>      513       32        1
```

## Where to go next

- [`vignette("pmxmodel-demo")`](https://iamstein.github.io/synpmx/articles/pmxmodel-demo.md)
  — this generator run on a study end to end.
- [`vignette("pmxmodel-algorithm")`](https://iamstein.github.io/synpmx/articles/pmxmodel-algorithm.md)
  — how each quantity above is arrived at.
- [`vignette("pca-fingerprint")`](https://iamstein.github.io/synpmx/articles/pca-fingerprint.md)
  — the same walk for the PCA generator, which carries the same
  apparatus and a completely different description of the profiles.
- [`vignette("avatar-scorecard")`](https://iamstein.github.io/synpmx/articles/avatar-scorecard.md)
  — the checks that read a generated dataset against its source.
