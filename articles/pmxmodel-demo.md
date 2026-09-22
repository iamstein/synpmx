# Demo: Using synpmx_model

Fit a population model to a study, look at what the fit carries, and
generate a synthetic dataset by simulating from it. Generation uses
fitted population parameters and study summaries: structural curves,
between-subject variability, residual error, endpoint frequencies, and
dosing and visit models per arm. It does not read the original patient
rows.

It makes no formal privacy claim, and **it is not for estimation** — the
fitted parameters exist to make simulated profiles look like the source
study, and the object prints that warning with itself. The full
specification is in
[`vignette("pmxmodel-algorithm")`](https://iamstein.github.io/synpmx/articles/pmxmodel-algorithm.md).

The dataset is [`xgxr::mad`](https://rdrr.io/pkg/xgxr/man/mad.html): 60
patients in six treatment arms, with a pharmacokinetic (PK)
concentration endpoint and continuous, count, ordinal and binary
pharmacodynamic (PD) endpoints. The two-compartment fit passed the
convergence and generation checks. The public-data survey also shows a
fit accepted with a false-convergence warning.

## Configuration

To run this on your own study, edit the two chunks below that:

1.  Read in the dataset.
2.  Define the column roles.

The rest of the code can be kept as is.

``` r

library(dplyr)
raw <- as.data.frame(get(utils::data(list = "mad", package = "xgxr")))
SEED <- 808
```

The column meanings. Columns not specified here are dropped from the
synthetic dataset.

``` r

roles <- pmx_roles(
  id           = "ID",
  time         = "TIME",
  dv           = "LIDV",
  mdv          = "MDV",
  amt          = "AMT",
  evid         = "EVID",
  cmt          = "CMT",
  dvid         = "NAME",        # which endpoint each observation row is
  nominal_time = "NOMTIME",     # the grid the visit model is built on
  strata       = c("TRTACT", "DOSE"), # treatment arms, fitted one at a time
  covariates   = c("WEIGHTB", "SEX") # drawn independently of the profiles
)
```

## Fitting

[`synpmx_model_estimate()`](https://iamstein.github.io/synpmx/reference/synpmx_model_estimate.md)
is the only stage that reads patient data, and the only one that needs
`nlmixr2`. It works out which endpoint is the drug, what design produced
it, and tries two compartments with checked fallback to one.

``` r

fit <- synpmx_model_estimate(raw, roles, seed = 1)
```

The chunk above is shown rather than run. Fitting compiles a model, so
this document reads a stored
[`synpmx_model_estimate()`](https://iamstein.github.io/synpmx/reference/synpmx_model_estimate.md)
fit and `R CMD check` never needs a compiler.

The fit report gives the elapsed time and the number of subjects used.
Each likelihood evaluation reads the dosing history, so repeated-dose
studies cost more than single-dose studies. `max_fit_subjects` caps the
population fit; the dosing, visit and covariate models still read every
subject.

`pk = "1cmt_oral"` requests one compartment directly. The requested fit
must pass the same convergence, parameter and generation checks, and is
never replaced by another model automatically.

A clearance of 6.25 L/h and a volume of 48.5 L. Whether those are the
right numbers for this compound is not the question the generator asks:
they exist to put the simulated profiles where the source’s are, and the
object prints that warning with itself.

## What the fit carries

Two halves.
[`model_report()`](https://iamstein.github.io/synpmx/reference/model_report.md)
separates them, because they answer to different things — one half is an
estimate with all an estimate’s caveats, the other is a summary of the
study’s apparatus. Everything the generator simulates from is in here,
which is why printing the fitted object shows the same account.

``` r

model_report(fit)
#> Summarized from the source, not estimated
#>   cohort             60 patients in 6 arm(s)
#>                      Placebo / 0 (10)
#>                      100 mg / 100 (10)
#>                      200 mg / 200 (10)
#>                      400 mg / 400 (10)
#>                      800 mg / 800 (10)
#>                      1600 mg / 1600 (10)
#>   dose changes       none
#>   visit grid         5 endpoint(s) at 28 nominal time(s), 66 slot(s) in all
#>   visit attendance   median 100% of an arm attends a slot (0% to 100%)
#>   covariates         WEIGHTB lognormal, SEX categorical, drawn once for the
#>                      whole study, independently of the profiles
#>   discrete endpoints PD - Binary, PD - Ordinal: drawn from each arm's
#>                      recorded frequencies at each visit, not simulated
#>   columns emitted    ID, TIME, NOMTIME, LIDV, AMT, EVID, CMT, NAME, MDV,
#>                      WEIGHTB, SEX, TRTACT, DOSE
#> 
#> Values at the lower limit of what was observed
#>     PK Concentration   0.025, half the smallest value seen, no assay limit
#>     PD - Continuous    0.0825, half the smallest value seen, no assay limit
#> 
#> Each non-PK continuous endpoint, fitted as constant, linear, or exponential
#>   PD - Continuous    exponential
#>                        plateau          31.47
#>                        baseline         1.637
#>                        rate             0.01344
#>                        between-subject  1.33 (SD on the log baseline)
#>                        residual         additive 8.13
#>                        chosen on AIC from constant, linear, exponential
#>   PD - Count         exponential
#>                        plateau          2.888
#>                        baseline         10.38
#>                        rate             0.01484
#>                        between-subject  0.248 (SD on the log baseline)
#>                        residual         additive 2.78
#>                        chosen on AIC from constant, linear, exponential
#> 
#> PK endpoint for the PopPK model
#>   PK Concentration   inferred from the following data characteristics:
#>                        absent before the first dose
#>                        dose-proportional: the peak scales with the dose
#>                        recorded in the compartment the doses go into
#>                        rises to one peak and comes back down within one
#>                          dose interval
#>   route              oral: the median profile rises to a peak at 2 before
#>                      declining, and 100% of subjects do too
#>   sampling           median 13 distinct times after a dose, 9 after the
#>                      peak: richer within-subject sampling
#> 
#> The PopPK model
#> 
#> Estimated by nlmixr2
#>   structural model   2cmt_oral 
#>   selected by        two-compartment model passed acceptance checks
#>   fitted on          all 50 patients with a concentration
#>   fixed effects      cl 6.254, v 48.55, q 4.871, v2 147, ka 1.208 
#>   between-subject    cl 0.42, v 0.445, ka 0.372, q 0.534, v2 0.424 (as SD on the log scale)
#>   residual error     proportional 0.372 
#>   time to fit        5 min 42 s
#>   whole call         5 min 43 s, against 5 min 42 s in the fitter
```

The concentration has a proportional residual error. The continuous and
count PD endpoints have fitted time courses with no exposure term.
Binary and ordinal endpoints are drawn from arm-by-visit frequencies
rather than from the compartment model. The signals table explains why
`PK Concentration` was identified as the concentration endpoint.

The correlation block is the other thing to read. A covariate that moves
with a random effect and is not in the model above is generated
independently of the profiles, so the synthetic data carries no
relationship between them.
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
keeps those relationships without modelling them.

## The candidates

``` r

show(model_candidates(fit), "Models attempted and acceptance results")
```

The table contains only models actually attempted. An accepted
two-compartment fit needs no one-compartment comparison. A rejected fit
retains its convergence status and rejection reason, and the fallback
has its own row.
[`model_report()`](https://iamstein.github.io/synpmx/reference/model_report.md)
explains the selection. A vector of model names in `pk` requests an
Akaike information criterion (AIC) comparison of those models after
acceptance checks.

A study containing only pharmacodynamic observations can instead declare
`endpoint_roles = c(pd = "response")`, using its endpoint name. That
skips compartment fitting and uses the existing PD time courses. A
positive baseline can also identify the study as PD-only automatically,
as the public-data survey demonstrates on `wbcSim`.

## Generating

[`synpmx_model_generate()`](https://iamstein.github.io/synpmx/reference/synpmx_model_generate.md)
reads no patient data. Its arguments are the fit and a subject count, so
everything about the source that reaches the output has already passed
through the fit.

``` r

synthetic <- synpmx_model_generate(fit, seed = SEED)
```

The generated table follows the declared schema: retained columns and
classes, the same compartment numbers, and new identifiers.

``` r

show(head(synthetic, 10), "The first ten rows")
```

## Plot synthetic data and original data

Compare each endpoint on the recorded time scale, with source and
synthetic subjects shown separately within each treatment arm.

``` r

library(ggplot2)
library(xgxr)
xgx_theme_set()
comparison_colours <- c(source = "#1B6CA8", synthetic = "#D95F02")

both <- rbind(
  cbind(raw[, names(synthetic)], DATA = "source"),
  cbind(synthetic, DATA = "synthetic")
)
obs <- both[both$EVID == 0 & !is.na(both$LIDV), ]
obs$TRTACT <- factor(obs$TRTACT,
                     levels = c("Placebo", "100 mg", "200 mg", "400 mg",
                                "800 mg", "1600 mg"))
```

``` r

ggplot(obs[obs$NAME == "PK Concentration" & obs$TIME <= 24, ],
       aes(TIME, LIDV, group = ID, colour = DATA)) +
  geom_line(alpha = 0.4) +
  facet_grid(DATA~TRTACT) +
  xgx_scale_y_log10() +
  xgx_scale_x_time_units("hours", breaks = seq(0, 24, by = 6)) +
  scale_colour_manual(values = comparison_colours) +
  labs(x = "Time (hours)", y = "PK concentration", colour = NULL) +
  theme(legend.position = "top") +
  ggtitle("Day 1 Conc. Profile")
```

![](pmxmodel-demo_files/figure-html/overlay-pk-1.png)

``` r

last_dose <- obs[obs$NAME == "PK Concentration" &
                   obs$TIME >= 120 & obs$TIME <= 144, ]
last_dose$TAD <- last_dose$TIME - 120

ggplot(last_dose, aes(TAD, LIDV, group = ID, colour = DATA)) +
  geom_line(alpha = 0.4) +
  facet_grid(DATA~TRTACT) +
  xgx_scale_y_log10() +
  xgx_scale_x_time_units("hours", breaks = seq(0, 24, by = 6)) +
  scale_colour_manual(values = comparison_colours) +
  labs(x = "Time after dose (hours)", y = "PK concentration", colour = NULL) +
  theme(legend.position = "top") +
  ggtitle("Last Dose (120 h) Conc. Profile")
```

![](pmxmodel-demo_files/figure-html/overlay-pk-last-1.png)

Every synthetic concentration profile uses the selected structural curve
with new parameter draws and residual error. Compare the rise, decline
and spread within each arm: pooled generation checks can miss
discrepancies confined to one dose group or part of the sampling window.

``` r

ggplot(obs[obs$NAME == "PD - Continuous", ],
       aes(TIME, LIDV, group = ID, colour = DATA)) +
  geom_line(alpha = 0.35) +
  xgx_scale_x_time_units(units_dataset = "hours", units_plot = "days") +
  facet_grid(DATA~TRTACT) +
  scale_colour_manual(values = comparison_colours) +
  labs(x = "Time (days)", y = "PD (continuous)", colour = NULL) +
  theme(legend.position = "top") +
  ggtitle("PD Response")
```

![](pmxmodel-demo_files/figure-html/overlay-pd-1.png)

The continuous PD endpoint has an exponential time course. Its
parameters are fitted to the pooled study, and subjects differ in their
baseline draws and residual error.

**What is missing is the dose ordering.** One shape is fitted to the
pooled observations and each subject gets a draw on its baseline, so the
arm a subject was assigned to does not reach its response.

``` r

late_pd <- function(data, label) {
  rows <- data$NAME == "PD - Continuous" & data$EVID == 0 &
    data$TIME > 144 & !is.na(data$LIDV)
  out <- aggregate(list(mean_pd = data$LIDV[rows]),
                   list(arm = data$TRTACT[rows]), function(x) round(mean(x)))
  stats::setNames(out, c("arm", paste0("mean_pd_", label)))
}
pd_table <- merge(late_pd(raw, "source"), late_pd(synthetic, "synthetic"),
                  by = "arm")
show(pd_table, "Mean PD response after 144 h, by arm")
```

This table compares late responses by arm. The continuous PD model
carries no exposure term, so any source relationship between dose and
response is not preserved by construction. A dataset whose purpose is
exposure-response needs that limitation considered separately from
numerical fit acceptance.

## Distributions of Synthetic and Original Data

One panel per endpoint and per baseline covariate, source against
synthetic.
[`compare_pmx_distributions_height()`](https://iamstein.github.io/synpmx/reference/compare_pmx_distributions_height.md)
sizes the figure, since it grows a row of panels at a time;
`output = "tables"` gives the numbers behind it.

``` r

compare_pmx_distributions(raw, synthetic, roles)
```

![](pmxmodel-demo_files/figure-html/distributions-1.png)

## Scoring it

[`synpmx_scorecard()`](https://iamstein.github.io/synpmx/reference/synpmx_scorecard.md)
does not know what produced the dataset it scores, so it reads this
output the same way it reads an AVATAR or PCA one.

``` r

card <- synpmx_scorecard(raw, synthetic, roles)
synpmx_scorecard_datatable(card)
```

Read the fidelity checks alongside the plots. A passing acceptance
screen excludes gross fitting and generation failures; it does not
guarantee that the simulated spread, attendance or endpoint
relationships match the source. Checks requiring an AVATAR run record
are not applicable to this generator.

[`vignette("scorecard")`](https://iamstein.github.io/synpmx/articles/scorecard.md)
documents what each row asks and what its pass criterion is.

## One call

[`synpmx_model()`](https://iamstein.github.io/synpmx/reference/synpmx_model.md)
is the two stages together, for when you do not need to look at the fit
first. It is on the result either way, as an attribute.

``` r

synthetic <- synpmx_model(raw, roles, seed = SEED)
attr(synthetic, "pmx_fitted_model")
```

## See also

- [`vignette("pmxmodel-algorithm")`](https://iamstein.github.io/synpmx/articles/pmxmodel-algorithm.md)
  — what each step does and why.
- [`vignette("pmxmodel-fingerprint")`](https://iamstein.github.io/synpmx/articles/pmxmodel-fingerprint.md)
  — every quantity the fit carries out of the study, in detail.
- [Evaluating the model generator on public
  data](https://iamstein.github.io/synpmx/articles/pmxmodel-public-data-examples.html)
  — the same generator over the public studies the other two surveys
  use.
- [`vignette("pca-demo")`](https://iamstein.github.io/synpmx/articles/pca-demo.md)
  — a worked study through the principal-component generator.
- [`vignette("avatar-demo")`](https://iamstein.github.io/synpmx/articles/avatar-demo.md)
  — a worked study through blending.
- [`vignette("scorecard")`](https://iamstein.github.io/synpmx/articles/scorecard.md)
  — the checks above, in detail.
