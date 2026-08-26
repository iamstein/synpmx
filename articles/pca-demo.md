# Demo: Using synpmx_pca

Generate a synthetic dataset from a fitted model rather than from real
patients’ values, look at it, and read the checks of the data.

[`synpmx_pca_summarize()`](https://iamstein.github.io/synpmx/reference/synpmx_pca_summarize.md)
reduces the study to summaries and
[`synpmx_pca_generate()`](https://iamstein.github.io/synpmx/reference/synpmx_pca_generate.md)
builds a dataset from those summaries alone. No number a patient
measured reaches the output. What it carries out of the source is a
mean, a scale, a set of principal-component loadings, one mean score
vector per arm, and a residual covariance.
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
blends real values from neighbouring patients and is the subject of
[`vignette("avatar-demo")`](https://iamstein.github.io/synpmx/articles/avatar-demo.md),
which runs the same study through that algorithm.

It makes no formal privacy claim, and it is not for estimation. The full
specification is in
[`vignette("pca-algorithm")`](https://iamstein.github.io/synpmx/articles/pca-algorithm.md)
and the checks in
[`vignette("avatar-scorecard")`](https://iamstein.github.io/synpmx/articles/avatar-scorecard.md).

The dataset is
[`xgxr::case1_pkpd`](https://rdrr.io/pkg/xgxr/man/case1_pkpd.html): 180
patients, six arms from placebo to 300 mg, with a pharmacokinetic (PK)
concentration endpoint and a continuous pharmacodynamic (PD) endpoint.
It is the dataset
[`vignette("avatar-demo")`](https://iamstein.github.io/synpmx/articles/avatar-demo.md)
uses, so the two runs can be read against each other.

## Configuration

To run this on your own study, edit the two chunks below that:

1.  Read in the dataset.
2.  Define the column roles.

The rest of the code can be kept as is.

``` r

library(dplyr)
raw <- as.data.frame(get(utils::data(list = "case1_pkpd", package = "xgxr"))) |>
  mutate(CENS = ifelse(NAME == "PD - Continuous", 0, CENS))
SEED <- 808
```

The column meanings. Columns not specified here are dropped from the
synthetic dataset.

``` r

roles <- pmx_roles(
  id           = "ID",
  time         = "TIME",
  dv           = "LIDV",
  cens         = "CENS",
  amt          = "AMT",
  evid         = "EVID",
  cmt          = "CMT",
  dvid         = "NAME",        # which endpoint each observation row is
  nominal_time = "NOMTIME",     # required: the grid the model is built on
  strata       = c("TRTACT", "DOSE"), # treatment arms, and the score model's groups
  covariates   = "WEIGHTB",     # drawn jointly with the trajectories
  keep         = "STUDY"        # carried through verbatim
)
```

## Summarize the Study

Generation happens in two steps, and they are separate so that the first
can be looked at before the second runs.
[`synpmx_pca_summarize()`](https://iamstein.github.io/synpmx/reference/synpmx_pca_summarize.md)
is the only step that reads patient data.

``` r

model <- synpmx_pca_summarize(raw, roles)
model
#> A synpmx PCA model
#> 
#>   fitted on    180 patients, 6 arm(s): Placebo / 0 (30), 3 mg / 3 (30), 10 mg / 10 (30), 30 mg / 30 (30), 100 mg / 100 (30), 300 mg / 300 (30) 
#>   endpoints    PD - Continuous (9 visits modelled), PK Concentration (24 visits modelled) 
#>   covariates   WEIGHTB 
#>   components   16 (90% of variance) 
#>   dose term    factor 
#>   dosing       85 dose(s) per arm | shared by 100%-100% of each arm 
#> 
#> Generation reads this object and nothing else. To look inside it:
#>   pmx_pca_report(model)      what it read out of the source data
#>   pmx_pca_dosing(model)      the dose schedule each arm shares
#>   pmx_pca_visits(model)      the probability of a visit, per arm
#>   pmx_pca_components(model)  the loadings, over time
```

**`nominal_time` is required.** The grid it names is the axis every
feature sits on, and dose rows and observation rows are placed on it
together, so the timing comes through at its source values: doses per
patient, the last dose, the end of follow-up, the baseline sample before
the first dose.
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
will derive a grid when none is declared;
[`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md)
refuses, because choosing where a visit belongs is a statement about the
protocol rather than something to infer from recorded times.

`strata` does more work here than it does in
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md).
It names the arms, and the arms are the groups the score model is fitted
against: each one gets its own mean score vector and its own residual
covariance. An arm of fewer than three patients has no spread of its own
to model, so the function refuses rather than generating from one or two
people. `min_arm_patients` is where that floor is set.

## What the summary read out of the source

``` r

pmx_pca_report(model)
#> What the PCA fit read out of the source data
#> 
#>   subjects: 180  components retained: 16 
#> 
#>             quantity                                               what numbers
#>           visit grid               Nominal times modelled, per endpoint      33
#>      feature centers               Mean of each grid cell and covariate      34
#>       feature scales                     Standard deviation of the same      34
#>             loadings                 Component loadings on each feature     544
#>          score means                         Mean score vector, per arm      96
#>     score covariance             Residual covariance between components     256
#>  endpoint transforms                      Log or identity, per endpoint       2
#>         assay limits                   Censoring boundary, per endpoint       1
#>         dosing model             Dose times and amounts each arm shares    1020
#>          visit model Probability of a visit, per arm, endpoint and time     198
#>        arm constants         Strata and kept columns, one value per arm      18
#>  min_patients
#>           150
#>           150
#>           150
#>           180
#>            30
#>            30
#>           180
#>           180
#>            30
#>            30
#>            30
```

`min_patients` is the column to read. A grid cell or a covariate mean is
backed by most of the cohort; a per-arm quantity is backed by one arm.
Nothing in this table is a patient record, but a quantity standing on
few patients is where disclosure risk sits, not in the loadings.

## The dosing model

Dosing is not copied from anybody. Each arm contributes the schedule its
patients hold in common, and that schedule is what every generated
subject in the arm receives.

``` r

dosing <- pmx_pca_dosing(model)
aggregate(cbind(doses = dose) ~ arm + amt + share + distinct, dosing, length)
#>            arm amt share distinct doses
#> 1  Placebo / 0   0     1        1    85
#> 2     3 mg / 3   3     1        1    85
#> 3   10 mg / 10  10     1        1    85
#> 4   30 mg / 30  30     1        1    85
#> 5 100 mg / 100 100     1        1    85
#> 6 300 mg / 300 300     1        1    85
```

`share` is the fraction of the arm holding that schedule and `distinct`
is how many schedules the arm actually contained. Both at their best
here: one schedule per arm, held by everyone. A study recording actual
dose times rather than planned ones would show `distinct` close to the
arm size and `share` close to `1/n`, and the generated data would hold
one schedule where the source held thirty. The timing would still be
right; the variety would be gone.

## The visit model

``` r

library(ggplot2)
library(xgxr)
xgx_theme_set()
comparison_colours <- c(source = "#1B6CA8", synthetic = "#D95F02")
```

``` r

ggplot(pmx_pca_visits(model), aes(time, probability, colour = arm)) +
  geom_line() +
  facet_wrap(~endpoint, scales = "free_x") +
  labs(x = "Nominal time (hours)", y = "P(observation)", colour = NULL) +
  theme(legend.position = "top")
```

![](pca-demo_files/figure-html/visits-1.png)

Attendance is drawn per visit from these probabilities, so no real
patient’s set of attended visits is reused.

## Generate Synthetic Data

Everything above is what the synthetic data will be built from.
[`synpmx_pca_generate()`](https://iamstein.github.io/synpmx/reference/synpmx_pca_generate.md)
reads the model and a subject count, and no patient row is in scope
while it runs.

``` r

synthetic <- synpmx_pca_generate(model, seed = SEED)
```

`synpmx_pca(raw, roles, seed = SEED)` is the same two calls in one, for
when there is no reason to stop and look.

## Plot synthetic data and original data

``` r

both <- rbind(
  cbind(raw[, names(synthetic)], DATA = "source"),
  cbind(synthetic, DATA = "synthetic")
)
obs <- both[both$EVID == 0 & !is.na(both$LIDV), ]
obs$TRTACT <- factor(obs$TRTACT,
                     levels = c("Placebo", "3 mg", "10 mg", "30 mg",
                                "100 mg", "300 mg"))
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

![](pca-demo_files/figure-html/overlay-pk-1.png)

``` r

ggplot(obs[obs$NAME == "PD - Continuous", ],
       aes(TIME, LIDV, group = ID, colour = DATA)) +
  geom_line(alpha = 0.35) +
  xgx_scale_x_time_units(units_dataset = "hours", units_plot = "weeks") +
  facet_grid(DATA~TRTACT) +
  scale_colour_manual(values = comparison_colours) +
  labs(x = "Time (hours)", y = "PD (continuous)", colour = NULL) +
  theme(legend.position = "top") +
  ggtitle("PD Response")
```

![](pca-demo_files/figure-html/overlay-pd-1.png)

The synthetic profiles are smoother than the source ones. The model
holds a handful of components, and everything outside them, including
measurement noise visit to visit, is not reproduced.

The low-dose arms sit flat on the assay limit in both panels. That is
the censoring being put back: the drawn value is the latent one, and
where the source declared a lower limit of quantification the reported
value, `CENS` and `LIMIT` are rebuilt from it together. The basis itself
is fitted on values drawn inside the censoring region rather than on the
stack of identical limits, so it describes the patients rather than the
assay.

## Reading the components

A principal component is a direction in the space of whole trajectories,
and its loadings say how much each visit contributes to it. Plotted
against time rather than tabulated, the components become curves: a
component that is flat and positive across a profile is overall
magnitude, and one that crosses zero separates early visits from late
ones.

``` r

components <- pmx_pca_components(model)
ggplot(subset(components, component %in% c("PC1", "PC2", "PC3")),
       aes(time, loading, colour = component)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_line() + geom_point(size = 1) +
  facet_wrap(~endpoint, scales = "free_x") +
  labs(x = "Nominal time (hours)", y = "Loading", colour = NULL) +
  theme(legend.position = "top")
```

![](pca-demo_files/figure-html/components-1.png)

``` r

attr(pmx_pca_components(model), "variance_explained")
#>    component variance_explained cumulative
#> 1        PC1         0.50953784  0.5095378
#> 2        PC2         0.06276067  0.5722985
#> 3        PC3         0.03853618  0.6108347
#> 4        PC4         0.03471444  0.6455491
#> 5        PC5         0.03029726  0.6758464
#> 6        PC6         0.02788634  0.7037327
#> 7        PC7         0.02738581  0.7311185
#> 8        PC8         0.02552563  0.7566442
#> 9        PC9         0.02416383  0.7808080
#> 10      PC10         0.02273215  0.8035402
#> 11      PC11         0.02039767  0.8239378
#> 12      PC12         0.01777112  0.8417089
#> 13      PC13         0.01687197  0.8585809
#> 14      PC14         0.01509856  0.8736795
#> 15      PC15         0.01428485  0.8879643
#> 16      PC16         0.01274791  0.9007122
```

Describe what the curve shows rather than naming a mechanism for it. A
component that separates early from late looks like a difference in
decline rate, and calling it a half-life would be a claim the fit cannot
support: no clearance was estimated, and a score standard deviation is a
variance decomposition rather than a pharmacokinetic parameter. That is
the price of this basis, and
[`vignette("pca-algorithm")`](https://iamstein.github.io/synpmx/articles/pca-algorithm.md)
states what the alternative would cost.

## Distributions of Synthetic and Original Data

``` r

compare_pmx_distributions(raw, synthetic, roles)
```

![](pca-demo_files/figure-html/distributions-1.png)

## The scorecard

The scorecard computes checks of the synthetic dataset. It is described
further in
[`vignette("avatar-scorecard")`](https://iamstein.github.io/synpmx/articles/avatar-scorecard.md).

``` r

scorecard <- synpmx_scorecard(raw, synthetic, roles)
synpmx_scorecard_datatable(scorecard)
```

**Fifteen of the eighteen rows mean the same thing here as they do for
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md).**
They read the source and synthetic tables and compute the same way
whatever produced the data, which is what the scorecard was for.

**Three read `unavailable`, and that is not the same as unanswered.**
B1a, B1b and C2 read the AVATAR run record rather than the tables, and
this table does not carry one. For B1a and B1b the guarantee is
structural instead: no individual’s visit set and no individual’s dose
schedule exists in the model to be copied. B4a and B4b ask the copy
question directly, on the finished tables, and both read zero.

**C2 is the one to be careful about.** It asks how many distinct
dose-time schedules survived, and this is the row where the loss is
largest — one schedule per arm, against one per patient in a study
recording actual dose times. Read
[`pmx_pca_dosing()`](https://iamstein.github.io/synpmx/reference/pmx_pca_dosing.md)
in its place, where `distinct` and `share` say the same thing per arm.

## Where to go next

- [`vignette("pca-algorithm")`](https://iamstein.github.io/synpmx/articles/pca-algorithm.md)
  — every step of this algorithm, and what each one reads out of the
  source.
- [`vignette("avatar-demo")`](https://iamstein.github.io/synpmx/articles/avatar-demo.md)
  — the same study through
  [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md),
  which blends real values instead.
- [`vignette("avatar-scorecard")`](https://iamstein.github.io/synpmx/articles/avatar-scorecard.md)
  — every check, what it asks, and what passing means.
