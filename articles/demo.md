# Demo: Using synpmx_avatar

Generate a synthetic dataset, look at it, and read the checks of the
data, sumamrized in a scorecard.

[`vignette("scorecard-synthetic-data-checks")`](https://iamstein.github.io/synpmx/articles/scorecard-synthetic-data-checks.md)
is the reference for better understanding the synthetic data checks.

The dataset is
[`xgxr::case1_pkpd`](https://rdrr.io/pkg/xgxr/man/case1_pkpd.html): 180
patients, six arms from placebo to 300 mg, with a pharmacokinetic (PK)
concentration endpoint and a continuous pharmacodynamic (PD) endpoint.

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
  nominal_time = "NOMTIME",     # the protocol grid
  strata       = c("TRTACT", "DOSE"), #treatment arms.
  covariates   = "WEIGHTB",     # blended like everything else
  keep         = "STUDY"        # carried through verbatim
)
```

## Generate Synthetic Data

``` r

synthetic <- synpmx_avatar(raw, roles, seed = SEED)
#> Warning: Synthetic generation used documented small-group/profile fallbacks:
#> - Endpoint `PD - Continuous` has generated observations at times no donor
#>   was measured at; the cohort's median trajectory was used for those rows.
#>   Expected wherever subjects were sampled on different schedules.
#> - Endpoint `PK Concentration` has generated observations at times no donor
#>   was measured at; the cohort's median trajectory was used for those rows.
#>   Expected wherever subjects were sampled on different schedules.
```

## Plot synthetic data and original data

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

![](demo_files/figure-html/overlay-pk-1.png)

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

![](demo_files/figure-html/overlay-pd-1.png)

## Summary Tables of Synthetic and Original Data

``` r

compare_pmx_distributions(raw, synthetic, roles)
```

| variable | dataset | n | n_subjects | mean | sd | min | q25 | median | q75 | max |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| PD - Continuous | source | 1620 | 180 | 135 | 238 | -621 | -21.6 | 123 | 283 | 936 |
| PD - Continuous | synthetic | 1620 | 180 | 130 | 153 | -297 | 27.9 | 114 | 221 | 657 |
| PK Concentration | source | 3600 | 150 | 0.36 | 0.737 | 0.05 | 0.05 | 0.0634 | 0.263 | 7 |
| PK Concentration | synthetic | 3600 | 150 | 0.313 | 0.649 | 0.05 | 0.05 | 0.0517 | 0.227 | 5.84 |

Endpoints (dependent variable on observation rows) {.table
style="width:100%;"}

| variable | dataset   |   n | mean |   sd |  min | q25 | median | q75 | max |
|:---------|:----------|----:|-----:|-----:|-----:|----:|-------:|----:|----:|
| WEIGHTB  | source    | 180 |  116 | 20.5 |   80 |  98 |    117 | 133 | 150 |
| WEIGHTB  | synthetic | 180 |  115 | 14.5 | 85.8 | 104 |    114 | 127 | 148 |

Continuous covariates (baseline, per patient) {.table}

Medians track closely on both endpoints. The standard deviations are
visibly smaller in the synthetic cohort — for the PD endpoint and for
baseline weight alike. This behavior is expected from the AVATAR
blending algorithm.

## The scorecard

The scorecard computes checks of the synthetic dataset. It is described
further in
[`vignette("scorecard-synthetic-data-checks")`](https://iamstein.github.io/synpmx/articles/scorecard-synthetic-data-checks.md).

``` r

scorecard <- synpmx_scorecard(raw, synthetic, roles)
if (requireNamespace("DT", quietly = TRUE)) {
  DT::datatable(as.data.frame(scorecard), options = list(paging = FALSE))
} else {
  scorecard
}
```

| check | question | reads | result | verdict | explore |
|:---|:---|:---|:---|:---|:---|
| A1 | Synthetic table is a legal PMX dataset | synthetic | TRUE | pass | validate_pmx(synthetic, roles) |
| A2 | Source is legal under the declared roles | source | TRUE | pass | validate_pmx(source, roles, strict = FALSE) |
| A3 | Every endpoint survived | both | 2 of 2 | pass | compare_pmx_distributions(source, synthetic, roles) |
| A4 | Cohort size survived | both | 180 -\> 180 | pass | pmx_masking_report(synthetic, source, roles, section = “anchors”) |
| A5a | Observations per patient | both | 30.7 -\> 30.7 | review | compare_pmx_distributions(source, synthetic, roles) |
| A5b | Doses per patient | both | 70.8 -\> 70.8 | review | pmx_masking_report(synthetic, source, roles, section = “dose_schedules”) |
| A6 | Discrete endpoints keeping their source scale | both | no discrete endpoint | pass | pmx_endpoint_types(source, roles) |
| B1a | Avatars with a visit set nobody else shares | run settings | 0 | pass | unmaskable_strata(source, roles) |
| B1b | Avatars with a dose schedule nobody else shares | run settings | 0 | pass | unmaskable_strata(source, roles) |
| B2 | Synthetic patients unusual within their stratum | synthetic | 0 of 180 | review | flag_identifiable_subjects(synthetic, roles) |
| B3 | Adversarial accuracy inside its null interval | both | 0.617 in \[0.415, 0.578\] | review | compare_pmx_proximity(source, synthetic, roles) |
| B4a | Generated time vectors copying an exposed real one | both | 0 | pass | skeleton_uniqueness(source, roles, coarsen_time = TRUE) |
| B4b | Generated DV vectors copying an exposed real one | both | 0 | pass | compare_pmx_proximity(source, synthetic, roles) |
| B5a | Patients holding the least-held categorical level | synthetic | TRTACT = 10 mg: 30 | pass | table(synthetic\$TRTACT) |
| B5b | Rare source levels copied into the output | both | 0 of 0 exposed | pass | compare_pmx_rare_levels(source, synthetic, roles) |
| C1 | Strata keeping their source size | both | 6 of 6 | pass | compare_pmx_strata_sizes(source, synthetic, roles) |
| C2 | Distinct dose-time schedules represented | run settings | 2 of 2 | pass | pmx_masking_report(synthetic, source, roles, section = “dose_schedules”) |
| D1 | Values landing in the same range | both | sd x0.64 on PD - Continuous (furthest of 3) | review | compare_pmx_distributions(source, synthetic, roles) |

*D1 reports numbers, not shapes. Plot source and synthetic on the same
axes – `DV` against time, and each covariate – with whatever you
normally use.*

**`review` is not a soft `pass`.** Five rows here have no pass mark
because no threshold would be honest. Doses per patient is the clearest:
on this study it is unchanged, but on a study with individualised dosing
it can halve.

**Every row names the call that explains it.** When a row reads oddly,
run what its `explore` column says: the numbers are a summary, and the
call behind each one is where the answer is.

## Where to go next

- [`vignette("scorecard-synthetic-data-checks")`](https://iamstein.github.io/synpmx/articles/scorecard-synthetic-data-checks.md)
  — every check, what it asks, and what passing means.
- [`vignette("public-data-examples")`](https://iamstein.github.io/synpmx/articles/public-data-examples.md) -
  AVATAR algorithm applied to 8 public datasets.
