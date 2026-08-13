# Demo: Using synpmx_avatar

Generate a synthetic dataset, look at it, and read the checks of the
data, sumamrized in a scorecard.

[`vignette("synthetic-data-checks")`](https://iamstein.github.io/synpmx/articles/synthetic-data-checks.md)
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
```

## Plot synthetic data and original data

``` r

library(ggplot2)
library(xgxr)
xgx_theme_set()

both <- rbind(
  cbind(raw[, names(synthetic)], DATA = "source"),
  cbind(synthetic, DATA = "synthetic")
)
obs <- both[both$EVID == 0 & !is.na(both$LIDV), ]
obs$TRTACT <- factor(obs$TRTACT,
                     levels = c("Placebo", "3 mg", "10 mg", "30 mg",
                                "100 mg", "300 mg"))
arms <- c(source = "red4", synthetic = "blue4")
```

``` r

ggplot(obs[obs$NAME == "PK Concentration" & obs$TIME <= 24, ],
       aes(TIME, LIDV, group = ID, colour = DATA)) +
  geom_line(alpha = 0.4) +
  facet_grid(DATA~TRTACT) +
  xgx_scale_y_log10() +
  xgx_scale_x_time_units("hours", breaks = seq(0, 24, by = 6)) +
  scale_colour_manual(values = arms) +
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
  scale_colour_manual(values = arms) +
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
| PD - Continuous | synthetic | 1620 | 180 | 130 | 161 | -297 | 17.4 | 116 | 231 | 657 |
| PK Concentration | source | 3600 | 150 | 0.36 | 0.737 | 0.05 | 0.05 | 0.0634 | 0.263 | 7 |
| PK Concentration | synthetic | 3600 | 150 | 0.316 | 0.65 | 0.05 | 0.05 | 0.0539 | 0.233 | 5.84 |

RESTRICTED – endpoints (dependent variable on observation rows) {.table
style="width:100%;"}

| variable | dataset   |   n | mean |   sd |  min | q25 | median | q75 | max |
|:---------|:----------|----:|-----:|-----:|-----:|----:|-------:|----:|----:|
| WEIGHTB  | source    | 180 |  116 | 20.5 |   80 |  98 |    117 | 133 | 150 |
| WEIGHTB  | synthetic | 180 |  115 | 14.5 | 85.8 | 104 |    114 | 127 | 148 |

RESTRICTED – continuous covariates (baseline, per patient) {.table}

Medians track closely on both endpoints. The standard deviations are
visibly smaller in the synthetic cohort — for the PD endpoint and for
baseline weight alike. This behavior is expected from the AVATAR
blending algorithm.

## The scorecard

The scorecard computes checks of the synthetic dataset. It is described
further in
[`vignette("synthetic-data-checks")`](https://iamstein.github.io/synpmx/articles/synthetic-data-checks.md).

``` r

pmx_scorecard(raw, synthetic, roles)
```

| check | question | reads | result | verdict |
|:---|:---|:---|:---|:---|
| A1 | Synthetic table is a legal PMX dataset | synthetic | TRUE | pass |
| A2 | Source is legal under the declared roles | source | TRUE | pass |
| A3 | Every endpoint survived | both | 2 of 2 | pass |
| A4 | Cohort size survived | both | 180 -\> 180 | pass |
| A5 | Observations per patient | both | 30.7 -\> 30.7 | review |
| A5 | Doses per patient | both | 70.8 -\> 70.8 | review |
| B1a | Avatars wearing one real patient’s visit set | run report | 0 | pass |
| B1b | Avatars wearing one real patient’s dose schedule | run report | 0 | pass |
| B2 | Synthetic patients unusual within their stratum | synthetic | 0 of 180 | review |
| B3 | Adversarial accuracy inside its null interval | both | 0.544 in \[0.422, 0.561\] | pass |
| B4a | Generated time vectors copying an exposed real one | both | 0 | pass |
| B4b | Generated DV vectors copying an exposed real one | both | 0 | pass |
| B5 | Synthetic patients holding the least-held level | synthetic | 30 | pass |
| C3 | Strata keeping their source size | both | 6 of 6 | pass |
| C4 | Dose regimens represented | both | 2 of 2 | pass |

**`review` is not a soft `pass`.** Three rows have no pass mark because
no threshold would be honest. Doses per patient is the clearest: on this
study it is unchanged, but on a study with individualised dosing it can
change.

## Where to go next

- [`vignette("synthetic-data-checks")`](https://iamstein.github.io/synpmx/articles/synthetic-data-checks.md)
  — every check, what it asks, and what passing means.
- [`vignette("avatar-evaluation-public-data")`](https://iamstein.github.io/synpmx/articles/avatar-evaluation-public-data.md) -
  AVATAR algorithm applied to 8 public datasets.
