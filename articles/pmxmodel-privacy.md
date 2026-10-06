# Privacy Protections in the PMX Model Generator

## Introduction

The pharmacometric (PMX) model generator,
[`synpmx_model()`](https://iamstein.github.io/synpmx/reference/synpmx_model.md),
fits a set of models to a study and then simulates from these models a
synthetic dataset, for the purpose of transferring the data outside the
GxP environment, but staying within the secure system ecosystem. The
models that are fit are:

1.  PopPK model fit to the PK data
2.  PD model fit to each PD observation type, pooled across cohorts
3.  An attendance model, to allow for simulation of missed visits
4.  A dosing model to simulate dose reductions and missed doses

This article clarifies what information from a study (via the model
parameters) can leave the GxP computing environment, and what risk
remains. While this approach offers many privacy protections which
reduce how much any one patient can impact the synthetic data, the
algorithm does not offer a formal differential privacy (DP) guarantee.

## Why a data generation method that does not formally guarantee DP was chosen

DP limits how much any one person’s data can change what is released,
whatever else an attacker knows \[1\]. A DP mechanism achieves this
limit by adding calibrated noise to every released number.

For a Phase 1 study (12-60 patients), a DP guarantee would require
adding so much noise to the synthetic dataset that it would no longer
meaningfully resemble the clinical data. Moreover, the DP guarantee
would be protecting the data from an attacker that the use does not
face. That use is synthetic data for developing analysis code, shared
with people who work under the study’s own access controls or an
agreement like them. The case rests on three points:

1.  **The release holds what the study’s own reports ultimately make
    public, more coarsely.** A population PK report gives each typical
    value, its between-subject variability and the residual error, with
    their standard errors, and regulators ask for exactly that \[11,
    12\]. The study’s main publication gives each arm’s baseline
    characteristics and how many patients discontinued \[13\]. This
    release holds the same kinds of number at one significant figure,
    pooled over the arms and the visits, without standard errors, and a
    few dozen at most ([Released Parameters at a
    Glance](#released-parameters-at-a-glance)).
2.  **Any attacker already holds the answer they are seeking.**
    Membership inference needs the candidate’s own measurements to
    assess whether that candidate was in the study; measurements such as
    concentrations, visit attendance, and dose changes. Outside the
    company sponsoring the study, this information exists only at the
    sites and in routine care of the patient. Thus whoever holds a
    candidate’s trial record knows the candidate was in the trial, and
    learning it again from the release discloses nothing.
3.  **DP at small cohort sizes replaces the study with noise**, as [the
    appendix](#dp-noise-small-study) computes, and for this use noise is
    worse than no data. The [calibrated generator’s
    evaluation](https://iamstein.github.io/synpmx/articles/calibrated-public-data-examples.html)
    finds a DP release at phase 1 sizes worse than generating from a
    public prior, which reads nothing.

## Privacy Assessment Summary

### What the Generator Releases

A population PK/PD model is a short list of numbers: a **typical value**
for each parameter, such as clearance; a **between-subject variance**
for each, which says how far patients spread around the typical value;
and a **residual error**, the scatter of single measurements around a
patient’s own curve.
[`synpmx_model()`](https://iamstein.github.io/synpmx/reference/synpmx_model.md)
estimates these with `nlmixr2`.

The generator also summarizes the study as it was carried out: the
number of patients in each treatment group (**arm**), each arm’s planned
dose schedule, how often doses were reduced, skipped or stopped, the
share of scheduled visits attended for each kind of measurement, the
distribution of each baseline covariate such as body weight, age or sex,
and how often each level of a yes/no or graded endpoint was recorded.
Each of these is one set of numbers for the whole study, pooled over the
arms and the visits.

These numbers are the **release**, and they are all that leaves the
environment that holds the study. Synthetic data is simulated from the
release and random numbers alone, so it discloses nothing the release
does not. The release is the study’s fingerprint, and [The PMX model
study
fingerprint](https://iamstein.github.io/synpmx/articles/pmxmodel-fingerprint.html)
walks through one field by field. [Released Parameters at a
Glance](#released-parameters-at-a-glance) counts the numbers in each
public study’s release.

### What Is Protected, and What Is Not

**✅ No patient record, identifier or per-patient estimate leaves the
study.** The release holds numbers about the population and about arms.
Each patient’s estimated random effects are used during estimation and
dropped from the release, and synthetic subject IDs are new numbers that
cannot collide with real ones.

**✅ No released number is one patient’s value, and every released
frequency rests on at least three patients.** A few settings the
generator needs (the lowest value it generates, the reference weight for
body-weight scaling, a planned dose where doses varied) are computed so
that no single patient decides them, rather than taken as half the
smallest measurement, the median weight or the most common dose. Each
frequency counts the patients who had something and those who did not:
an attendance rate of 30 in 32 also describes the 2 who missed. A
frequency is released only where both groups hold no patients or at
least three.

**✅ A patient who moves a released estimate by 15% or more is left out
of the estimates.** Estimation reads how far each patient moves each
model parameter estimate, leaves out anyone who moves one by 15% or more
(for a between-subject spread, by 0.15 in its standard deviation (SD) on
the log scale), and estimates again. A patient left out of one
estimation step is still used in the others. The reading is approximate
and does not cover every estimate; its limits are among the outstanding
risks.

**⚠️ Membership inference remains possible for an attacker who holds a
candidate’s own trial data.** Such an attacker can check whether the
released model fits the candidate slightly better than it would fit
someone outside the study. In simulation the advantage is small: given
one member and one non-member, the attacker picks the member 53% of the
time in a 60-patient study when its population parameters come from
another study (20% off), and 70% of the time in a 12-patient study when
it knows the population’s true parameters, where 50% is a coin flip.
[Appendix: Membership Inference Against the
Release](#membership-inference-against-the-release) provides further
details on these calculations.

**❌ No formal guarantee.** Nothing bounds what an attacker who knows
every other patient learns about the remaining one.

## The Privacy Protections

### 1. No Single Patient’s Value is Used

Every quantity that would naturally have been one patient’s value is
computed another way:

| Quantity | How it’s computed |
|----|----|
| Lowest value generated, used for LOQ where no assay limit is declared | half the lowest value at least three patients reached, rounded down to 1, 2 or 5 times a power of ten |
| Reference weight for body-weight scaling | the median to one significant figure |
| Planned dose, where patients’ doses differ | the arm’s mean to two significant figures |
| Continuous covariate | a trimmed mean and SD after the highest and lowest 5% of patients are set aside |

### 2. At Least Three Patients Behind Every Frequency

If only one or two patients in a study missed a visit or had a dose
reduced, any frequency describing that event reveals something about
them. The same is true when all but one or two patients had it: 30 of 32
patients attending also describes the 2 who did not. Statistical
disclosure control (SDC), the discipline that protects published
statistical tables, handles this with a threshold rule: publish a count
only where at least *k* people stand behind it \[2\]. Here *k* is 3, and
the rule is applied to both groups, the patients who had the event and
those who did not:

- an arm, visit, dose level or planned cycle is kept only where at least
  three patients reached it;
- an attendance rate is set to 1, and a dose-change rate to 0, where
  fewer than three patients missed a visit or changed a dose;
- a category level, of a covariate or of a graded endpoint, held by
  fewer than three patients is folded into the most common level.

The threshold stops a frequency from singling out one or two patients.
It does not stop membership inference, which adds up small signals
across many frequencies. Pooling is what weakens that: attendance is one
rate per endpoint, the dose-change rates are three, and a graded
endpoint is one set of level frequencies, each pooled over the arms and
the visits.

### 3. Covariate Summaries Without Their Extremes

A continuous covariate is summarized by a mean and SD taken after the
highest and lowest 5% of patients are set aside, with the SD corrected
for the trimming. No single extreme patient, such as one weighing 250
kg, can then move the summary. This is the trimmed mean of robust
statistics \[3\], and it is also how DP mean estimators bound one
person’s influence \[4\].

### 4. One Significant Figure

Every released PK and PD estimate is rounded to one significant figure;
covariate summaries are rounded to two, since one would put a mean
height of 172 cm at 200. One figure moves an estimate by 9% on average
and 33% at worst. In exchange, the released value is usually the same
whether a given patient is in the study or not, because a one-figure
step is about as wide as an estimate’s sampling variation. Rounding does
not hide an outlier’s effect on a variance, which protection 5 measures.
The full-precision estimates stay on the fit, outside the release.

### 5. How Far One Patient Moves Each Estimate

For each released estimate, estimation approximates how far it would
move if each patient in turn were left out, from that patient’s random
effects and without refitting. This is *local sensitivity*, the
study-specific version of DP’s sensitivity \[1\], and in regression the
same idea is Cook’s distance \[5\]. A move of 15 (percent for a typical
value; points of between-subject SD on the log scale, close to CV%, for
a spread) triggers a review, and 30 fails.

A patient who moves an estimate by 15 or more is left out of the PK fit,
the PD fits and the covariate summaries, and estimation runs again. The
patient stays in the dosing, visit and arm models, so cohort and arm
sizes do not change. Rounds repeat until nobody is flagged, because
leaving out the most influential patients makes the next ones the most
influential \[6\], and stop at a tenth of the cohort. Patients left out
are named on the console only; the fit records how many.
`drop_influential = FALSE` keeps every patient.

## The Five Checks

[`model_privacy_checks()`](https://iamstein.github.io/synpmx/reference/model_privacy_checks.md)
verifies the protections on a fit; a release carries the verdicts of P4
and P5.

| Check | Question | Passes when |
|----|----|----|
| P1 | No per-patient table is released | the release holds only population- and arm-level fields |
| P2 | No source identifier is released | no source subject ID appears in the release |
| P3 | No single patient’s value is released | protection 1 holds |
| P4 | Every released frequency rests on at least three patients | protection 2 holds for every frequency |
| P5 | No single patient moves a released estimate far | every estimate moves less than 15 (protection 5) |

``` r

checks <- model_privacy_checks(fit)
as.data.frame(checks)[, c("check", "verdict", "result")]
#>   check verdict
#> 1    P1    pass
#> 2    P2    pass
#> 3    P3    pass
#> 4    P4    pass
#> 5    P5    pass
#>                                                                                                                                                             result
#> 1                                                                                                                                                             none
#> 2                                                                                                                                                             none
#> 3 none: floors on the 1-2-5 series, covariates summarized without their extremes, no median, minimum or maximum stored, no factor level beyond what the arms carry
#> 4                                                                                                smallest group: 5, patients holding a categorical covariate level
#> 5                                                                                                                  largest: cp: ka between-subject SD, 9.84 points
```

The reading behind P5, one row per released estimate: the largest move
any one patient causes, the number of patients it was read over, the
shrinkage of the random effect it came from, and the verdict. No row is
about a patient.

``` r

digits3(attr(checks, "influence"))
#>        group   name            quantity released change   unit patients
#> 1         PK cp: cl       typical value      0.1   2.62      %       32
#> 2         PK cp: cl  between-subject SD    0.265   4.13 points       32
#> 3         PK  cp: v       typical value        7   1.01      %       32
#> 4         PK  cp: v  between-subject SD      0.2    1.9 points       32
#> 5         PK cp: ka       typical value      0.4    3.2      %       32
#> 6         PK cp: ka  between-subject SD    0.548   9.84 points       32
#> 7         PK  cp: q       typical value      0.1 0.0146      %       32
#> 8         PK  cp: q  between-subject SD   0.0775  0.829 points       32
#> 9         PK cp: v2       typical value        2   2.31      %       32
#> 10        PK cp: v2  between-subject SD    0.632   3.63 points       32
#> 11        PD    pca    typical baseline      100   2.11      %       32
#> 12        PD    pca  between-subject SD      0.1   5.86 points       32
#> 13 covariate     wt      geometric mean       69   1.22      %       32
#> 14 covariate     wt SD on the log scale     0.18   1.71 points       32
#> 15 covariate    age      geometric mean       29   2.05      %       32
#> 16 covariate    age SD on the log scale     0.34    2.3 points       32
#>    shrinkage verdict
#> 1     0.0153    pass
#> 2     0.0153    pass
#> 3      0.359    pass
#> 4      0.359    pass
#> 5      0.447    pass
#> 6      0.447    pass
#> 7      0.977    pass
#> 8      0.977    pass
#> 9      0.474    pass
#> 10     0.474    pass
#> 11        NA    pass
#> 12        NA    pass
#> 13        NA    pass
#> 14        NA    pass
#> 15        NA    pass
#> 16        NA    pass
```

The recount behind P4: for each kind of released frequency, the smallest
group of patients behind any one of them, and how many values the rules
changed to meet the floor.

``` r

attr(checks, "frequencies")[, c("quantity", "smallest", "threshold",
                                "adjusted")]
#>                                                   quantity smallest threshold
#> 1                                       patients in an arm       32         3
#> 2 patients on either side of an endpoint's attendance rate       21         3
#> 3              patients with the dose change behind a rate       NA         3
#> 4           patients holding a categorical covariate level        5         3
#> 5              patients holding a level of a factor column       32         3
#>   adjusted
#> 1        0
#> 2        0
#> 3        0
#> 4        0
#> 5        0
```

## Three Constructed Cases

Each case is `warfarin` with one or more patients changed, fitted once
and stored with the package. None is a real study.

``` r

warfarin <- as.data.frame(nlmixr2data::warfarin)
warfarin$ntime <- warfarin$time
roles <- pmx_roles(
  id = "id", time = "time", nominal_time = "ntime", dv = "dv", amt = "amt",
  evid = "evid", dvid = "dvid", covariates = c("wt", "age", "sex")
)
patients <- sort(unique(warfarin$id))

misdosed <- warfarin
on_cp <- misdosed$id == patients[[5L]] & misdosed$evid == 0 &
  misdosed$dvid == "cp"
misdosed$dv[on_cp] <- misdosed$dv[on_cp] * 1000

extreme <- warfarin
extreme$wt[extreme$id == patients[[3L]]] <- 250
extreme$age[extreme$id == patients[[8L]]] <- 95
extreme$sex <- as.character(extreme$sex)
extreme$sex[extreme$id == patients[[12L]]] <- "not recorded"
extreme$sex <- factor(extreme$sex)
```

``` r

misdosed_fit <- synpmx_model_estimate(misdosed, roles, seed = 1)
misdosed_kept_fit <- synpmx_model_estimate(misdosed, roles, seed = 1,
                                           drop_influential = FALSE)
extreme_fit <- synpmx_model_estimate(extreme, roles, seed = 1)
```

### A Patient Given a Thousand Times the Recorded Dose

One patient’s concentrations are a thousand times what the recorded dose
implies. Kept in with `drop_influential = FALSE`, the patient fails P5:

``` r

kept_checks <- model_privacy_checks(misdosed_kept_fit)
as.data.frame(kept_checks)[, c("check", "verdict", "result")]
#>   check verdict
#> 1    P1    pass
#> 2    P2    pass
#> 3    P3    pass
#> 4    P4    pass
#> 5    P5    FAIL
#>                                                                                                                                                             result
#> 1                                                                                                                                                             none
#> 2                                                                                                                                                             none
#> 3 none: floors on the 1-2-5 series, covariates summarized without their extremes, no median, minimum or maximum stored, no factor level beyond what the arms carry
#> 4                                                                                                smallest group: 5, patients holding a categorical covariate level
#> 5                                                                                                                  largest: cp: v2 between-subject SD, 70.7 points
```

By default the patient is left out, and the checks pass:

``` r

misdosed_fit$privacy$left_out$patients
#> [1] 1
as.data.frame(model_privacy_checks(misdosed_fit))[, c("check", "verdict",
                                                      "result")]
#>   check verdict
#> 1    P1    pass
#> 2    P2    pass
#> 3    P3    pass
#> 4    P4    pass
#> 5    P5    pass
#>                                                                                                                                                             result
#> 1                                                                                                                                                             none
#> 2                                                                                                                                                             none
#> 3 none: floors on the 1-2-5 series, covariates summarized without their extremes, no median, minimum or maximum stored, no factor level beyond what the arms carry
#> 4                                                                                                smallest group: 4, patients holding a categorical covariate level
#> 5                                                                 largest: cp: ka between-subject SD, 9.02 points; after leaving 1 patient(s) out of the estimates
```

Kept, the patient’s exposure becomes between-subject variability and the
synthetic patients spread wider than the source. Left out, they do not:

``` r

concentration <- function(data, label) {
  rows <- data$evid == 0 & data$dvid == "cp" & data$dv > 0
  data.frame(dataset = label, id = as.character(data$id[rows]),
             time = data$time[rows], dv = data$dv[rows])
}
kept_synthetic <- suppressWarnings(
  synpmx_model_generate(misdosed_kept_fit, seed = 11))
synthetic <- suppressWarnings(synpmx_model_generate(misdosed_fit, seed = 11))
plotted <- rbind(concentration(misdosed, "Source"),
                 concentration(kept_synthetic, "Synthetic, patient kept"),
                 concentration(synthetic, "Synthetic, patient left out"))
plotted$dataset <- factor(plotted$dataset, levels = unique(plotted$dataset))
ggplot2::ggplot(plotted, ggplot2::aes(time, dv, group = id,
                                      colour = dataset)) +
  ggplot2::geom_line(alpha = 0.4) +
  ggplot2::facet_wrap(~dataset) +
  ggplot2::scale_y_log10() +
  ggplot2::scale_colour_manual(values = c(unname(comparison_colours),
                                          "#7570B3")) +
  ggplot2::labs(x = "Time (hours)", y = "Warfarin concentration") +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "none")
```

![](pmxmodel-privacy_files/figure-html/misdosed-plot-1.png)

### Extreme Covariates

One patient weighs 250 kg and another is 95 years old. Trimming sets
both aside, so the covariate summaries barely change:

``` r

summary_of <- function(model) {
  data.frame(covariate = c("wt", "age"),
             geometric_mean = exp(c(model$wt$meanlog, model$age$meanlog)),
             sd_log = c(model$wt$sdlog, model$age$sdlog))
}
rbind(cbind(study = "warfarin", summary_of(fit$covariates)),
      cbind(study = "with the two extreme patients",
            summary_of(extreme_fit$covariates)))
#>                           study covariate geometric_mean sd_log
#> 1                      warfarin        wt             69   0.18
#> 2                      warfarin       age             29   0.34
#> 3 with the two extreme patients        wt             70   0.19
#> 4 with the two extreme patients       age             29   0.36
```

### A Category Nobody Else Holds

A third patient has a sex recorded as a level no other patient holds.
The level is dropped, so no synthetic patient can carry it:

``` r

extreme_fit$covariates$sex$levels
#> [1] "female" "male"
levels(extreme_fit$schema$prototypes$sex)
#> [1] "female" "male"
```

## The Checks Across the Public Studies

The checks on every stored fit of the [public-data
evaluation](https://iamstein.github.io/synpmx/articles/pmxmodel-public-data-examples.html).
`left_out` counts the patients left out for moving an estimate by 15 or
more, and `largest` is the largest move among those that remain; the
same patient moves a twelve-patient fit about five times as far as a
sixty-patient one.

``` r

files <- c(warfarin = "warfarin-model-fit.rds",
           theo_md = "theo-md-model-fit.rds",
           mad = "mad-model-fit.rds",
           case1_pkpd = "case1-pkpd-model-fit.rds",
           wbcSim = "wbcsim-model-fit.rds",
           mavoglurant = "mavoglurant-model-fit.rds",
           nimoData = "nimo-model-fit.rds",
           pheno_sd = "pheno-model-fit.rds",
           mixroute_sim = "mixroute-sim-model-fit.rds",
           onc_sim = "onc-sim-model-fit.rds")
rows <- lapply(names(files), function(study) {
  one <- stored_fit(files[[study]])
  if (!inherits(one, "pmx_fitted_model")) return(NULL)
  checks <- model_privacy_checks(one)
  verdicts <- stats::setNames(checks$verdict, checks$check)
  data.frame(study = study, patients = one$n_source,
             left_out = if (is.null(one$privacy$left_out)) 0L else
               one$privacy$left_out$patients, t(verdicts),
             largest = sub("; after leaving.*$", "",
                           checks$result[checks$check == "P5"]),
             check.names = FALSE)
})
do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
#>           study patients left_out   P1   P2   P3   P4   P5
#> 1      warfarin       32        0 pass pass pass pass pass
#> 2       theo_md       12        0 pass pass pass pass pass
#> 3           mad       60        1 pass pass pass pass pass
#> 4    case1_pkpd      180        0 pass pass pass pass pass
#> 5        wbcSim       45        0 pass pass pass pass pass
#> 6   mavoglurant      120        0 pass pass pass pass pass
#> 7      nimoData       12        0 pass pass pass pass pass
#> 8      pheno_sd       59        0 pass pass pass pass pass
#> 9  mixroute_sim       90        0 pass pass pass pass pass
#> 10      onc_sim      200        0 pass pass pass pass pass
#>                                                           largest
#> 1                 largest: cp: ka between-subject SD, 9.84 points
#> 2                           largest: DV: ka typical value, 9.92 %
#> 3        largest: PD - Continuous between-subject SD, 10.3 points
#> 4        largest: PD - Continuous between-subject SD, 6.44 points
#> 5                     largest: DV between-subject SD, 2.25 points
#> 6                 largest: DV: v2 between-subject SD, 2.09 points
#> 7                 largest: DV: v2 between-subject SD, 12.7 points
#> 8                  largest: APGR SD on the log scale, 4.32 points
#> 9                 largest: DV: ka between-subject SD, 4.25 points
#> 10 largest: Everolimus trough: ka between-subject SD, 3.39 points
```

## Checks on the Synthetic Table

[`synpmx_scorecard()`](https://iamstein.github.io/synpmx/reference/synpmx_scorecard.md)
compares the synthetic table with its source. Three of its rows ask
whether a synthetic patient reproduces a real one: B3 compares how close
synthetic patients sit to real ones with how close real patients sit to
each other, B4b asks whether any synthetic patient copies a real
patient’s values, and B5 asks whether a level fewer than three real
patients held reached the output.

``` r

synthetic <- synpmx_model_generate(fit, seed = 11)
card <- synpmx_scorecard(warfarin, synthetic, roles)
as.data.frame(card)[card$check %in% c("B3", "B4b", "B5"),
                    c("check", "question", "result", "verdict")]
#>    check                                         question
#> 11    B3    Adversarial accuracy inside its null interval
#> 13   B4b Generated DV vectors copying an exposed real one
#> 14    B5        Rare source levels copied into the output
#>                     result verdict
#> 11 0.656 in [0.344, 0.712]    pass
#> 13                       0    pass
#> 14          0 of 0 exposed    pass
```

## Outstanding Risks

The ways a release can still disclose something about a patient, most
serious first:

- **An attacker who holds a candidate’s own trial data can infer
  membership**, with the small advantage measured in the
  [appendix](#membership-inference-against-the-release). Only noise
  would reduce it.
- **An attacker who knows every other patient is not stopped.** They can
  compute the exact effect of the one patient they do not know. Only a
  DP mechanism stops them.
- **Releases compose.** Two releases from overlapping data, such as an
  interim and a final analysis, differ by exactly what the patients
  between them contributed. Nothing here accounts for that.
- **Three patients is the smallest threshold in common use.** A stricter
  rule needs `min_arm_patients` and `min_category_patients` raised.
- **The influence reading is approximate, and its thresholds are
  empirical.** They were set so that ordinary patients pass and a gross
  data error fails, not derived from a privacy target. Shrinkage makes
  it understate a patient’s effect, and the residual error, a parameter
  without between-subject variability, and PD shape parameters other
  than the baseline are not read.
- **The full fit and the estimation log are not safe to share.** The fit
  holds diagnostics the release does not, and the console names patients
  left out. Share the release or the synthetic data only.
- **Rarity in the world is not measured.** A covariate level many
  patients in this study hold can still identify someone if few people
  hold it.

## Released Parameters at a Glance

Every number in a release that was estimated from patients, counted
once, for the stored fit of each public study. A rate pooled over the
arms is listed in the release once per arm and counts once here, and the
planned dose schedule, which is the protocol, is counted apart.

``` r

# An arm's dosing model, or one per drug where the study declared which doses
# drive which endpoint.
per_drug <- function(entry) if (!is.null(entry$planned)) list(entry) else entry
parameter_count <- function(release) {
  arms <- names(release$arms$sizes)
  shape <- function(s) length(s$typical) + 2
  rates <- unique(unlist(lapply(release$dosing, function(entry) {
    lapply(per_drug(entry), function(d) {
      paste(d$reduction, d$interruption, d$discontinuation,
            paste(d$levels, collapse = "/"))
    })
  })))
  # One rate per endpoint, listed at each of its visits in each arm.
  attendance <- unique(unlist(lapply(release$visits, function(v) {
    p <- as.numeric(v$probability)
    paste(release$cells$endpoint, p)[p > 0 & p < 1]
  })))
  marginals <- unique(Filter(Negate(is.null),
                             unlist(release$discrete, recursive = FALSE)))
  counts <- c(
    `PK model` = sum(vapply(release$pk_models, function(m) {
      length(m$parameters$fixed) + nrow(m$parameters$omega) + 1
    }, numeric(1))),
    `PD time courses` = sum(vapply(release$pd, function(s) {
      if (length(s$arms)) sum(vapply(s$arms, shape, numeric(1))) else shape(s)
    }, numeric(1))),
    covariates = sum(vapply(release$covariates, function(c) {
      switch(c$kind, lognormal = 2, normal = 2,
             categorical = length(c$levels) - 1, 0)
    }, numeric(1))),
    `dose changes` = sum(vapply(strsplit(rates, " "), function(r) {
      3 + length(strsplit(r[[4]], "/")[[1]]) - 1
    }, numeric(1))),
    attendance = length(attendance),
    `discrete endpoints` = sum(vapply(marginals, function(m) {
      max(0, length(m$levels) - 1)
    }, numeric(1))),
    `counts and floors` = length(arms) + 1 +
      length(release$quantification_floor) +
      length(unique(unlist(lapply(release$covariate_effects,
                                  function(e) e$reference)))))
  c(counts, total = sum(counts),
    `planned cycles (protocol)` = sum(vapply(release$dosing, function(entry) {
      sum(vapply(per_drug(entry), function(d) nrow(d$planned), numeric(1)))
    }, numeric(1))))
}
studies <- c(warfarin = "warfarin-model-fit.rds",
             theo_md = "theo-md-model-fit.rds",
             mad = "mad-model-fit.rds",
             case1_pkpd = "case1-pkpd-model-fit.rds",
             wbcSim = "wbcsim-model-fit.rds",
             mavoglurant = "mavoglurant-model-fit.rds",
             nimoData = "nimo-model-fit.rds",
             pheno_sd = "pheno-model-fit.rds",
             mixroute_sim = "mixroute-sim-model-fit.rds",
             onc_sim = "onc-sim-model-fit.rds")
counted <- do.call(rbind, lapply(names(studies), function(study) {
  release <- model_release(stored_fit(studies[[study]]))
  data.frame(study = study, patients = release$n_source,
             t(parameter_count(release)), check.names = FALSE)
}))
knitr::kable(counted, row.names = FALSE)
```

| study | patients | PK model | PD time courses | covariates | dose changes | attendance | discrete endpoints | counts and floors | total | planned cycles (protocol) |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| warfarin | 32 | 11 | 5 | 5 | 3 | 2 | 0 | 4 | 30 | 1 |
| theo_md | 12 | 11 | 0 | 2 | 3 | 1 | 0 | 3 | 20 | 7 |
| mad | 60 | 11 | 10 | 3 | 3 | 0 | 3 | 9 | 39 | 36 |
| case1_pkpd | 180 | 11 | 5 | 2 | 3 | 0 | 0 | 7 | 28 | 510 |
| wbcSim | 45 | 0 | 4 | 0 | 3 | 1 | 0 | 3 | 11 | 1 |
| mavoglurant | 120 | 9 | 0 | 8 | 3 | 1 | 0 | 3 | 24 | 1 |
| nimoData | 12 | 9 | 0 | 6 | 3 | 1 | 0 | 3 | 22 | 10 |
| pheno_sd | 59 | 9 | 0 | 4 | 6 | 1 | 0 | 3 | 23 | 14 |
| mixroute_sim | 90 | 12 | 0 | 2 | 3 | 0 | 0 | 4 | 21 | 9 |
| onc_sim | 200 | 11 | 4 | 5 | 4 | 2 | 0 | 4 | 30 | 842 |

What each column holds:

- **PK model**: a typical value and a between-subject variance for each
  PK parameter, and one residual error: 11 numbers for the five
  parameters of a two-compartment oral model. A parameter fitted without
  between-subject variability, such as a bioavailability, adds its
  typical value alone.
- **PD time courses**: for each other continuous endpoint, one to three
  shape parameters, the spread of the patients’ baselines, and a
  residual error.
- **Covariates**: a trimmed mean and SD for each continuous covariate,
  and one fewer frequency than levels for each categorical one.
- **Dose changes**: the three pooled rates of reducing, skipping and
  stopping, and each level of the dose ladder below the full dose, a
  fraction of the starting dose to one significant figure.
- **Attendance**: one rate per endpoint, and none for an endpoint whose
  rate is 1, as where every patient was observed at every visit of it
  their arm had.
- **Discrete endpoints**: one fewer level frequency than levels, for
  each binary or ordinal endpoint.
- **Counts and floors**: the size of each arm and of the cohort, an
  emission floor for each endpoint without a declared assay limit, and
  the reference weight where body-weight scaling was requested.
- **Planned cycles**: each arm’s planned dose times and amounts, which
  are the protocol rather than estimates. The nominal visit grid and the
  description of the table to generate are the protocol too, and are not
  counted.

## Appendix: Membership Inference Against the Release

A membership-inference attack scores each candidate by how far the
release has moved, away from what the population would give, toward the
candidate’s own values. Members score higher on average than
non-members, and the attack’s strength is the area under the curve
(AUC): the probability that a randomly chosen member scores above a
randomly chosen non-member, where 0.5 is chance and 1 is certain
identification. The attacker holds the candidate’s own measurements and
reference values for the population, and does not know the other
patients. [Membership inference against the PMX model
release](https://iamstein.github.io/synpmx/articles/pmxmodel-privacy-simulation.html)
works through the calculation. Three findings:

- **Cohort size decides the population model’s exposure.** With exact
  reference values the AUC falls from 0.70 at 12 patients to 0.56 at 60;
  with reference values off by 20%, to 0.59 and 0.53.
- **One significant figure lowers it; two would not.** Rounding to one
  figure lowers the AUC by up to 0.05, while two figures move it by at
  most 0.001.
- **A table per visit would be the most exposed part of a release, so
  the release holds none.** A share of patients sampled at each of 20
  visits gives 0.83 with 10 patients behind each share and 0.62 with
  100, and the three-patient floor moves it by at most 0.01. One
  attendance rate per endpoint is about as exposed as one visit’s share:
  0.67 for five endpoints behind 10 patients.

An AUC averages over members, and the most exposed member is the one far
from the typical values, which is what protection 5 reads. The
membership-inference literature reports success against such outliers
separately for this reason \[7, 8\], and synthetic data from generative
models leaves them exposed \[9\].

## Appendix: The Noise a DP Guarantee Needs in a Small Study

Two routes give a population model a DP guarantee. The per-patient route
estimates each patient’s parameters from that patient’s own
measurements, clips them to a public range and releases their mean with
noise, which is how
[`synpmx_calibrated()`](https://iamstein.github.io/synpmx/reference/synpmx_calibrated.md)
releases its one number. Subsample-and-aggregate splits the patients
into disjoint groups of about twenty, fits the model in each, and
releases the mean of the groups’ clipped estimates with noise \[10,
14\].

With $`\varepsilon = 1`$ split evenly over *d* released numbers, each
clipped to a public range spanning sixteenfold on the log scale
(fourfold either side of a prior, as
[`synpmx_calibrated()`](https://iamstein.github.io/synpmx/reference/synpmx_calibrated.md)
uses), each number gets Laplace noise of scale $`b = dS/(n\varepsilon)`$
on the log scale by the per-patient route, and $`b = dS/(k\varepsilon)`$
with $`k = n/20`$ groups by subsample-and-aggregate. The table gives the
median fold error of one released number, $`e^{b \ln 2}`$, and the
noise’s spread against the uncertainty an estimate already has from
sampling, for a between-subject SD of 0.3:

``` r

span <- log(16)
dp_scale <- function(patients, released) released * span / patients
dp_fold <- function(patients, released) {
  exp(log(2) * dp_scale(patients, released))
}
dp_over_sampling <- function(patients, released) {
  sqrt(2) * dp_scale(patients, released) / (0.3 / sqrt(patients))
}
grid <- expand.grid(patients = c(30, 60, 200, 1000), released = c(11, 30))
fold_text <- function(x) {
  ifelse(x > 100, "over 100-fold", sprintf("%.2f-fold", x))
}
knitr::kable(data.frame(
  patients = grid$patients,
  released_numbers = grid$released,
  per_patient_median_error = fold_text(dp_fold(grid$patients, grid$released)),
  noise_over_sampling_error = sprintf("%.0f times", dp_over_sampling(
    grid$patients, grid$released)),
  subsample_median_error = fold_text(dp_fold(grid$patients / 20,
                                             grid$released))))
```

| patients | released_numbers | per_patient_median_error | noise_over_sampling_error | subsample_median_error |
|---:|---:|:---|:---|:---|
| 30 | 11 | 2.02-fold | 26 times | over 100-fold |
| 60 | 11 | 1.42-fold | 19 times | over 100-fold |
| 200 | 11 | 1.11-fold | 10 times | 8.28-fold |
| 1000 | 11 | 1.02-fold | 5 times | 1.53-fold |
| 30 | 30 | 6.83-fold | 72 times | over 100-fold |
| 60 | 30 | 2.61-fold | 51 times | over 100-fold |
| 200 | 30 | 1.33-fold | 28 times | over 100-fold |
| 1000 | 30 | 1.06-fold | 12 times | 3.17-fold |

Eleven numbers are a two-compartment population model alone, and thirty
a whole release. At 60 patients the per-patient route leaves half of the
population model’s numbers off by more than 1.42-fold, with noise 19
times the sampling error each estimate already carries, and a whole
release off by 2.6-fold. At 30 patients the population model alone is
off by 2.0-fold. Even at 1,000 patients the noise is 4.5 times the
sampling error. Subsample-and-aggregate, which keeps the fit as it is,
is unusable below about a thousand patients. Gaussian noise with tighter
accounting grows with the square root of the count rather than the
count, and does not change this picture at these cohort sizes. At the 12
to 200 patients this generator sees, a DP release would say less about
the study than a public prior does. Neither route is implemented beyond
the one number
[`synpmx_calibrated()`](https://iamstein.github.io/synpmx/reference/synpmx_calibrated.md)
releases.

## References

1.  Dwork C, McSherry F, Nissim K, Smith A. Calibrating noise to
    sensitivity in private data analysis. *Theory of Cryptography
    Conference (TCC).* 2006.
2.  Hundepool A, Domingo-Ferrer J, Franconi L, et al. *Statistical
    Disclosure Control.* Wiley; 2012.
3.  Huber PJ, Ronchetti EM. *Robust Statistics.* 2nd ed. Wiley; 2009.
4.  Bun M, Steinke T. Average-case averages: private algorithms for
    smooth sensitivity and mean estimation. *Advances in Neural
    Information Processing Systems (NeurIPS).* 2019.
5.  Cook RD. Detection of influential observation in linear regression.
    *Technometrics.* 1977;19(1):15–18.
6.  Carlini N, Jagielski M, Zhang C, Papernot N, Terzis A, Tramèr F. The
    privacy onion effect: memorization is relative. *Advances in Neural
    Information Processing Systems (NeurIPS).* 2022.
7.  Shokri R, Stronati M, Song C, Shmatikov V. Membership inference
    attacks against machine learning models. *IEEE Symposium on Security
    and Privacy.* 2017.
8.  Carlini N, Chien S, Nasr M, Song S, Terzis A, Tramèr F. Membership
    inference attacks from first principles. *IEEE Symposium on Security
    and Privacy.* 2022.
9.  Stadler T, Oprisanu B, Troncoso C. Synthetic data – anonymisation
    groundhog day. *USENIX Security Symposium.* 2022.
10. Nissim K, Raskhodnikova S, Smith A. Smooth sensitivity and sampling
    in private data analysis. *ACM Symposium on Theory of Computing
    (STOC).* 2007.
11. U.S. Food and Drug Administration. *Population Pharmacokinetics:
    Guidance for Industry.* 2022.
12. European Medicines Agency. *Guideline on Reporting the Results of
    Population Pharmacokinetic Analyses.* CHMP/EWP/185990/06. 2007.
13. Schulz KF, Altman DG, Moher D. CONSORT 2010 Statement: updated
    guidelines for reporting parallel group randomised trials. *BMJ.*
    2010;340:c332.
14. Smith A. Privacy-preserving statistical estimation with optimal
    convergence rates. *ACM Symposium on Theory of Computing (STOC).*
    2011.
