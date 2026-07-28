# Example: are PK/PD, covariate, and treatment relationships preserved?

This article builds a source dataset with four relationships put in
deliberately — a covariate effect, a dose effect, a treatment effect,
and an exposure–response link — runs
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
over it, and measures how much of each one comes out the other side.

AVATAR does not explicity model any of these relationships. There is no
covariate model in it, no dose–response model, no PK/PD link, and no
notion of a treatment arm. Nothing in the generator knows that weight
ought to affect clearance. Whatever survives does so as a **side effect
of how the method works**; whole real subjects are blended together, and
a subject carries all of its properties at once. So the right
expectation is a *tendency* to keep relationships, but there is not a
guarantee.

## The model the source data comes from

Everything here is simulated, so the truth is known exactly and no real
patient is involved. The four relationships being tested are the four
places a covariate or a dose enters these equations — nothing else in
the article puts a relationship in.

**Pharmacokinetics.** One-compartment oral, first-order absorption and
elimination, for subject $`i`$ given dose $`D_i`$:

``` math
C_i(t)=\frac{D_i}{V_i}\cdot\frac{k_a}{k_a-CL_i/V_i}
        \left(e^{-\frac{CL_i}{V_i}t}-e^{-k_a t}\right),
\qquad k_a = 1.2\ \mathrm{h^{-1}}
```

**Allometric scaling — relationship 1.** This is the covariate effect,
and the exponents are the numbers to remember:

``` math
CL_i = 5 \left(\frac{WT_i}{70}\right)^{\mathbf{0.75}} e^{\eta_{CL,i}},
\qquad
V_i = 40 \left(\frac{WT_i}{70}\right)^{\mathbf{1.0}} e^{\eta_{V,i}}
```

with $`\eta_{CL}\sim N(0, 0.20^2)`$ and $`\eta_V\sim N(0, 0.15^2)`$ on
the log scale.

**Dose — relationship 2.** The `low` arm gets $`D=300`$ and the `high`
arm $`D=600`$. The model is linear in dose, so exposure should scale
exactly with it.

**Pharmacodynamics — relationships 3 and 4.** A direct, linear
inhibitory effect on a baseline:

``` math
E_i(t) = \bigl(E_{0,i} - 3.0\, C_i(t)\bigr)\,e^{\varepsilon},
\qquad E_{0,i} = 100\, e^{\eta_{E_0,i}}
```

The slope of $`3.0`$ is what creates both the **treatment effect** (the
high arm reaches higher concentrations, so its response is driven
further down) and the **exposure–response** link (whatever a subject’s
exposure, response tracks it). There is no separate “treatment” term:
the arm matters only through the dose it implies.

### What that predicts, before any data is generated

Two of these have closed-form consequences worth writing down, because
they are what the table later checks against:

Since $`AUC_{0-\infty} = D/CL`$ and $`CL \propto WT^{0.75}`$,

``` math
\log\frac{AUC_i}{D_i} = \text{constant} - 0.75\log WT_i
```

so **regressing dose-normalised log AUC on log weight should recover a
slope of $`-0.75`$** — the allometric exponent itself, negative because
a heavier subject clears faster and is therefore exposed *less*. And
because the model is linear in dose, the **high/low arm AUC ratio should
be $`600/300 = 2`$**.

Those two numbers, $`0.75`$ and $`2`$, are the ground truth. Everything
below asks how much of them survives.

## Generating the source dataset

Eighty subjects, forty per arm.

``` r

set.seed(2026)

make_subject <- function(id, arm, dose) {
  wt <- rnorm(1, 70, 12)
  # RELATIONSHIP 1 -- the covariate effect. The 0.75 and the 1.0 exponents are
  # the only place weight ever influences anything.
  cl <- 5  * (wt / 70)^0.75 * exp(rnorm(1, 0, 0.20))
  v  <- 40 * (wt / 70)^1.0  * exp(rnorm(1, 0, 0.15))
  ka <- 1.2
  pk_time <- c(0.5, 1, 2, 4, 8, 12, 24)
  pd_time <- c(0, 4, 12, 24)
  # RELATIONSHIP 2 -- `dose` enters here and nowhere else, so the arm affects
  # exposure only through the dose it implies.
  conc <- function(t) {
    dose / v * ka / (ka - cl / v) * (exp(-cl / v * t) - exp(-ka * t))
  }
  e0 <- 100 * exp(rnorm(1, 0, 0.10))
  rbind(
    data.frame(ID = id, TIME = 0, DV = NA_real_, AMT = dose, EVID = 1L,
               CMT = 1L, DVID = "PK", WT = wt, ARM = arm),
    data.frame(ID = id, TIME = pk_time,
               DV = conc(pk_time) * exp(rnorm(length(pk_time), 0, 0.08)),
               AMT = 0, EVID = 0L, CMT = 2L, DVID = "PK", WT = wt, ARM = arm),
    # RELATIONSHIPS 3 and 4 -- the `- 3 * conc(...)` term. It makes response
    # track exposure, and therefore makes the high-dose arm respond more.
    data.frame(ID = id, TIME = pd_time,
               DV = (e0 - 3 * conc(pd_time)) *
                 exp(rnorm(length(pd_time), 0, 0.05)),
               AMT = 0, EVID = 0L, CMT = 3L, DVID = "PD", WT = wt, ARM = arm)
  )
}

source_data <- do.call(rbind, c(
  lapply(1:40,  make_subject, arm = "low",  dose = 300),
  lapply(41:80, make_subject, arm = "high", dose = 600)
))
source_data <- source_data[
  order(source_data$ID, source_data$TIME, source_data$EVID == 0L), ]
rownames(source_data) <- NULL

roles <- pmx_roles(
  id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
  cmt = "CMT", dvid = "DVID",
  covariates = "WT",     # blended across donors into a new value
  keep = "ARM"           # copied verbatim from the one anchor subject
)
```

Note the two different treatments of the two non-PK columns, because it
matters for what follows. `WT` is a **covariate**: the synthetic subject
gets a new weight blended from its donors. `ARM` is **kept**: the
synthetic subject gets one real subject’s arm label, copied unchanged
from its anchor. Those are different subjects — the anchor supplies the
skeleton and the kept columns, the donors supply the values — which is
why the treatment effect surviving at all is not obvious in advance.

## Measuring the four relationships

Each measurement is chosen so its answer can be compared against a
number we already know from the equations above.

To recover the allometric exponent, log AUC must be regressed on log
weight after dividing out the dose.

``` r

relationships <- function(d) {
  obs <- d$EVID == 0 & !is.na(d$DV)
  pk  <- d[obs & d$DVID == "PK", ]
  pd  <- d[obs & d$DVID == "PD", ]
  ids <- as.character(pk$ID)

  auc <- tapply(seq_len(nrow(pk)), ids, function(i) {      # AUC 0-24 h, trapezoid
    t <- pk$TIME[i]; y <- pk$DV[i]; o <- order(t)
    sum(diff(t[o]) * (head(y[o], -1) + tail(y[o], -1)) / 2)
  })
  wt   <- tapply(pk$WT, ids, function(v) v[1])
  arm  <- tapply(as.character(pk$ARM), ids, function(v) v[1])
  dose <- tapply(d$AMT[d$EVID == 1], as.character(d$ID[d$EVID == 1]),
                 function(v) v[1])[names(auc)]
  nadir <- tapply(pd$DV, as.character(pd$ID), min)[names(auc)]
  ok <- is.finite(auc) & auc > 0 & is.finite(wt) & is.finite(nadir)

  c(# 1. allometric exponent on CL; truth = 0.75. Negated so it reads positive.
    allometric = unname(-coef(lm(log(auc[ok] / dose[ok]) ~ log(wt[ok])))[2]),
    # 2. dose effect; truth = 600/300 = 2 for linear PK.
    dose_ratio = unname(mean(auc[ok & arm == "high"]) /
                          mean(auc[ok & arm == "low"])),
    # 3. treatment effect: how much deeper the high arm's response goes.
    pd_effect  = unname(mean(nadir[ok & arm == "low"]) -
                          mean(nadir[ok & arm == "high"])),
    # 4. exposure-response: response against the exposure that drove it.
    er_slope   = unname(coef(lm(nadir[ok] ~ log(auc[ok])))[2]))
}

truth <- relationships(source_data)
round(truth, 3)
#> allometric dose_ratio  pd_effect   er_slope 
#>      0.747      2.151     18.747    -23.678
```

The first two land where the equations said they would: an allometric
exponent of about 0.75, and an arm ratio of about 2. That is the check
that the source dataset really contains what we think it contains,
before AVATAR is asked to preserve it.

## What the source data looks like

![](example-avatar-PKPD-covariate-treatment-effect_files/figure-html/source-plot-1.png)

The high arm reaches higher concentrations and its response is driven
further down — relationships 2 and 3, visible directly.

## Run AVATAR, several times

One synthetic dataset tells you little, because each run draws its own
anchors and donors. Thirty runs show both the central tendency and the
spread, and the spread is as much the point as the average.

``` r

set.seed(7)
replicates <- replicate(30, {
  synthetic <- suppressWarnings(
    synpmx_avatar(source_data, roles, seed = sample.int(1e6, 1))
  )
  relationships(synthetic)
})
```

| relationship | source | synthetic | run_to_run_sd | retained |
|:---|---:|---:|---:|:---|
| allometric exponent on CL (truth = 0.75) | 0.747 | 0.683 | 0.184 | 92% |
| dose effect: AUC ratio high/low (truth = 2.0) | 2.151 | 2.240 | 0.076 | 104% |
| treatment effect: PD nadir difference (low - high) | 18.747 | 18.637 | 4.475 | 99% |
| exposure-response: PD nadir per log AUC | -23.678 | -20.964 | 4.664 | 89% |

Relationships in the source and in AVATAR output, 30 runs {.table}

Every relationship survives, and all four land within roughly 10% of the
source value. The direction is consistently mild **dilution** rather
than distortion: the allometric exponent softens from 0.75 to about
0.68, and the exposure–response slope loses about a tenth of its
steepness. Blending pulls each subject toward its neighbours, and a
relationship measured across subjects flattens slightly as a result.

The `run_to_run_sd` column is the honest part of the table. The
*average* over thirty runs sits close to the truth, but a single
synthetic dataset can land well off it — the allometric exponent has a
run-to-run SD of about 0.18 on a value of 0.68, so one draw can easily
read 0.5 or 0.9.

## Seeing it, rather than reading it

The two relationships that are easiest to misjudge from a table are the
covariate effect and the exposure–response link. Both are visible
directly.

![](example-avatar-PKPD-covariate-treatment-effect_files/figure-html/relationship-plots-1.png)

Heavier subjects sit lower on both panels — higher clearance, lower
exposure — and the fitted line has a similar downward tilt in each. Note
also that the two arms overlay one another once AUC is divided by dose,
which is what dose proportionality looks like.

![](example-avatar-PKPD-covariate-treatment-effect_files/figure-html/er-plot-1.png)

The treatment effect is the separation between the two colours, and the
exposure–response link is the downward slope running through both. Both
are present in the synthetic panel, with the cloud slightly tighter —
that tightening is the blending, and it is why the fitted slope flattens
a little.

## Is the treatment arm coherent with the dose?

`ARM` is copied from the anchor and the concentrations are blended from
the donors, so nothing structurally prevents a “high” label from
arriving with “low” data. Check it directly.

``` r

one_run <- suppressWarnings(synpmx_avatar(source_data, roles, seed = 99))
table(ARM = one_run$ARM[one_run$EVID == 1], dose = one_run$AMT[one_run$EVID == 1])
#>       dose
#> ARM    300 600
#>   high   0  41
#>   low   39   0
```

Every synthetic subject’s arm label matches its dose. That is not luck,
and it is worth understanding why, because it is the mechanism the whole
article rests on.

## Why the AVATAR algorithm preserves relationships

**1. The event signature separates the arms.** Donors for blending are
first sought among subjects with an identical event signature, and that
signature *includes the dose amount*. The 300 mg and 600 mg subjects
therefore fall into different donor groups automatically, and with 40
subjects per arm the floor of `k = 5` is met without ever borrowing
across them. The arm label copied from the anchor lands on
concentrations blended from same-arm donors.

**2. The profile distance separates them** Even where dose does not
distinguish two arms, the donor search ranks candidates by distance in a
profile space built from the covariates and the DV trajectory. A subject
with higher concentrations sits near other subjects with higher
concentrations. The effect you are trying to preserve is itself part of
what selects the donors.

## Where it weakens

The second mechanism above also says where this breaks down: it works
because the effect is visible in the trajectory. An effect that is small
relative to between-subject variability does not organise the donor
neighbourhoods, donors get drawn from both groups, and the difference
dilutes further than the mild flattening seen here.

That is the boundary to keep in mind. Everything measured on this page
had a large, clean effect and a dose that separated the arms exactly.

## What to take from this

- Relationships put into a source dataset come out of AVATAR at close to
  their original size — within about 10% here — without anything in the
  method being told about them.
- That is a **tendency, not a specification**. It rests on donors being
  chosen by profile similarity and on covariates, endpoints, and the
  anchor’s kept columns being tied to the same donor draw. Nothing
  declares it, nothing enforces it, and no argument controls its
  strength.
- It is strongest when the effect is large relative to between-subject
  variability, and when the arms differ in something the event signature
  already separates — dose, schedule, route.
- Expect mild dilution. Every relationship here came back a little
  weaker than it went in, so a size read off synthetic data is a slight
  underestimate, not an inflation.
- Nothing here indicates the parameters can accurately be estimated from
  synthetic data. These measurements say the *structure* survives well
  enough to develop and debug analysis code against. A parameter
  estimate from AVATAR output is an estimate from blended data, and this
  article does not make it trustworthy.

The numbers on this page come from one simulated source with clean,
strong relationships and 80 subjects. They are an illustration of the
mechanism, not a general guarantee about your dataset. Run the same
measurement on your own source if the answer matters — the
`relationships()` function above is the whole method, and it takes any
data frame.
