# synpmx

📖 **Website and documentation: <https://iamstein.github.io/synpmx/>**

`synpmx` builds **synthetic pharmacometric datasets** from actual
datasets.

## Which job do you want synthetic data to do?

There are many reasons to generate “synthetic data.” It is important to
be be explicit about your use case because the use case determines the
appropriate data generation algorithm.

| Use case | Served here? |
|----|----|
| **Develop code.** You need data with the right *shape* — schema, event grammar, covariates, dosing and sampling pattern. You’ll use this data to develop code for data processing, diagnostics, and model building outside the environment that holds the real study. | **Yes — this is what the package is for.** |
| **Send data past a trust boundary.** The output will reach people who cannot see the real data: a partner, a publication, a public repository. | **Only with care.** This needs a formal guarantee; see the privacy modes below, which are illustrative rather than audited. |
| **Answer the scientific question.** Estimate parameters, select a model, quantify a covariate effect, choose a dose, or stand in for real patients as a synthetic control arm. | **No.** Use the real data for this. |
| **Teach and compare.** Show what the different synthetic data generation methods do. | **Yes**, secondarily; that is why the non-default modes ship. |

The use case this package was built for: sharing realistic-looking study
data outside the GxP (Good Practice regulated) computing environment but
still within the organization, so that code development can happen
without the real data. In some cases the GxP environment does not permit
the most advanced agentic coding tools, because of the risk of
misalignment or unintended agent behavior. Working with synthetic data
lets those tools be used without exposing them to patient data.

## The main deliverable: AVATAR method

The main deliverable of the package is an implementation of the AVATAR
method \[1, 2\], which offers some masking but no formal privacy
guarantee. AVATAR works by blending together patient profiles, and
requires no model to be specified. “AVATAR” is a method name rather than
an initialism, from the patient-centric *avatarization* literature: the
original method was developed in Guillaudeux and colleagues \[2\], and
Destere and colleagues benchmark a modified AVATAR for population PK
\[1\]. Destere et al \[1\], did not test their method using
pharmacometrics dataset from actual clinical trials. We have done here,
and in the process we have added features such as handling BLOQ data and
masking patients with unique dosing schedules. The key function is
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md).

The package also provides code for other data masking methods: trial
simulation from prior knowledge, and two differential privacy (DP)
methods that give more formal protection. These methods are provided
mainly to illustrate the tradeoffs between ways of generating synthetic
data, and are not actively maintained.

## Installation

`synpmx` is not on CRAN; install it from GitHub, then load it as usual:

``` r

# install.packages("remotes")
remotes::install_github("iamstein/synpmx")

library(synpmx)
```

`remotes` is the smallest thing that does the job, and it is already
present in most environments. `pak::pak("iamstein/synpmx")` and
`devtools::install_github("iamstein/synpmx")` do the same, if you
already have either.

To pin a branch or commit, append it:
`remotes::install_github("iamstein/synpmx@main")`. Add
`build_vignettes = TRUE` if you want
[`vignette("synpmx-4-methods")`](https://iamstein.github.io/synpmx/articles/synpmx-4-methods.md)
offline; otherwise the same material is on the
[website](https://iamstein.github.io/synpmx/).

**If your environment blocks installing from GitHub** — likely in the
validated environments this package is meant for — download the source
archive from
`https://github.com/iamstein/synpmx/archive/refs/heads/main.tar.gz`,
move it across, and install from the file. This needs nothing but base
R, and no compiler, since AVATAR is pure R:

``` r

install.packages("synpmx-main.tar.gz", repos = NULL, type = "source")
```

If you were given the repository as a **ZIP** rather than a tarball,
unzip it first and install the resulting directory —
[`install.packages()`](https://rdrr.io/r/utils/install.packages.html)
cannot read a GitHub ZIP directly:

``` r

unzip("synpmx-main.zip")
install.packages("synpmx-main", repos = NULL, type = "source")
```

AVATAR needs nothing beyond base R. The other methods additionally
require the official [OpenDP R
package](https://docs.opendp.org/en/stable/api/r/):

``` r

install.packages("opendp", repos = "https://opendp.r-universe.dev")
```

## Running it on your own study

AVATAR needs two things: the data, and a declaration of what its columns
mean. There is no model to specify and nothing to fit.

The declaration is also the **manifest of what survives**.
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
drops every column no role names, so a column you forget is dropped
rather than quietly copied out of a real patient. Only `id`, `time`,
`dv`, and `evid` are required; everything else is optional.

The block below stands in for your study — replace the first dozen lines
with your own data frame and edit the role names to match your columns.
It is built to carry every declarable column at once, so most studies
will use a subset:

``` r

library(synpmx)

study <- pmx_simulated_fixture(24)                  # 24 subjects, 2 endpoints
study$YTYPE <- ifelse(study$DVID == "cp", 1L, 2L)   # endpoint key, numeric
study$NAME  <- as.character(study$DVID)             # same endpoint, as text
study$DVID  <- NULL
study$TRTN  <- ifelse(study$ID %% 2L == 1L, 1L, 2L) # assigned arm, numeric
study$TRT   <- ifelse(study$TRTN == 1L, "100 mg QD", "200 mg QD")
study$AMT[study$EVID != 0] <- 100 * study$TRTN[study$EVID != 0]
study$STUDYID <- "EXAMPLE-001"
bloq <- study$EVID == 0 & study$NAME == "cp" & study$DV < 1.2
study$DV[bloq]   <- 1.2                             # DV reports the limit
study$CENS[bloq] <- 1L                              # 1 = left-censored (BLOQ)
study$LIMIT <- ifelse(bloq, 0, NA_real_)            # the other boundary

roles <- pmx_roles(
  id                 = "ID",                    # subject identifier
  time               = "TIME",                  # actual elapsed time, numeric
  dv                 = "DV",                    # dependent variable
  evid               = "EVID",                  # event identifier
  amt                = "AMT",                   # dose amount
  rate               = "RATE",                  # infusion rate
  cmt                = "CMT",                   # compartment
  dvid               = c("YTYPE", "NAME"),      # endpoint key; several columns
                                                #   may label the same endpoint,
                                                #   first is authoritative. If
                                                #   CMT is your only endpoint
                                                #   key, name it twice:
                                                #   cmt = "CMT", dvid = "CMT"
  mdv                = "MDV",                   # missing-DV indicator
  nominal_time       = "NTIME",                 # protocol visit time
  tad                = "TAD",                   # time after dose; recomputed
  occasion           = "OCC",                   # dosing occasion
  cens               = "CENS",                  # 1 = BLOQ, -1 = above, 0 = not
  limit              = "LIMIT",                 # the other interval boundary
  covariates         = c("WT", "AGE", "SEX"),   # measured; blended across donors
  subject_properties = c("TRT", "TRTN"),        # assigned stratum; groups the
                                                #   dose rule and visit patterns
  keep               = "STUDYID"                # carried through verbatim
  # addl, ii          -- accepted and carried, but not expanded; expand
  #                      ADDL doses into explicit rows before synthesis
  # assigned_dose, exclude -- differential-privacy engines only
)

validate_pmx(study, roles)$valid
#> [1] TRUE
```

**Which role does a column want?** The three that are easy to confuse:

- `covariates` are *measured* characteristics. They are **blended**
  across the donors, so a synthetic subject’s weight is a new number
  nobody had.
- `subject_properties` are *assigned* strata — arm, dose group, cohort.
  They are copied from the anchor, and the stratum is what groups the
  dose rule and the pool of visit patterns an avatar may be given.
- `keep` is the escape hatch, for anything else you want carried through
  untouched: a study identifier, a randomization sequence, a units
  column, a redundant label. A kept value is **one real subject’s real
  value**, so keep only what the source data’s own access controls
  already permit.

## The masking options, and their defaults

Every masking argument below is shown at its default, so this call
behaves exactly like `synpmx_avatar(study, roles, seed = 2026)`:

``` r

synthetic <- synpmx_avatar(
  study, roles,
  n_subjects         = NULL,   # cohort size; NULL matches the source
  seed               = 2026,
  # --- how much of one real patient can reach one synthetic patient ---
  k                  = 5,      # donors blended into each avatar (the floor)
  max_donor_weight   = 0.50,   # no one donor is more than half an avatar
  on_donor_shortfall = "drop", # a route arm below k + 1 subjects is dropped
  # --- what structure is masked ---
  screen             = TRUE,   # never anchor on a subject whose follow-up or
                               #   dose count exceeds 2x the cohort's 90th pct
  coarsen_time       = TRUE,   # snap times onto a shared visit grid, then
                               #   resample pooled deviations back, so no avatar
                               #   carries one real visit schedule. Uses
                               #   `nominal_time` when declared; K-means centres
                               #   of the pooled times otherwise
  min_pattern_share  = 2L,     # an avatar's attended-visit pattern must be one
                               #   at least 2 real subjects share, so no
                               #   synthetic schedule is unique to one patient
  # --- how much noise is added ---
  subject_noise_sd   = 0.15,   # per-subject perturbation
  residual_noise_sd  = 0.05,   # within-trajectory noise
  residual_phi       = 0.6,    # AR(1) correlation in observation order
  time_jitter        = 0,      # realism only, NOT privacy: jitter is clamped
                               #   within half a gap of the source visit
  pca_variance       = 0.90    # profile variance retained for donor distances
)

validate_pmx(synthetic, roles)$valid
#> [1] TRUE
dim(synthetic)
#> [1] 384  21
head(synthetic[, c("ID", "TIME", "NTIME", "OCC", "NAME", "DV", "CENS", "TRT")], 6)
#>   ID TIME NTIME OCC NAME        DV CENS       TRT
#> 1 25 0.00  0.00   1   cp  0.000000    0 100 mg QD
#> 2 25 0.00  0.00   1   pd 90.823217    0 100 mg QD
#> 3 25 0.25  0.25   1   cp  1.200000    1 100 mg QD
#> 4 25 1.00  1.00   1   cp  8.249096    0 100 mg QD
#> 5 25 2.00  2.00   1   cp  3.766242    0 100 mg QD
#> 6 25 4.00  4.00   1   pd 78.636255    0 100 mg QD
```

Two things the defaults do that are worth knowing about. Where dosing is
proportional to a covariate (mg/kg, mg/m²), each avatar’s `AMT` is
**recomputed from its own blended covariate** rather than copied, so the
synthetic patient’s dose matches the synthetic patient’s weight. And a
source attendance pattern held by fewer than `min_pattern_share`
subjects is **lost, not approximated** — that loss is the mechanism
working, and every run reports how much of it happened.

Measure what the masking achieved with
[`skeleton_uniqueness()`](https://iamstein.github.io/synpmx/reference/skeleton_uniqueness.md)
on the source and
[`compare_pmx_proximity()`](https://iamstein.github.io/synpmx/reference/compare_pmx_proximity.md)
on the pair.

## The four modes

| Mode | Function | Output built from | Guarantee | Works at |
|----|----|----|----|----|
| **1. AVATAR blending** | [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md) | Real subject templates and blended real trajectories | None; governance only | At least 5 subjects |
| **2. Prior only** | [`synpmx_prior()`](https://iamstein.github.io/synpmx/reference/synpmx_prior.md) | A public model and protocol only | `epsilon = 0` (no data read) | Any (data-independent) |
| **3. Calibration** | [`synpmx_calibrated()`](https://iamstein.github.io/synpmx/reference/synpmx_calibrated.md) | A public model, magnitude corrected by 2 private releases | `(epsilon, delta)` DP | At least 20 subjects |
| **4. Empirical** | [`synpmx_empirical()`](https://iamstein.github.io/synpmx/reference/synpmx_empirical.md) | Dozens of noised population summaries | `(epsilon, delta)` DP | At least hundreds |

**The trust boundary decides whether you need differential privacy.**

If the generated data will not be accessible to anyone who cannot access
the original data, formal privacy guarantees are not needed and the
AVATAR blending approach is the recommended method because it is the
simplest method to use as it doesn’t require the specification of a
model.

On the other hand, if the synthetic data will reach those who do not
have access to the original data, then more formal methods with
mathematical trust guarantees are the appropriate methods of choice.

[`vignette("synpmx-4-methods")`](https://iamstein.github.io/synpmx/articles/synpmx-4-methods.md)
runs all four methods on the same dataset and shows the results side by
side.

## Maintenance status

**AVATAR blending is the primary, maintained code.** It has no
dependencies beyond base R, and is what to reach for when the output
stays within the source data’s own access controls and obligations.
However, AVATAR does not offer any formal, mathematical guarantees
around privacy.

The three other modes (**prior**, **calibration**, **empirical**) are
secondary; but, they are present in this repository because they cover
scenarios where data crosses a trust boundary and formal privacy
conditions must be met. The methods are provided as-is: not under active
development and **not independently privacy-audited**. Treat them as a
principled demonstration of the privacy/utility tradeoff, not as a
production ready. That status is enforced in that
[`synpmx_calibrated()`](https://iamstein.github.io/synpmx/reference/synpmx_calibrated.md)
and
[`synpmx_empirical()`](https://iamstein.github.io/synpmx/reference/synpmx_empirical.md)
refuse to run until
[`synpmx_enable_dp_engines()`](https://iamstein.github.io/synpmx/reference/synpmx_enable_dp_engines.md)
has been called once in the session.

## Key Documentation

| Document | Question it answers |
|----|----|
| [The four synthetic generation modes](https://iamstein.github.io/synpmx/articles/synpmx-4-methods.html) | What are the modes, and which one do I want? **Start here.** |
| [Demo; Using synpmx AVATAR with 5 datasets](https://iamstein.github.io/synpmx/articles/synpmx-demo.html) | How do I run this on my own study? |
| [The AVATAR Algorithm](https://iamstein.github.io/synpmx/articles/avatar-algorithm.html) | How does the default generator work, step by step? |

## References

1.  Destere A, Lombardi R, Labriffe M, et al. *Can synthetic data
    overcome the privacy and fidelity bottleneck in Pharmacometrics? A
    comparative benchmark using a daptomycin population pharmacokinetic
    model.* medRxiv preprint, posted June 2, 2026. doi:
    [10.64898/2026.05.30.26354512](https://doi.org/10.64898/2026.05.30.26354512).

2.  Guillaudeux M, Rousseau O, Petot J, et al. Patient-centric synthetic
    data generation, no reason to risk re-identification in biomedical
    data analysis. *npj Digital Medicine.* 2023;6. doi:
    [10.1038/s41746-023-00771-5](https://doi.org/10.1038/s41746-023-00771-5).

## License

[MIT](https://iamstein.github.io/synpmx/LICENSE.md) © 2026 Andrew Stein.
