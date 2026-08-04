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
| **Teach and compare.** Show what the different synthetic data generation methods do. | **Yes**, secondarily. |

The main use caes of this package is for sharing realistic-looking study
data outside the GxP computing environment but still within the
organization, so that code development can occur without the real data.
In some cases the GxP environment does not permit the most advanced
agentic coding tools, because of the risk of misalignment or unintended
agent behavior. Working with synthetic data lets those tools be used
without exposing them to patient data.

## The main deliverable: AVATAR method

The main deliverable of this package is an implementation of the AVATAR
method \[1, 2\], which offers some data masking but no formal privacy
guarantee. AVATAR works by blending together patient profiles, and
requires no model to be specified. The original method was developed in
Guillaudeux and colleagues \[2\], and Destere and colleagues benchmark a
modified AVATAR for population PK datasets \[1\]. However, they did not
test their method using pharmacometrics dataset from actual clinical
trials. We have done here, and in the process we have added features
such as handling BLOQ data and masking patients with unique dosing
schedules and observation times.

The key function is
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md),
which builds artificial profiles from real patient profiles, It masks
identifiable characteristics. but does not offer formal privacy
guarantees. For didactic purposes, the package also provides code for
other data masking methods, using trial simulation from prior knowledge,
and differential privacy methods. These methods give more formal privacy
protection and they are included here to illustrate the tradeoffs
between ways of generating synthetic data. These methods are not
actively maintained.

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

**If your environment blocks installing from GitHub** , then download
the source archive from
`https://github.com/iamstein/synpmx/archive/refs/heads/main.tar.gz` and
install from the file. This package needs nothing but base R, and no
compiler, since AVATAR is pure R:

``` r

install.packages("synpmx-main.tar.gz", repos = NULL, type = "source")
```

AVATAR needs nothing beyond base R. The DP methods additionally require
the official [OpenDP R
package](https://docs.opendp.org/en/stable/api/r/):

``` r

install.packages("opendp", repos = "https://opendp.r-universe.dev")
```

## Running it on your own study

AVATAR, called by
`synpmx_avatar() needs two things: the data, and a declaration of what its columns mean. A model is not needed. The function drops every column that is not described, so a column you forget is describe is dropped rather than quietly copied out of a real patient. Only`id`,`time`,`dv`, and`evid\`
are required; everything else is optional.

The block below stands in for your study — replace the first dozen lines
with your own data frame and edit the column descriptions to match your
dataset. The example below shows every declarable column, so most
studies will use a subset.

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
  dvid               = c("YTYPE", "NAME"),      # endpoint key; several columns may label the same endpoint,
                                                # If CMT is your only endpoint key, name it in both
                                                # cmt = "CMT", dvid = "CMT"
  mdv                = "MDV",                   # missing-DV indicator
  nominal_time       = "NTIME",                 # protocol visit time
  tad                = "TAD",                   # time after dose; recomputed
  occasion           = "OCC",                   # dosing occasion
  cens               = "CENS",                  # 1 = BLOQ, -1 = above, 0 = not
  limit              = "LIMIT",                 # the other interval boundary
  covariates         = c("WT", "AGE", "SEX"),   # measured; blended across donors
  dose_covariate     = "WT",                    # in this case, dose is body-weight based and it should be declared
  subject_properties = c("TRT", "TRTN"),        # assigned patient stratification variable
  keep               = "STUDYID"                # columns to be carried through verbatim
  # addl, ii          -- accepted and carried, but not expanded; expand
  #                      ADDL doses into explicit rows before synthesis
)

validate_pmx(study, roles)$valid
#> [1] TRUE
```

**Covariates and Subject Properties** are easy roles to confuse:

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
  study,             #study data
  roles,             #column desrciption
  n_subjects = NULL, # cohort size; NULL matches the source
  seed       = 2026)

validate_pmx(synthetic, roles)$valid
#> [1] TRUE
```

## The four modes for generating synthetic data

| Mode | Function | Output built from | Guarantee | Works at |
|----|----|----|----|----|
| **1. AVATAR blending** | [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md) | Real subject templates and blended real trajectories | None; governance only | At least 6 subjects |
| **2. Prior only** | [`synpmx_prior()`](https://iamstein.github.io/synpmx/reference/synpmx_prior.md) | A public model and protocol only | `epsilon = 0` (no data read) | Any (data-independent) |
| **3. Calibration** | [`synpmx_calibrated()`](https://iamstein.github.io/synpmx/reference/synpmx_calibrated.md) | A public model, magnitude corrected by 2 private releases | `(epsilon, delta)` DP | At least 20 subjects |
| **4. Empirical** | [`synpmx_empirical()`](https://iamstein.github.io/synpmx/reference/synpmx_empirical.md) | Dozens of noised population summaries | `(epsilon, delta)` DP | At least hundreds |

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
