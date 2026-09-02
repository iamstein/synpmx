# synpmx

📖 **Website and documentation: <https://iamstein.github.io/synpmx/>**

`synpmx` builds **synthetic pharmacometric datasets** from actual
datasets.

## Will the `synpmx` package support your use case?

There are many reasons to generate “synthetic data.” It is important to
be be explicit about your use case because the use case determines
whether `synpmx` can fully support you.

**✅ Develop code (Intended Use Case)** You need synthetic data that
resembles the true data — schema, event grammar, covariates, dosing,
sampling, censoring, and drop-out pattern. You’ll use this data to
develop code for data processing, diagnostics, and model building
outside the environment that holds the real study.

**✅ Teaching tool for comparing synthetic data methods (Yes).**
Illustrate the difference between synthetic data generation methods.

**⚠️ Send data past a trust boundary (Use Caution).** If the output will
reach people who cannot see the real data: a partner, a publication, a
public repository, this package should be used with caution. The formal
privacy-protecting methods provided with this package are illustrative,
but not audited. Carefully assess what level of privacy protection is
needed.

**❌ Answer scientific questions about the the data (No).** Use the real
data for estimating parameters, selecting a model, quantifying a
covariate effect or choosing a dose.

The main use caes of this package is for sharing realistic-looking study
data outside the GxP computing environment but still within the
organization, so that code development can occur without the real data.
In some cases the GxP environment does not permit the most advanced
agentic coding tools, because of the risk of misalignment or unintended
agent behavior. Working with synthetic data lets those tools be used
without exposing them to patient data.

## Three ways to build a dataset from a real one

Three generators read a study and build a synthetic one from it. They
are peers rather than a default and its alternatives: each carries
something different out of the source, and which one suits a given study
is still an open question. None offers a formal privacy guarantee.

| Function | Synthetic data based on |  |  |
|----|----|----|----|
| `1. synpmx_model()` | Fixed effects, a covariance matrix, a residual error, and model for dose changes and missed visits |  |  |
| `2. synpmx_pca()` | Principle component analysis from vector of all observations and covariates |  |  |
| `3. synpmx_avatar()` | Blended values from real neighbouring patients |  |  |

All three take a declaration of what the columns mean. Only the roles of
`id`, `time`, `nominal_time`, `dv`, `evid` are generally required.

Three additional functions
([`synpmx_prior()`](https://iamstein.github.io/synpmx/reference/synpmx_prior.md),
`synpmx_calibration()`,
[`synpmx_empirical()`](https://iamstein.github.io/synpmx/reference/synpmx_empirical.md))
take a public structural model and either base the data off purely a
public model and design information about the study, or they
differential-privacy budget to correct and use data summaries. These
methods cover the case where data crosses a trust boundary and formal
privacy conditions must be met. Treat them as a principled demonstration
of the privacy/utility tradeoff rather than as production ready — a
status enforced in that
[`synpmx_calibrated()`](https://iamstein.github.io/synpmx/reference/synpmx_calibrated.md)
and
[`synpmx_empirical()`](https://iamstein.github.io/synpmx/reference/synpmx_empirical.md)
refuse to run until
[`synpmx_enable_dp_engines()`](https://iamstein.github.io/synpmx/reference/synpmx_enable_dp_engines.md)
has been called once in the session.

## Generating Synthetic Data

[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
is the easiest place to start. Every column that is not described is
dropped.

``` r

library(synpmx)

study <- as.data.frame(get(utils::data(list = "case1_pkpd", package = "xgxr")))
study$CENS[study$NAME == "PD - Continuous"] <- 0  # CENS here flags the PK assay limit only

# ?pmx_roles` describes the options here
roles <- pmx_roles(
  id             = "ID",                 # subject identifier - REQUIRED
  time           = "TIME",               # actual elapsed time, numeric - REQUIRED
  dv             = "LIDV",               # dependent variable - REQUIRED
  evid           = "EVID",               # event identifier - REQUIRED
  amt            = "AMT",                # dose amount
  cmt            = "CMT",                # compartment
  dvid           = "NAME",               # endpoint key: which endpoint the row reports
  mdv            = NULL,                 # missing-dependent-variable flag
  rate           = NULL,                 # infusion rate
  nominal_time   = "NOMTIME",            # protocol visit time
  tad            = NULL,                 # time after dose; this is not used by AVATAR, it is recomputed
  occasion       = NULL,                 # set if TIME resets by occasion
  cens           = "CENS",               # 1 = BLOQ, -1 = above, 0 = not
  limit          = NULL,                 # other end of the censoring interval
  addl           = NULL,                 # additional doses
  ii             = NULL,                 # interdose interval
  covariates     = "WEIGHTB",            # patient baseline covariates; blended across donors
  strata         = c("TRTACT", "DOSE"),  # assigned arm / dose group / cohort (default is to balance synthetic data by strata)
  dose_covariate = NULL,                 # covariate the dose is a fixed multiple of (e.g. WEIGHTB for weight based dosing)
  endpoint_types = NULL,                 # value kind of each DV variable (continuous, binary, ordinal) per endpoint; inferred when NULL
  keep           = "STUDY",              # columns carried through verbatim
)

synthetic <- synpmx_avatar(
  study,             #study data
  roles,             #column desrciption
  n_subjects = NULL, # cohort size; NULL matches the source
  seed       = 2026)
```

## Maintenance status

All three data-reading generators are maintained and none is retired.
Which of them to reach for is not settled, and the documents are written
to be read by comparison for that reason: each has an algorithm document
and a demo over the same ground, and
[`synpmx_scorecard()`](https://iamstein.github.io/synpmx/reference/synpmx_scorecard.md)
scores any of their outputs the same way. None of them offers a formal,
mathematical privacy guarantee.

## Installation

`synpmx` is not on CRAN; install it from GitHub, then load it as usual:

``` r

# install.packages("remotes")
remotes::install_github("iamstein/synpmx")

library(synpmx)
```

If `remotes` is not available in your environment, you can also try:
`pak::pak("iamstein/synpmx")` or
`devtools::install_github("iamstein/synpmx")`. either.

If your environment blocks installing from GitHub, then download the
source archive from
`https://github.com/iamstein/synpmx/archive/refs/heads/main.tar.gz` and
install from the file. This package needs nothing but base R.

``` r

install.packages("synpmx-main.tar.gz", repos = NULL, type = "source")
```

The Data Privacy methods additionally require the official [OpenDP R
package](https://docs.opendp.org/en/stable/api/r/):

``` r

install.packages("opendp", repos = "https://opendp.r-universe.dev")
```

## Key Documentation

| Document | Question it answers |
|----|----|
| [Demo: one dataset, end to end](https://iamstein.github.io/synpmx/articles/avatar-demo.html) | What does a whole run look like, from raw event table to checked synthetic one? |
| [Evaluating AVATAR on public data](https://iamstein.github.io/synpmx/articles/avatar-public-data-examples.html) | How well does it work, and what did the masking cost, on eight public datasets? |
| [The AVATAR Algorithm](https://iamstein.github.io/synpmx/articles/avatar-algorithm.html) | How does blending work, step by step? |
| [The PCA Algorithm](https://iamstein.github.io/synpmx/articles/pca-algorithm.html) and its [demo](https://iamstein.github.io/synpmx/articles/pca-demo.html) | How does the component-basis generator work, and what does a run look like? |
| [Evaluating PCA on public data](https://iamstein.github.io/synpmx/articles/pca-public-data-examples.html) | How well does it work across eight public datasets? |
| [The PMX Model Algorithm](https://iamstein.github.io/synpmx/articles/pmxmodel-algorithm.html) and its [demo](https://iamstein.github.io/synpmx/articles/pmxmodel-demo.html) | How does the population-model generator work, and what does a run look like? |
| [Evaluating the PMX model generator on public data](https://iamstein.github.io/synpmx/articles/pmxmodel-public-data-examples.html) | Which studies can it fit at all, and what does it lose on the ones it can? |
| [Scorecard: Checks of the synthetic data](https://iamstein.github.io/synpmx/articles/avatar-scorecard.html) | I have a synthetic dataset. Should I use it? |
| [The synthetic generation modes](https://iamstein.github.io/synpmx/articles/synpmx-methods.html) | What are the six modes, all six on one study, and which one do I want? |

## License

[MIT](https://iamstein.github.io/synpmx/LICENSE.md) © 2026 Andrew Stein.
