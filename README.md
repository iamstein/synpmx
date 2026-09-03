# synpmx
📖 **Website and documentation: <https://iamstein.github.io/synpmx/>**

`synpmx` builds **synthetic pharmacometric datasets** from actual datasets.

## Will the `synpmx` package support your use case?

There are many reasons to generate "synthetic data."  It is important to be be explicit about your use case because the use case determines whether `synpmx` can fully support you.

**✅ Develop code (Intended Use Case)** You need synthetic data that resembles the true data — schema, event grammar, covariates, dosing, sampling, censoring, and drop-out pattern.  You'll use this data to develop code for data processing, diagnostics, and model building outside the environment that holds the real study.  

**✅ Teaching tool for comparing synthetic data methods (Yes).** Illustrate the difference between synthetic data generation methods.

**⚠️ Send data past a trust boundary (Use Caution).** If the output will reach people who cannot see the real data: a partner, a publication, a public repository, this package should be used with caution. The formal privacy-protecting methods provided with this package are illustrative, but not audited.  Carefully assess what level of privacy protection is needed.

**❌ Answer scientific questions about the the data (No).** Use the real data for estimating parameters, selecting a model, quantifying a covariate effect or choosing a dose. 

The main use caes of this package is for sharing realistic-looking study data
outside the GxP computing environment but still within
the organization, so that code development can occur without the real data. In
some cases the GxP environment does not permit the most advanced agentic coding
tools, because of the risk of misalignment or unintended agent behavior. Working
with synthetic data lets those tools be used without exposing them to patient
data.

## Three ways to build a dataset from a real one

Three generators read a study and build a synthetic one from it. They are peers
rather than a default and its alternatives: each carries something different out
of the source, and which one suits a given study is still an open question.
None offers a formal privacy guarantee.

Synthetic data based on:

1. `synpmx_model()` — Fits model to the data, and bases synthetic data off of 
Fixed effects, a residual error, and a model for dose changes and missed visits
2. `synpmx_pca()` — Principal component analysis from vector of all observations
   and covariates, with a model for dose changes and missed visits.
3. `synpmx_avatar()` — Blended values from real patients

All the above algorithms impute assay LOQ. All three take a declaration of what the columns mean. 
Only the roles of `id`, `time`, `nominal_time`, `dv`, `evid` are generally required.

There are three additional generation algorithms provided that cover the case where data crosses a trust boundary
and more formal privacy protections are needed.  Treat them
as a principled demonstration of the privacy/utility tradeoff rather than as
production ready — a status enforced in that `synpmx_calibrated()` and
`synpmx_empirical()` refuse to run until `synpmx_enable_dp_engines()` has been
called once in the session.

4. `synpmx_prior()` - take a public structural model and prespecified design, not using the real data all.
5. `synpmx_calibration()` - use data summaries of a few parameters (for smaller trials), and uses a differential-privacy budget to protect privacy.
6. `synpmx_empirical()`  - use data summaries of many parameters (for larger trials), and uses a differential-privacy budget to protect privacy.

## Generating Synthetic Data

`synpmx_avatar()` is the easiest place to start. Every column that is not described is dropped.

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

## Installation

`synpmx` is not on CRAN; install it from GitHub, then load it as usual:

``` r
# install.packages("remotes")
remotes::install_github("iamstein/synpmx")

library(synpmx)
```

If `remotes` is not available in your environment, you can also try:
`pak::pak("iamstein/synpmx")` or `devtools::install_github("iamstein/synpmx")`.
either.

If your environment blocks installing from GitHub, then download the source archive from
`https://github.com/iamstein/synpmx/archive/refs/heads/main.tar.gz` and install from the file. 
This package needs nothing but base R.

``` r
install.packages("synpmx-main.tar.gz", repos = NULL, type = "source")
```

The Data Privacy methods additionally require the official [OpenDP R package](https://docs.opendp.org/en/stable/api/r/):

``` r
install.packages("opendp", repos = "https://opendp.r-universe.dev")
```

## Documentation

Use the navigation bar at the top of the published website for further documentation.
<https://iamstein.github.io/synpmx/>
## License

[MIT](LICENSE.md) © 2026 Andrew Stein.
