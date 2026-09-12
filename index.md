# synpmx

`synpmx` builds synthetic pharmacometric datasets from actual datasets.

Documentation is available in the navigation bar at the top of the
website: <https://iamstein.github.io/synpmx/>

## Will the `synpmx` package support your use case?

There are many reasons to generate “synthetic data” and `synpmx` is
designed primarily to support certain use cases.

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

**❌ Answer scientific questions about the data (No).** Use the real
data for estimating parameters, selecting a model, quantifying a
covariate effect or choosing a dose.

The main use case of this package is for sharing realistic-looking study
data outside the GxP computing environment but still within the
organization, so that code development can occur without the real data.
In some cases the GxP environment does not permit the most advanced
agentic coding tools, because of the risk of misalignment or unintended
agent behavior. Working with synthetic data lets those tools be used
without exposing them to patient data.

## Synthetic data generation methods

Five generators build a synthetic dataset. The first three read a study
and build from it:

1.  [`synpmx_model()`](https://iamstein.github.io/synpmx/reference/synpmx_model.md)
    — Fits simple PK and PD models to the observation data, and
    statistical models to the dosing and missed visit data. Builds
    synthetic data by simulating from the models.
2.  [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
    — Blended values from real patients
3.  [`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md)
    — Principal component analysis from vector of all observations and
    covariates. Uses model to simulate any dose changes and missed
    visits.

These three all impute assay LOQ, and each takes a declaration of what
the dataset columns mean. The roles of `id`, `time`, `nominal_time`,
`dv`, `evid` are generally required.

The other two assert a public structural model instead, and cover the
case where data crosses a trust boundary and more formal privacy
protections are needed:

4.  [`synpmx_prior()`](https://iamstein.github.io/synpmx/reference/synpmx_prior.md)
    — Simulates a public model and a public protocol. Reads no data at
    all.
5.  [`synpmx_calibrated()`](https://iamstein.github.io/synpmx/reference/synpmx_calibrated.md)
    — Keeps the public model’s shape, and spends a small privacy budget
    correcting its magnitude.

All five are reviewed in
[synpmx-methods](https://iamstein.github.io/synpmx/articles/synpmx-methods.html),
which numbers them the same way the navigation bar does.

A sixth generator,
[`synpmx_empirical()`](https://iamstein.github.io/synpmx/reference/synpmx_empirical.md),
measures the trajectory shape rather than asserting it, and releases far
more numbers to do so — which splits the same privacy budget many ways
and makes it unusable below a few hundred subjects. It is exported so
the tradeoff can be checked, but it is not developed further and has no
documents of its own; one section of
[synpmx-methods](https://iamstein.github.io/synpmx/articles/synpmx-methods.html)
says why.

## Example (with AVATAR)

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
  roles,             #column description
  n_subjects = NULL, # cohort size; NULL matches the source
  seed       = 2026)
```

## Installation

Install `synpmx` from GitHub, then load it as usual:

``` r

# install.packages("remotes")
remotes::install_github("iamstein/synpmx")

library(synpmx)
```

If `remotes` is not available in your environment, you can also try:
`pak::pak("iamstein/synpmx")` or
`devtools::install_github("iamstein/synpmx")`.

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

## License

[MIT](https://iamstein.github.io/synpmx/LICENSE.md) © 2026 Andrew Stein.
