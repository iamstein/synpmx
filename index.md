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

## The main deliverable: AVATAR method

The main deliverable of this package is an implementation of the AVATAR
method \[1, 2\], which offers some data masking but no formal privacy
guarantee. AVATAR works by blending together patient profiles, and
requires no model to be specified. The original method was developed in
Guillaudeux and colleagues \[2\], and Destere and colleagues benchmark a
modified AVATAR for population PK datasets \[1\]. However, they did not
test their method using pharmacometrics dataset from actual clinical
trials. We have developed the synthetic data method here using actual
clinical data and thus it handles realistic aspects of trials such as
BLOQ data, weight-based dosing, multiple regimens, and it also has
masking methods added so that patients cannot be identified by a unique
dosing schedule or set of observation times. The key function is
`synpmx_avatar`, which builds artificial profiles from real patient
profiles, It masks identifiable characteristics, but does not offer
formal privacy guarantees.

For teaching purposes, the package also provides code for other data
masking methods, using trial simulation from prior knowledge, and
differential privacy methods. These methods give more formal privacy
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

While the AVATAR method needs nothing beyond base R. The DP methods
additionally require the official [OpenDP R
package](https://docs.opendp.org/en/stable/api/r/):

``` r

install.packages("opendp", repos = "https://opendp.r-universe.dev")
```

## Generating Synthetic Data with AVATAR algorithm

AVATAR is called by
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md).
It needs two things: the data, and a declaration of what its columns
mean. A model is not needed. Every column that is not described is
dropped. Only `id`, `time`, `dv`, and `evid` are required; everything
else is optional.

The example below runs on
[`xgxr::case1_pkpd`](https://rdrr.io/pkg/xgxr/man/case1_pkpd.html), the
same dataset as the
[demo](https://iamstein.github.io/synpmx/articles/demo.html): 180
subjects, six arms from placebo to 300 mg, a PK and a PD endpoint. Swap
in your own data frame and edit the column descriptions to match it.

``` r

library(synpmx)

study <- as.data.frame(get(utils::data(list = "case1_pkpd", package = "xgxr")))
study$CENS[study$NAME == "PD - Continuous"] <- 0  # CENS here flags the PK assay limit only

roles <- pmx_roles(
  id           = "ID",                 # subject identifier
  time         = "TIME",               # actual elapsed time, numeric
  dv           = "LIDV",               # dependent variable
  evid         = "EVID",               # event identifier
  amt          = "AMT",                # dose amount
  cmt          = "CMT",                # compartment
  dvid         = "NAME",               # endpoint key: which endpoint the row reports
  nominal_time = "NOMTIME",            # protocol visit time
  cens         = "CENS",               # 1 = BLOQ, -1 = above, 0 = not
  covariates   = "WEIGHTB",            # measured; blended across donors
  strata       = c("TRTACT", "DOSE"),  # assigned arm / dose group / cohort
  keep         = "STUDY"               # carried through verbatim
)
```

This study has no column for the other roles, so they are left out:
`rate`, `mdv`, `occasion`, `tad`, `limit` (the other boundary when
`cens` is set), `dose_covariate` (the covariate a weight-based dose is
computed from), and `addl`/`ii` (carried but not expanded — expand ADDL
doses into explicit rows before synthesis).
[`?pmx_roles`](https://iamstein.github.io/synpmx/reference/pmx_roles.md)
describes each one.

``` r

synthetic <- synpmx_avatar(
  study,             #study data
  roles,             #column desrciption
  n_subjects = NULL, # cohort size; NULL matches the source
  seed       = 2026)
```

## Maintenance status

**AVATAR blending is the primary, maintained code.** It has no
dependencies beyond base R, and is what to reach for when the output
stays within the source data’s own access controls and obligations.
However, AVATAR does not offer any formal, mathematical guarantees
around privacy.

The three other modes for generating synthetic data (**prior**,
**calibration**, **empirical**) are secondary; but, they are present in
this repository because they cover scenarios where data crosses a trust
boundary and formal privacy conditions must be met. Treat them as a
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
| [Demo: one dataset, end to end](https://iamstein.github.io/synpmx/articles/demo.html) | What does a whole run look like, from raw event table to checked synthetic one? |
| [Evaluating AVATAR on public data](https://iamstein.github.io/synpmx/articles/public-data-examples.html) | How well does it work, and what did the masking cost, on eight public datasets? |
| [The AVATAR Algorithm](https://iamstein.github.io/synpmx/articles/avatar-algorithm.html) | How does the default generator work, step by step? |
| [Scorecard: Checks of the synthetic data](https://iamstein.github.io/synpmx/articles/scorecard-synthetic-data-checks.html) | I have a synthetic dataset. Should I use it? |
| [The four synthetic generation modes](https://iamstein.github.io/synpmx/articles/synpmx-4-methods.html) | What are the modes, and which one do I want? |

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
