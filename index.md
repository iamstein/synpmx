# synpmx

`synpmx` builds **synthetic pharmacometric datasets** from actual
datasets. The intended use case is to provide realistic looking data of
a study that can be shared outside the GxP computing environment (but
stay within the organization) so that code for data-assembly code,
diagnostic plots, and modeling can be developed outside the restricted
computing environment that holds the real data. in some cases the GxP
xomputijg environment does not permit isebof the most advanced agentic
coding tools due to possible risks from misalignment and unintended
behavior. However the most advanced coding agents can be useful for code
development.

In general, when developing synthetic data, it’s important to think
carefully about **how much information about the real data is allowed to
survive into the synthetic data**. This is a privacy question and
`synpmx` offers four options. The main deliverable of the package uses
an implementation of the AVATAR method which offers some blinding, but
not formal privacy guarantees.

## Installation

``` r

# install.packages("pak")
pak::pak("iamstein/synpmx")
```

AVATAR needs nothing beyond base R. Other additionally require the
official [OpenDP R package](https://docs.opendp.org/en/stable/api/r/):

``` r

install.packages("opendp", repos = "https://opendp.r-universe.dev")
```

## A first synthetic dataset

The default AVATAR synthetic data algorithm needs only the data and a
declaration of what the columns mean.

``` r

library(synpmx)
data("theo_md", package = "nlmixr2data")

roles <- pmx_roles(
  id = "ID", time = "TIME", dv = "DV", amt = "AMT",
  evid = "EVID", cmt = "CMT", covariates = "WT"
)

synthetic <- suppressWarnings(synpmx_avatar(theo_md, roles, seed = 101))
validate_pmx(synthetic, roles)$valid
#> [1] TRUE
head(synthetic, 4)
#>   ID TIME          DV    AMT EVID CMT       WT
#> 1 13 0.00  0.00000000 267.84  101   1 85.25496
#> 2 13 0.00  0.02403989   0.00    0   2 85.25496
#> 3 13 0.30 10.00272146   0.00    0   2 85.25496
#> 4 13 0.63 11.89216329   0.00    0   2 85.25496
```

The output dataset keeps the same structure.

The other three modes for generating synthetic data cover the scenarios
where formal data privacy requirements are needed. They are included in
the package as part of a conceptual framework to help the modeler think
about various options that are available, but these methods have not
been audited for their ability to protect data privacy.

| Mode | Function | Output built from | Guarantee | Works at |
|----|----|----|----|----|
| **1. AVATAR blending** | [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md) | Real subject templates and blended real trajectories | None; governance only | ~12 subjects up |
| **2. Prior only** | [`synpmx_prior()`](https://iamstein.github.io/synpmx/reference/synpmx_prior.md) | A public model and protocol only | `epsilon = 0` (no data read) | Any (data-independent) |
| **3. Calibration** | [`synpmx_calibrated()`](https://iamstein.github.io/synpmx/reference/synpmx_calibrated.md) | A public model, magnitude corrected by 2 private releases | `(epsilon, delta)` DP | ~20 subjects up |
| **4. Empirical** | [`synpmx_empirical()`](https://iamstein.github.io/synpmx/reference/synpmx_empirical.md) | Dozens of noised population summaries | `(epsilon, delta)` DP | A few hundred up |

Two rules of thumb decide between them:

- **The trust boundary decides whether you need differential privacy.**
  Ask whether the generated data can reach anyone the source data could
  not. Moving the synthetic data to a new cmouting environment where the
  same people can access it is different from sending the data to
  external vendors or making it public.  
  If no one new can access the data, AVATAR is much easier to use.

[`vignette("synpmx-method")`](https://iamstein.github.io/synpmx/articles/synpmx-method.md)
runs all four methods on the same dataset and shows the results side by
side.

## Maintenance status

**AVATAR blending is the primary, maintained code.** It has no
dependencies beyond base R, and is what to reach for when the output
stays within the source data’s own access controls and obligations.

The three differentially private modes (**prior**, **calibration**,
**empirical**) are secondary. They are present in this repository
because formal privacy is the right answer when data crosses a trust
boundary. They are provided as-is: not under active development and and
**not independently privacy-audited**. Treat them as a principled
demonstration of the privacy/utility tradeoff, not as a production
release mechanism. A real regulated release needs the specialist review.

That status is enforced in that
[`synpmx_calibrated()`](https://iamstein.github.io/synpmx/reference/synpmx_calibrated.md)
and
[`synpmx_empirical()`](https://iamstein.github.io/synpmx/reference/synpmx_empirical.md)
refuse to run until
[`synpmx_enable_dp_engines()`](https://iamstein.github.io/synpmx/reference/synpmx_enable_dp_engines.md)
has been called once in the session.

## Documentation

| Document | Question it answers |
|----|----|
| [The four generation modes](https://iamstein.github.io/synpmx/articles/synpmx-method.html) | What are the modes, and which one do I want? **Start here.** |
| [Using synpmx](https://iamstein.github.io/synpmx/articles/synpmx-demo.html) | How do I run this on my own study? |
| [Privacy in synpmx](https://iamstein.github.io/synpmx/articles/synpmx-privacy.html) | What does differential privacy guarantee, does my release need it, and what epsilon? |
| [AVATAR mathematics](https://iamstein.github.io/synpmx/articles/avatar-mathematics.html) | How does the default generator work, step by step? |
| [Privacy background](https://iamstein.github.io/synpmx/articles/privacy-background.html) | Where do `d`, `f`, and the error law come from? |
| [Feasibility by cohort size](https://iamstein.github.io/synpmx/articles/feasibility.html) | What can actually be released, and from how many patients? |
| [Mechanism-level privacy argument](https://iamstein.github.io/synpmx/articles/privacy-argument.html) | The formal argument, for a reviewer. |
| [Model elicitation](https://iamstein.github.io/synpmx/articles/model-elicitation.html) / [data elicitation](https://iamstein.github.io/synpmx/articles/data-elicitation.html) | How do I produce the public inputs modes 2–4 need? |

The full function reference is at
<https://iamstein.github.io/synpmx/reference/>.

## License

[MIT](https://iamstein.github.io/synpmx/LICENSE.md) © 2026 Andrew Stein.
