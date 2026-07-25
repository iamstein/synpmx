<!-- Edit this file directly; it is the only README. The example below is not
     executed when the file is built, so it is pinned by a regression test:
     tests/testthat/test-readme.R runs the same code and asserts the same
     output. If that test fails, the example here is stale — update both. -->

# synpmx

`synpmx` builds **synthetic pharmacometric datasets** by using dosing and
measurement event tables with the same setup as a real study.  The intended use case
is to provide realistic looking data of a study so that can be shared outside the 
GxP computing environment (but stay within the organization) so that data-assembly code, diagnostic plots, and model-run plumbing can be developed outside the restricted
computing environment that holds the real data, where more advanced coding agents
can more often be used.

In general, when developing synthetic data, it's important to think carefully about 
**how much information about the real data is allowed to survive into the synthetic data**.
This is a privacy question and `synpmx` offers four generation modes at different points on that scale, and helps you pick one, though the main deliverable of the package uses an
implementation of the AVATAR method which offers some blinding, but not formal privacy guarantees.

## Installation

``` r
# install.packages("pak")
pak::pak("iamstein/synpmx")
```

AVATAR needs nothing beyond base R. The two differentially private modes
additionally require the official [OpenDP R package](https://docs.opendp.org/en/stable/api/r/):

``` r
install.packages("opendp", repos = "https://opendp.r-universe.dev")
```

## A first synthetic dataset

The default AVATAR mode needs only the data and a declaration of what
the columns mean.

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

The output keeps the same structure.  

The other three modes for generating synthetic data,  cover the scenarios where formal data privacy requirements are needed.  They are included in the package 
as part of a conceptual framework to help the modeler think about various
options that are available, but these methods have not been yet been audited
for their ability to protect data privacy.  

| Mode | Function | Output built from | Guarantee | Works at |
|----|----|----|----|----|
| **1. AVATAR blending** | `synpmx_avatar()` | Real subject templates and blended real trajectories | None; governance only | ~12 subjects up |
| **2. Prior only** | `synpmx_prior()` | A public model and protocol only | `epsilon = 0` (no data read) | Any (data-independent) |
| **3. Calibration** | `synpmx_calibrated()` | A public model, magnitude corrected by 2 private releases | `(epsilon, delta)` DP | ~20 subjects up |
| **4. Empirical** | `synpmx_empirical()` | Dozens of noised population summaries | `(epsilon, delta)` DP | A few hundred up |

Two rules of thumb decide between them:

- **The trust boundary decides whether you need differential privacy.**
  Ask whether the generated data can reach anyone the source data could
  not.  Moving the synthetic data to an environment where the same people can access it
  is different from sending it to external vendors or making it public.  
  If no one new can access the data, AVATAR is more useful and easy to use.  

`vignette("synpmx-method")` runs all four methods on the same dataset and shows
the results side by side.

## Maintenance status

**AVATAR blending is the primary, maintained code.** It has no
dependencies beyond base R, and is what to reach for when the output
stays within the source data’s own access controls and obligations.

The three differentially private modes (**prior**, **calibration**,
**empirical**) are secondary.  They are present in this repository
because formal privacy is the right answer when data crosses a trust
boundary. They are provided as-is: not under active development,
 and **not independently privacy-audited**.
Treat them as a principled demonstration of the privacy/utility
tradeoff, not as a production release mechanism. A real regulated
release needs the specialist review.

That status is enforced in that `synpmx_calibrated()` and
`synpmx_empirical()` refuse to run until `synpmx_enable_dp_engines()`
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

## What this is not

Generated data exercises cleaning, joins, reshaping, plots, control-file
plumbing, repeated-dose PK code, longitudinal PD/biomarker code,
infusion events, and censoring conventions. It aims for broad magnitude
and shape.

It is **not** appropriate for parameter estimation, inference, model
selection, dose selection, or clinical conclusions, and it does not
reproduce source distributions, parameter estimates, or
covariate-response relationships.

Differential privacy, where used, is mathematically bounded rather than
absolute. It does not guarantee impossibility of linkage or
re-identification, establish legal anonymity, authorize release, secure
a compromised environment, or validate public-input claims. Independent
privacy, legal, information-security, and data-governance review remains
required.

## License

[MIT](LICENSE.md) © 2026 Andrew Stein.
