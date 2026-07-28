# synpmx

📖 **Website and documentation: <https://iamstein.github.io/synpmx/>**

`synpmx` builds **synthetic pharmacometric datasets** from actual datasets. The intended use case
is to provide realistic looking data of a study that can be shared outside the 
GxP computing environment (but stay within the organization) so that data-assembly code, diagnostic plots, and modeling can be developed outside the restricted
computing environment that holds the real data. 

In some cases the GxP computing environment does not allow use of the most advanced agentic coding tools due to risks from misalignment and unintended behavior of the agents.  Using synthetic data allows the most advanced agentic coding tools to be used without exposing them to patient data.  

In general, when developing synthetic data, it's important to think carefully about 
**how much information about the real data is allowed to survive into the synthetic data**.
This is a privacy question and `synpmx` offers four options.  The main deliverable of the package uses an
implementation of the AVATAR method which offers some blinding, but not formal privacy guarantees.

## Installation

``` r
# install.packages("pak")
pak::pak("iamstein/synpmx")
```

AVATAR needs nothing beyond base R. The other methods 
additionally require the official [OpenDP R package](https://docs.opendp.org/en/stable/api/r/):

``` r
install.packages("opendp", repos = "https://opendp.r-universe.dev")
```

## A first synthetic dataset

The default AVATAR synthetic data algorithm needs only the data and a declaration of what the columns mean.

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
#>   ID TIME       DV    AMT EVID CMT       WT
#> 1 13 0.00 0.000000 267.84  101   1 70.30935
#> 2 13 0.00 0.000000   0.00    0   2 70.30935
#> 3 13 0.30 2.429472   0.00    0   2 70.30935
#> 4 13 0.63 3.706058   0.00    0   2 70.30935
```

The output dataset keeps the same structure.  

The other three modes for generating synthetic data cover the scenarios where formal data privacy requirements are needed.  They are included in the package 
as part of a conceptual framework to help the modeler think about various
options that are available, but these methods have not been audited
for their ability to protect data privacy.  

| Mode | Function | Output built from | Guarantee | Works at |
|----|----|----|----|----|
| **1. AVATAR blending** | `synpmx_avatar()` | Real subject templates and blended real trajectories | None; governance only | At least 5 subjects |
| **2. Prior only** | `synpmx_prior()` | A public model and protocol only | `epsilon = 0` (no data read) | Any (data-independent) |
| **3. Calibration** | `synpmx_calibrated()` | A public model, magnitude corrected by 2 private releases | `(epsilon, delta)` DP | At least 20 subjects |
| **4. Empirical** | `synpmx_empirical()` | Dozens of noised population summaries | `(epsilon, delta)` DP | At least hundreds |

**The trust boundary decides whether you need differential privacy.**

If the generated data will not be accessible to anyone who cannot access the original data, formal privacy guarantees are not needed and the AVATAR blending approach is the recommended method because it is the simplest method to use.  

On the other hand, if the synthetic data will reach those who do not have access to the original data, then  more formal methods with mathematical trust guarantees are the appropriate methods of choice.  

`vignette("synpmx-method")` runs all four methods on the same dataset and shows the results side by side.

## Maintenance status

**AVATAR blending is the primary, maintained code.** It has no
dependencies beyond base R, and is what to reach for when the output
stays within the source data’s own access controls and obligations.

The three differentially private modes (**prior**, **calibration**,
**empirical**) are secondary.  They are present in this repository
because formal privacy is needed when data crosses a trust
boundary. They are provided as-is: not under active development and **not independently privacy-audited**.
Treat them as a principled demonstration of the privacy/utility
tradeoff, not as a production ready. That status is enforced in that `synpmx_calibrated()` and `synpmx_empirical()` refuse to run until `synpmx_enable_dp_engines()` has been called once in the session.

## Key Documentation

| Document | Question it answers |
|----|----|
| [The four generation modes](https://iamstein.github.io/synpmx/articles/synpmx-method.html) | What are the modes, and which one do I want? **Start here.** |
| [Using synpmx AVATAR with 5 datasets](https://iamstein.github.io/synpmx/articles/synpmx-demo.html) | How do I run this on my own study? |
| [The AVATAR Algorithm](https://iamstein.github.io/synpmx/articles/avatar-algorithm.html) | How does the default generator work, step by step? |

## License

[MIT](LICENSE.md) © 2026 Andrew Stein.
