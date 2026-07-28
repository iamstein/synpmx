# synpmx

📖 **Website and documentation: <https://iamstein.github.io/synpmx/>**

`synpmx` builds **synthetic pharmacometric datasets** from actual datasets.

## Which job do you want synthetic data to do?

There are many reasons to generate "synthetic data."  It is worth being explicit about which one you are after, because the goal decides which method is appropriat.

| The goal | Served here? |
|----|----|
| **Build the code.** You need data with the right *shape* — schema, event grammar, covariates, dosing and sampling pattern — to develop data-assembly code, diagnostic plots, and model-run plumbing outside the environment that holds the real study. Only the structure has to be right. | **Yes — this is what the package is for.** |
| **Use better tools.** The restricted environment forbids, or has not yet approved, the tooling you want to work with. Synthetic data moves the development work somewhere those tools are allowed. | **Yes**, same job as above. |
| **Send data past a trust boundary.** The output will reach people who cannot see the real data: a partner, a publication, a public repository. | **Only with care.** This needs a formal guarantee; see the privacy modes below, which are illustrative rather than audited. |
| **Answer the scientific question.** Estimate parameters, select a model, quantify a covariate effect, choose a dose, or stand in for real patients as a synthetic control arm. | **No.** No method here — see below. |
| **Teach and compare.** Show what the different generation methods actually do and what each one costs. | **Yes**, secondarily; that is why the non-default modes ship. |

The "answer the scientific question" row is the one worth stating plainly:
**synthetic data contains no
information about the drug that was not already in its source.** Analyze it and
you recover, at best, a noisier version of what the real data would have told
you, and at worst your own modeling assumptions handed back to you. A generated
dataset is also not anonymous data, and is not evidence that source subjects
cannot be re-identified. This package aims at structural usefulness, not
scientific equivalence.

The concrete use case it was built for: sharing realistic-looking study data
outside the GxP (Good Practice regulated) computing environment but still within
the organization, so that development work can happen without the real data. In
some cases that environment does not permit the most advanced agentic coding
tools, because of the risk of misalignment or unintended agent behavior. Working
against synthetic data lets those tools be used without exposing them to patient
data.

## How much of the real data survives?

Whatever the goal, the design question is
**how much information about the real data is allowed to survive into the synthetic data**.
That is a privacy question before it is a technical one, and `synpmx` offers
four answers.

The main deliverable of the package is an implementation of the AVATAR method,
which offers some blinding but no formal privacy guarantee. AVATAR works by
blending together actual patient profiles, and requires no model to be
specified. The package also provides code for other methods: trial simulation
from prior knowledge, and two differential privacy (DP) methods that give more formal
protection. Those are provided mainly to illustrate the tradeoffs between ways
of generating synthetic data, and are not actively maintained.

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
