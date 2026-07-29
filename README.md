# synpmx
📖 **Website and documentation: <https://iamstein.github.io/synpmx/>**

`synpmx` builds **synthetic pharmacometric datasets** from actual datasets.

## Which job do you want synthetic data to do?

There are many reasons to generate "synthetic data."  It is important to be be explicit about your use case because the use case determines the appropriate data generation algorithm.  

| Use case | Served here? |
|----|----|
| **Develop code.** You need data with the right *shape* — schema, event grammar, covariates, dosing and sampling pattern.  You'll use this data to develop code for data processing, diagnostics, and model building outside the environment that holds the real study.  | **Yes — this is what the package is for.** |
| **Send data past a trust boundary.** The output will reach people who cannot see the real data: a partner, a publication, a public repository. | **Only with care.** This needs a formal guarantee; see the privacy modes below, which are illustrative rather than audited. |
| **Answer the scientific question.** Estimate parameters, select a model, quantify a covariate effect, choose a dose, or stand in for real patients as a synthetic control arm. | **No.**   Use the real data for this. |
| **Teach and compare.** Show what the different synthetic data generation methods do. | **Yes**, secondarily; that is why the non-default modes ship. |

The use case this package was built for: sharing realistic-looking study data
outside the GxP (Good Practice regulated) computing environment but still within
the organization, so that code development can happen without the real data. In
some cases the GxP environment does not permit the most advanced agentic coding
tools, because of the risk of misalignment or unintended agent behavior. Working
with synthetic data lets those tools be used without exposing them to patient
data.

## The main deliverable: AVATAR method

The main deliverable of the package is an implementation of the AVATAR method
[1, 2], which offers some masking but no formal privacy guarantee. AVATAR works
by blending together patient profiles, and requires no model to be
specified. "AVATAR" is a method name rather than an initialism, from the
patient-centric *avatarization* literature: the original method was developed in
Guillaudeux and colleagues [2], and Destere and colleagues benchmark a modified
AVATAR for population PK [1]. Destere et al [1], did not test their method using
pharmacometrics dataset from actual clinical trials.  We have done here, and in the process
we have added features such as handling BLOQ data and masking patients with unique dosing schedules.  The key function is `synpmx_avatar()`.

The package also provides code for other data masking methods: trial simulation
from prior knowledge, and two differential privacy (DP) methods that give more formal
protection. These methods are provided mainly to illustrate the tradeoffs between ways
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
#> 1 13 0.0000000 0.000000 267.84  101   1 69.80615
#> 2 13 0.0000000 0.000000   0.00    0   2 69.80615
#> 3 13 0.2816667 2.861275   0.00    0   2 69.80615
#> 4 13 0.5416667 3.600730   0.00    0   2 69.80615
```

The output dataset keeps the same structure.  

| Mode | Function | Output built from | Guarantee | Works at |
|----|----|----|----|----|
| **1. AVATAR blending** | `synpmx_avatar()` | Real subject templates and blended real trajectories | None; governance only | At least 5 subjects |
| **2. Prior only** | `synpmx_prior()` | A public model and protocol only | `epsilon = 0` (no data read) | Any (data-independent) |
| **3. Calibration** | `synpmx_calibrated()` | A public model, magnitude corrected by 2 private releases | `(epsilon, delta)` DP | At least 20 subjects |
| **4. Empirical** | `synpmx_empirical()` | Dozens of noised population summaries | `(epsilon, delta)` DP | At least hundreds |

**The trust boundary decides whether you need differential privacy.**

If the generated data will not be accessible to anyone who cannot access the original data, formal privacy guarantees are not needed and the AVATAR blending approach is the recommended method because it is the simplest method to use as it doesn't require the specification of a model.  

On the other hand, if the synthetic data will reach those who do not have access to the original data, then  more formal methods with mathematical trust guarantees are the appropriate methods of choice.  

`vignette("synpmx-method")` runs all four methods on the same dataset and shows the results side by side.

## Maintenance status

**AVATAR blending is the primary, maintained code.** It has no
dependencies beyond base R, and is what to reach for when the output
stays within the source data’s own access controls and obligations.  However,
AVATAR does not offer any formal, mathematical guarantees around privacy.  

The three other modes (**prior**, **calibration**,
**empirical**) are secondary; but, they are present in this repository
because they cover scenarios where data crosses a trust
boundary and formal privacy conditions must be met. 
The methods are provided as-is: not under active development and **not independently privacy-audited**.
Treat them as a principled demonstration of the privacy/utility
tradeoff, not as a production ready. That status is enforced in that `synpmx_calibrated()` and `synpmx_empirical()` refuse to run until `synpmx_enable_dp_engines()` has been called once in the session.

## Key Documentation

| Document | Question it answers |
|----|----|
| [The four synthetic generation modes](https://iamstein.github.io/synpmx/articles/synpmx-method.html) | What are the modes, and which one do I want? **Start here.** |
| [Demo; Using synpmx AVATAR with 5 datasets](https://iamstein.github.io/synpmx/articles/synpmx-demo.html) | How do I run this on my own study? |
| [The AVATAR Algorithm](https://iamstein.github.io/synpmx/articles/avatar-algorithm.html) | How does the default generator work, step by step? |

## References

1. Destere A, Lombardi R, Labriffe M, et al. *Can synthetic data overcome the
   privacy and fidelity bottleneck in Pharmacometrics? A comparative benchmark
   using a daptomycin population pharmacokinetic model.* medRxiv preprint,
   posted June 2, 2026. doi:
   [10.64898/2026.05.30.26354512](https://doi.org/10.64898/2026.05.30.26354512).

2. Guillaudeux M, Rousseau O, Petot J, et al. Patient-centric synthetic data
   generation, no reason to risk re-identification in biomedical data analysis.
   *npj Digital Medicine.* 2023;6.
   doi: [10.1038/s41746-023-00771-5](https://doi.org/10.1038/s41746-023-00771-5).

## License

[MIT](LICENSE.md) © 2026 Andrew Stein.
