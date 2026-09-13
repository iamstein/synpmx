# Enabling AI-Assisted Pharmacometric Workflows Using Synthetic Data

Andrew Stein, Alex Pogodaev

## Poster Layout

Three columns: workflow and available methods on the left; the fitted-model
generator in the centre; the worked example and evaluation on the right.
Give the example the most space. Target approximately 650–800 words of visible
copy, including tables and captions. Figure instructions and the preparation
notes below are working material, not poster copy.

## Objective

`synpmx` is an R package for generating synthetic pharmacometric event tables
to support artificial intelligence (AI)-assisted code development when clinical
data must remain in a secure computing environment. Synthetic data provides
the columns, dosing records and observations needed to develop exploration,
model-fitting and diagnostic code. Scientific estimation and interpretation
use the real data inside the secure environment.

## Workflow

**Figure 1: Draw two labelled environments, with arrows for synthetic data
leaving and code returning.**

1. **Secure environment:** Declare column roles with `pmx_roles()`, generate a
   synthetic dataset, and compare its structure and distributions with the source.
2. **Release review:** Assess the synthetic dataset and attached information
   under the organization's requirements before transfer.
3. **AI-enabled environment:** Develop and debug analysis code using the
   approved synthetic dataset.
4. **Secure environment:** Run the returned code on real data, validate the
   analysis and interpret the results.

The fitted-model generator has no formal privacy guarantee. A passing
scorecard does not establish that an output is safe to release.

## Methods Available in `synpmx`

Differential privacy (DP) bounds how much an individual's participation can
affect a release. Its privacy budget controls the strength of that bound.

| Function | Generation approach | Required specification beyond column roles | Formal privacy |
|---|---|---|---|
| `synpmx_model()` | Fit population pharmacokinetic (PK) models and pharmacodynamic (PD) time courses; simulate new subjects | Nominal times; endpoint and route declarations where inference is ambiguous | No |
| `synpmx_avatar()` | AVATAR-style blending of neighbouring patient profiles | Role and masking choices | No |
| `synpmx_pca()` | Principal component analysis (PCA): fit shared profile patterns and draw new coefficients | Nominal times | No |
| `synpmx_prior()` | Simulate a public model and trial design | Public model, parameters and design | Reads no protected study data, provided inputs are independent of it |
| `synpmx_calibrated()` | Privately correct the magnitude of a public model | Public model, design, correction bounds and budget | DP mechanism; experimental implementation* |
| `synpmx_empirical()` | Reconstruct from noisy study summaries | Public bounds, contribution limits and budget allocation | DP mechanism; experimental implementation* |

*The DP engines require the OpenDP privacy backend for their DP claim, have
known open findings and have not been independently privacy-audited. The public
fixture backend adds no privacy noise. The empirical engine is not under active
development.*

## Fitted-Model Generation

**Figure 2: Four model boxes feeding a synthetic event table. Draw the dose
schedule feeding the PK simulation explicitly.**

| Component | What the generator does |
|---|---|
| Dosing | Summarizes planned amounts and times by arm, including scheduled escalation; draws reductions, skipped cycles and discontinuation. |
| Observations | Estimates attendance by endpoint and nominal visit; draws which observations each synthetic subject contributes. |
| Population PK | Fits a one-compartment model by default using `nlmixr2`, an R population-model fitting tool. Draws subject parameters and simulates concentrations against each generated dosing history. Additional model candidates can be requested. |
| PD | Fits constant, linear or exponential time courses to continuous endpoints, with subject variability and residual error. The default pools arms; `pd_by_arm = TRUE` allows arm-specific curves. Discrete endpoints are drawn from arm-and-visit frequencies. |

`synpmx_model_estimate()` builds a reusable fitted object;
`synpmx_model_generate()` draws datasets from that object without rereading
patient records. `synpmx_model()` combines both stages.

Declared below-limit PK observations enter the fit as censored observations,
and generated values are censored at the output boundary. Output uses declared
column names and classes, with new subject identifiers. Columns outside the
role declaration are dropped.

The PD time courses have no exposure term. Dose reductions affect simulated
PK, but the generator does not model a patient's response causing a subsequent
dose change. Sparse sampling and small cohorts can produce weak fits.

## Worked Example: `xgxr::mad`

The publicly available multiple-ascending-dose example in `xgxr` combines a PK
concentration with continuous, ordinal, count and binary PD endpoints. It
exercises a shared event-table schema across different observation types.

**Figure 3: Source versus synthetic data, with matching axes.** Use PK profiles
over one dosing interval and continuous PD over study time, grouped by arm.
Place a small distribution comparison beneath them. Give each panel a caption
describing the observed agreement or discrepancy after inspecting the run.

**Results strip:** Display the scorecard verdict tally, including unanswered
checks, plus endpoint retention, observation and dose counts per subject, and
the largest change in standard deviation. Obtain these from
`synpmx_scorecard()` and `compare_pmx_distributions()` for the plotted run.

**Workflow demonstration:** Develop one concentration-time plotting function
against the synthetic table, then execute the same function on the public source
table. Record any changes needed. This demonstration still needs to be performed;
successful data generation alone does not establish successful code transfer.

## Evaluation Across Study Designs

**Compact results table:** One row per dataset, with a design feature, counts
of failed/review/unanswered checks, and the main observed limitation. Populate
the outcomes from a fresh run of the package evaluation.

The existing model survey covers ten publicly available examples:
`case1_pkpd`, `mad`, `warfarin`, `wbcSim`, `mavoglurant`, `theo_md`, `nimoData`,
`pheno_sd`, `mixroute_sim` and `onc_sim`. These include simulated datasets;
public availability and simulation origin are overlapping categories.

Select rows illustrating mixed administration routes, sparse sampling,
censoring and dose changes. Link to the complete survey. Add an aggregate
internal-study evaluation only after its outcomes and permissible wording have
been confirmed.

## Availability

R package prototype, source code, worked examples and evaluation methods:
[iamstein.github.io/synpmx](https://iamstein.github.io/synpmx/).

**Footer:** Add a scannable link to the package website and the approved
poster materials. Internal discussions on use are ongoing.

---

## Preparation Notes — Not Poster Copy

### Example Generation and Verification

Use the existing stored `mad` fit for the first figures. The following code
uses the same role declaration as the public-data survey. The fit was estimated
previously; this code runs generation and evaluation only.

```r
library(synpmx)
source_data <- as.data.frame(get(utils::data(list = "mad", package = "xgxr")))
roles <- pmx_roles(
  id = "ID", time = "TIME", dv = "LIDV", amt = "AMT", evid = "EVID",
  cmt = "CMT", dvid = "NAME", mdv = "MDV", nominal_time = "NOMTIME",
  strata = c("TRTACT", "DOSE"), covariates = c("WEIGHTB", "SEX")
)
fit <- readRDS(system.file("extdata", "mad-model-fit.rds", package = "synpmx"))
synthetic <- synpmx_model_generate(fit, seed = 909)
card <- synpmx_scorecard(source_data, synthetic, roles)
verdicts <- table(factor(
  card$verdict, levels = c("pass", "review", "FAIL", "not applicable")
))
verdicts
card[card$check %in% c("A3", "A5a", "A5b", "D1"), ]
compare_pmx_distributions(source_data, synthetic, roles)
```

Generation and the scorecard were rerun against the working-tree R code while
preparing this outline. Endpoint retention passed; the spread comparison needs
review, and several checks are not applicable to this generator. Use the live
output for numbers and name the unanswered checks in the linked evaluation.
The population fit was not re-estimated in this preparation run.

### Feature Claims

| Topic from the initial outline | Wording supported by the implementation |
|---|---|
| Below the limit of quantification (BLOQ) | Declared PK censoring is passed to the fitting likelihood; output censoring is applied after simulation. Do not promise identical censored fractions by arm. |
| Weight-based dosing | Do not claim automatic preservation of each subject's dose–weight relationship. The model generator draws arm-level schedules independently of generated covariates. Optional allometric effects on PK parameters are a separate feature and default to off. |
| Event handling | Describe declared event-table fields and supported routes concretely. Regression coverage includes routes and dosing histories; avoid an unrestricted claim that every input event convention survives. |
| Extreme outliers | Do not advertise general outlier removal for this generator. It filters sparsely supported arms/visits and applies measurement boundaries; that is not a general filter on extreme times, doses and observations. |
| Multiple drugs and routes | Endpoint-specific models and dosing assignments exist. State which configuration is demonstrated in the evaluation, rather than implying a joint parent–metabolite or interacting-drug model. |
| Internal datasets | The original outline reports three complex internal studies. Their outcomes were not verified here. Keep source-derived plots and tables out of public poster materials. |

### Package Sources for Figures and Technical Review

- [Model algorithm](../vignettes/pmxmodel-algorithm.Rmd): model choices,
  declarations, dosing and attendance, and limitations.
- [Estimation implementation](../R/model-estimate.R),
  [generation implementation](../R/model-generate.R) and
  [dose/visit implementation](../R/dose-visit-models.R): checked for the
  behavioural descriptions above.
- [Public-data evaluation](../vignettes/articles/pmxmodel-public-data-examples.Rmd):
  source preparation, stored fits and cross-dataset reporting. Recompute its
  results rather than copying numerical claims from its prose.
- [Scorecard implementation](../R/scorecard.R): checks and verdict meanings.
- [Model-generation regression tests](../tests/testthat/test-model-generate.R):
  focused generation checks; their existence is not a fresh test-suite result.

Before producing the poster, run the example and selected evaluation rows with
the same package revision and seeds, then write captions from those outputs.
Retain full tables and the method comparison online. Clear the public artifact
through the organization's internal review before submission.
