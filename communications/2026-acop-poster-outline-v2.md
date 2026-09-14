# Enabling AI-Assisted Pharmacometric Workflows Using Synthetic Data

Andrew Stein, Alex Pogodaev

## Objective

`synpmx` is an R package for generating synthetic pharmacometric datasets
to support AI-assisted code development when clinical
data must remain in a secure computing environment that does not allow modern agentic AI coding tools. Synthetic data provides
the columns, dosing records and observations needed to develop exploration,
model-fitting and diagnostic code. Scientific estimation and interpretation
use the real data inside the secure environment.

## Workflow

```mermaid
%%{init: {"themeVariables": {"fontSize": "24px"}}}%%
flowchart LR
  subgraph secure["Secure clinical-data environment"]
    real[("Real clinical data")]
    generate["Generate synthetic data<br/>with synpmx"]
    review["Check synthetic data<br/>and review for release"]
    validate["Run code on real data<br/>Validate and interpret"]
    real --> generate --> review
    real --> validate
  end

  subgraph development["AI-enabled development environment"]
    synthetic[("Approved synthetic data")]
    coding["Develop and debug code<br/>with AI assistance"]
    synthetic --> coding
  end

  review -->|"Synthetic data out"| synthetic
  coding -->|"Analysis code back"| validate

  classDef data fill:#DBEAFE,stroke:#2563EB,color:#172554
  classDef process fill:#FFFFFF,stroke:#64748B,color:#0F172A
  classDef result fill:#DCFCE7,stroke:#15803D,color:#14532D
  class real,synthetic data
  class generate,review,coding process
  class validate result
  style secure fill:#F1F5F9,stroke:#64748B,color:#0F172A
  style development fill:#FFF7ED,stroke:#C2410C,color:#7C2D12
```

**Figure 1.** Workflow

1. **Within secure environment:** Declare column roles with `pmx_roles()`, generate a synthetic dataset, and compare its structure and distributions with the source.
2. **Review synthetic data:** Assess the synthetic dataset and attached information
   under the organization's requirements before transfer.
3. **Develop code in AI-enabled environment:** Develop and debug analysis code using the approved synthetic dataset.
4. **Return code to secure environment:** Run the code on real data, validate the analysis and interpret the results.

## Selected methods available in `synpmx`

| Function | Generation approach | Required specification  | Formal privacy | Source patients | Reportable fingerprint |
|---|---|---|---|---|---|
| `synpmx_prior()` | Simulate a public model and trial design |      Public model,       Study design | Yes  | 0 | Yes |
| `synpmx_model()` | Fit longitudinal PK and PD models; simulate new subjects    | Nominal times;    | No | 10s | Yes |
| `synpmx_avatar()` | blending of similar patient profiles | Nominal times        | No | 10s–100s | No |
|  |

Additional methods available, including some methods with formal differential privacy guarantees, but for datasets with less than 100 patients, the noise needed for formal privacy guarantees is large.

## Preferred generation algorithm `synpmx_model()`

**Figure 2:** Overview of `synpmx_model`

| Component | What the generator does |
|---|---|
| Privacy Protection  | Generate new IDs; Drop all columns without defined roles; Drop cohorts with fewer than 3 patients; Resample categorical covariates that fewer than 3 patients take |
Dosing | For each cohort, fits hazard model for missed doses, reduced doses, and discontinuation.  Covariate-based dosing (e.g. weight-based) can be specified.
| Observations | Estimate rates of missed visits for each endpoint and at at each visit; simulate missingness. |
| Population PK | Fit a one-compartment model using `nlmixr2`.  Simulate PK from this model.  Alternative model candidates can be specified. |
| PD | Fit constant, linear or exponential time courses to continuous endpoints, with subject variability and residual error. The default pools arms.  If capturing dose-response is desired, `pd_by_arm = TRUE` allows arm-specific curves. Discrete endpoints are sampled from distribution of all values taken |
| LOQ | Uses or estimates LOQ for all continuous observations |

## Worked Example: `xgxr::mad`

The publicly available multiple-ascending-dose example in `xgxr` combines a PK concentration with continuous, ordinal, count and binary PD endpoints.

![Source and synthetic PK profiles by treatment arm](2026-acop-poster-figures/mad-pk-profiles.png)

**Figure 3a. Pharmacokinetic (PK) profiles over the first dosing interval.**
Source data vs synthetic data.  Because one compartment model used by default, there is some mismatch between source and synthetic.  If greater fidelity to the data is desired, alterative models can be fit.

![Source and synthetic continuous PD profiles by treatment arm](2026-acop-poster-figures/mad-pd-profiles.png)

**Figure 3b. Continuous pharmacodynamic (PD) profiles over study time.**
Source data vs synthetic data.  Because by default all cohorts are pooled together, the dose-response relationship is lost.  If better fidelity is desired, one can fit a PD profile to each cohort separately.  

![Source and synthetic distributions of continuous PD, PK, baseline weight and sex](2026-acop-poster-figures/mad-distributions.png)

**Figure 3c. Pooled distributions.** Distributions of observations and covariates provide one graphical tool for assessing the synthetic dataset.

`synpmx` also includes a scorecard that checks dataset structure, distribution
changes and potential patient copying, helping identify synthetic outputs
that need review. [Scorecard documentation](https://iamstein.github.io/synpmx/articles/avatar-scorecard.html).

## Evaluation across study designs

Ten publicly available dataset examples span controlled dosing, routine care and simulated
studies. The datasets below are covered in the
[model evaluation article](https://iamstein.github.io/synpmx/articles/pmxmodel-public-data-examples.html).

| Dataset | Patients | Design feature |
|---|---:|---|
| `case1_pkpd` | 180  (simulated)  | Multiple arms; censored concentrations |
| `mad` | 60 (simulated) | Multiple doses; mixed endpoint types (continuous, categorical) |
| `warfarin` | 32 (actual) | Single oral dose; delayed response |
| `wbcSim` | 45 (actual) | Infusions; white-cell nadir and recovery |
| `mavoglurant` | 120 (actual)| Repeated occasions; resetting clock |
| `theo_md` | 12 (actual)|  Repeated oral doses; small cohort |
| `nimoData` | 12 (actual) | Weekly infusions; small cohort |
| `pheno_sd` | 59 (actual) | Neonatal care; sparse, irregular sampling |
| `mixroute_sim` | 90 (simulated) | Intravenous and subcutaneous dosing |
| `onc_sim` | 200 (simulated) | Trough sampling; tumour response; dose changes |

## Availability

R package prototype, source code, worked examples and evaluation methods:
[iamstein.github.io/synpmx](https://iamstein.github.io/synpmx/).

**Footer:** Add a scannable link to the package website and the approved
poster materials. Internal discussions on use of this package are ongoing.

---

## Preparation Notes — Not Poster Copy

### Poster Layout

Three columns: workflow and available methods on the left; the fitted-model
generator in the centre; the worked example and evaluation on the right.
Give the example the most space. Target approximately 650–800 words of visible copy, including tables and captions. Figure instructions and the preparation notes are working material, not poster copy.

### Workflow Figure Typography

Use 24–28 pt labels at the final printed poster size, including arrow labels,
with larger environment headings. Keep labels short and enlarge the diagram's
allocated space if needed; do not shrink the lettering to fit the column.
The Mermaid preview uses a 24 px base font; check the physical text size again
after placing the figure on the poster.

### Example Generation and Verification

Run [2026-acop-poster-figures.R](2026-acop-poster-figures.R) from the
repository root:

```sh
Rscript communications/2026-acop-poster-figures.R
```

The script loads the working-tree package, reads public `xgxr::mad` data and
its stored fit, and generates synthetic subjects with seed 909. It does not
refit the population model. Required R packages are `devtools`, `xgxr`,
`ggplot2`, `patchwork`, and the package's own dependencies.

Figures are saved beside the script in `2026-acop-poster-figures/` as PNG
previews and vector PDF files. Each main figure is 18 inches wide, with large
labels for poster placement. Use the PDF for scaling; check label sizes if
reducing its width.

- [PK profiles, PDF](2026-acop-poster-figures/mad-pk-profiles.pdf)
- [Continuous PD profiles, PDF](2026-acop-poster-figures/mad-pd-profiles.pdf)
- [Compact distributions, PDF](2026-acop-poster-figures/mad-distributions.pdf)
- [All endpoint and covariate distributions, PDF](2026-acop-poster-figures/mad-distributions-full.pdf)
  — supporting figure, outside the main poster.

The same directory holds aggregate distribution tables, the scorecard and
`run-info.txt` with the seed, fit and source-code checksums, and R session
information. Scorecard output stays in the supporting evaluation.
An optional output directory and seed can be passed as the first and second
arguments. Inputs are fixed to the public example; this script accepts no
internal-study data.

### Evaluation Table Generation

Run [2026-acop-poster-evaluation.R](2026-acop-poster-evaluation.R) to regenerate
the public-data table:

```sh
Rscript communications/2026-acop-poster-evaluation.R
```

The script uses the evaluation article's dataset preparation, roles, stored
fits and generation seeds. It saves a Markdown table and aggregate CSV,
supporting scorecards and run provenance in `2026-acop-poster-evaluation/`.
Patient counts are computed from the public source datasets.

Add internal-study rows using the same three columns: **Dataset**, **Patients**,
**Design feature**. Describe the design briefly. The script regenerates only the public table
and does not overwrite the outline or any internal rows added to it.

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
