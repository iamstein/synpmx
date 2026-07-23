# synpmx

`synpmx` generates structurally faithful **synthetic data for model-workflow
exploration**. Its primary `synthesize_pmx()` API uses AVATAR-style profile
blending inside a trusted environment; the calibrated structural/DP workflow
remains available when formal privacy is required.

The primary architecture is deliberately simple: blend complete longitudinal
profiles while preserving event-table structure. For formal privacy, the
package also supports a fit-once, generate-many workflow:

1. `fit_private_pmx()` runs inside the restricted source environment, bounds
   each complete subject contribution, and releases only OpenDP-noised
   population summaries plus a privacy ledger.
2. `generate_pmx()` reads only that fitted model. Repeated datasets are
   post-processing and consume no additional privacy budget.

Generated data can exercise cleaning, joins, reshaping, plots, control-file
plumbing, repeated-dose PK code, longitudinal PD/biomarker code, infusion
events, and censoring conventions. It aims for broad magnitude and shape—not
source distributions, parameter estimates, covariate-response relationships,
scientific fidelity, inference, model selection, or clinical conclusions.

## Documentation map

Start with the introduction vignette; it applies all four generation modes to
one dataset and says which is appropriate where.

| Vignette | Question it answers | Read it when |
|---|---|---|
| [`synpmx-intro`](vignettes/synpmx-intro.Rmd) | What are the four generation modes, and which one do I want? | **Start here.** |
| [`synpmx-demo`](vignettes/synpmx-demo.Rmd) | How do I actually run this on my data? | You have picked a mode and want the worked workflow across five public datasets. |
| [`synpmx-simulation-method`](vignettes/synpmx-simulation-method.Rmd) | How does the default AVATAR generator work, step by step? | You need to defend or debug what the generator did. |
| [`synpmx-privacy-intro`](vignettes/synpmx-privacy-intro.Rmd) | What does differential privacy guarantee, and does my release need it? | The generated data might cross a trust boundary. |
| [`synpmx-epsilon-exploration`](vignettes/synpmx-epsilon-exploration.Rmd) | What does a given epsilon cost me in accuracy? | You are choosing an epsilon and cohort size. |

The design documents in `design/` are the internal record behind those
vignettes: `TODO.md` is the working queue and indexes the rest. `AGENTS.md`
holds the repository conventions.

## Privacy contract

The production claim is:

> Generated from a subject-level `(epsilon, delta)`-differentially private model.

One subject's complete longitudinal record is the privacy unit. Neighboring
datasets differ by adding or removing that complete subject. Epsilon is the
one-person influence limit: smaller is stronger, and no universal default is
chosen. Delta is a very small additive allowance in the probability bound; it
is not a re-identification probability or fraction of unprotected patients.

Privacy is mathematically bounded, not absolute. Differential privacy does not
guarantee impossibility of linkage or re-identification, establish legal
anonymity, authorize release, secure a compromised environment, or validate
public-input claims. Independent privacy, legal, information-security, and
data-governance review remains required.

## Production dependency

Fitting confidential data requires the official [OpenDP R
package](https://docs.opendp.org/en/stable/api/r/). The package fails closed if
OpenDP is unavailable and never falls back to hand-written or ordinary R noise.

```r
install.packages("opendp", repos = "https://opendp.r-universe.dev")
dp_backend_status()
run_dp_backend_tests()  # canonical adapter checks; requires OpenDP
```

An explicitly noiseless `backend = "public"` exists only for fully public
fixtures and requires `public_source = TRUE`. Its privacy report makes no DP
claim; it must never process confidential data.

## Public API

- `pmx_roles()` declares PMX semantics, including nominal time, TAD, occasion,
  CENS/LIMIT, ADDL/II, baseline covariates, treatment-like subject properties,
  occasion-assigned dose, and explicit exclusions.
- `pmx_endpoint()` declares each DVID's dose-relative, study-time, occasion, or
  hybrid scientific clock.
- `pmx_bounds()`, `pmx_schema()`, `pmx_public_design()`,
  `pmx_contribution_limits()`, and `pmx_budget_allocation()` make proof-relevant
  public inputs explicit.
- `fit_private_pmx()` is the only confidential-data stage.
- `generate_pmx()` constructs new event tables from the fitted model and, by
  default, uses its privacy-accounted subject-count release.
- `subject_property_summary()` reports released treatment/property strata and
  their jointly fitted generalized regimens.
- `privacy_report()` and `validate_private_model()` expose accounting,
  assumptions, the release ledger, and leakage guards.
- `validate_pmx()` checks generated PMX structure and censoring coherence.
- `compare_pmx()` is a restricted diagnostic and marks every source-derived
  component `restricted_not_releasable`.

## Minimal shape of a confidential workflow

```r
roles <- pmx_roles(
  id = "ID", time = "TIME", nominal_time = "NTIME", tad = "TAD",
  occasion = "OCC", dv = "DV", amt = "AMT", evid = "EVID",
  cmt = "CMT", dvid = "DVID", mdv = "MDV", rate = "RATE",
  cens = "CENS", limit = "LIMIT", covariates = "WT"
)

endpoints <- list(
  cp = pmx_endpoint(
    dvid = "cp", alignment = "dose_relative",
    transform = "log", shape = "occasion", cmt = 2
  ),
  response = pmx_endpoint(
    dvid = "response", alignment = "study_time",
    transform = "auto", shape = "global", cmt = 3
  )
)

# Bounds, schema, levels, grids, and protocol values must be justified public
# inputs, not exact extrema silently derived from confidential patients.
private_model <- fit_private_pmx(
  confidential_data, roles, endpoints,
  epsilon = approved_epsilon,
  delta = approved_delta,
  delta_justification = approved_delta_justification,
  bounds = approved_bounds,
  public_design = approved_public_design,
  contribution_limits = approved_contribution_limits,
  budget_allocation = approved_budget,
  backend = "opendp"
)

privacy_report(private_model)
synthetic_1 <- generate_pmx(private_model, seed = 101)
synthetic_2 <- generate_pmx(private_model, seed = 202)
validate_pmx(synthetic_1, roles, endpoints, strict = TRUE)
```

There is intentionally no fitting seed. OpenDP controls private mechanism
randomness; ordinary seeds control only generation from an already released
model. Omitted `n_subjects` uses the fitted noisy count release (or the exact
count for an explicitly public fixture); pass it explicitly only when a
different public workflow cohort size is intended.

## Endpoint and event behavior

- Dose-relative endpoints use a small private TAD curve and create a new
  excursion after every generated dose.
- Study-time endpoints use one global curve and do not restart after a dose.
- Occasion endpoints use related within-occasion profiles.
- Hybrid endpoints combine a global baseline with a small dose-relative
  excursion.
- Multiple DVIDs are learned and generated separately.
- Actual-like times are generated around generalized nominal cells; tied
  collection blocks share jitter, TAD/occasion are re-derived, and the released
  total observation count limits how densely endpoint grids are instantiated.
- For dose-relative endpoints, the private timing release separately learns
  each occasion's sampling probability and its bounded observation count
  conditional on being sampled. Generation can therefore distinguish an
  uncommon dense visit from a sparse visit in every subject, without copying a
  source visit vector. `sampling_summary()` exposes these fitted quantities as
  releasable post-processing.
- By default, dose count, interval, amount, infusion behavior, occasion
  activation, conditional sample counts, and timing-cell occupancy are inferred
  from the input through budgeted summaries. A generic grid is only a
  discretization basis, not a supplied sampling schedule.
- A public `endpoint_occasion_grids` schedule remains available only as an
  exceptional override when the protocol is independently public; the package
  demonstrations do not use it.
- Factor-valued ID columns retain the factor class but never retain source ID
  levels; generated IDs receive a fresh synthetic-only level set.
- Dose and infusion fields are created coherently. A generated infusion start
  and negative stop share the generated amount/rate and duration.
- A declared `subject_properties` field such as ACTARM, TRT, or nominal dose
  group is generated jointly with its property-conditioned regimen. Its finite
  category domain must be independently public. A declared `assigned_dose`
  field is reconstructed from positive AMT and held constant within
  subject/occasion.
- Numeric-coded covariates are modeled categorically when their public levels
  are supplied in `category_levels`.
- Censoring is applied to a generated latent value, then DV, CENS, and LIMIT are
  reconstructed together under Monolix-style conventions.
- Source IDs, raw rows, complete profiles, schedules, residuals, and unnoised
  aggregates are absent from the fitted model.

## Public examples

The practical vignette and `scripts/demo_nlmixr2data.R` exercise:

- `nlmixr2data::theo_md`: the privacy-accounted event/timing fit discovers the
  seven-dose Q24H regimen, dense first/final profiles, sparse occasion-2
  sampling, and no observations after doses 3--6;
- `nlmixr2data::warfarin`: lower-case schema, factor preservation, separate
  dose-relative `cp` and global study-time `pca`;
- `nlmixr2data::wbcSim`: coherent infusion starts/stops, generalized follow-up,
  and delayed decline/nadir/recovery without reproducing the singleton
  multi-thousand-hour regimen;
- `nlmixr2data::nimoData`: four nominal-dose subject properties linked to ten
  inferred weekly infusions, declared OCC/TAD, and terminal washout follow-up;
  and
- `nlmixr2data::mavoglurant`: reset occasion clocks, numeric categorical SEX,
  and an occasion-assigned DOSE kept coherent with generated AMT.

`pmx_censoring_fixture()` supplies a fully simulated public example with
uncensored, left-censored, right-censored, and interval-censored records.
`pmx_simulated_fixture()` supplies a deterministic, two-endpoint repeated-dose
study with 60 subjects by default for broader privacy-utility evaluation.

## Important limitations

- Six- or twelve-subject studies can satisfy the same formal privacy definition
  but may yield very noisy summaries. The package warns when private count and
  requested dimensionality imply weak utility.
- Broad generated variability is intentionally public rather than precisely
  estimated from a small source study.
- Baseline covariates are marginal and subject-constant. Time-varying
  covariates, covariate correlations, and covariate-response relationships are
  not modeled; declare unsupported longitudinal fields in `exclude`.
- `assigned_dose` guarantees event-table consistency but does not by itself
  reconstruct a crossover-sequence distribution. Supply a genuine
  subject-level ACTARM/TRT/sequence property when one exists.
- Dose-relative AR(1) perturbations restart at each generated occasion. If the
  released coarse curve is already approximately unimodal, source-free
  post-processing prevents residual noise from adding a second PK peak.
- The current mechanism uses conservative L1 sensitivity proportional to each
  released vector's dimension and basic sequential composition.
- OpenDP's mechanism, the R/OpenDP boundary, serialization, public-input
  assertions, floating-point behavior, and side channels still require an
  independent specialist audit.
- Empirical attacks can discover bugs but cannot prove differential privacy.
- Generated data remain inappropriate for scientific analysis.

See the documentation map above for which vignette answers which question;
`vignette("synpmx-intro")` is the entry point.

## Development

```r
testthat::test_local()
roxygen2::roxygenise()
```

Run a source build and `R CMD check` after behavioral changes. The repository
keeps package functions in `R/`, tests in `tests/testthat/`, and runnable public
demonstrations in `scripts/`.

## License

[MIT](LICENSE.md) © 2026 Andrew Stein.
