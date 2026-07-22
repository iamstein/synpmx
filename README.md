# pmxmock

`pmxmock` is an R package prototype for generating structurally faithful mock
pharmacometric event datasets from an existing dataset without fitting or
calling a PK, PD, or nonlinear mixed-effects model.

It is designed to create **mock data for model-workflow exploration**: run the
generator where restricted source data are available, then use the result to
develop data manipulation, diagnostics, plotting, and model-run plumbing in a
less restricted environment.

The output is not anonymous data. The package provides no formal anonymization,
differential privacy, re-identification risk bound, or guarantee of scientific
fidelity. Do not use mock output for estimation, inference, model selection, or
clinical decisions.

## How it works

`mock_pmx()` uses an AVATAR-inspired method adapted to longitudinal PMX data. It
is an AVATAR-like variant, not an exact reproduction of published AVATAR
software.

For each generated subject, the prototype:

1. separates baseline covariates, longitudinal measurements, and event-control
   fields using an explicit `pmx_roles()` mapping;
2. forms source-subject profiles from baseline covariates and endpoint-specific
   trajectories on common grids;
3. median-imputes profile features only for distance calculation, standardizes
   them, and uses a rank-safe PCA representation when possible;
4. samples one source subject's complete event template, preserving row order,
   tied times, doses, infusion start/stop records, compartments, and endpoints;
5. finds compatible non-anchor neighbors and assigns capped randomized weights
   based on profile distance and randomized rank;
6. uses the same weights for baseline covariates and all endpoint trajectories;
7. interpolates and blends each endpoint separately, then adds a subject shift
   and modest AR(1) noise; and
8. restores source names, column order, practical classes, factor levels, and
   deterministic MDV logic when present.

For compatible donors, the raw weighting rule is
`Exp(1) / max(distance, epsilon) * 2^(-randomized_rank)`. Weights are normalized,
and a dominant donor is capped at 0.80 with the excess redistributed. If a
donor lacks an interpolated value, the available weights are renormalized for
that value only. Positive-like endpoints use an offset log transformation and
are constrained to be nonnegative after back-transformation; transform choices
are endpoint-specific and recorded in `attr(mock, "pmx_settings")`.

Event fields such as `EVID`, `AMT`, `RATE`, `CMT`, `DVID`, `ADDL`, and `II` are
not averaged, PCA-transformed, or generated independently. They remain a
coherent sampled template. Profile and donor times are aligned relative to each
subject's first positive dose within compatible schedule groups. When this
aligned-time interpolation cannot cover an irregular subject's observation
window, the prototype uses a normalized within-window fallback rather than
unbounded extrapolation.

## Public API

- `pmx_roles()` declares all critical and optional column roles explicitly.
- `mock_pmx()` generates an ordinary data frame or tibble with a lightweight
  `pmx_settings` attribute.
- `validate_pmx()` returns a structured report and can fail in strict mode.
- `compare_pmx()` returns structural summaries and exploratory plots without
  claiming statistical equivalence.

## Installation and development

```r
install.packages(c("testthat", "roxygen2"))
# From this repository:
install.packages(".", repos = NULL, type = "source")
```

Open `pmxmock.Rproj` for RStudio development. Run tests with:

```r
testthat::test_local()
```

## Theophylline repeated dosing

```r
library(pmxmock)

data("theo_md", package = "nlmixr2data")
theo_roles <- pmx_roles(
  id = "ID", time = "TIME", dv = "DV", amt = "AMT",
  evid = "EVID", cmt = "CMT", covariates = "WT"
)

theo_mock <- mock_pmx(theo_md, theo_roles, seed = 101)
theo_check <- validate_pmx(theo_mock, theo_roles)
theo_comparison <- compare_pmx(theo_md, theo_mock, theo_roles)
theo_comparison$plots$faceted
```

The complete seven-dose subject templates are retained, including dose and
observation rows tied at the same time. Mock IDs are new and concentrations are
generated as coherent nonnegative longitudinal blends.

## Warfarin PK and PD endpoints

```r
data("warfarin", package = "nlmixr2data")
warfarin_roles <- pmx_roles(
  id = "id", time = "time", dv = "dv", amt = "amt",
  evid = "evid", dvid = "dvid",
  covariates = c("wt", "age", "sex")
)

warfarin_mock <- mock_pmx(warfarin, warfarin_roles, seed = 202)
warfarin_comparison <- compare_pmx(warfarin, warfarin_mock, warfarin_roles)
warfarin_comparison$plots$overlay
```

The lower-case schema and factor levels are retained. The `cp` concentration
and `pca` response endpoints are transformed, interpolated, blended, and
constrained separately.

## WBC infusion and delayed response

```r
data("wbcSim", package = "nlmixr2data")
wbc_roles <- pmx_roles(
  id = "ID", time = "TIME", dv = "DV", amt = "AMT",
  evid = "EVID", cmt = "CMT", rate = "RATE",
  covariates = c("V2I", "V1I", "CLI")
)

wbc_mock <- mock_pmx(wbcSim, wbc_roles, seed = 303)
wbc_comparison <- compare_pmx(wbcSim, wbc_mock, wbc_roles)
wbc_comparison$plots$faceted
```

Positive infusion-start and negative infusion-stop `AMT`/`RATE` records are
copied together from a subject template. Delayed nadir/recovery shapes come
from blended observed trajectories; the known Friberg model is never called.

A runnable version of all three examples is in
[`scripts/demo_nlmixr2data.R`](scripts/demo_nlmixr2data.R).

## Important limitations

- Mock records can retain structural information from the source, including a
  complete event template. Treat the result as potentially sensitive until an
  appropriate privacy review says otherwise.
- The algorithm intentionally does not preserve parameter estimates,
  covariate-response relationships, or the source distribution exactly.
- Sparse endpoints and small compatible event-pattern groups may require a
  documented fallback and emit a warning. Inspect `attr(x, "pmx_settings")`.
- Normalized observation-window interpolation is pragmatic for irregular
  follow-up; it is not an occasion model or a biological time-warping model.
- Nonnegative constraints are applied only to endpoints classified as
  positive-like from their observed support.
- The prototype is intended for workflow development, not scientific analysis.

## Repository layout

- `R/` — package implementation
- `tests/testthat/` — fixture-based unit tests and optional integration tests
- `scripts/` — runnable demonstrations
- `references/` — papers, specifications, and notes excluded from package builds
- `data/` and `output/` — local/generated artifacts ignored by Git

## License

[MIT](LICENSE.md) © 2026 Andrew Stein.
