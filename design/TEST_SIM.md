# Simulation evaluation test plan

## Purpose

This document defines the continuously maintained evaluation suite for
`pmxSynthData`. It turns failures found while developing the demonstrations
into reproducible checks. The suite must answer three separate questions:

1. Is the generated PMX table structurally valid?
2. Does it retain the coarse regimen, sampling, and trajectory behavior needed
   for workflow testing?
3. Do the demonstrations display that behavior honestly and clearly?

The tests are regression and utility tests, not evidence that generated data
are scientifically interchangeable with the source. Privacy-mechanism tests
and the proof argument remain separate. Repeated fitting for evaluation must
use only public or fully simulated fixtures; it must never repeatedly consume
a confidential dataset without explicit privacy accounting and authorization.

This is a living document. Every newly discovered simulation or demonstration
failure should add:

- a row to the issue registry;
- an objective metric or invariant when one is possible;
- a minimal automated regression test; and
- a multi-seed or visual diagnostic when a single deterministic assertion is
  insufficient.

## Dataset registry

The first evaluation set contains every dataset used in the practical demo.
Configuration values below are public domains, contribution limits, schema
semantics, and endpoint alignments. Source-derived regimen and sampling
schedules must not be supplied to the fit.

| ID | Dataset | Origin | Main behavior exercised | Required checks |
|---|---|---|---|---|
| `censoring` | Eight-subject expansion of `pmx_censoring_fixture()` | Package-owned fully simulated fixture | Minimal end-to-end workflow, study-time endpoint, schema restoration, left/right/interval censoring | The source covers all conventions; generated censoring states and limits are coherent; source and synthetic tables are displayed; validation passes |
| `theo_md` | `nlmixr2data::theo_md` | Public package data | Seven Q24H doses, dose-relative log-PK, dense profiles after occasions 1 and 7, occasional trough after occasion 2, no samples after occasions 3--6 | Regimen and sampling are inferred; dose rows are not counted as samples; intensive profiles have one directional peak; inactive occasions stay inactive |
| `warfarin` | `nlmixr2data::warfarin` | Public package data | Lower-case schema, one dose, dose-relative `cp`, study-time `pca`, factors and multiple DVIDs, late PK follow-up | Both endpoints and all subjects remain; CP extends beyond 24 hours; endpoint-specific point counts and follow-up are retained broadly |
| `wbcSim` | `nlmixr2data::wbcSim` | Public package data | Infusion start/stop pairs, study-time log-WBC response, delayed nadir and recovery, numeric covariates | Infusion rows are coherent; no source sentinel schedule is copied; WBC declines and recovers; cohort and follow-up remain comparable |
| `nimoData` | `nlmixr2data::nimoData` | Public package data | Four nominal dose groups, ten approximately weekly infusions, declared OCC/TAD, long terminal follow-up, and a time-varying weight column | DOS is a subject property that conditions amount/rate/duration; every subject has ten coherent infusions; WGT is explicitly excluded; dose-relative sample count and terminal coverage remain broad |
| `mavoglurant` | `nlmixr2data::mavoglurant` | Public package data | One- and two-period profiles, TIME reset within OCC, occasion-varying assigned DOSE, numeric-coded SEX, infusion rows | Reset clocks validate within ID/OCC; DOSE equals positive AMT and is constant within ID/OCC; SEX is categorical; cohort and two-occasion event structure remain |

The next expansion dataset should be `pmx_simulated_fixture(60)`. It is already
used by package tests and provides a larger repeated-dose, two-endpoint study
for seed sweeps and privacy-utility experiments. Later additions should cover
IV bolus decay, multiple infusions, irregular dosing, below-quantification
patterns, missing covariates, more than two endpoints, and genuinely
multi-phasic dose-relative profiles.

Public datasets must be loaded from their installed packages during a run. Do
not copy their records into this repository.

## Issue registry

The registry below captures the problems found during prototype and demo
review so far. “Gate” describes the check that should prevent recurrence.

| ID | Dataset/area | Failure observed | Cause or interpretation | Required gate |
|---|---|---|---|---|
| `SIM-001` | All PMX data | Dose/event rows were visually interpreted as samples | Plots and summaries did not consistently separate `EVID != 0` from observations | Observation metrics and plots must exclude event rows; report event and observation counts separately |
| `SIM-002` | All named demos | Regimen and sampling times were specified from disclosed datasets | Demo configuration supplied dose times, counts, amounts, rates, endpoint grids, or occasion schedules | Inspect the fitted public configuration and assert those overrides are absent; infer them inside `fit_private_pmx()` |
| `SIM-003` | `theo_md` | Q24H dosing incorrectly appeared to imply Q24H sampling | One local endpoint grid was repeated after every generated dose | Test occasion activation separately from conditional sample count; occasions 3--6 must receive no observations |
| `SIM-004` | `theo_md` | The sparse occasion-2 trough was confused with a full profile | Presence probability and sample density were represented by one unconditional count | Test the two-part release: occasion presence probability and observations conditional on being sampled |
| `SIM-005` | `theo_md` | Automatic local grid extended to the next dose and reassigned a late sample to the following occasion | Generic horizon and timing jitter allowed equality with the next dose boundary | Every dose-relative time must remain strictly before the next dose and retain its generating occasion |
| `SIM-006` | `theo_md` | Synthetic PK profiles differed markedly and developed an artificial trough followed by another peak | Unoccupied cells on a wide log domain were decoded as domain-midpoint measurements instead of missing support | Interpolate unoccupied cells from released occupied neighbors; never treat absence as a DV measurement |
| `SIM-007` | `theo_md` | Individual profiles still showed secondary peaks after the mean curve was repaired | AR(1) noise continued across distant dose occasions and unconstrained residuals reversed the post-peak direction | Restart serial noise by occasion; if the released curve is approximately unimodal, require at most one directional peak per generated occasion |
| `SIM-008` | `warfarin` | CP observations disappeared after about 24 hours | Count matching retained the first grid cells and systematically trimmed late cells | Select timing cells with released presence probabilities; require late CP coverage and a median subject maximum time of at least 72 hours |
| `SIM-009` | `warfarin` | CP appeared to contain fewer patients or missing data | Endpoint allocation and combined plotting obscured subject/endpoint coverage | Assert cohort size, patients per endpoint, endpoint set, and mean points per patient; facet endpoints separately |
| `SIM-010` | All demos | Synthetic datasets sometimes contained fewer subjects than their sources | Generation size was an independent default rather than the fitted privacy-accounted count | Omitted `n_subjects` must use the fitted count; public-fixture demos require exact cohort equality |
| `SIM-011` | All demos | Number of generated time points could differ materially from source | Full grids were repeated or trimmed without respecting released observation totals | Compare endpoint-specific observations per subject with explicit tolerances and test the released total observation count |
| `SIM-012` | `wbcSim` | Infusion/event behavior could be incoherent or reproduce an exceptional source schedule | Start/stop construction and generalized-regimen logic were incomplete | Require paired positive/negative amount and rate rows, coherent duration, bounded values, and absence of the 4580-hour source schedule |
| `SIM-013` | `wbcSim` | Longitudinal response could lose delayed decline, nadir, or recovery | Study-time behavior was at risk of being restarted at a dose or shortened by grid selection | Require a value below baseline followed by recovery and broad late follow-up |
| `SIM-014` | All endpoints | Exact source timing vectors could be copied | Earlier designs considered source anchors or exact schedules | Compare complete source and generated time vectors; no generated vector may be identical to a source vector |
| `SIM-015` | All data | IDs, schema classes, factor levels, or endpoint columns could be lost or reused | Schema restoration and ID generation were incomplete | Require new IDs, original column order/classes, declared public factor levels, and all endpoints |
| `SIM-016` | Demo plots | Source and synthetic observations were not overlaid or comparably displayed | Separate plotting code and incorrect grouping | Use one comparison data frame, consistent colors, and source-above/synthetic-below facets |
| `SIM-017` | Demo plots | Lines stopped between profile segments | Lines were grouped by event/occasion rather than by subject for study-time displays | Study-time plots group all observed points by dataset, subject, and endpoint; event rows remain excluded |
| `SIM-018` | Demo plots | Log scaling obscured Warfarin and other comparisons | `xgx_scale_log10()` was applied when a linear comparison was more interpretable | Demo comparison plots use linear DV axes unless a dataset-specific design decision explicitly changes this |
| `SIM-019` | Demo plots | Multiple DVIDs and source/synthetic differences were difficult to see | Overlay-only panels were too dense | Facet source above synthetic and endpoint columns side by side while retaining dataset colors |
| `SIM-020` | Decoding, all noisy fits | Grid cells with no source support decoded to the bottom of the endpoint working domain, producing spurious deep troughs on log endpoints at every cohort size | Decoders compared released presence fields against the bare constant `0.25`, but those fields are unnormalized subject counts of order N, not probabilities. Laplace noise around zero passed the gate, then the separately noised value sum clamped to zero and `0 / noise` decoded to the domain floor. `.fill_unoccupied_curve()` could not repair it because the cell was already labelled occupied | Gate released presence on `.support_threshold(count, noise_scale)`, derived from the release's own `sensitivity / epsilon`. A cell with no support must never decode to a working-domain endpoint; see `tests/testthat/test-decode-support-threshold.R` |
| `SIM-020` | Demo output | The actual source and generated tables were not shown | Vignette focused only on plots and summaries | Every demo prints a small source preview and a small synthetic preview, with public-data labeling |
| `SIM-021` | Documentation | Readers expected an ODE/NLME or spline model | The implemented fixed-grid, smoothing, interpolation, timing, and variability model was not stated early or precisely | Vignette test/check confirms the method introduction names fixed-grid summaries, 1--2--1 smoothing, linear interpolation, and the absence of ODE/NLME/splines |
| `SIM-022` | Vignette artifacts | A repaired source vignette could still appear broken | An ignored, previously rendered HTML file or already loaded installed namespace was stale | Evaluation rebuilds vignettes from the current source package in a clean library and records package/source identity in its manifest |
| `SIM-023` | Demo plotting | A narrow-bin median line itself could look jagged or multi-peaked | Sparse subjects at slightly different actual times populated different 0.25-hour bins | Treat thick summary lines as diagnostics only; compute peak tests from complete ordered profiles and include the released mean curve in the report |
| `SIM-024` | Cohort assignments | `ACTARM`, `TRT`, or nominal dose group could be sampled independently of the generated regimen | Ordinary categorical covariates retained only marginal proportions and had no link to event summaries | Declare subject-level assignment fields as `subject_properties`; release stratum count and regimen jointly; require generated property-to-event coherence |
| `SIM-025` | Numeric category codes | A numeric-coded category such as mavoglurant `SEX` could be modeled as a continuous covariate | Numeric storage class was treated as sufficient evidence of continuity | A declared public category domain forces categorical modeling regardless of storage type |
| `SIM-026` | Occasion-assigned dose | Mavoglurant `DOSE` could be missing, vary within an occasion, or disagree with generated AMT | Nominal dose was restored as an unmodeled schema default | An `assigned_dose` role is reconstructed from each generated positive event AMT and validated on every ID/OCC profile |
| `SIM-027` | Reset clocks | Mavoglurant source TIME restarts at zero in the second OCC and was rejected as globally decreasing | Validation and plotting assumed one monotone clock per subject | When OCC is declared, validate time order within ID/OCC and group dose-relative display lines by subject/occasion |
| `SIM-028` | Positive-rate infusions | NimoData and mavoglurant encode a positive rate without explicit negative source stop rows, so duration could decode as zero | Duration inference looked only for a later negative event row | If no stop row exists, infer bounded duration as AMT/RATE; generated starts and stops must remain coherent |
| `SIM-029` | Terminal washout | A final profile longer than one dose interval was compressed into the ordinary interdose window | The automatic local basis was rescaled to one interval for every occasion | Keep nonterminal observations before the next dose, but let the final occasion use its released occupied horizon |
| `SIM-030` | Time-varying covariates | NimoData WGT changes within subject and could be mislabeled as a baseline covariate | The prototype generates only subject-constant baseline covariates | Validation rejects varying baseline covariates; explicitly exclude WGT until a longitudinal-covariate model exists |

When an issue is fixed, keep its row. The registry is a permanent record of
the behavior that must continue to work.

### Related setup and documentation issues

Several problems encountered during the same review are important preflight
checks but are not simulation-fidelity metrics:

- The OpenDP production backend is an R package dependency. The evaluation
  runner should print `dp_backend_status()` and clearly distinguish an
  unavailable production backend from the explicitly nonprivate public-fixture
  backend.
- Vignette PDF/manual builds require a working `pdflatex`; TinyTeX installation
  and PATH discovery are environment checks rather than generator tests.
- Privacy concepts and the simulation algorithm were initially mixed in one
  method vignette. The package now requires separate privacy-introduction and
  simulation-method vignettes, plus the practical demo.
- The practical vignette must explain near the top that the implementation is
  a fixed-grid population generator, not an ODE, NLME likelihood, or spline.
- Evaluation should report missing optional packages and rendering tools as
  explicit preflight failures or skips, never as simulation successes.

### Existing coverage reused

The evaluation layer builds on, rather than silently replacing, these existing
checks:

- `tests/testthat/test-integration-nlmixr2data.R` covers inferred Theophylline
  dosing/sampling and single peaks, late Warfarin CP, WBC infusion/recovery,
  cohort size, schema, and timing non-copying.
- `tests/testthat/test-private-fitting.R` covers automatic grid bases,
  unoccupied-cell interpolation, and the separate occasion-presence and
  conditional-count release.
- `tests/testthat/test-generation.R` covers cohort-size defaults, timing-cell
  selection, schema restoration, repeated-dose versus study-time behavior,
  chronological coherence, reproducibility, subject-property/regimen
  coherence, and assigned-dose reconstruction.
- `tests/testthat/test-censoring.R` covers generated censoring conventions.
- `scripts/demo_nlmixr2data.R` currently performs cohort, point-count, and
  follow-up checks and constructs the visual comparison panels.

The shared evaluator described below closes the earlier gaps in metric reuse,
systematic multi-seed evaluation, machine-readable results, plot-semantic
assertions, and one report spanning all datasets and known gates. Future
coverage should extend this layer rather than create a competing metric path.

## Implemented evaluation

Use one shared metric layer with two callers.

### `tests/testthat/helper-simulation-evaluation.R`

This helper contains no test expectations and writes no files. It provides:

- a dataset registry returning source data, roles, endpoints, public bounds,
  public schema semantics, contribution limits, and budget allocation;
- row classifiers that distinguish dose/event, infusion-stop, and observed
  endpoint rows;
- regimen summaries by subject: dose count, interval, amount, rate, and
  infusion duration;
- sampling summaries by endpoint and occasion: activation, conditional count,
  total count, timing-cell coverage, first/last time, and late-follow-up
  coverage;
- trajectory summaries on study-time or TAD bins: bounded quantiles, peak or
  nadir time, directional peak count, decline/recovery indicators, and broad
  range;
- schema, ID, endpoint, timing-vector-copy, and PMX-validity checks; and
- plot-data checks for facet order, group identifiers, event exclusion, axis
  scale, and consistent dataset labels/colors.

Metrics should be ordinary data frames with stable column names so the same
code can feed `testthat` expectations and the longer report script.

### `tests/testthat/test-simulation-evaluation.R`

This is the fast deterministic regression gate run by `devtools::test()` and
`R CMD check`. It:

1. Runs package-owned fixtures unconditionally.
2. Runs `nlmixr2data` cases with `skip_if_not_installed("nlmixr2data")`.
3. Uses the guarded noiseless public-fixture backend and fixed generation
   seeds.
4. Fits without source-derived regimen or sampling overrides.
5. Asserts all hard invariants and the dataset-specific gates below.
6. Avoids writing plots or evaluation artifacts.

Some overlapping assertions remain in
`tests/testthat/test-integration-nlmixr2data.R` as focused regression tests.
Consolidate them only after equivalent coverage is demonstrated, and avoid
duplicating two independent implementations of the same metric.

### `scripts/evaluate_simulations.R`

This is the longer evaluation runner. A typical invocation is:

```sh
Rscript scripts/evaluate_simulations.R \
  --datasets=all \
  --seeds=101:200 \
  --backend=public \
  --output=output/simulation-evaluation
```

The script uses base-R argument parsing to avoid a new dependency. It:

1. Loads the registry from package/test helper code without duplicating dataset
   definitions.
2. Fits each public dataset once for an ordinary generation-seed sweep.
3. Generates every requested seed and computes the shared metrics.
4. Records hard failures immediately but finishes the run so all problems are
   visible together.
5. Produces source-above/synthetic-below study-time and scientific-clock plots,
   endpoint sampling panels, regimen summaries, and released-curve panels.
6. Writes only derived metrics and figures by default, not source or synthetic
   row-level datasets.
7. Exits nonzero when a hard gate fails.

Outputs under the ignored `output/` directory are:

```text
output/simulation-evaluation/
  run-manifest.txt
  metrics-by-seed.csv
  gate-results.csv
  failures.csv
  regimen-by-fit.csv
  sampling-by-fit.csv
  subject-properties-by-fit.csv
  summary.html
  figures/
    censoring-study-time.png
    theo_md-study-time.png
    theo_md-tad.png
    warfarin-study-time.png
    wbcSim-study-time.png
    nimoData-tad.png
    mavoglurant-tad.png
```

The manifest should record package version, Git commit when available, dirty
worktree status, R and dependency versions, dataset package versions, backend,
bounds, contribution limits, generation seeds, timestamp, and platform. It
must also state whether the run used the public-fixture or formal DP backend.

### Optional privacy-utility sweep

A separate mode may repeatedly fit the public/simulated datasets across
epsilon values and production-backend noise draws. Because fitting intentionally
does not accept a user-controlled privacy-noise seed, this mode should report
distributions rather than exact snapshots. It must never be pointed at a
confidential dataset merely for convenient tuning.

## Initial gates and tolerances

### Hard gates for every dataset

- `validate_pmx(..., strict = TRUE)` passes.
- Generated IDs are disjoint from source IDs.
- Column names, order, practical classes, declared factor levels, and endpoint
  set are restored.
- Every observation has finite DV, a declared endpoint, and zero or missing
  amount; every event has coherent event fields.
- Times are finite and bounded. They are nondecreasing within subject, or
  within subject and declared occasion when the source clock resets by OCC.
- No complete generated timing vector is identical to a source timing vector.
- The fitted model contains no raw source rows, identifiers, or unnoised
  aggregates.
- Named demo configurations contain no `dose_times`, `dose_interval`,
  `n_doses`, `dose_amount`, `dose_rate`, `infusion_duration`,
  `endpoint_grids`, or `endpoint_occasion_grids` derived from the source.
- With the public-fixture backend and omitted `n_subjects`, generated and
  source cohort sizes are equal.
- Every expected endpoint is represented in every dataset-level result.

### General utility gates

- Mean observations per subject and endpoint differ by no more than
  `max(1, 25% of the source mean)` in the deterministic public-fixture run.
- First and last endpoint observation times differ by no more than 20% of the
  declared public time span, unless a stricter dataset-specific gate applies.
- Event rows and observation rows are counted and reported separately.
- Multi-seed evaluation reports median, 5th, and 95th percentiles rather than
  relying on one favorable seed. Hard-invariant failures are never tolerated;
  soft utility gates should pass for at least 95% of the planned seeds.

These are initial engineering tolerances, not scientific equivalence margins.
Any change must be justified in this file and must not be made only to excuse a
new regression.

### `theo_md`

- Fitted generalized dose count is 7 and interval is 24 hours within 0.1 hour.
- Generated subjects retain all seven dose events.
- Released sampling probabilities identify occasions 1 and 7 as active,
  occasion 2 as sparse, and occasions 3--6 as inactive.
- Intensive occasions contain approximately 10--11 observations and the
  sparse occasion contains zero or one according to its fitted activation
  draw.
- Dose-relative observations stay strictly within their generating occasion.
- Every intensive generated profile has at most one directional peak after
  zero-length/flat differences are removed.
- Peak occurs after the first observation and is followed by an overall
  decline to the last observation.

### `warfarin`

- Both `cp` and `pca` are present with the lower-case source schema.
- Generated cohort and patients represented per endpoint equal the source
  cohort for the public fixture.
- CP observations per subject differ from source by no more than one in the
  deterministic run.
- At least one generated CP observation is later than 24 hours and the median
  subject-specific final CP time is at least 72 hours.
- Pca remains a study-time trajectory and is not restarted as a dose-relative
  excursion.

### `wbcSim`

- Positive infusion-start and negative infusion-stop rows both exist.
- Paired generated amount and rate fields are coherent and bounded.
- The exceptional 4580-hour source schedule is not reproduced.
- Every subject has a post-baseline value below baseline and a later recovery
  above its nadir when enough observations are generated.
- Late WBC follow-up remains within the general coverage tolerance.

### `nimoData`

- `subject_property_summary()` contains public DOS strata 50, 100, 200, and
  400 with ten inferred doses in every stratum.
- Every generated subject has one constant DOS value, ten positive dose events,
  and positive AMT equal to DOS.
- Positive infusion starts and negative generated stop rows both exist; fitted
  duration can be inferred from source AMT/RATE when no source stop row exists.
- The final occasion may extend beyond one weekly interval, while observations
  on occasions 1--9 stay before the following dose.
- WGT is absent from generated schema because it is explicitly excluded; BSA,
  AGE, and HGT remain constant within subject.

### `mavoglurant`

- The 120-subject public cohort is preserved by default.
- TIME ordering is valid within ID/OCC even though the source clock resets.
- Every generated subject has two positive dose occasions in the deterministic
  evaluation case.
- Generated DOSE is finite on every row, constant within ID/OCC, and equals
  positive event AMT.
- Numeric-coded SEX uses only its declared categorical levels.
- The current generalized-regimen case does not promise preservation of the
  source crossover-sequence distribution when no ACTARM/TRT/sequence property
  is present; that remains an explicit utility limitation, not a passed gate.

### Censoring fixture

- The source fixture contains uncensored, left-, right-, and interval-censored
  examples; generated states use only valid codes and include censoring.
- DV, CENS, and LIMIT combinations satisfy the declared convention.
- Censoring is applied to generated latent values, not copied source rows.

## Plot and vignette evaluation

Pixel snapshots are likely to be brittle across R, graphics-device, and font
versions. Automated tests should inspect plot semantics instead:

- plot data contain observations only;
- grouping connects a subject's chronological observations across study-time
  gaps;
- facet levels place Source above Synthetic and endpoints in columns;
- colors map consistently to Source and Synthetic;
- demo DV scales are linear; and
- every plot contains all expected subjects and endpoints.

The long evaluation script should still render PNGs for human review. Reviewers
should look for disappearing patients, lost late samples, unintended repeated
profiles, false peaks or nadirs, axis compression, connections through event
rows, and misleading summary lines. A visual concern should be converted into
an objective regression metric whenever possible.

Vignettes must be rendered from a clean temporary library containing the
current source build. Do not validate a stale HTML file or an older namespace
already loaded in an interactive R session.

## Where the tests should live

The core regression checks belong in `tests/testthat/`. Your instinct is right:
patient counts, endpoint coverage, inferred occasion sampling, late Warfarin
CP, infusion coherence, one-peak Theophylline profiles, schema restoration,
and timing-vector non-copying are package behavior. They should fail during
ordinary development and `R CMD check`, close to the code change that caused
the regression.

The complete evaluation simulation should not live only in the test suite.
Multi-seed runs, repeated privacy-noise fits, figure rendering, HTML reports,
and artifact writing are slower and partly judgment-based. Put that runner in
`scripts/evaluate_simulations.R` and keep its outputs under ignored `output/`.

The recommended division is therefore:

| Location | Responsibility |
|---|---|
| `tests/testthat/` | Fast deterministic invariants and minimal regressions for every issue that can be asserted objectively |
| `scripts/evaluate_simulations.R` | Multi-seed, multi-backend, distributional, and visual evaluation with persisted reports |
| `design/TEST_SIM.md` | Dataset registry, issue history, gates, tolerances, and rationale |

If only one location were allowed, choose `tests/testthat/`, because a test
that is not run automatically will eventually be forgotten. In practice the
hybrid is stronger: every issue first receives the smallest reliable automated
test, while the evaluation script answers broader questions that are too slow
or visual for package checks.

## Maintenance workflow

1. Reproduce a newly reported problem with a public or simulated dataset.
2. Add it to the issue registry before changing thresholds or implementation.
3. Add the smallest deterministic failing test.
4. Fix the implementation and run the focused test.
5. Run the complete `testthat` suite and `R CMD check`.
6. Run the multi-seed evaluation for affected datasets and inspect its plots.
7. Record intentional threshold or dataset changes in this document and NEWS.

No issue is considered closed solely because one plot or one seed looks good.
