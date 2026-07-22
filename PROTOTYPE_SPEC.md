You are a senior R package developer with pharmacometrics experience. Work directly in the current directory and implement a working prototype—not merely a design document or code snippets.

First inspect the repository, including AGENTS.md if present, git status, DESCRIPTION, README, and existing R/package structure. Preserve existing work and conventions. If this is not yet an R package, scaffold one named `pmxmock`. Do not commit changes unless asked.

## Objective

Build an R package prototype that generates structurally faithful mock pharmacometric datasets from an existing PMX event dataset without requiring a PK, PD, or NLME model.

The intended use is:

1. Run the generator inside a restricted environment that can access the real data.
2. Export the mock data to a less restricted development environment.
3. Develop and debug data-processing, plotting, and model-workflow code there.
4. Port that code back to the restricted environment.

This prototype is not intended to:

- preserve parameter estimates or covariate relationships;
- support scientific model selection;
- reproduce the source distribution precisely;
- provide formal anonymization or privacy guarantees.

Use the phrase “mock data for model-workflow exploration” in the documentation. Do not market the output as anonymous data.

## Development sequence

Implement and verify the datasets in this order:

1. `nlmixr2data::theo_md`
2. `nlmixr2data::warfarin`
3. `nlmixr2data::wbcSim`

Finish a working theophylline implementation and its tests before generalizing to warfarin and WBC. Nevertheless, the final prototype should demonstrate all three datasets if feasible. Do not leave three partially implemented paths. If an unexpected limitation prevents completion of all three, deliver a completely working theophylline implementation with explicit extension points and explain the blocker.

Inspect the actual objects at runtime. Do not rely blindly on documentation or assume capitalization, factor levels, missing-value conventions, column names, or the meaning of zero versus `NA`.

## Core design principle

Separate the source dataset into:

1. Subject-level data:
   - ID
   - weight, age, sex, and other baseline covariates

2. Event skeleton:
   - TIME
   - EVID
   - AMT
   - RATE
   - CMT
   - DVID
   - dosing sequence
   - observation rows
   - row ordering, including ordering of tied times

3. Deterministic row logic:
   - which rows are doses versus observations;
   - where DV is allowed;
   - MDV, if present;
   - missing-observation conventions;
   - factor levels and column types.

4. Measurement trajectories:
   - PK concentrations;
   - PD measurements;
   - WBC values.

Never average, PCA-transform, or independently generate event-control variables such as EVID, MDV, CMT, DVID, ADDL, II, or RATE. Construct those as a coherent event skeleton.

## Public API

Implement an API approximately like this, adjusting only where a clearly better design is justified:

    roles <- pmx_roles(
      id = "ID",
      time = "TIME",
      dv = "DV",
      amt = "AMT",
      evid = "EVID",
      cmt = "CMT",
      dvid = NULL,
      mdv = NULL,
      rate = NULL,
      covariates = "WT"
    )

    mock <- mock_pmx(
      data,
      roles,
      n_subjects = NULL,
      seed = 123,
      event_method = "template",
      dv_method = "avatar_blend",
      k = 5,
      pca_variance = 0.90,
      subject_noise_sd = 0.15,
      residual_noise_sd = 0.05,
      residual_phi = 0.6,
      time_jitter = 0
    )

Requirements:

- `n_subjects = NULL` means the same number of subjects as the source.
- Preserve source column names, column order, and practical column classes.
- Generate new subject IDs. Preserve numeric IDs when the source ID is numeric and character IDs when it is character.
- Return an ordinary data frame or tibble usable by standard PMX and plotting code.
- Store generator settings as a lightweight attribute, but do not require downstream users to understand a custom class.
- The seed must make the complete result reproducible.
- Avoid unnecessarily changing the caller’s global random-number state.
- Column-role matching must be explicit. Do not guess critical roles silently.
- Optional roles must work when absent. For example, `theo_md` has no DVID or MDV.

Also implement:

    validate_pmx(data, roles)
    compare_pmx(source, mock, roles)

`validate_pmx()` should return a useful structured report and optionally fail in strict mode. `compare_pmx()` should produce concise structural summaries and useful exploratory plots, not claim statistical equivalence.

## Mock-data algorithm

Implement an AVATAR-inspired method adapted to longitudinal PMX data. Clearly document that it is an AVATAR-like variant, not necessarily an exact reproduction of the published AVATAR software.

### A. Build subject profiles

Construct one feature vector per source subject using:

- one representative value for each baseline covariate;
- transformed DV values evaluated on a common grid separately for each endpoint;
- no event-control variables in the numerical profile.

For aligned data, the common grid may use the observed nominal times. For irregular data, construct a reasonable dataset-level grid and interpolate within endpoint. Avoid excessive extrapolation.

For repeated-dose data, first assess whether absolute TIME adequately aligns subjects. Where dosing schedules differ, group subjects by compatible dosing/event patterns or use occasion/TAD-style alignment. Do not mix clearly incompatible dosing schedules without warning.

For DV transformation:

- use a log transform with a defensible small offset for positive PK-like or biomarker data;
- use identity transformation if values can legitimately be zero or negative;
- make the choice endpoint-specific;
- record the transformation used.

Before PCA:

- median-impute feature values only for the purpose of distance calculation;
- standardize continuous features;
- drop zero-variance and unusable columns;
- handle the small-subject, many-feature case safely.

Use PCA and retain the smallest number of components explaining `pca_variance`, subject to the available rank. Provide a sensible fallback if PCA cannot be performed.

### B. Choose an event skeleton

For each mock subject:

- sample a source subject’s complete event skeleton;
- assign a new synthetic ID;
- preserve event ordering and tied-time ordering;
- preserve or coherently perturb the skeleton according to options;
- default to no time jitter in the prototype;
- preserve the source dataset’s conventions for AMT and DV on event rows.

The sampled skeleton subject is the anchor used to identify a compatible neighborhood. The anchor itself should not be used as a measurement donor when enough other subjects exist.

### C. Select neighbors

Find up to `k` nearest compatible subjects in PCA space, excluding the anchor.

Handle small compatible groups explicitly:

- reduce `k` when needed;
- require at least two donors when possible;
- issue an informative warning or use a documented fallback when there are too few compatible donors;
- never divide by a zero distance.

### D. Generate randomized weights

Use a documented randomized weighting rule similar to:

    raw_weight_j =
      Exp(1) / max(distance_j, epsilon) *
      2^(-randomized_rank_j)

Then normalize the weights to sum to one.

Use the same subject-level donor weights across covariates and the longitudinal trajectory so that the generated subject remains internally coherent. If a donor lacks a value for a particular endpoint/time, renormalize the available weights for that value.

Prevent a single donor from silently receiving essentially all the weight. Either cap the maximum normalized weight or redraw weights above a documented threshold.

### E. Generate covariates

- Continuous covariates: blend on an appropriate scale and add modest subject-level perturbation.
- Positive skewed covariates: consider log-scale blending.
- Categorical covariates: sample a donor category using the donor weights; do not numerically average factor codes.
- Preserve factor levels and source column types.
- Preserve subject-level constancy across rows.

Scientific preservation of covariate–response relationships is not an objective.

### F. Generate DV trajectories

For every observation row in the selected skeleton:

1. Identify its endpoint.
2. Interpolate each donor’s transformed trajectory to the target observation time or appropriate aligned time.
3. Blend donor values using the subject-level weights.
4. Add:
   - a subject-level random shift; and
   - modest correlated within-trajectory noise, such as AR(1) noise in observation order.
5. Transform back to the original scale.
6. Enforce only defensible domain constraints, such as nonnegative concentrations or WBC counts.

Process each DVID separately. Do not blend warfarin PK concentrations and PCA response measurements together.

Only replace DV on rows that are genuine observations according to EVID, MDV if present, and the source convention. Preserve missing-observation patterns from the chosen event skeleton unless an explicit option says otherwise.

### G. Reconstruct and validate

After synthesis:

- restore the original schema and column order;
- restore factor and other practical column classes;
- sort by synthetic ID, TIME, and original within-time order;
- derive row-logic fields when present rather than generating them stochastically;
- remove internal helper columns;
- run structural validation.

## Dataset-specific demonstrations

Create a runnable script, vignette, or README section demonstrating each dataset.

### 1. Theophylline: `theo_md`

Use an explicit role mapping based on the actual data, expected to resemble:

    pmx_roles(
      id = "ID",
      time = "TIME",
      dv = "DV",
      amt = "AMT",
      evid = "EVID",
      cmt = "CMT",
      covariates = "WT"
    )

This is the first release gate. Demonstrate:

- repeated dosing;
- preserved dose and observation rows;
- coherent positive concentration-time profiles;
- new IDs;
- reproducibility;
- individual concentration-time plots;
- an overlay or faceted comparison of source and mock data.

### 2. Warfarin: `warfarin`

Inspect the actual lower-case schema and factor levels. Provide an explicit role mapping for:

- ID;
- TIME;
- DV;
- AMT;
- EVID;
- DVID;
- weight, age, and sex.

Handle PK and PD endpoints separately. Preserve DVID categories and their event semantics. Demonstrate endpoint-faceted plots.

### 3. WBC: `wbcSim`

Inspect the actual object because its documented covariate names may be inconsistent. Preserve:

- infusion or dosing rows;
- AMT;
- RATE;
- EVID;
- CMT;
- delayed longitudinal WBC profiles;
- subject-level input covariates.

Demonstrate that the mock data retain plausible-looking delayed/nadir/recovery trajectories without calling the known Friberg model.

## Validation and tests

Use `testthat`. Include small internal test fixtures so core tests do not require `nlmixr2data` to be installed. Add integration tests for the three public datasets and skip them gracefully if the package is unavailable.

At minimum test:

1. Identical results with the same seed.
2. Different results with different seeds.
3. Requested number of subjects is produced.
4. New IDs have the appropriate type and do not preserve a source ID linkage.
5. Original columns and order are preserved.
6. Relevant column classes and factor levels are preserved.
7. Covariates are constant within subject.
8. Dose rows retain coherent EVID, AMT, RATE, and CMT values.
9. Observation rows have appropriate DV behavior.
10. DVs are finite and respect endpoint constraints.
11. Multiple DVIDs remain separated.
12. Subjects are correctly ordered, including tied TIME rows.
13. Small `n`, `k >= n`, zero-variance features, missing DVs, and duplicated profiles are handled.
14. No generated subject has an entire DV trajectory exactly equal to one source donor, except where generation is mathematically impossible; such fallbacks must warn.
15. `validate_pmx()` passes for generated examples.
16. The three `nlmixr2data` demonstrations run end to end.

Run:

- the full test suite;
- package documentation generation;
- `R CMD check` or `devtools::check()` with no errors or warnings attributable to the package.

Do not suppress legitimate warnings merely to make the check appear clean.

## Package structure and dependencies

Prefer small, testable internal functions, for example:

- role validation;
- observation/event classification;
- event signatures;
- subject-profile construction;
- PCA/neighborhood selection;
- randomized weights;
- endpoint interpolation;
- trajectory blending;
- schema restoration;
- PMX validation.

Use minimal mainstream R dependencies. Do not depend on `rxode2`, `nlmixr2`, `mrgsolve`, NONMEM, Monolix, Python, a GAN framework, or a fitted PK/PD model. `nlmixr2data` may be a suggested development/example dependency rather than a runtime dependency.

Add roxygen documentation and a README that explains:

- what problem the package solves;
- the event-template plus trajectory-blend algorithm;
- why event fields are not synthesized numerically;
- intended and prohibited uses;
- a short example for all three datasets;
- known limitations, especially the lack of formal privacy guarantees and scientific fidelity.

## Scope control

Prioritize a clear, working prototype over speculative abstraction.

Do not implement:

- formal differential privacy;
- CTGAN or TVAE;
- a general synthetic-data benchmarking framework;
- PK parameter recovery;
- automated NLME model selection;
- a Shiny application.

Leave clean extension points for later work.

## Final response

Do not stop after giving a plan. Implement the code and verify it.

At the end, report concisely:

- the public API;
- files created or modified;
- what works for each of the three datasets;
- exact test and package-check results;
- important limitations;
- any decisions where the implementation differs from this specification and why.