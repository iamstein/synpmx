# Version 2 implementation specification: private PMX generator

You are a senior R package developer with pharmacometrics, longitudinal-data, and differential-privacy experience. Work directly in the existing package repository and implement a complete working package—not merely a design, plan, or collection of code snippets. Do not commit changes unless asked.  Th

## Version 2 is a fresh synthesis engine

Version 1 used an AVATAR-like anchor-and-donor blending algorithm. Version 2 is a fresh implementation centered on a subject-level differentially private population generator.

Do not preserve the Version 1 synthesis architecture merely for backward compatibility. In particular, Version 2 must not:

- select a source subject as an event-template anchor;
- copy a complete source event skeleton;
- choose raw source subjects as nearest-neighbor donors;
- blend raw subject trajectories;
- expose a `blend`, `avatar`, `template`, or exact-source-timing mode in the public API;
- describe jitter, generalization, k-anonymity, or disclosure testing as a substitute for differential privacy.

Git preserves the Version 1 implementation. Treat the main Version 2 code as a new engine. Reuse existing low-level utilities only after verifying that they do not retain raw subjects or conflict with the private design. Potentially reusable pieces include role validation, PMX row classification, schema reconstruction, plotting helpers, and structural validation.

Use the package name already present in `DESCRIPTION`. Do not rename the package as part of this work.

## Objective

Build a source-calibrated generator of **mock data for model-workflow exploration**. It must read an actual pharmacometric dataset inside a restricted environment, learn only differentially private population summaries, and generate structurally coherent mock PMX event datasets outside that environment.

The generated data should be useful for:

- data cleaning and transformation code;
- joins, reshaping, filtering, and derivation code;
- exploratory visualizations;
- control-file and model-run plumbing;
- testing whether PK, PD, biomarker, and repeated-dose workflows execute;
- preliminary model-code debugging before returning to the restricted environment.

The generator should capture broad features such as:

- the correct order of magnitude for each endpoint;
- whether an endpoint rises and falls after each dose;
- whether an endpoint changes over a longer study-level timescale;
- approximate dosing and sampling-time scales;
- coarse subject-to-subject and observation-level variability;
- endpoint, censoring, and event-table conventions.

It is not intended to preserve:

- exact source distributions;
- detailed covariance or covariate-response relationships;
- rare individual trajectories;
- parameter estimates;
- scientific conclusions;
- inferential validity;
- model-selection results.

The deliberate lack of statistical fidelity is an advantage for privacy and is acceptable for the intended workflow-development use.

## Start by inspecting the repository

Before editing:

1. Read `AGENTS.md`, `DESCRIPTION`, `NAMESPACE`, `README.md`, `NEWS.md`, the existing prototype specification, all current vignettes, and the test suite.
2. Inspect Git status and preserve unrelated work.
3. Run the current tests and package check to establish a baseline.
4. Inspect the actual `nlmixr2data::theo_md`, `nlmixr2data::warfarin`,
   `nlmixr2data::wbcSim`, `nlmixr2data::nimoData`, and
   `nlmixr2data::mavoglurant` objects at runtime.
5. Identify and remove or isolate all code paths that use raw anchors, donor subjects, nearest neighbors, nonprivate PCA, or exact source timing.
6. Inspect any existing `METHODS.md`. Its Version 1 AVATAR/template instructions are superseded by this specification.

# Privacy contract

## Formal guarantee

The private fitting procedure must provide subject-level, add-or-remove `(epsilon, delta)` differential privacy. One subject's complete longitudinal contribution is the privacy unit, including:

- all rows and visits;
- dosing and infusion records;
- actual and nominal timing;
- baseline covariates;
- all endpoint observations;
- DVID, censoring, and missingness information;
- rare schedules or protocol deviations.

For every pair of neighboring source datasets `D` and `D'` differing by one complete subject, and every possible set of released outputs `S`, the mechanism `M` must satisfy:

```text
Pr[M(D) in S] <= exp(epsilon) * Pr[M(D') in S] + delta
```

This definition limits how strongly one subject can influence the distribution of the released private model and everything generated solely from it.

The implementation and documentation must not claim zero privacy risk, guaranteed impossibility of re-identification, or automatic legal anonymity. The accurate claim is:

> Generated from a subject-level `(epsilon, delta)`-differentially private model.

## Identity, attributes, and participation

Explain clearly in the documentation:

- A release can expose that an unusual record existed without directly identifying its owner.
- Identification generally requires linkage to outside information, but absence of direct identifiers does not prevent linkage or attribute disclosure.
- Membership is one protected question: did a person's complete record influence the study release?
- Attribute disclosure is another: did a person's unusual measurements, schedule, response, or censoring pattern materially affect the output?
- Differential privacy limits the additional information attributable to one person's participation, even when an attacker has auxiliary information.

## Epsilon in plain language

Describe epsilon as the **one-person influence limit**:

- `epsilon = 0` means the output distribution cannot depend on the source patients at all;
- smaller positive epsilon means stronger protection and less allowed influence;
- larger epsilon permits more influence and therefore weaker protection;
- a mechanism can technically satisfy differential privacy with a very large epsilon while providing weak practical protection.

Require explicit epsilon. Do not choose or market a universal default.

## Delta in plain language

Describe delta as a very small additive allowance in the probability bound. It permits a small amount of output probability that is not controlled solely by the multiplicative `exp(epsilon)` term.

Do not describe delta as:

- the probability that a patient is re-identified;
- the percentage of patients who are unprotected;
- a direct breach probability.

Smaller delta is stronger. Permit `delta = 0` where supported. Require explicit justification and validation of any positive delta.

## Is privacy protected?

The documentation must give a direct answer:

> Privacy is formally protected only when output is generated entirely from a valid subject-level differentially private fitted model, with the reported epsilon, delta, contribution bounds, public inputs, cumulative accounting, and implementation assumptions. The protection is mathematically bounded rather than absolute.

Do not use `private`, `anonymous`, `safe`, or `de-identified` as unqualified binary labels.

# Public API

## Explicit PMX roles

Implement or revise an API resembling:

```r
roles <- pmx_roles(
  id = "ID",
  time = "TIME",
  nominal_time = NULL,
  tad = NULL,
  occasion = NULL,
  dv = "DV",
  amt = "AMT",
  evid = "EVID",
  cmt = "CMT",
  dvid = NULL,
  mdv = NULL,
  rate = NULL,
  cens = NULL,
  limit = NULL,
  covariates = "WT",
  subject_properties = "ACTARM",
  assigned_dose = NULL,
  exclude = NULL
)
```

Critical roles must be explicit. Do not silently infer column meaning from names.
Ordinary covariates are baseline attributes generated from marginal summaries.
Subject properties such as `ACTARM`, `TRT`, or a nominal dose group are
categorical subject-level assignments whose released strata must condition the
generated regimen. An occasion-varying nominal assigned-dose column must be
derived from the generated event amount, not sampled independently.

## Endpoint behavior specification

Allow the user to describe the scientific clock of each endpoint without specifying a PK or PD model:

```r
endpoints <- list(
  cp = pmx_endpoint(
    dvid = "cp",
    alignment = "dose_relative",
    transform = "log",
    shape = "occasion"
  ),
  response = pmx_endpoint(
    dvid = "response",
    alignment = "study_time",
    transform = "auto",
    shape = "global"
  )
)
```

Support at least:

- `alignment = "dose_relative"`: behavior organized by TAD after each dose;
- `alignment = "study_time"`: one trajectory over time from study or treatment start;
- `alignment = "occasion"`: separate but related occasion profiles;
- `alignment = "hybrid"`: dose-related excursions around a longer-term baseline.

Endpoint alignment is preferably user-declared scientific metadata. If the package attempts to infer endpoint behavior from confidential data, that decision must be differentially private and consume privacy budget.

## Fit once, generate many times

Use a two-stage API resembling:

```r
private_model <- fit_private_pmx(
  data,
  roles,
  endpoints,
  epsilon,
  delta,
  bounds,
  public_design = NULL,
  contribution_limits,
  budget_allocation
)

mock <- generate_pmx(
  private_model,
  n_subjects = NULL,
  seed = 123
)

privacy_report(private_model)
```

Names may change if a clearer API results, but preserve the separation:

- `fit_private_pmx()` is the only stage that reads confidential source data and consumes privacy budget.
- `generate_pmx()` reads only the already-private model.
- Generating additional datasets from the same fitted model is post-processing and consumes no additional budget.
- Refitting against the confidential data consumes additional privacy budget.
- An ordinary seed may control generation from the private model.
- By default, generation uses the fitted privacy-accounted subject-count
  release. This equals the source cohort size for an explicitly public fixture;
  for confidential data it is the differentially private count, not an exact
  unbudgeted disclosure. An explicit `n_subjects` remains an override.
- Privacy noise during fitting must not be user-seeded, logged, or exposed.

Do not provide a one-step convenience function that encourages users to unknowingly refit repeatedly unless it returns the private model, emits a release-ledger record, and explains the additional budget consumption.

# Private fitting algorithm

## Public configuration and bounds

Strongly prefer configuration established without inspecting confidential patient values:

- endpoint definitions and units;
- plausible continuous bounds;
- maximum rows, doses, occasions, and observations contributed per subject;
- possible factor levels;
- expected schema;
- allowed transformations;
- any protocol information already considered public.

The actual data must still calibrate approximate magnitudes and behavior through private mechanisms. Public configuration defines safe domains; it does not replace source calibration.

If a bound, category set, nominal grid, dose class, censoring threshold, endpoint type, schema element, or subject count is learned from the confidential data, include that computation in the privacy budget.

Do not automatically use exact source minima and maxima as bounds. Exact extrema are particularly sensitive in small studies.

## Contribution bounding

Before computing any private statistic, deterministically bound each subject's total contribution:

- maximum number of rows;
- maximum number of doses and occasions;
- maximum observations per endpoint and occasion;
- maximum number of timing cells;
- bounded dose, rate, time, covariate, DV, and censoring-limit domains;
- bounded norms for any subject-level feature or trajectory-coefficient vectors.

All sensitivity calculations must protect the subject's complete bounded contribution, not a single row.

Do not release which subjects were clipped or truncated. Release counts of affected subjects only if those counts are public or privately computed.

## Fixed-dimensional subject representation

Convert each subject into a bounded representation on endpoint-appropriate grids. This representation exists only inside the restricted fitting process and is never released unprivatized.

Possible components include:

- generalized regimen class;
- dose and infusion summaries;
- endpoint/occasion observation-presence indicators;
- bounded observation counts conditional on an occasion being sampled;
- categorical and continuous baseline covariates;
- categorical subject-property strata jointly represented with regimen
  summaries;
- censoring indicators and public/private limits;
- low-dimensional endpoint trajectory coefficients.

Use a deliberately low-dimensional representation. The package does not need to reproduce detailed source statistics.

## Privately learned information

Privately estimate only the coarse information required for plausible workflow data:

- approximate log-scale endpoint magnitude;
- a small number of endpoint shape coefficients;
- coarse between-subject and residual scales;
- generalized dose levels, intervals, rates, and infusion durations;
- generalized sampling-time cells, occasion-presence probabilities, and
  conditional sample counts;
- broad covariate centers, scales, or category proportions;
- coarse frequencies of declared treatment/property strata and their
  associated generalized regimens;
- censoring frequency and applicable limit classes.

Do not privately learn high-dimensional joint distributions merely because they are available. Prefer low-order marginals, bounded moments, or a small public basis.

## Population variability

It is acceptable—and often preferable—to generate broad variability from public generator settings rather than estimate it precisely from a six- or twelve-subject study.

For example:

- privately learn the approximate log-DV center;
- use a public broad distribution for subject variability;
- privately learn a coarse mean curve;
- use public residual and smoothness settings;
- privately learn whether sampling is hourly, daily, or weekly;
- generate new actual times around generalized nominal cells.

This preserves visual plausibility without spending privacy budget on detailed statistical fidelity.

# Endpoint trajectory generation

## Dose-relative PK-like endpoints

For `alignment = "dose_relative"`:

1. Assign source observations internally to the most recent qualifying dose.
2. Construct or privately learn a generalized TAD grid.
3. Represent each subject's within-dose behavior using a small fixed public basis or a few bounded shape features.
4. Privately estimate an approximate population magnitude and coarse rise-and-fall shape.
5. Generate a new within-dose excursion after every generated dose.
6. Permit broad generated subject variability, residual variability, predose baseline, and accumulation behavior.
7. Restart serial residual perturbations at each generated occasion. When the
   released coarse curve is already approximately unimodal, use source-free
   post-processing to keep each perturbed profile single-peaked; residual noise
   must not manufacture a second absorption peak.

The resulting mock PK data should plausibly rise and fall after each dose when the privately learned coarse behavior supports that pattern. It need not preserve the source half-life, exposure, accumulation ratio, or parameter distribution.

## Study-time PD or biomarker endpoints

For `alignment = "study_time"`:

1. Align source observations internally to study or treatment time.
2. Represent each subject with the fixed public low-dimensional grid basis.
3. Privately estimate approximate magnitude and a few coarse trajectory coefficients.
4. Generate one continuous subject trajectory over the study.
5. Permit delayed onset, peak/nadir, plateau, and recovery when supported by the private coarse model.

Do not restart a study-time endpoint after each dose.

## Hybrid endpoints

For `alignment = "hybrid"`, generate a smooth study-time baseline plus dose- or occasion-related excursions. Keep the representation low-dimensional and budgeted.

## Multiple DVIDs

Process endpoints separately according to their declared alignment and transformation. Do not place PK and PD measurements into one undifferentiated trajectory model.

Cross-endpoint relationships are optional and should generally come from public generator assumptions unless a small, explicitly budgeted private relationship is necessary for workflow usefulness.

# Event, dose, and timing generation

## Sampling-design inference

Dose events and observation rows are different parts of the PMX design. A
Q24H regimen must not, by itself, cause an endpoint to be sampled after every
dose. When an occasion-specific sampling schedule is not independently public,
infer it through the timing privacy budget rather than repeating one endpoint
grid at every generated occasion.

For each dose-relative or occasion-aligned endpoint, construct fixed-dimensional
bounded subject features for every allowed occasion:

- a binary indicator that the endpoint was sampled on that occasion;
- the observation count divided by the public per-endpoint/per-occasion limit;
- coarse within-occasion timing-cell presence on a fixed public basis.

Privately release the vector sums with subject-level sensitivity accounting.
Decode the first feature as the population probability that the occasion is
sampled. Decode the second, using the released presence denominator, as the
mean number of observations conditional on that occasion being sampled. Do not
use an unconditional mean count as a presence probability: an occasion sampled
densely in a minority of subjects is not equivalent to a sparsely sampled
occasion in every subject.

During generation, first draw whether each endpoint/occasion is active, then
allocate observations according to the released conditional count and coarse
timing cells. Preserve the separately generated dose regimen even when an
occasion has no observations. A public occasion schedule may override this
model only as an exceptional protocol override when its public status is
independently justified and recorded. The default workflow and every
named-data demonstration must omit source-derived dose times, dose counts,
intervals, amounts, rates, infusion durations, endpoint grids, and occasion
schedules. The fit must infer regimen quantities and timing occupancy through
their allocated private summaries. An automatic source-independent fixed grid
is a discretization basis required for a finite-dimensional release; it is not
a claimed or disclosed visit schedule.
When the allocated count is shorter than the timing grid, select cells using
the released cell-presence probabilities. Never implement count matching by
taking the first N grid cells, because that systematically erases late PK
follow-up. Public demonstrations must fail if cohort counts differ or if mean
endpoint time points per patient or bounded time-range coverage differ
materially between source and synthetic data.

Diagnostics and demonstrations must label dose/event rows separately from
samples. Study-time plots may connect a subject's observed points across gaps,
but must not add event rows as observations or imply that intermediate samples
exist. Dose-relative plots should retain occasion grouping because their time
axis resets after each dose. Provide a releasable post-processing summary of
fitted occasion probabilities and conditional sample counts.

## Event structure

Construct coherent PMX records from the private population model and public structural rules. Do not copy source row blocks.

Generate or derive together:

- EVID;
- AMT;
- RATE or infusion duration;
- CMT;
- DVID;
- ADDL and II when supported;
- observation rows;
- MDV;
- censoring fields;
- tied-time ordering.

Do not numerically average event-control fields.

## Nominal and actual time

Support both actual and nominal time when available. If no nominal column exists, privately learn a coarse timing grid using prespecified bins or another budgeted mechanism.

Generate chronological actual-like times around the private/generalized nominal design while preserving:

- nondecreasing within-subject time, or within subject and declared occasion
  when the input clock resets by occasion;
- correct dose occasion;
- postdose observations after the qualifying dose;
- tied collection-time blocks;
- dose/predose ordering;
- nonnegative time where required;
- consistent TAD and occasion derivations.

Never copy a subject's complete time vector, missing-visit pattern, or infusion duration.

## Doses and infusions

Learn only generalized, private distributions of dose amounts, intervals, rates, and durations. Generate new coherent regimens.

When a declared categorical subject property describes treatment assignment,
release its stratum count and bounded event/regimen vector together so one
subject contributes to only one stratum. Generate the property and regimen as
one post-processed draw. Public category levels define the finite domain; they
must not be extracted from confidential records outside the private fit.

When a nominal dose column varies by occasion, regenerate it from that
occasion's positive generated AMT and require it to be complete and constant
within subject/occasion. This guarantees internal consistency but does not, by
itself, recover a crossover-sequence distribution when no sequence/arm
property is declared.

Rare source regimens must not be reproduced merely because they occurred. A regimen may be generated when it is:

- part of declared public protocol information; or
- represented through a valid private mechanism with adequate reliability.

Derive infusion stop times from generated starts and generated durations. Preserve internally coherent AMT/RATE/CMT/EVID conventions.

# Censoring support

## Monolix-style roles

Support:

- `CENS = 0`: uncensored;
- `CENS = 1`: left censored, with DV containing the upper boundary such as LLOQ;
- `CENS = -1`: right censored, with DV containing the lower boundary such as ULOQ;
- optional `LIMIT`: the other interval boundary.

Limits may vary by endpoint, occasion, or row. Treat public assay limits as public only when their status is justified; otherwise learn or classify them privately.

## Private censoring model

Privately estimate coarse censoring frequency or applicable limit class by endpoint and broad time region. Generate a latent mock trajectory first, then derive DV, CENS, and LIMIT coherently from the generated value and generated/public limits.

Never blend, average, or independently perturb CENS values.

## Censoring validation

Validate:

- allowed CENS values;
- boundary direction;
- finite limits where required;
- lower boundary not exceeding upper boundary;
- DV equals the reported censoring boundary under the Monolix convention;
- consistency among EVID, MDV, DV, DVID, CENS, and LIMIT.

# Identifiers, schema, and unreleased source information

- Generate entirely new IDs of a practical type.
- Never include source IDs in the fitted private model.
- Reject or explicitly remove direct identifiers.
- Do not export exact calendar dates or datetimes.
- Stop on unmodeled Date/POSIXct columns unless explicitly excluded or privately handled.
- Treat factor levels, schema, and requested output size as public only when justified.
- Preserve the declared public schema, names, order, practical classes, and factor levels in generated output.

The fitted private model must not contain:

- raw rows;
- raw subject profiles;
- source IDs;
- exact event templates;
- raw residuals;
- unnoised aggregates;
- confidential category sets or bounds;
- debugging caches that depend nonprivately on source data.

# Privacy accounting and implementation safety

## Budget allocation

Account for every source-dependent released computation, including:

- subject count if not public;
- category or regimen selection;
- timing grids;
- visit frequencies;
- dose and covariate distributions;
- endpoint magnitude and shape summaries;
- censoring summaries;
- data-dependent bounds or hyperparameter choices;
- released source-versus-synthetic diagnostics.

The composed realized privacy loss must not exceed the requested epsilon and delta.

## Repeated releases

Repeated calls to `fit_private_pmx()` against the same confidential dataset compose and consume additional privacy budget. Produce a machine-readable release-ledger entry for every fit. Organizational governance must combine this release with prior fits and other private releases.

Generating additional datasets from one fitted private model does not consume additional privacy budget.

## Validated backend

Do not hand-code production Laplace or Gaussian mechanisms using ordinary R random-number functions. Use a reputable, tested differential-privacy backend with an R interface where feasible, such as OpenDP, isolated behind a small adapter.

Document:

- backend and version;
- adjacency definition;
- mechanism and accountant;
- epsilon and delta;
- contribution bounds;
- proof assumptions;
- handling of floating-point and random-number issues.

Fail closed when the validated backend is unavailable. Never silently fall back to ordinary random noise.

## Post-processing discipline

Once a valid private model has been fitted, generation, constraint repair, plotting, and validation may use that model without additional privacy cost, provided they do not consult the confidential source data.

`compare_pmx(source, mock)` is not automatically releasable merely because `mock` is private. Any source-derived diagnostic that leaves the restricted environment must itself be public or privately computed and budgeted.

# Small-study behavior

The method must operate mathematically with six or twelve subjects, but it must not promise detailed fidelity. One person represents a large portion of a small interim dataset, so strong subject-level privacy can substantially perturb learned summaries.

Design the generator to remain useful by learning only:

- broad order of magnitude;
- coarse endpoint shape;
- generalized sampling and dosing scale;
- minimal variability and censoring information.

Allow the generator to warn or stop when dimensionality, bounds, privacy budget, or effective sample size cannot support even these goals. Do not weaken the privacy unit or silently increase epsilon to make examples look better.

Use the three small public `nlmixr2data` examples to verify structure and algorithm behavior. Add a larger fully simulated public fixture for privacy-utility evaluation.

# Validation and comparison

Implement or revise:

```r
validate_pmx(data, roles)
validate_private_model(private_model)
privacy_report(private_model)
compare_pmx(source, mock, roles)
```

`compare_pmx()` is a restricted-environment diagnostic unless its source-derived output is privatized. Mark every returned component with its release status.

Generated-data validation should cover:

- schema and classes;
- PMX event coherence;
- chronological and tied-row ordering;
- endpoint alignment;
- repeated-dose PK-like behavior where configured;
- global PD-like behavior where configured;
- censoring consistency;
- finite and defensible domain values;
- new IDs;
- absence of direct identifiers.

# Dataset-specific demonstrations

## `nlmixr2data::theo_md`

Demonstrate:

- repeated Q24H dosing;
- dose-relative PK alignment;
- concentrations that rise and fall following generated doses;
- new sampling-time vectors;
- correct order of magnitude without claiming distributional fidelity.

## `nlmixr2data::warfarin`

Demonstrate:

- lower-case schema;
- separate `cp` PK and `pca` PD endpoints;
- dose-relative PK behavior;
- study-time PD behavior;
- generated endpoint-specific observation presence;
- factor and class preservation.

## `nlmixr2data::wbcSim`

Demonstrate:

- coherent infusion starts and stops;
- a study-time delayed WBC decline, nadir, and recovery;
- no exact reproduction of singleton source regimens;
- generated follow-up times and observation counts.

## `nlmixr2data::nimoData`

Demonstrate:

- four public nominal-dose property levels that condition inferred amount,
  rate, duration, and ten-dose regimen summaries;
- declared TAD/OCC sampling and a final washout profile that may extend beyond
  one dosing interval;
- coherent generated infusion start/stop records even though source duration
  is represented by positive AMT/RATE rather than explicit stop rows; and
- explicit exclusion of longitudinal WGT because the prototype supports only
  baseline covariates.

## `nlmixr2data::mavoglurant`

Demonstrate:

- TIME reset and validation within subject/occasion;
- one- and two-occasion dose-relative profiles;
- an occasion-specific assigned DOSE reconstructed from generated positive
  AMT and constant within subject/occasion; and
- numeric-coded SEX treated categorically through a public category domain.

## Censoring fixture

Add a small fully public fixture demonstrating uncensored, left-censored, right-censored, and interval-censored rows with endpoint-specific limits.

# Tests and acceptance criteria

Use `testthat`. Include internal fixtures so core tests do not require `nlmixr2data`.

At minimum test:

1. The only synthesis engine is the private generator; no raw blend/template mode remains public.
2. Neighboring datasets differ by one complete subject.
3. Contribution limits and clipping occur before every private aggregate.
4. Composed accounting does not exceed requested epsilon and delta.
5. Missing, invalid, or incomplete privacy configuration fails closed.
6. Missing DP backend fails closed.
7. The private fitting noise is not user-seeded or returned.
8. Generation from a fitted private model is reproducible by seed.
9. Repeated generation does not change privacy accounting.
10. Every fit creates a release-ledger entry.
11. The fitted model contains no raw IDs, rows, subject profiles, templates, residuals, or unnoised aggregates.
12. Public versus private configuration is explicit and recorded.
13. Source-derived diagnostics are marked nonreleasable unless privately budgeted.
14. Direct identifiers and unmodeled datetimes fail clearly.
15. New IDs have the expected practical type.
16. Dose-relative endpoints rise/fall after generated doses in suitable fixtures.
17. Study-time endpoints form a single long-timescale trajectory rather than restarting after each dose.
18. Hybrid and occasion alignment behave as documented.
19. Multiple DVIDs remain separate.
20. Generated times, TAD, occasions, dosing, and infusion records are coherent.
21. Censoring states and limits are coherent.
22. Output schema, classes, column order, and factor levels are restored from declared public schema.
23. Same generation seed gives identical output from one private model; different seeds change output.
24. Repeated-dose sampling inference distinguishes occasion activation from
    conditional sample density and does not create observations merely because
    a dose event exists.
25. A sparse-design integration fixture recovers inactive occasions (for
    example, intensive first/final profiles with no samples after intervening
    doses), while retaining all required dose events.
26. Small-study and high-dimensional limitations warn or stop as documented.
27. `theo_md`, `warfarin`, `wbcSim`, `nimoData`, and `mavoglurant`
    examples run end to end.
28. A larger simulated fixture supports privacy-utility evaluation.
29. Empirical privacy audits detect an intentionally broken test mechanism while documentation states that auditing is not a proof.
30. The selected DP backend's canonical mechanism tests pass.
31. Omitted `n_subjects` generates the privacy-accounted fitted cohort size,
    and public demonstrations therefore use the same number of source and
    synthetic subjects.
32. Endpoint sampling checks compare observations per subject and time-range
    coverage; trimming a dense public grid must not systematically discard
    late-time cells.
33. Every demonstration comparison facets source data above synthetic data,
    while retaining consistent dataset colors and endpoint separation.
34. Named-data demonstrations supply no source-derived regimen or sampling
    schedule. With automatic generic grid bases, fitting infers dose count,
    interval, amount, infusion behavior, occasion activation, conditional
    counts, and timing-cell occupancy from the input through budgeted releases.
35. If a released dose-relative curve is approximately unimodal, generated
    intensive profiles contain at most one directional peak per occasion;
    unoccupied cells and serial noise must not create secondary PK peaks.
36. Declared subject properties are complete and constant within subject, and
    their generated regimen is drawn from the same privacy-accounted stratum.
37. An assigned-dose column is finite, constant within subject/occasion, and
    equal to the positive generated event amount.
38. Numeric covariates with declared public category levels are generated as
    categories rather than continuous values.
39. Positive-rate infusions without an explicit source stop row use bounded
    AMT/RATE duration inference and generate coherent start/stop pairs.
40. Declared occasion clocks may reset; row-order validation remains strict
    within each subject/occasion.

# Required vignettes

Create exactly four complementary package vignettes. Their purposes must
remain distinct: practical use, privacy, simulation mechanics, and epsilon
exploration.

## 1. `vignettes/pmxSynthData-demo.Rmd`

Title: **“Using pmxSynthData”**

This is a practical, example-driven vignette. It should help a pharmacometrician generate and inspect mock data without first understanding every algorithmic detail.

Include:

1. Installation and setup.
2. A minimal end-to-end workflow:
   - declare roles;
   - declare endpoint alignment;
   - provide bounds and privacy parameters;
   - fit the private model once;
   - inspect `privacy_report()`;
   - generate one or more mock datasets;
   - validate and plot them.
3. A clear warning that illustrative epsilon/delta values are examples, not universal recommendations.
4. Examples for `theo_md`, `warfarin`, `wbcSim`, `nimoData`, and
   `mavoglurant`.
5. Endpoint-faceted plots showing dose-relative PK versus study-time
   PD/biomarker behavior, with source panels above synthetic panels.
6. Repeated generation from the same private model with different ordinary generation seeds.
7. The censoring fixture.
8. Common errors and how to fix them:
   - missing bounds;
   - undeclared time roles;
   - unmodeled datetimes;
   - too many requested private summaries;
   - insufficient privacy budget;
   - attempting to refit unnecessarily.
9. A concise checklist before exporting generated data.

Keep mathematical detail light. Link separately to the privacy introduction
and simulation-method vignette.

## 2. `vignettes/pmxSynthData-privacy-intro.Rmd`

Title: **“Privacy in pmxSynthData”**

This is a rigorous but beginner-accessible introduction to privacy. Write for
pharmacometricians who may know nothing about differential privacy. Keep the
statistical-generation mechanics at overview level and link to the dedicated
simulation-method vignette.

The explanation must proceed in layers: intuition first, then equations and implementation details.

### Required beginner privacy section

Explain:

1. Removing names is not sufficient because unusual values, schedules, or external information may permit linkage or attribute inference.
2. Revealing that an unusual record existed is not necessarily the same as identifying its owner, but it can become identifying when combined with outside information.
3. The one-person-influence interpretation of subject-level differential privacy.
4. The two-world comparison: the source with one person's complete record versus without it.
5. Why `epsilon = 0` would prevent using patient data at all.
6. Epsilon as the allowed one-person influence, with a small table of `epsilon` and `exp(epsilon)` values.
7. Delta as a small additive allowance in the probability bound—not the probability of re-identification and not the fraction of unprotected patients.
8. A simple numerical probability example containing both epsilon and delta.
9. Why privacy is genuinely protected but not absolute.
10. What membership, attribute disclosure, linkage, and re-identification mean and how they differ.
11. Why small Phase 1 studies have a harder privacy-utility tradeoff without losing the mathematical guarantee.
12. Composition: why refitting consumes budget and generating again from the same private model does not.
13. What DP does and does not protect.
14. Why legal anonymity and release authorization remain separate questions.

Include this direct statement:

> Output generated entirely from a valid private model is formally subject-level `(epsilon, delta)`-differentially private. This means one person's complete record has a mathematically bounded influence on the release. It does not mean zero disclosure risk or guaranteed impossibility of re-identification.

### Required technical privacy section

Explain and document:

1. The privacy unit and neighboring-dataset definition.
2. The formal `(epsilon, delta)` inequality and every symbol.
3. Public configuration versus confidential data-dependent quantities.
4. Contribution bounding and sensitivity.
5. Privacy mechanisms and budget allocation.
6. Composition and post-processing.
7. The validated DP backend and implementation assumptions.
8. The fixed-dimensional subject representation.
9. The categories of privately estimated population summaries.
10. Small-sample privacy-utility limitations.
11. Privacy ledger, repeat fits, and releasable versus restricted diagnostics.

Use equations where they materially clarify the implementation, but define every symbol immediately. Include:

- a full pipeline diagram;
- a six-patient worked privacy example;
- a table distinguishing population information, one-person information, identification, and linkage;
- a table distinguishing releasable private outputs from restricted diagnostics;
- a final practical interpretation checklist.

Do not describe the Version 1 AVATAR algorithm except, at most, in a brief historical note explaining why raw anchors and donors were removed.

## 3. `vignettes/pmxSynthData-simulation-method.Rmd`

Title: **“How pmxSynthData Simulates Patients”**

This document describes the implemented statistical generator in enough detail
for a pharmacometrician to understand exactly where event rows, sample times,
endpoint values, variability, covariates, and censoring come from. State near
the top that the prototype fits neither an ODE/NLME model nor splines.

Explain and document:

1. The fit-once/generate-many boundary, with a link to the privacy vignette.
2. Public roles, bounds, transforms, alignments, and automatic generic grid
   bases; distinguish those bases from inferred sampling schedules.
3. Contribution-bounded fixed-dimensional event, timing, trajectory,
   covariate, and censoring representations.
4. Nearest-grid assignment and subject-level cell presence/value summaries.
5. The exact one-pass 1--2--1 smoother and piecewise-linear
   `stats::approx()` interpolation; explicitly explain why this is not a
   spline.
6. Dose-relative, occasion, study-time, and hybrid clocks.
7. The separate sampling model:
   - timing-cell presence probability;
   - dose-occasion activation probability;
   - sample count conditional on an active occasion;
   - probability-weighted cell selection when count matching;
   - exceptional independently public occasion-schedule overrides.
8. Default cohort size from the privacy-accounted subject-count release and
   explicit public overrides.
9. Generalized regimen, dose/event, and infusion-pair generation.
10. Nominal-time jitter, tied collections, chronological constraints, TAD, and
    occasion reconstruction.
11. Exact public subject, occasion, and AR(1) variability multipliers.
12. Numeric/categorical covariate generation and the absence of a learned
    covariate-response model.
13. Subject-property/regimen strata, occasion-assigned dose reconstruction,
    and their public-category requirements and limitations.
14. CENS/LIMIT reconstruction, schema restoration, new IDs, and validation.
15. Expected structural fidelity versus quantities that are not preserved.
16. Compact implementation pseudocode and current limitations.

Include a pipeline diagram, a fixed-grid/interpolation figure, a sampling-time
example with early and late cells, equations that define all symbols, and a
table distinguishing expected versus unsupported fidelity.

## 4. `vignettes/pmxSynthData-epsilon-exploration.Rmd`

Title: **“Exploring Epsilon in pmxSynthData”**

Use the production OpenDP adapter to show `theo_md`, `warfarin`, and `wbcSim`
at three clearly labeled illustrative epsilon values. For every dataset,
display one epsilon per panel with source facets above synthetic facets,
endpoint separation, connected observed profiles, consistent colors, and
linear DV axes. Hold the displayed cohort size fixed as explicit public
post-processing and separately report the noisy subject-count release.

Explain that:

- the values are relative engineering examples, not recommendations;
- smaller epsilon gives stronger privacy and generally more distortion;
- one noise draw need not improve monotonically in every coordinate;
- privacy-noise fits are deliberately not seedable;
- all three source datasets are already public;
- performing and releasing all three fits on one confidential source would
  compose their privacy losses; and
- the public-fixture backend cannot be used for epsilon exploration because it
  is intentionally noiseless and makes no DP claim.

The vignette must remain buildable when optional dependencies are unavailable,
must fail closed rather than substitute ordinary random noise, and must show an
explicit empty/failure panel if a valid high-noise fit cannot generate a table.

## Vignette verification

For all four vignettes:

- include valid vignette metadata;
- ensure examples match the implemented API exactly;
- render them;
- visually inspect the HTML output;
- check tables, equations, figures, captions, links, and code output;
- eliminate stale Version 1 claims;
- avoid unnecessary dependencies;
- ensure package check builds all four vignettes successfully.

# Other documentation

Update README, roxygen documentation, examples, and `NEWS.md` to match Version 2.

The README should lead with:

- source-calibrated mock data for workflow exploration;
- subject-level differential privacy;
- fit-once/generate-many architecture;
- broad PK/PD behavior rather than scientific fidelity;
- explicit limitations.

Retire or rewrite old AVATAR-focused `METHODS.md` instructions so future developers are not directed to recreate the obsolete algorithm.

# Independent review

Before claiming production-ready differential privacy:

1. Inventory every confidential-data access.
2. Write a mechanism-level privacy argument and budget composition.
3. Audit serialization for raw-data leakage.
4. Obtain review from a differential-privacy specialist independent of the implementation.
5. Use empirical privacy attacks only as bug-finding supplements, never as substitutes for the formal argument.
6. Complete the organization's privacy, legal, information-security, and data-governance review.

The package may implement a mathematical guarantee. It cannot by itself declare a release legally anonymous or organizationally authorized.

# Scope control

Do not add:

- an AVATAR/blend/template synthesis mode;
- a PK, PD, Friberg, or NLME model;
- CTGAN, TVAE, or a deep generative framework;
- scientific model fitting or selection;
- high-dimensional private learning merely to improve resemblance;
- claims of k-anonymity as a formal guarantee;
- hand-written production DP noise mechanisms;
- a Shiny application.

Prefer a small number of transparent, testable, low-dimensional private summaries.

# Verification and final response

Run:

- documentation generation;
- all unit and integration tests;
- canonical DP-backend tests;
- private-model serialization and leakage tests;
- accounting and release-ledger tests;
- all four vignette renders and visual inspections;
- `R CMD check` or `devtools::check()`.

Do not suppress legitimate warnings merely to obtain a clean check.

At completion, report concisely:

- the implemented public API;
- files created, rewritten, or removed;
- the endpoint-alignment and trajectory-generation approach;
- the privacy unit, epsilon, delta, contribution bounds, budget accounting, backend, and proof assumptions;
- behavior on six- and twelve-subject fixtures;
- what works for each public example dataset;
- censoring support;
- exact test, vignette-render, and package-check results;
- independent-review status;
- remaining privacy and scientific-fidelity limitations;
- any deviations from this specification and why.
