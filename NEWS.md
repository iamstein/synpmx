# pmxSynthData 0.0.0.9000

## Version 3 calibrated structural generator

* Added a low-dimensional generator built around a *public* structural model.
  `pmx_structural_model()` declares one-compartment IV, oral, or infusion PK
  with optional direct-effect or indirect-response PD, using built-in analytic
  solutions that need no compiler. `pmx_trial_design()` declares dose levels,
  cohorts, and the protocol sampling schedule. Both require a `source` string
  recording data-independent provenance.
* Added `pmx_generate()`, which works in two modes from the same code path.
  Given a structural model it generates entirely from public inputs, reading no
  confidential data and making no privacy claim. Given a calibrated model it
  uses the privately corrected parameters, which is post-processing and
  consumes no further budget.
* Added `fit_calibrated_pmx()`, which releases a small multiplicative
  correction to the model's prediction rather than absolute parameters. Each
  subject is summarized by non-compartmental analysis on that subject's own
  rows, so clipping to a public prior bounds the per-subject sensitivity at one.
  A prior on how wrong a preclinical prediction is spans roughly eight-fold and
  costs nothing, where a prior on absolute clearance for a new compound spans a
  hundred-fold.
* Added `pmx_prior()` and `pmx_priors()`. A prior range is now the dominant
  sensitivity term, so `source` provenance is mandatory.
* Added `pmx_preflight()`, which reports `f = d / (epsilon * N)` and a verdict
  before any budget is spent. `f` is the fraction of the prior's width that
  survives as noise, so it answers whether a release would beat the prior at
  all. `fit_calibrated_pmx()` warns when it would not.
* A released correction pressed against its prior boundary is reported as a
  diagnostic: the prior was probably wrong and the release is censored.
* `typical` parameter values are the median of the lognormal population, the
  usual population-PK convention. They previously centred the arithmetic mean
  when generating while calibration estimated the geometric mean, leaving a
  systematic `exp(sigma^2 / 2)` gap between fitting and generation that did not
  shrink with cohort size or privacy budget.
* `pmx_preflight()` caps its reported fold-error at the prior's half-width.
  Clipping prevents a release from landing outside the prior, so the error
  saturates rather than diverging; the uncapped formula reported absurd values
  in exactly the regime where the release is worthless anyway.
* Added two-compartment PK, `"2cmt_iv"` and `"2cmt_oral"`, as analytic
  solutions verified to conserve `Dose/CL`.
* The PD correction now solves for the `emax` reproducing the released effect
  rather than multiplying `emax` by it. Response magnitude saturates in `emax`,
  so multiplying overshoots badly. The per-subject statistic is a signed area
  between response and baseline; a peak is biased upward by residual noise, and
  an absolute-deviation area accumulates noise on the observed side while the
  prediction carries none.
* The PD correction is experimental. It is exact without residual error but
  biased low when the response deviation is small relative to residual error,
  because a geometric mean of noisy per-subject ratios sits below the ratio of
  their means. PK is largely unaffected. See `design/FEASIBILITY.md`.
* Covariate columns are explicitly out of scope. Generated output exercises
  event-table, timing, dosing, and censoring code, but not covariate handling.

## Fixes

* Fixed released presence fields being decoded against a fixed constant as
  though they were probabilities. They are unnormalized subject counts, so
  privacy noise around zero passed the "this cell had support" test and the
  cell then decoded to the bottom of the endpoint working domain, producing
  spurious deep troughs on log-scale endpoints. Decoding now gates on a
  threshold derived from the release's own mechanism scale
  (`sensitivity / epsilon`). At epsilon 5 this removes the artifact entirely
  for cohorts of about 600 subjects or more and halves population-curve error
  at 2000 subjects. Noiseless public-fixture fits are unchanged by
  construction. Tracked as `SIM-020` in `design/TEST_SIM.md`.

## Version 2 private population generator

* Replaced the source-subject synthesis architecture with a fit-once,
  generate-many subject-level differential-privacy design.
* Added explicit endpoint clocks, public schema/design declarations, numeric
  bounds, contribution limits, and budget allocation.
* Added an OpenDP adapter that fails closed when unavailable; privacy noise is
  neither user-seeded nor returned. A guarded public-fixture backend supports
  only data explicitly asserted to be public.
* Made the OpenDP Laplace scale conservatively robust to floating-point
  rounding for fractional budget allocations and added a canonical fractional
  budget backend check.
* Added fixed-dimensional bounded subject summaries, basic composed accounting,
  machine-readable release ledgers, privacy reports, and private-model leakage
  validation.
* Added new event-table generation with nominal/actual time, TAD, occasion,
  repeated doses, infusion start/stop pairs, multiple DVIDs, and schema/class
  restoration without source event-row copying.
* Fixed repeated-dose generation to honor the released per-subject observation
  total instead of repeating a complete endpoint grid after every dose.
* Interpolate across unoccupied trajectory cells during release
  post-processing instead of treating them as midpoint-scale measurements;
  this removes artificial troughs and secondary peaks on log-scale PK curves.
* Restart serial perturbations at each generated dose occasion. When a
  released dose-relative population curve is already approximately unimodal,
  project each perturbed profile back to a single rise-and-decline shape so
  individual noise cannot introduce spurious secondary absorption peaks.
* Added a two-part, privacy-accounted sampling model by dose occasion: sampling
  probability plus conditional observation count. `sampling_summary()` exposes
  the fitted design. All named-data demos now omit regimen and sampling
  schedules, use source-independent automatic grid bases, and infer dose and
  visit behavior inside the private fit. The Theophylline example reports its
  inferred visit probabilities explicitly.
* Updated demonstration figures to connect each subject's chronological
  observations and use endpoint-specific linear DV axes.
* Made `generate_pmx()` default to the fitted privacy-accounted cohort size.
  Timing-count trimming now respects fitted timing-cell probabilities rather
  than systematically keeping early cells and deleting late PK follow-up.
* Added Monolix-style uncensored, left-, right-, and interval-censoring support.
* Added a deterministic 60-subject public repeated-dose fixture for
  privacy-utility evaluation.
* Marked every source-derived comparison component as restricted unless
  separately privatized.
* Added a shared simulation-evaluation registry, deterministic hard gates, and
  a multi-seed `scripts/evaluate_simulations.R` report runner covering every
  demo dataset and the accumulated regression history in `design/TEST_SIM.md`.
* Added explicit subject-property/regimen strata for treatment-arm fields such
  as ACTARM, TRT, and nominal dose group, plus occasion-assigned dose
  reconstruction from generated AMT. Numeric-coded public categories are now
  modeled categorically when levels are declared.
* Added `nlmixr2data::nimoData` and `nlmixr2data::mavoglurant` demonstrations
  and regression gates, including reset occasion clocks, positive-rate
  duration inference, treatment-group coherence, and terminal washout timing.
* Rewrote tests, README, and scripts for Version 2. Split the documentation
  into practical, privacy-introduction, detailed simulation-method, and formal
  epsilon-exploration vignettes, including five practical nlmixr2data demos, a
  public censoring fixture, and empirical bug-finding tests that are explicitly
  not privacy proofs.
