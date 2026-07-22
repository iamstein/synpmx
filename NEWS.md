# pmxSynthData 0.0.0.9000

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
