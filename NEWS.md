# pmxSynthData 0.0.0.9000

## Version 2 private population generator

* Replaced the source-subject synthesis architecture with a fit-once,
  generate-many subject-level differential-privacy design.
* Added explicit endpoint clocks, public schema/design declarations, numeric
  bounds, contribution limits, and budget allocation.
* Added an OpenDP adapter that fails closed when unavailable; privacy noise is
  neither user-seeded nor returned. A guarded public-fixture backend supports
  only data explicitly asserted to be public.
* Added fixed-dimensional bounded subject summaries, basic composed accounting,
  machine-readable release ledgers, privacy reports, and private-model leakage
  validation.
* Added new event-table generation with nominal/actual time, TAD, occasion,
  repeated doses, infusion start/stop pairs, multiple DVIDs, and schema/class
  restoration without source event-row copying.
* Fixed repeated-dose generation to honor the released per-subject observation
  total instead of repeating a complete endpoint grid after every dose.
* Added a two-part, privacy-accounted sampling model by dose occasion: sampling
  probability plus conditional observation count. `sampling_summary()` exposes
  the fitted design. The theophylline demo declares neither its regimen nor its
  occasion schedule and reports its inferred visit probabilities explicitly.
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
* Rewrote tests, README, scripts, and exactly two vignettes for Version 2,
  including `theo_md`, `warfarin`, `wbcSim`, a public censoring fixture, and
  empirical bug-finding tests that are explicitly not privacy proofs.
