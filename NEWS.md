# synpmx 0.0.0.9000

* Development version. There has not been a released version yet, so there is
  nothing here to upgrade from; changes are tracked in the git history until
  the first release, at which point this file starts recording user-visible
  changes by version.

* `synpmx_avatar()` gains `coarsen_time`, defaulting to `TRUE`. Source times are
  collapsed onto a shared visit grid before generation and per-visit deviations
  are pooled across the cohort and resampled onto each avatar. This closes
  `SIM-014` ("no generated vector may be identical to a source vector") on the
  AVATAR engine, where the gate had only ever been enforced against the
  structural/DP path even though AVATAR is the mode that copies an event
  skeleton verbatim. **Generated output changes** for any source with actual
  recorded times; a source already on its nominal grid is byte-identical. Set
  `coarsen_time = FALSE` to keep exact source timing.

* New `skeleton_uniqueness()` reports, per source subject, how many others share
  its observation time vector, its observation count, and its event signature —
  the "is this subject alone?" complement to `flag_identifiable_subjects()`,
  which finds subjects that are extreme. `scripts/measure_skeleton_uniqueness.R`
  runs the before-and-after over the public datasets.

* `synpmx_avatar()` now recomputes the dose from the avatar's own blended
  covariate when dosing is proportional to one (mg/kg, mg/m^2) within each
  assigned stratum. **Generated `AMT` changes** for such studies. This closes
  `REV-027`: the amount was previously copied verbatim from the anchor while
  covariates were blended, so it both disclosed the anchor's weight exactly and
  left every avatar violating its own protocol — a cohort dosed at exactly
  5 mg/kg produced avatars from 4.41 to 5.25. Several dose levels are found by
  clustering the observed ratios, so a 1/2/3 mg/kg escalation is recognised
  without declaring the arm, as is intra-patient escalation. Detection fails
  closed and is recorded as `dose_basis` / `dose_levels` in the settings.

* `synpmx_avatar()` gains `min_pattern_share`, default 2. Each avatar's set of
  attended visits is drawn from patterns at least that many source subjects
  share, so no one-of-a-kind attendance pattern is reproduced — the residual
  `coarsen_time` cannot reach, since it is which visits a subject attended rather
  than when they occurred. Nobody is dropped to achieve it: a subject with a rare
  pattern still contributes measurements as a donor. Dose events are never
  sampled. `1` restores copying the anchor's pattern. This partly closes
  `REV-026`, which coarsening made far cheaper by reducing a schedule to a bitmap
  over a shared grid.

  The draw is two-stage, because matching exact patterns alone discards nearly
  everything: two patients who each missed one visit count as different patterns
  if they missed different visits. A **shape** is drawn first — how many visits
  were missed and whether the misses were terminal, contiguous, or scattered —
  then a real pattern of that shape if one clears the floor, and only otherwise a
  generated arrangement, rejected and redrawn if it lands on a pattern too rare
  to reuse. On `warfarin` this takes the loss from 12 patterns to 2, keeping 5 of
  the source's 6 distinct sample counts. What is lost is resolution: how much
  missingness and of what kind survive, which specific visits does not.

  The default of 2 states exactly one thing: no synthetic patient carries a
  schedule unique to a real patient. Patterns below the floor are **discarded**,
  not approximated, so real dropout and dose-interruption patterns are lost. That
  cost is dataset-dependent and can be large — on `warfarin` the default excludes
  12 of 14 patterns held by 12 of 32 patients — so every run now reports it as a
  loud alert and as `patterns_total`, `patterns_dropped` and
  `subjects_with_dropped_pattern` in the settings. Note that at the default floor
  the last two are necessarily equal, since a pattern is discarded exactly when
  one patient holds it.

* New `compare_pmx_proximity()` measures whether synthetic subjects landed too
  close to real ones — the measurement for donor blending, the one masking
  mechanism that acts on the values rather than the structure and previously had
  none. It reports a nearest-neighbour adversarial accuracy against a null built
  by splitting the source cohort in half and running the identical statistic, so
  small-sample artefacts cancel. Raw distance to the closest record is
  deliberately not the headline: it has no natural scale and mostly tracks cohort
  size. A regression test hands it a verbatim copy and requires it to object.

* `subject_properties` is no longer rejected by `synpmx_avatar()`. It now names
  the assigned stratum — treatment arm, dose group, cohort — carried verbatim and
  used to group the dose basis and the attendance-pattern pools. Declaring a dose
  group is what lets a multi-level study be recognised as weight-based within
  each level. It is **not** a blending barrier; only route of administration is.

* `time_jitter` is documented as a realism control rather than a privacy one.
  Every jittered time is clamped inside its own Voronoi cell, so no value of
  `time_jitter` moves a visit more than half a gap from the source subject's
  visit and the source schedule stays recoverable.

* The last observation's jitter cell is now bounded above by half the final gap
  instead of being left open, so a large `time_jitter` can no longer stretch
  follow-up arbitrarily past what `screen = TRUE` permits.

* `tad`, when declared, is recomputed from each avatar's own generated dose times
  rather than carried over from the anchor, where it described a schedule the
  avatar no longer has.
