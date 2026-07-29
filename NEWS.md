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
