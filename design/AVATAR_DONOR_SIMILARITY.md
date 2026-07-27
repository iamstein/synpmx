# AVATAR donor similarity: from exact event-signature matching to continuous blending

A scoping document, not an implementation plan. It exists so the redesign below
can be discussed and refined before touching `synpmx_avatar()`'s core
mechanism — the primary, maintained generation path (`design/PROTOTYPE_SPEC.md`
section 0), and the one with the widest blast radius of anything in the
package. Nothing here is implemented. See `design/TODO.md` and `REV-024` in
`design/REVIEW_BACKLOG.md` for status.

## 1. The problem, in one sentence

`.select_donors()` (`R/synthesis.R:363-398`) requires an anchor and a candidate
donor to have an **exactly equal** `.event_signature()` string
(`R/profiles.R:1-52`) before any distance-based neighbor selection runs at all.
For daily-dosing regimens where patients differ in exactly when they held or
resumed dosing, or how long they were followed, that string is close to unique
per patient — so the "compatible group" `.select_donors()` searches within
collapses to the anchor alone for nearly everyone, and blending silently stops
doing anything.

## 2. Where this actually bites

### 2a. Phase 3, long daily dosing, toxicity-driven interruptions

Everyone is on the same regimen — say, one tablet daily until progression or
unacceptable toxicity. The drug turns out to be poorly tolerated, so most
patients accumulate dose holds: a few days off for a grade 2/3 adverse event,
then resume, sometimes at the same dose and sometimes reduced. No two patients
hold on the same days.

### 2b. Phase 1 dose-escalation, variable follow-up

Patients enter at different dose levels as the study escalates, are followed
for different total durations (early cohorts have had more time to progress or
come off study; late cohorts are still on treatment at the data cut), and daily
dosing again means the event table is one row per administered dose.

In both shapes, the event table is exactly the kind of "structurally rich,
individually distinguishing" record `synpmx_avatar()` is supposed to be good
at blending — and exactly the kind the current signature gate excludes from
blending entirely.

## 3. Why the signature diverges: a worked example

`.event_signature()`'s schedule component (`R/profiles.R:29-46`) is built from
the dose-start times:

```r
schedule_token <- if (length(start_time) <= 1L) {
  paste0("starts=", length(start_time))
} else {
  interval <- signif(diff(start_time), 2L)
  paste0("starts=", length(start_time), ":gaps=",
         paste(format(interval, trim = TRUE, scientific = FALSE), collapse = ","))
}
```

Take two patients on an identical once-daily regimen over 40 days. Patient A
never interrupts: `start_time` has 40 entries, all 24 h apart.

```
starts=40:gaps=24,24,24,...,24        (39 values, all 24)
```

Patient B is the closest possible match in every clinical sense — same
regimen, same dose, same 40-day course — but held the drug for 3 days around
day 12 for a grade 3 event, then resumed at the same dose:

```
starts=37:gaps=24,24,...,24,96,24,...,24     (36 values, one 96 instead of 24)
```

(37 doses over 40 days — 3 missed — so 12 doses before the hold, a 4-day gap
spanning the 3 skipped days, then 25 more doses after resuming.)

Even a single held day makes the `starts=` count differ, and `paste(...,
collapse = ",")` makes the whole `gaps=` string differ (it is not decomposed
back into "one long inter-dose interval among otherwise-regular ones" anywhere
in the comparison — `==` on the full string is all `.select_donors()` does).
Patient A and Patient B are **not compatible donors for each other**, despite
being about as similar as two oncology patients on the same protocol can be.
Scale this to a cohort where dose interruption is common, and most patients
end up in a signature group of size one — themselves.

## 4. This is a known tradeoff, documented for a different reason

The failure mode itself is not new. It was written up twice already, both
times as an argument about **privacy**, not blending quality:

- `design/PROTOTYPE_SPEC.md:1075-1086` (why Version 1 AVATAR was originally
  abandoned): *"signature grouping partitions before blending, so an unusual
  regimen lands in a near-singleton group and is blended with itself."*
- `vignettes/articles/feasibility.Rmd:73-79`: *"A patient with an unusual
  regimen lands in a small — possibly singleton — signature group, and then
  gets 'blended' with themselves."*

Version 4 (`design/PROTOTYPE_SPEC.md:995-1030`) restored AVATAR anyway,
reasoning that a resampling method is acceptable for **trusted-environment**
use the same way `synadam`'s column resampling is — but it explicitly recorded
"the asymmetry": *"a resampled subject trajectory is strongly identifying,
close to a fingerprint"* (`design/PROTOTYPE_SPEC.md:1020-1023`). Singleton
donor groups push output precisely toward that fingerprint end of the
spectrum, and they do it hardest exactly where an event table is richest and
most individually shaped — the daily-dosing-with-interruptions case above is
close to a worst case on both axes at once: the least blending happens
precisely where the source record is most identifying.

Separately from that governance concern, singleton groups also undercut
AVATAR's plain value proposition for workflow realism: if the "blend" is one
donor (the anchor), the output for that subject is the anchor's own trajectory
plus subject/residual noise, not a blend of several real patients. For the
richest, most clinically interesting datasets, AVATAR currently does the least
work.

## 5. Proposed direction

Two coupled changes, matching the shape you described: rank the *k* most
similar event tables continuously rather than gating on an exact match, and
handle within-donor time misalignment by borrowing from the nearest actual
event rather than only interpolating.

### 5a. Split the signature into a coarse hard gate and a soft distance feature

Today's signature conflates two different kinds of compatibility:

1. **Coarse regimen compatibility** — same route/dose sign, same `CMT`/`DVID`
   pattern on event rows, same endpoint set. This does *not* degenerate under
   interruptions: a hold doesn't change the drug, the route, or which
   endpoints are measured.
2. **Fine schedule shape** — exact dose count and exact inter-dose gap
   sequence. This is what degenerates, because it demands the two schedules
   agree almost to the day.

Proposal: keep (1) as a hard gate — it is cheap, does not degenerate under
interruptions, and still rules out blending a QD patient with a BID patient
just because a distance metric found them numerically close. Drop (2) from the
gate entirely, and instead turn it into a handful of **numeric features** added
to the same profile vector `.build_profiles()` already builds for PCA distance
(`R/profiles.R:111-224`) — for example:

- `n_doses` — count of administration events;
- `total_followup` — last recorded time minus first dose time;
- `dose_days_fraction` — `n_doses` divided by the number of nominally scheduled
  dosing days over that follow-up, i.e. the fraction of scheduled doses
  actually taken, which directly encodes interruption burden;
- `n_interruption_episodes` — count of gaps exceeding some multiple of the
  nominal interval;
- `max_gap` — the single longest hold;
- `median_gap` — the typical inter-dose interval when on drug.

These features are standardized and centered exactly like the existing
covariate/trajectory features (`R/profiles.R:181-189`) and folded into the same
PCA space, so `.neighbor_distances()` and `.randomized_weights()`
(`R/profiles.R:226-253`) need no change — only the feature list going into
`.build_profiles()` changes, plus dropping the schedule component out of
`.event_signature()`'s equality gate. Patient A and Patient B from section 3
would now be compatible (same coarse regimen) and close in distance (nearly
identical `dose_days_fraction`, `n_doses`, `median_gap`; only `max_gap` and
`n_interruption_episodes` differ a little) — exactly "the five patients with
the most similar event tables" you described, ranked continuously rather than
filtered by exact match.

An alternative worth naming and rejecting for now: drop the hard gate
entirely and let distance do all the work, including route/dose-sign/endpoint
compatibility. This is simpler to implement but risks blending across
regimens that should never mix (a QD arm with a BID arm, an IV cohort with an
oral cohort) whenever noise in the other features happens to put them close.
Keeping a coarse hard gate and only softening the schedule-shape part is the
more conservative starting point.

### 5b. Borrowing at a target time when donor schedules don't line up

Today, `.interpolate_trajectory()` (`R/synthesis.R:218-238`) linearly
interpolates each donor's trajectory onto the anchor's own observation times,
and falls back to a proportional-position remapping when a target time falls
outside a donor's observed window (`vignettes/articles/avatar-mathematics.Rmd`
Step 3, "Different observation times across subjects").

The risk once donors no longer share the anchor's exact schedule: if a
donor's own interruption happened on a different day than the anchor's, a
target time that sits inside a real drug-holiday/rebound excursion for the
donor can get interpolated *across* that excursion from two more-distant real
points, quietly inventing a smoothed value the donor never actually measured
nearby. Three options, in increasing order of how much they change today's
mechanism:

1. **Nearest-event-time borrowing.** Instead of interpolating, take the
   donor's single closest recorded observation to the target time. Simplest,
   never invents an intermediate value, but is noisier and can flip
   discontinuously between two donor points as the target time crosses their
   midpoint.
2. **Bounded interpolation with a nearest-time fallback.** Interpolate only
   when the two bracketing donor observations are within some tolerance of the
   target time; beyond that tolerance, fall back to nearest-time borrowing
   rather than interpolating across a gap the donor's own record doesn't
   support. This most directly matches "borrow from the nearest event time
   when the schedules don't line up, interpolate when they roughly do" —
   recommended as the first thing to try, since it changes today's behavior
   only when the existing interpolation would otherwise be least trustworthy.
3. **Align on exposure rather than calendar time** — e.g. cumulative doses
   received, or a nominal protocol day that accounts for holds — instead of
   "time since first dose" as the interpolation axis, so two patients who are
   both "on their 40th administered dose" line up even if one took 55 calendar
   days to get there and the other took 40. This is the most faithful answer
   to "the interruption is what's messing up my time axis," but it is a
   materially bigger change: it touches `.aligned_time()` (`R/utils.R:244-263`),
   which other parts of the package also rely on, and probably deserves to be
   scoped separately once (5a) and a choice between 1/2 above are settled.

## 6. Open engineering questions

- Exact feature list for schedule shape, and whether it generalizes across
  trial designs (dose escalation, multiple arms, non-daily regimens) without
  per-study hand-tuning.
- Precisely where the coarse/soft split falls: is dropping only the schedule
  token enough, or should the endpoint-set token also move to a soft feature
  for studies that add endpoints partway through follow-up?
- Relative weighting between the new schedule-shape features and the existing
  covariate/trajectory features in the combined PCA distance — same footing as
  everything else, or a tunable knob (analogous to `pca_variance`)?
- Whether `k` (default 5) still means the same thing once there is no hard
  partition first — every subject in the coarse-compatible pool is now a
  candidate, so the neighbor pool can be much larger than today, and a warning
  when the 5th-nearest neighbor is still far away may be worth adding.
- Which of the three time-borrowing options (5b) to implement first, and
  whether it is a new opt-in argument to `synpmx_avatar()` or a change to
  today's default.
- Test/fixture impact: none of the current fixtures or the five `nlmixr2data`
  demo datasets exercise variable dose-interruption schedules, so validating
  this change needs new synthetic fixtures built specifically to have that
  shape — this is a real regression-coverage gap today, independent of this
  proposal.
- Documentation impact once a design is chosen and implemented:
  `vignettes/articles/avatar-mathematics.Rmd` steps on time alignment,
  donor compatibility, and trajectory synthesis all need rewriting to describe
  the new mechanism precisely — see `AGENTS.md`'s documentation-sync rule,
  which now names this vignette explicitly for exactly this kind of change.
