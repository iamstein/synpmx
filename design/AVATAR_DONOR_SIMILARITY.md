# AVATAR donor similarity: from exact event-signature matching to continuous blending

> **Status (2026-07-27): the direction is decided, nothing is implemented
> yet.** Andy decided (by voice, away from a machine) to drop
> `.select_donors()`'s hard event-signature gate **completely** — see section
> 5a. What is still open is *how* to implement that without regressing the
> thing the gate protected for free (never blending across a genuinely
> incompatible regimen); section 6 lists exactly what needs resolving before
> writing code. Pick this up by reading section 5a's decision note, then
> section 6.

This is a scoping document, not (yet) an implementation plan, for changes to
`synpmx_avatar()`'s core mechanism — the primary, maintained generation path
(`design/PROTOTYPE_SPEC.md` section 0), and the one with the widest blast
radius of anything in the package. See `design/TODO.md` and `REV-024` in
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

### 5a. Decision: drop the hard gate entirely

**Decided 2026-07-27 (Andy, by voice while away from the machine this needs to
be implemented on): drop `.select_donors()`'s hard signature-equality gate
completely.** Not narrowed, not split into a coarse/fine version — removed. No
part of event-table compatibility remains a boolean filter; everything about
an event table that currently distinguishes one subject from another becomes
a **numeric feature** feeding the same continuous distance that already ranks
covariate and trajectory similarity, and the `k` nearest subjects by that
combined distance are the donors, full stop.

This was raised and set aside as an alternative in the first draft of this
document (keep a coarse route/dose-sign/endpoint-set gate, only soften the
fine schedule-shape matching) — that alternative is now explicitly rejected in
favor of the simpler, fully continuous design. The risk that motivated the
coarse gate — blending a QD regimen with a BID regimen, or an IV cohort with
an oral cohort, because noise in unrelated features happened to put them
close — still needs an answer, but the answer is now "make that risk small
through the feature design and distance weighting," not "rule it out with a
filter." Concretely: route, dose sign, `CMT`/`DVID` pattern, and endpoint set
all become features in the same profile vector rather than disappearing —
e.g. one-hot/indicator features analogous to how categorical covariates are
already encoded (`R/profiles.R:132-146`) — so a QD/BID mismatch still shows up
as a large distance, just one computed continuously rather than enforced by
`==`. Getting the relative scaling right so that these become naturally
dominant contributors to distance when they differ (rather than getting
diluted by many small schedule-shape/covariate features) is now the central
open question for this section — see the list at the end of this document.

Alongside the coarse-compatibility signals, fold in the schedule-shape
features from the original proposal — a handful of **numeric features** added
to the same profile vector `.build_profiles()` already builds for PCA distance
(`R/profiles.R:111-224`):

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
`.build_profiles()` changes, and `.event_signature()`/the signature-equality
check in `.select_donors()` (`R/synthesis.R:363-398`) is deleted rather than
narrowed. Patient A and Patient B from section 3 would now both be candidates
and close in distance (nearly identical `dose_days_fraction`, `n_doses`,
`median_gap`; only `max_gap` and `n_interruption_episodes` differ a little) —
exactly "the five patients with the most similar event tables" you described,
ranked continuously with no filtering step at all.

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

Section 5a's core decision is made (no hard gate). What remains is how to
implement it without silently regressing the thing the old gate protected for
free — never blending across a genuinely incompatible regimen.

- **Feature encoding and weighting for what used to be the hard gate.** Route,
  dose sign, `CMT`/`DVID` pattern, and endpoint set need to become features
  that dominate the combined distance whenever they differ, not features that
  get averaged in alongside dozens of others and diluted. This may need an
  explicit weighting mechanism beyond today's uniform standardization
  (`R/profiles.R:181-189`) — e.g. a large fixed multiplier on structural-
  mismatch features, or computing distance in two blocks (structural,
  then everything else) and only using the second block to rank within ties
  on the first. Getting this wrong is the main way "drop the gate" could go
  wrong in practice, so it deserves the most scrutiny before writing code.
- Exact feature list for schedule shape, and whether it generalizes across
  trial designs (dose escalation, multiple arms, non-daily regimens) without
  per-study hand-tuning.
- Relative weighting between the new schedule-shape features and the existing
  covariate/trajectory features in the combined PCA distance — same footing as
  everything else, or a tunable knob (analogous to `pca_variance`)?
- Whether `k` (default 5) still means the same thing once there is no
  partition at all — every subject in the dataset is now a candidate for every
  anchor, so the neighbor pool can be much larger than today, and a warning
  when the 5th-nearest neighbor is still far away (in the structural-mismatch
  sense above, not just Euclidean distance overall) may be worth adding.
- Which of the three time-borrowing options (5b) to implement first, and
  whether it is a new opt-in argument to `synpmx_avatar()` or a change to
  today's default.
- Test/fixture impact: none of the current fixtures or the five `nlmixr2data`
  demo datasets exercise variable dose-interruption schedules, so validating
  this change needs new synthetic fixtures built specifically to have that
  shape. It also needs a fixture that actively tests the removed-gate risk —
  e.g. a mixed QD/BID dataset — to confirm the new distance still keeps them
  apart in practice, not just in theory. This is a real regression-coverage
  gap today, independent of this proposal.
- Documentation impact once a design is chosen and implemented:
  `vignettes/articles/avatar-mathematics.Rmd` steps on time alignment,
  donor compatibility, and trajectory synthesis all need rewriting to describe
  the new mechanism precisely — see `AGENTS.md`'s documentation-sync rule,
  which now names this vignette explicitly for exactly this kind of change.
