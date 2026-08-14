# Changelog

## synpmx 0.0.0.9000

- Development version. There has not been a released version yet, so
  there is nothing here to upgrade from; changes are tracked in the git
  history until the first release, at which point this file starts
  recording user-visible changes by version.

- **[`pmx_masking_report()`](https://iamstein.github.io/synpmx/reference/pmx_masking_report.md)
  takes a `section` argument.** The full report is thirty-odd rows,
  which is the right size to read once after a run and the wrong size
  under a paragraph making one point: `section = "dose_schedules"`
  prints the five dose rows instead. Sections are `anchors`, `donors`,
  `blend`, `visits`, `visit_sets`, `dose_schedules`, `dose_amounts`. The
  dose rows were previously the tail of the visit-set block and are now
  their own section, and the old `Dose` header is
  `Dose amounts: HOW MUCH each patient received`. Every masking table in
  `public-data-examples.Rmd` now prints one section.

- **Captions no longer say RESTRICTED.** The word was printed on the
  output of every comparison, census, screen, and scorecard, including
  runs over public package data where nothing is restricted. The
  `"release_status"` attribute still records
  `"restricted_not_releasable"`, which is where that judgement belongs.

- **Scorecard rows B1a and B1b read “Avatars with a visit set / dose
  schedule nobody else shares”**, not “wearing”.

- **Donor trajectories are no longer extrapolated past the times a donor
  was actually observed at** (`SIM-046`). The old fallback rescaled the
  anchor’s time range onto the donor’s, which reads as “the same shape
  on a stretched clock” and is wrong whenever the missing region is a
  distinct kinetic phase. On `warfarin`, where 19 of 32 subjects have no
  `cp` sample before 24 h, a donor’s 24 h elimination concentration
  stood in for an anchor’s 0.5 h absorption concentration: the generated
  median at 0.5 h was 8.8 against a source median of 0.0, and the
  absorption limb came back as a flat plateau. Rows outside a donor’s
  range are now blended over the donors that do cover them; rows no
  donor covers fall to a per-time cohort median rather than a single
  dataset-wide one.

- **A generated visit set can no longer reach deeper into the study than
  `min_pattern_share` real patients reached** (`SIM-047`). A `scattered`
  attendance shape was placed over the whole union grid, so one subject
  followed far past the rest lent their tail to every avatar: on
  `wbcSim` an avatar came out with two observations, at 501 h and
  3910 h. Row B1a could not see it – it asks whether a real patient’s
  visit set was *copied*, and a fabricated schedule is exactly what
  passes.

- **Scorecard row A6 now prints in order.** It was appended after the
  main row list, so wherever it fired it appeared after B4b instead of
  after A5.

- **`vignettes/avatar-evaluation-public-data.Rmd` is renamed
  `public-data-examples.Rmd`**, its figures put source and synthetic
  side by side on a shared y axis with one row per endpoint (a PK
  concentration in single digits no longer flattens against a PD score
  in the hundreds), and its model-based differential-privacy walkthrough
  is cut – that vignette is about
  [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md),
  and the DP modes are documented elsewhere. The old pkgdown URL
  redirects.

- **`pmx_scorecard()` is now
  [`synpmx_scorecard()`](https://iamstein.github.io/synpmx/reference/synpmx_scorecard.md)**,
  and every row carries an `explore` column naming the call that
  explains it. Three rows changed with it: C4 is titled “distinct
  dose-time schedules represented”, since the key it counts is dose
  *times* and never amounts; B5 names the column and level it found
  rather than reporting a bare count; and C3 marks a stratum that
  changed size `review` rather than `FAIL`, because dropping a subject
  for want of donors is the generator working as designed. The `reads`
  value `run report` is now `run settings`.
  `vignettes/synthetic-data-checks.Rmd` is renamed
  `scorecard-synthetic-data-checks.Rmd`.

- **[`compare_pmx_rare_levels()`](https://iamstein.github.io/synpmx/reference/compare_pmx_rare_levels.md)
  censuses the categorical levels too few *source* patients held**, and
  the scorecard reports it as row **B5b**. Categorical covariates are
  not blended – they are sampled from the donors’ values, and `strata`
  are copied from the anchor by design – so a synthetic patient’s
  category is always some real patient’s actual category, copied. What
  matters is therefore how many *real* patients held it, not how many
  synthetic ones do: a level two real patients held, appearing in a
  released table, says that someone with that attribute was in this
  study. The floor is `min_pattern_share`, the same rule that already
  protects visit sets. The existing synthetic-side check is now **B5a**;
  both are `review`, and B5b reads the source, so on a study where B5a
  was the only categorical row the card becomes restricted output.

- **The printed scorecard is narrower.** The `explore` calls moved out
  of the table and under it, listed only for the rows that did not pass,
  because a sixth column of calls pushed every line past 150 characters
  and wrapped. The B5b levels print one per line beneath that, and
  knitting the card emits them as a second `kable()` table, so a report
  says which levels without the reader running the census again. The
  returned data frame is unchanged: it still carries `explore` on every
  row.

- **Only seven scorecard rows can say `FAIL`**: A1, A3, A6, B1a, B1b,
  B4a and B4b — the output is not a legal dataset, it is not the study
  that went in, or it reproduces one real patient’s structure verbatim.
  A2, A4, B3 and B5 now say `review` instead, joining A5, B2, C3 and C4.
  Each of them can move for a legitimate reason: a subject dropped for
  want of donors (A4), a proximity statistic wandering at a small sample
  size and meaning opposite things on either side of its interval (B3),
  the weak synthetic-side form of the rare- level check (B5), or a real
  source a validator objects to (A2).

- **Discrete endpoints stay discrete, and this changes generated output
  for any dataset with one.** Blending is a weighted mean, so a weighted
  mean of several patients’ zeros and ones is a number between them: on
  [`xgxr::mad`](https://rdrr.io/pkg/xgxr/man/mad.html) a 0/1 endpoint
  came back as 600 distinct values spanning -0.13 to 1.08, an ordinal
  1/2/3 endpoint reached 4.69, and an integer count came back in
  fractions. Nothing caught it, because class restoration works on the
  DV column and a PMX table keeps every endpoint in the same one.
  [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
  now decides per endpoint, from the source, whether it is `binary`,
  `ordinal`, `integer` or `continuous`, and snaps generated values back
  onto that scale.
  [`pmx_endpoint_types()`](https://iamstein.github.io/synpmx/reference/pmx_endpoint_types.md)
  reports the decision and the evidence for it,
  `pmx_roles(endpoint_types = )` overrides it, and scorecard row **A6**
  checks the finished table. One study changes that nobody would call
  discrete:
  [`nlmixr2data::warfarin`](https://nlmixr2.github.io/nlmixr2data/reference/warfarin.html)’s
  `pca` is a percentage recorded without decimals, so its generated
  values are now whole numbers too.

- **`subject_properties` is now `strata`** in
  [`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md),
  and `subject_property_summary()` is now
  [`strata_summary()`](https://iamstein.github.io/synpmx/reference/strata_summary.md).
  The old name described where the columns live; the new one says what
  they do — group the donors, condition the regimen, and now carry the
  cohort balance — and separates them from `covariates` (blended) and
  `keep` (copied and inert). Pre-release, so there is no deprecation
  shim: the old argument name is an error.

- **`synpmx_avatar(preserve_strata_balance = TRUE)`, new and on by
  default, changes generated output for any dataset declaring
  `strata`.** An avatar never left its anchor’s stratum, but anchors are
  sampled with replacement, so the *balance* was left to the draw:
  [`xgxr::case1_pkpd`](https://rdrr.io/pkg/xgxr/man/case1_pkpd.html)’s
  exactly-30-per-arm design came back between 21 and 39 per arm
  depending only on the seed, and any downstream summary grouped by arm
  inherited that. Each stratum now gets its source share, as a
  proportion, so it holds when `n_subjects` differs from the source
  size. Strata holding fewer than three source patients are deliberately
  left stochastic — reproducing such a stratum’s size exactly would
  disclose it, since a joint cell with one real patient would get
  exactly one avatar on every seed — and `strata_balanced` /
  `strata_stochastic` report how many fell on each side.

- **[`flag_identifiable_subjects()`](https://iamstein.github.io/synpmx/reference/flag_identifiable_subjects.md)
  scores within stratum.** It was scoring every axis against the whole
  cohort, which reports the protocol back as a privacy finding on any
  study that assigns a dose: on
  [`xgxr::case1_pkpd`](https://rdrr.io/pkg/xgxr/man/case1_pkpd.html), a
  six-arm design from 3 mg to 300 mg, the top arm sits about 6.4
  modified-z units from the cohort median and **59 of 180 avatars were
  flagged, 31 of them for receiving the dose their arm was assigned**.
  Scored within arm it flags 1, and a patient given twice their own
  arm’s dose is still flagged. Strata under five subjects fall back to
  cohort-wide scoring, since a scale estimated from four patients
  describes the four rather than the one being screened. With no
  `strata` declared nothing changes.

- New vignette, “Checks of the synthetic data”
  (`vignettes/scorecard-synthetic-data-checks.Rmd`): the six categories
  of check to run on generated data, worked on
  [`xgxr::case1_pkpd`](https://rdrr.io/pkg/xgxr/man/case1_pkpd.html)
  with
  [`nlmixr2data::pheno_sd`](https://nlmixr2.github.io/nlmixr2data/reference/pheno_sd.html)
  as the case where they fail. Written for two audiences — someone
  deciding whether to use generated data, and someone building their own
  generator who wants the category list rather than these function
  names. It ends with its own gaps rather than implying the list is
  complete. The two items above were both found by writing it.

- The demo vignette is now an evaluation, and is named like one:
  `synpmx-demo.Rmd` becomes `avatar-evaluation-public-data.Rmd`
  (“Evaluating AVATAR on public data”). It is a measurement of
  [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
  across public data rather than a tutorial, so it now runs all eight
  datasets rather than five —
  [`xgxr::case1_pkpd`](https://rdrr.io/pkg/xgxr/man/case1_pkpd.html),
  [`xgxr::mad`](https://rdrr.io/pkg/xgxr/man/mad.html) and
  [`nlmixr2data::pheno_sd`](https://nlmixr2.github.io/nlmixr2data/reference/pheno_sd.html)
  join the original five — and opens with a table describing every
  dataset and what each is there to exercise. The old pkgdown URL
  redirects. Per-dataset numbers quoted in the prose were re-measured
  against current behavior; several had drifted, and the
  masking-mechanism table had been missing `M6`.

- Three more public datasets are covered by the test suite.
  [`xgxr::case1_pkpd`](https://rdrr.io/pkg/xgxr/man/case1_pkpd.html)
  (180 patients, declared nominal time, six arms, two endpoints,
  baseline weight) is the first public dataset shaped like a real study
  report — none of the nlmixr2data five has a declared nominal time, a
  treatment arm to stratify on, or a censoring column, so every feature
  a study report exercises was previously untested on public data.
  [`xgxr::mad`](https://rdrr.io/pkg/xgxr/man/mad.html) has five
  observation endpoints, where two was the previous maximum.
  [`nlmixr2data::pheno_sd`](https://nlmixr2.github.io/nlmixr2data/reference/pheno_sd.html)
  is 59 real patients whose individualised neonatal dosing genuinely
  cannot be masked, and is the registry’s honest example of that.

- **Dosing that stopped early is now represented rather than dropped
  (`SIM-044`).** Protecting a patient who stopped dosing at a point
  nobody else did meant not building on them, which removed the regimen
  from the output. Dose *times* are still never moved or invented – that
  would emit a regimen no protocol permits – but truncating a schedule
  at one of its own dose times yields a regimen the study did give
  someone, so an avatar now stops at a depth several patients used, or
  one nobody used. On a 32-patient fixture with four patients stopping
  at depths 4, 6, 7 and 9, the output carries depths 5, 8 and 10 instead
  of pushing every avatar onto the full schedule.

- [`pmx_masking_report()`](https://iamstein.github.io/synpmx/reference/pmx_masking_report.md)
  reports how many of the source’s distinct dose regimens are
  represented in the synthetic cohort. Declining to build on a patient
  whose dose schedule nobody shares is the only safe answer, but it
  removes that regimen from the output entirely, and that was silent:
  nineteen patients on three doses, one on two and one on one came out
  as twenty-one on three, with nothing saying so.

- **The dose side of the guarantee (`SIM-043`).**
  [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
  now also guarantees that no avatar carries a set of dose times only
  one real patient has. Dose events are not rewritten – resampling them
  would emit regimens the protocol never permitted – so a patient whose
  dose schedule nobody shares is not built upon, and the avatar is
  re-anchored instead. Where every patient is in that position the run
  says so rather than silently doing nothing.
  `identifying_dose_schedules` records it: 0 for `theo_md`, `warfarin`,
  `wbcSim` and `mavoglurant`, and an honest 12 of 12 for `nimoData`,
  whose dosing is individualised.

- **Bug fix (`SIM-042`), changes generated output.** An avatar could
  still be emitted carrying a visit set exactly one real patient holds.
  A schedule group of one – a single patient measuring a different set
  of endpoints from everybody else – has no shared visit set to draw
  from and nothing to substitute, and the loop then copied the anchor’s
  own. Now: a group with no pool of its own borrows from any group
  measuring the same endpoints; the visit set is decided before the
  avatar is built, so an anchor that cannot be masked causes that
  **avatar** to be re-anchored rather than the **patient** to be
  dropped; and a final check on the finished table records
  `identifying_visit_sets`, which must be 0. It is 0 across every public
  dataset the demo uses.

- [`pmx_masking_report()`](https://iamstein.github.io/synpmx/reference/pmx_masking_report.md)
  shows a count and a share on every row that counts patients or
  avatars. Neither reads on its own: “5%” of 21 patients is one patient,
  and “15” means nothing without the cohort size beside it.

- **Bug fix (`SIM-041`).**
  [`flag_identifiable_subjects()`](https://iamstein.github.io/synpmx/reference/flag_identifiable_subjects.md)
  flagged nearly every patient who stopped early. Its robust scale is
  the median absolute deviation, which is zero whenever more than half a
  cohort shares one exact value – the ordinary case on trial data, since
  most patients complete the protocol and stop at the same visit – and
  the zero branch then scored every other value as infinitely extreme.
  It now falls back to the mean absolute deviation from the median,
  which is zero only when nothing varies at all. On a clustered
  21-patient cohort this goes from 10 flagged to 0 while a genuine
  extreme still scores far past the threshold; `mavoglurant` goes from
  63 of 120 to 17.

- “Avatars keeping their anchor’s own visit set” was reported as a
  number to drive to zero, and it is not one. Copying an anchor’s visit
  set discloses nothing when several real patients share that set; it is
  a problem only when the set is unique to one of them. The report now
  separates the two and
  [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
  **fixes the identifying case by default**: where the anchor’s own set
  is held by nobody else and no arrangement is free, the group’s most
  widely held set is substituted instead. A run alerts only when the
  group has nothing shareable to substitute, which is the one case the
  generator cannot resolve on its own.

- Generated patterns no longer invent missing visits in an endpoint that
  has none. The shape was placed over the pooled endpoint-by-time grid
  in time order, so a biomarker drawn at all ten visits for all
  twenty-one patients came out of generation with seven to ten of them.
  Placements are now confined to endpoints whose visit times actually
  vary across the cohort.

- “Dropout” is no longer used to explain why visit sets differ. Visits
  can be missing because a patient discontinued, because a visit was
  missed, or because follow-up has not reached them yet, and the
  reports, alerts and vignettes said dropout throughout.

- [`plot_pmx_schedule()`](https://iamstein.github.io/synpmx/reference/plot_pmx_schedule.md)
  no longer overloads red. The endpoint palette was red-free-adjacent at
  best – the second endpoint was orange, which at screen distance is the
  same colour as the red used for a unique schedule and for a singleton
  visit time, so a red dot could mean an endpoint or a warning. Endpoint
  colours are now red-free, red means only “this identifies somebody”,
  unique-schedule patients are a shaded band rather than coloured label
  text, and endpoints are offset within each row so two measured at the
  same visits no longer hide each other.

- [`skeleton_uniqueness()`](https://iamstein.github.io/synpmx/reference/skeleton_uniqueness.md)
  gains `nearest_set_diff` and `n_visits`, and reports which endpoint
  drives the count. `n_share_schedule == 1` is exact-set equality, which
  on a real study is harsh: two patients differing by one missed sample
  score as unique exactly like two with nothing in common. A cohort
  reading 15 of 21 unique can have every one of those 15 a single
  missing sample from somebody else, and the count alone cannot say so.
  The printed output also states that the count is a property of the
  *source* and that what generation controls is the run report’s
  “avatars keeping their anchor’s own visit set”.

- **Bug fix (`SIM-039`), changes generated output.** Avatars kept their
  anchor’s own visit set far too often – 86% on a real 21-patient study
  – which copies one real patient’s exact pattern of absences and is
  what `min_pattern_share` exists to prevent. Two causes, both triggered
  by ordinary dropout: a `trailing` placement is deterministic, so the
  24 retries all re-proposed the one arrangement that had just been
  rejected; and under staggered discontinuation every (kind, count)
  shape has a single holder, so the group got no pool at all. Placements
  are now enumerated, the miss count walks outward when nothing at the
  wanted count is free, and a third abstraction – the kind of
  missingness alone – is reached when the finer one clears nothing.
  [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
  alerts when 10% or more of avatars still fall back, and
  `pattern_shifted_fraction` records how often the miss count moved.

- **[`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md)
  gains `dose_covariate`** (`SIM-040`). Name the covariate the dose is a
  fixed multiple of – weight, body surface area – and
  [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
  recomputes each avatar’s `amt` from the avatar’s own blended value
  instead of copying its anchor’s milligrams. This skips the
  conservative inference, which fails closed on studies that escalate
  within patient or dispense in vials, and holds each dose row’s own
  ratio, so intra-patient escalation is preserved exactly. The column
  must also be named in `covariates`. `theo_md` is the public
  demonstration: it is dosed by weight, its recorded mg/kg runs 3.1 to
  5.9, inference declines, and declaring `WT` takes avatars on a real
  dosing course from 0% to 100%
  (`tests/testthat/test-integration-nlmixr2data.R`).

- Alerts no longer print twice in a knitted report. The signalled
  condition carried the `warning` class, which `knitr` renders in
  addition to the message; it is now a plain `synpmx_alert` condition.
  Handle it with `withCallingHandlers(synpmx_alert = ...)`.

- Alerts that used to end “declare a `nominal_time` role” now say
  something useful to a caller who already declared one, naming the pool
  split from `strata` and the `min_pattern_share` floor instead.

- **Bug fix (`SIM-038`), changes generated output.** Times were keyed
  with `format(x, digits = 12)`, which fixes one layout for a whole
  vector, so the same visit keyed differently for a patient who also had
  a fractional sample. Two patients who attended exactly the same visits
  therefore produced different attendance keys: each pattern had a
  single holder, `min_pattern_share` discarded it, and the avatar fell
  back to its anchor’s own visit set — the outcome the mechanism exists
  to prevent. The same fault made
  [`skeleton_uniqueness()`](https://iamstein.github.io/synpmx/reference/skeleton_uniqueness.md)
  report one-off observation times that eighteen other patients in fact
  shared, firing the “unique observation times” alert and its
  `nominal_time` advice on cohorts already on a shared grid.

- Alerts from
  [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
  are wrapped to a fixed width with a headline, the count, why it
  matters, and the one thing to do about it, instead of a single
  unwrapped paragraph that needed horizontal scrolling in a knitted
  report. Each is emitted **once**: the condition is signalled rather
  than raised with [`warning()`](https://rdrr.io/r/base/warning.html),
  which used to print the identical text a second time.
  `tryCatch(warning = )` and
  [`suppressWarnings()`](https://rdrr.io/r/base/warning.html) behave as
  before; `options(warn = 2)` no longer promotes them, so rely on the
  printed banner. Alerts are also now split into `SYNPMX ALERT` (a
  property of the source that raises risk and that you can act on) and
  `SYNPMX NOTE` (a mechanism reporting what it cost).

- [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
  alerts when 10% or more of avatars were given no visit set from the
  pool and kept their anchor’s own, which is one real patient’s exact
  pattern of absences copied onto an avatar. Previously this could
  happen silently whenever every legal placement of a shape was already
  somebody’s, and only `pattern_sampled_fraction` recorded it.

- New
  [`pmx_masking_report()`](https://iamstein.github.io/synpmx/reference/pmx_masking_report.md)
  turns the `pmx_settings` attribute into the table to read after a run,
  with a sentence beside every number saying what it means. It states
  outright whether dose amounts were recomputed from a covariate and,
  when they were not, which covariates were tried and what failed —
  previously indistinguishable from detection never having run. The demo
  vignette and the study templates under `scripts_private/` all use it,
  so the labels are maintained in one place.

- [`compare_pmx_distributions()`](https://iamstein.github.io/synpmx/reference/compare_pmx_distributions.md)
  gains a `knit_print()` method, so its tables come out as tables when
  knitted. Each of the four study templates under `scripts_private/` had
  carried its own copy of that display code.

- New
  [`plot_pmx_schedule()`](https://iamstein.github.io/synpmx/reference/plot_pmx_schedule.md)
  draws a cohort’s dosing and observation schedule: one row per patient,
  one mark per event, with a per-visit patient count underneath and
  unique schedules marked in red. A uniqueness *count* cannot
  distinguish twelve one-off sampling times from twelve ordinary
  dropouts; the picture can.

- [`skeleton_uniqueness()`](https://iamstein.github.io/synpmx/reference/skeleton_uniqueness.md)
  gains `coarsen_time` (default `FALSE`) so the same cohort can be
  scored before and after the grid
  [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
  builds, and its printed output says which of the two you are looking
  at. Its columns are renamed to one consistent question — how many
  patients share this property, this patient included, so `1` means
  “nobody else”: `n_share_schedule`, `n_share_rarest_time`,
  `n_share_obs_count`, `n_share_dosing`, plus a `why_unique` column.
  [`print()`](https://rdrr.io/r/base/print.html) now leads with a
  verdict and two summary tables rather than twelve patient rows, and
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html) output is
  used automatically when knitting.

- [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
  gains `coarsen_time`, defaulting to `TRUE`. Source times are collapsed
  onto a shared visit grid before generation and per-visit deviations
  are pooled across the cohort and resampled onto each avatar. This
  closes `SIM-014` (“no generated vector may be identical to a source
  vector”) on the AVATAR engine, where the gate had only ever been
  enforced against the structural/DP path even though AVATAR is the mode
  that copies an event skeleton verbatim. **Generated output changes**
  for any source with actual recorded times; a source already on its
  nominal grid is byte-identical. Set `coarsen_time = FALSE` to keep
  exact source timing.

- New
  [`skeleton_uniqueness()`](https://iamstein.github.io/synpmx/reference/skeleton_uniqueness.md)
  reports, per source subject, how many others share its observation
  time vector, its observation count, and its event signature — the “is
  this subject alone?” complement to
  [`flag_identifiable_subjects()`](https://iamstein.github.io/synpmx/reference/flag_identifiable_subjects.md),
  which finds subjects that are extreme.
  `scripts/measure_skeleton_uniqueness.R` runs the before-and-after over
  the public datasets.

- [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
  now recomputes the dose from the avatar’s own blended covariate when
  dosing is proportional to one (mg/kg, mg/m^2) within each assigned
  stratum. **Generated `AMT` changes** for such studies. This closes
  `REV-027`: the amount was previously copied verbatim from the anchor
  while covariates were blended, so it both disclosed the anchor’s
  weight exactly and left every avatar violating its own protocol — a
  cohort dosed at exactly 5 mg/kg produced avatars from 4.41 to 5.25.
  Several dose levels are found by clustering the observed ratios, so a
  1/2/3 mg/kg escalation is recognised without declaring the arm, as is
  intra-patient escalation. Detection fails closed and is recorded as
  `dose_basis` / `dose_levels` in the settings.

- [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
  gains `min_pattern_share`, default 2. Each avatar’s set of attended
  visits is drawn from patterns at least that many source subjects
  share, so no one-of-a-kind attendance pattern is reproduced — the
  residual `coarsen_time` cannot reach, since it is which visits a
  subject attended rather than when they occurred. Nobody is dropped to
  achieve it: a subject with a rare pattern still contributes
  measurements as a donor. Dose events are never sampled. `1` restores
  copying the anchor’s pattern. This partly closes `REV-026`, which
  coarsening made far cheaper by reducing a schedule to a bitmap over a
  shared grid.

  The draw is two-stage, because matching exact patterns alone discards
  nearly everything: two patients who each missed one visit count as
  different patterns if they missed different visits. A **shape** is
  drawn first — how many visits were missed and whether the misses were
  terminal, contiguous, or scattered — then a real pattern of that shape
  if one clears the floor, and only otherwise a generated arrangement,
  rejected and redrawn if it lands on a pattern too rare to reuse. On
  `warfarin` this takes the loss from 12 patterns to 2, keeping 5 of the
  source’s 6 distinct sample counts. What is lost is resolution: how
  much missingness and of what kind survive, which specific visits does
  not.

  The default of 2 states exactly one thing: no synthetic patient
  carries a schedule unique to a real patient. Patterns below the floor
  are **discarded**, not approximated, so real dropout and
  dose-interruption patterns are lost. That cost is dataset-dependent
  and can be large — on `warfarin` the default excludes 12 of 14
  patterns held by 12 of 32 patients — so every run now reports it as a
  loud alert and as `patterns_total`, `patterns_dropped` and
  `subjects_with_dropped_pattern` in the settings. Note that at the
  default floor the last two are necessarily equal, since a pattern is
  discarded exactly when one patient holds it.

- New
  [`compare_pmx_proximity()`](https://iamstein.github.io/synpmx/reference/compare_pmx_proximity.md)
  measures whether synthetic subjects landed too close to real ones —
  the measurement for donor blending, the one masking mechanism that
  acts on the values rather than the structure and previously had none.
  It reports a nearest-neighbour adversarial accuracy against a null
  built by splitting the source cohort in half and running the identical
  statistic, so small-sample artefacts cancel. Raw distance to the
  closest record is deliberately not the headline: it has no natural
  scale and mostly tracks cohort size. A regression test hands it a
  verbatim copy and requires it to object.

- `strata` is no longer rejected by
  [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md).
  It now names the assigned stratum — treatment arm, dose group, cohort
  — carried verbatim and used to group the dose basis and the
  attendance-pattern pools. Declaring a dose group is what lets a
  multi-level study be recognised as weight-based within each level. It
  is **not** a blending barrier; only route of administration is.

- [`validate_pmx()`](https://iamstein.github.io/synpmx/reference/validate_pmx.md)
  no longer refuses a `nominal_time` column with gaps. Missing nominal
  times are ordinary — an unscheduled visit has no protocol slot — and
  [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
  already handles them row by row, snapping what has one and falling
  back to the inferred grid for the rest (reported as the `"mixed"`
  grid). A wholly missing column is still an error, since the role
  should just be left undeclared.

- `strata` may now be missing for some subjects. They are grouped as
  their own stratum and the check reports a warning instead of an error.
  A column that *varies* within a subject is still an error: it cannot
  be that subject’s assignment. A declared column that does not exist in
  the data remains fatal and names itself, because only role-named
  columns survive generation and silently skipping one would drop data
  on a typo.

- `time_jitter` is documented as a realism control rather than a privacy
  one. Every jittered time is clamped inside its own Voronoi cell, so no
  value of `time_jitter` moves a visit more than half a gap from the
  source subject’s visit and the source schedule stays recoverable.

- The last observation’s jitter cell is now bounded above by half the
  final gap instead of being left open, so a large `time_jitter` can no
  longer stretch follow-up arbitrarily past what `screen = TRUE`
  permits.

- `tad`, when declared, is recomputed from each avatar’s own generated
  dose times rather than carried over from the anchor, where it
  described a schedule the avatar no longer has.

- Attendance sampling now carries each visit’s own row metadata. It had
  cloned a single template row per endpoint — the anchor’s *first*
  observation of it — for every visit in the sampled pattern and
  rewritten only `time`, so every other row-varying column arrived from
  that one row: with `nominal_time` and `occasion` declared, the nominal
  column collapsed to one or two values while `time` was correct, and
  the second dose interval was labelled occasion 1. `tad` and `mdv`
  escaped only because they are rebuilt afterwards. The clone is now the
  anchor’s observation nearest the wanted visit, and where the source
  sits on its nominal grid the nominal column follows the sampled time
  exactly. **Generated output changes** for any source declaring
  `nominal_time` or `occasion`; `SIM-035`.
