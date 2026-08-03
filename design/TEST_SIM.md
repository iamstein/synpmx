# Simulation evaluation test plan

## Purpose

This document defines the continuously maintained evaluation suite for
`synpmx`. It turns failures found while developing the demonstrations
into reproducible checks. The suite must answer three separate questions:

1. Is the generated PMX table structurally valid?
2. Does it retain the coarse regimen, sampling, and trajectory behavior needed
   for workflow testing?
3. Do the demonstrations display that behavior honestly and clearly?

The tests are regression and utility tests, not evidence that generated data
are scientifically interchangeable with the source. Privacy-mechanism tests
and the proof argument remain separate. Repeated fitting for evaluation must
use only public or fully simulated fixtures; it must never repeatedly consume
a confidential dataset without explicit privacy accounting and authorization.

This is a living document. Every newly discovered simulation or demonstration
failure should add:

- a row to the issue registry;
- an objective metric or invariant when one is possible;
- a minimal automated regression test; and
- a multi-seed or visual diagnostic when a single deterministic assertion is
  insufficient.

## Dataset registry

The first evaluation set contains every dataset used in the practical demo.
Configuration values below are public domains, contribution limits, schema
semantics, and endpoint alignments. Source-derived regimen and sampling
schedules must not be supplied to the fit.

| ID | Dataset | Origin | Main behavior exercised | Required checks |
|---|---|---|---|---|
| `censoring` | Eight-subject expansion of `pmx_censoring_fixture()` | Package-owned fully simulated fixture | Minimal end-to-end workflow, study-time endpoint, schema restoration, left/right/interval censoring | The source covers all conventions; generated censoring states and limits are coherent; source and synthetic tables are displayed; validation passes |
| `theo_md` | `nlmixr2data::theo_md` | Public package data | Seven Q24H doses, dose-relative log-PK, dense profiles after occasions 1 and 7, occasional trough after occasion 2, no samples after occasions 3--6 | Regimen and sampling are inferred; dose rows are not counted as samples; intensive profiles have one directional peak; inactive occasions stay inactive |
| `warfarin` | `nlmixr2data::warfarin` | Public package data | Lower-case schema, one dose, dose-relative `cp`, study-time `pca`, factors and multiple DVIDs, late PK follow-up | Both endpoints and all subjects remain; CP extends beyond 24 hours; endpoint-specific point counts and follow-up are retained broadly |
| `wbcSim` | `nlmixr2data::wbcSim` | Public package data | Infusion start/stop pairs, study-time log-WBC response, delayed nadir and recovery, numeric covariates | Infusion rows are coherent; no source sentinel schedule is copied; WBC declines and recovers; cohort and follow-up remain comparable |
| `nimoData` | `nlmixr2data::nimoData` | Public package data | Four nominal dose groups, ten approximately weekly infusions, declared OCC/TAD, long terminal follow-up, and a time-varying weight column | DOS is a subject property that conditions amount/rate/duration; every subject has ten coherent infusions; WGT is explicitly excluded; dose-relative sample count and terminal coverage remain broad |
| `skeleton_uniqueness` | Public datasets via `scripts/measure_skeleton_uniqueness.R` | Package-owned report over public package data | Observation-time, observation-count, and event-signature equivalence classes, before and after `coarsen_time` | Coarsening must not increase any class exposure; a source already on its nominal grid must be generated unchanged; subjects still alone after coarsening must raise an alert |
| `mavoglurant` | `nlmixr2data::mavoglurant` | Public package data | One- and two-period profiles, TIME reset within OCC, occasion-varying assigned DOSE, numeric-coded SEX, infusion rows | Reset clocks validate within ID/OCC; DOSE equals positive AMT and is constant within ID/OCC; SEX is categorical; cohort and two-occasion event structure remain |
| `case1_pkpd` | `xgxr::case1_pkpd` | Public package data (`xgxr` is already a Suggests) | **The first public dataset shaped like a real study report.** 180 patients, a declared `NOMTIME`, six treatment arms as `subject_properties`, two endpoints keyed by a character `NAME`, and a baseline weight. Nothing in the nlmixr2data five has a declared nominal time, an arm to stratify on, or a `CENS` column | Grid is `nominal`; 0 patients with a unique observation schedule; both guarantees hold at 0; every arm survives into the output. Note `CENS` is meaningful only for the PK endpoint -- the PD effect is signed, and declaring `cens` is correctly refused by validation |
| `mad` | `xgxr::mad` | Public package data | **Six endpoints**, including ordinal, count and binary PD alongside continuous PD and PK. `warfarin`'s two endpoints are the most the rest of the registry offers, and `SIM-036`'s endpoint-loss failure mode is invisible below that | All six endpoints survive generation; schedule guarantee holds at 0 |
| `pheno_sd` | `nlmixr2data::pheno_sd` | Public package data | 59 **real** patients, the largest real cohort available. Neonatal phenobarbital: individualised dosing, sparse irregular sampling, and a time-varying weight | The observation-side guarantee holds at 0; the dose side does NOT and must say so -- 17 of 59 patients have a dose schedule nobody shares and there is nobody safe to anchor on instead. This is the registry's honest example of a study whose dosing cannot be masked |

Two more worth adding when there is a reason to. `nlmixr2data::nmtest` (54
subjects) is a NONMEM 7.4.3 event-grammar torture test -- steady state, `ADDL`,
`II`, lag time, bioavailability, duration and rate modes together -- which is
the hardest available exercise of role handling and validation, and the one
dataset that would say whether the `addl`/`ii` carry-through is sound. The ACOP
2016 simulated sets (`Oral_1CPT` and siblings, 120 subjects, 7920 rows each)
are the only public sources with `SS`/`ADDL`/`II` populated at scale and would
serve for performance work.

The next expansion dataset should be `pmx_simulated_fixture(60)`. It is already
used by package tests and provides a larger repeated-dose, two-endpoint study
for seed sweeps and privacy-utility experiments. Later additions should cover
IV bolus decay, multiple infusions, irregular dosing, below-quantification
patterns, missing covariates, more than two endpoints, and genuinely
multi-phasic dose-relative profiles.

Public datasets must be loaded from their installed packages during a run. Do
not copy their records into this repository.

## Issue registry

The registry below captures the problems found during prototype and demo
review so far. “Gate” describes the check that should prevent recurrence.

| ID | Dataset/area | Failure observed | Cause or interpretation | Required gate |
|---|---|---|---|---|
| `SIM-001` | All PMX data | Dose/event rows were visually interpreted as samples | Plots and summaries did not consistently separate `EVID != 0` from observations | Observation metrics and plots must exclude event rows; report event and observation counts separately |
| `SIM-002` | All named demos | Regimen and sampling times were specified from disclosed datasets | Demo configuration supplied dose times, counts, amounts, rates, endpoint grids, or occasion schedules | Inspect the fitted public configuration and assert those overrides are absent; infer them inside `synpmx_empirical()` |
| `SIM-003` | `theo_md` | Q24H dosing incorrectly appeared to imply Q24H sampling | One local endpoint grid was repeated after every generated dose | Test occasion activation separately from conditional sample count; occasions 3--6 must receive no observations |
| `SIM-004` | `theo_md` | The sparse occasion-2 trough was confused with a full profile | Presence probability and sample density were represented by one unconditional count | Test the two-part release: occasion presence probability and observations conditional on being sampled |
| `SIM-005` | `theo_md` | Automatic local grid extended to the next dose and reassigned a late sample to the following occasion | Generic horizon and timing jitter allowed equality with the next dose boundary | Every dose-relative time must remain strictly before the next dose and retain its generating occasion |
| `SIM-006` | `theo_md` | Synthetic PK profiles differed markedly and developed an artificial trough followed by another peak | Unoccupied cells on a wide log domain were decoded as domain-midpoint measurements instead of missing support | Interpolate unoccupied cells from released occupied neighbors; never treat absence as a DV measurement |
| `SIM-007` | `theo_md` | Individual profiles still showed secondary peaks after the mean curve was repaired | AR(1) noise continued across distant dose occasions and unconstrained residuals reversed the post-peak direction | Restart serial noise by occasion; if the released curve is approximately unimodal, require at most one directional peak per generated occasion |
| `SIM-008` | `warfarin` | CP observations disappeared after about 24 hours | Count matching retained the first grid cells and systematically trimmed late cells | Select timing cells with released presence probabilities; require late CP coverage and a median subject maximum time of at least 72 hours |
| `SIM-009` | `warfarin` | CP appeared to contain fewer patients or missing data | Endpoint allocation and combined plotting obscured subject/endpoint coverage | Assert cohort size, patients per endpoint, endpoint set, and mean points per patient; facet endpoints separately |
| `SIM-010` | All demos | Synthetic datasets sometimes contained fewer subjects than their sources | Generation size was an independent default rather than the fitted privacy-accounted count | Omitted `n_subjects` must use the fitted count; public-fixture demos require exact cohort equality |
| `SIM-011` | All demos | Number of generated time points could differ materially from source | Full grids were repeated or trimmed without respecting released observation totals | Compare endpoint-specific observations per subject with explicit tolerances and test the released total observation count |
| `SIM-012` | `wbcSim` | Infusion/event behavior could be incoherent or reproduce an exceptional source schedule | Start/stop construction and generalized-regimen logic were incomplete | Require paired positive/negative amount and rate rows, coherent duration, bounded values, and absence of the 4580-hour source schedule |
| `SIM-013` | `wbcSim` | Longitudinal response could lose delayed decline, nadir, or recovery | Study-time behavior was at risk of being restarted at a dose or shortened by grid selection | Require a value below baseline followed by recovery and broad late follow-up |
| `SIM-014` | All endpoints | Exact source timing vectors could be copied | Earlier designs considered source anchors or exact schedules | Compare complete source and generated time vectors; no generated vector may be identical to a source vector |
| `SIM-015` | All data | IDs, schema classes, factor levels, or endpoint columns could be lost or reused | Schema restoration and ID generation were incomplete | Require new IDs, original column order/classes, declared public factor levels, and all endpoints |
| `SIM-016` | Demo plots | Source and synthetic observations were not overlaid or comparably displayed | Separate plotting code and incorrect grouping | Use one comparison data frame, consistent colors, and source-above/synthetic-below facets |
| `SIM-017` | Demo plots | Lines stopped between profile segments | Lines were grouped by event/occasion rather than by subject for study-time displays | Study-time plots group all observed points by dataset, subject, and endpoint; event rows remain excluded |
| `SIM-018` | Demo plots | Log scaling obscured Warfarin and other comparisons | `xgx_scale_log10()` was applied when a linear comparison was more interpretable | Demo comparison plots use linear DV axes unless a dataset-specific design decision explicitly changes this |
| `SIM-019` | Demo plots | Multiple DVIDs and source/synthetic differences were difficult to see | Overlay-only panels were too dense | Facet source above synthetic and endpoint columns side by side while retaining dataset colors |
| `SIM-020` | Decoding, all noisy fits | Grid cells with no source support decoded to the bottom of the endpoint working domain, producing spurious deep troughs on log endpoints at every cohort size | Decoders compared released presence fields against the bare constant `0.25`, but those fields are unnormalized subject counts of order N, not probabilities. Laplace noise around zero passed the gate, then the separately noised value sum clamped to zero and `0 / noise` decoded to the domain floor. `.fill_unoccupied_curve()` could not repair it because the cell was already labelled occupied | Gate released presence on `.support_threshold(count, noise_scale)`, derived from the release's own `sensitivity / epsilon`. A cell with no support must never decode to a working-domain endpoint; see `tests/testthat/test-decode-support-threshold.R` |
| `SIM-020` | Demo output | The actual source and generated tables were not shown | Vignette focused only on plots and summaries | Every demo prints a small source preview and a small synthetic preview, with public-data labeling |
| `SIM-021` | Documentation | Readers expected an ODE/NLME or spline model | The implemented fixed-grid, smoothing, interpolation, timing, and variability model was not stated early or precisely | Vignette test/check confirms the method introduction names fixed-grid summaries, 1--2--1 smoothing, linear interpolation, and the absence of ODE/NLME/splines |
| `SIM-022` | Vignette artifacts | A repaired source vignette could still appear broken | An ignored, previously rendered HTML file or already loaded installed namespace was stale | Evaluation rebuilds vignettes from the current source package in a clean library and records package/source identity in its manifest |
| `SIM-023` | Demo plotting | A narrow-bin median line itself could look jagged or multi-peaked | Sparse subjects at slightly different actual times populated different 0.25-hour bins | Treat thick summary lines as diagnostics only; compute peak tests from complete ordered profiles and include the released mean curve in the report |
| `SIM-024` | Cohort assignments | `ACTARM`, `TRT`, or nominal dose group could be sampled independently of the generated regimen | Ordinary categorical covariates retained only marginal proportions and had no link to event summaries | Declare subject-level assignment fields as `subject_properties`; release stratum count and regimen jointly; require generated property-to-event coherence |
| `SIM-025` | Numeric category codes | A numeric-coded category such as mavoglurant `SEX` could be modeled as a continuous covariate | Numeric storage class was treated as sufficient evidence of continuity | A declared public category domain forces categorical modeling regardless of storage type |
| `SIM-026` | Occasion-assigned dose | Mavoglurant `DOSE` could be missing, vary within an occasion, or disagree with generated AMT | Nominal dose was restored as an unmodeled schema default | An `assigned_dose` role is reconstructed from each generated positive event AMT and validated on every ID/OCC profile |
| `SIM-027` | Reset clocks | Mavoglurant source TIME restarts at zero in the second OCC and was rejected as globally decreasing | Validation and plotting assumed one monotone clock per subject | When OCC is declared, validate time order within ID/OCC and group dose-relative display lines by subject/occasion |
| `SIM-028` | Positive-rate infusions | NimoData and mavoglurant encode a positive rate without explicit negative source stop rows, so duration could decode as zero | Duration inference looked only for a later negative event row | If no stop row exists, infer bounded duration as AMT/RATE; generated starts and stops must remain coherent |
| `SIM-029` | Terminal washout | A final profile longer than one dose interval was compressed into the ordinary interdose window | The automatic local basis was rescaled to one interval for every occasion | Keep nonterminal observations before the next dose, but let the final occasion use its released occupied horizon |
| `SIM-030` | Time-varying covariates | NimoData WGT changes within subject and could be mislabeled as a baseline covariate | The prototype generates only subject-constant baseline covariates | Validation rejects varying baseline covariates; explicitly exclude WGT until a longitudinal-covariate model exists |
| `SIM-031` | AVATAR schema | Excluded columns were copied back into the generated table and factor IDs became missing when assigned fresh labels | Synthesis operated on the full input and assigned unseen values into the source factor levels | Remove excluded columns before synthesis; extend factor ID levels; regression coverage in `test-avatar.R` |
| `SIM-033` | AVATAR timing | Every avatar reproduced one source subject's exact observation time vector, and nothing caught it | `SIM-014`'s gate ("no generated vector may be identical to a source vector") was implemented only against the structural/DP path, in `test-structural-v3.R` and `scripts/evaluate_simulations.R`, which never call `synpmx_avatar()`. AVATAR is the mode that actually copies an event skeleton verbatim (`synthesis.R` stage d), and `time_jitter` defaults to 0. Raising it does not help either: `.offset_unique_times()` clamps each time into its own Voronoi cell, so at any magnitude a visit stays within half a gap of the source value | `coarsen_time = TRUE` snaps the source onto a shared visit grid before profiling and resamples pooled deviations afterwards; `tests/testthat/test-avatar-coarsen.R` asserts SIM-014 on AVATAR output and asserts it *fails* with `coarsen_time = FALSE`, so the defect stays pinned |
| `SIM-034` | AVATAR timing, inferred grid | Coarsening could silently collapse nothing, and the residual it legitimately cannot reach looked like a failure | Two faults. (1) A single *global* merge threshold is hostage to the tightest pair of samples anywhere in the study: one subject with two draws minutes apart pinned the whole timeline, leaving warfarin and mavoglurant entirely uncoarsened. Merging agglomeratively, closest pair first, with the subject-disjointness constraint checked per merge, makes the decision local; a median-spacing guard stops an outlier's lone far-out observation from merging with a shared visit because their subject sets happen not to overlap. theo_md then goes 12/12 alone to 0/12 and mavoglurant's unshared times 35 to 1. (2) `obs_time_class` conflated two causes needing opposite remedies: a subject observed at a moment nobody else was (the grid's job) versus one whose times are all shared and only whose *attendance pattern* is unique (dropout, which no grid can touch) | `.derive_time_grid()` merges locally; `skeleton_uniqueness()` reports `min_time_share` so the two causes separate, and `synpmx_avatar()` raises a different alert for each, since only the first is fixed by declaring `nominal_time`. `tests/testthat/test-avatar-coarsen.R` pins the outlier guard |
| `SIM-037` | AVATAR endpoint identity, undeclared | The residual of `SIM-036`. With **no** `dvid` and no `cmt` to test against, a two-endpoint source is pooled into one endpoint and the values are not merely rescaled but destroyed: measured on a 12-subject PK/PD source, PD rows of 45.00-54.95 came back as 3.41-5.38 while PK rows of 7.20-8.79 came back as 14.88-23.19 — inverted, because the interleaved series `0, 50, 8, 40, 4, 30, 1, 45` is read as one endpoint's time course, smoothed, and blended. Row count is preserved, and an undeclared `CMT` is dropped from the output, so nothing looks wrong. Raised 2026-07-30 as "what if they declare neither?" | Same root cause as `SIM-036` — `.endpoint()` labels every row `"DV"` when `dvid` is absent — but undetectable in general: two endpoints measured at disjoint times (the `warfarin` shape) are indistinguishable from one endpoint measured at all of them. Repeated subject-and-time observations are the one visible hint, and they were **rejected as a trigger for refusal** by the owner, correctly: replicate measurements at one time are ordinary and refusing would be wrong | `.note_single_endpoint()` states the consequence whenever `dvid` is absent — one endpoint, one transform, one censoring boundary — and reports the repeated subject-and-time count when there is one, as information rather than a refusal. Accepted as a partial defense: a two-endpoint source sampled at disjoint times with no endpoint key declared is **not** detectable and will still be pooled. `tests/testthat/test-avatar-dose-attendance.R` asserts the note fires without `dvid`, names the repeats when present, does not error, and is silent once `dvid` is declared |
| `SIM-036` | AVATAR endpoint identity | A NONMEM-style source whose endpoints are distinguished only by `CMT`, with PK and PD drawn at the same visits, lost an entire endpoint: 108 rows in, **60 out**, every `CMT = 3` row gone, under default settings and with no warning. Found 2026-07-30 from the question "how is CMT used vs DVID, are both needed?" | `.endpoint()` returns the constant `"DV"` when no `dvid` is declared — it never falls back to `cmt`, which is read only on *dose* rows by `.route_key()`. Every observation therefore carried one endpoint label, so `.attendance_key()`'s `endpoint@time` cells collapsed a PK and a PD observation at the same time into one, and `.apply_attendance()` rebuilt one row where the source had two. Predates `SIM-035`: cloning `candidate[1L]` lost the same rows. Compounding it, `pmx_roles()` rejected `cmt = "CMT", dvid = "CMT"` as a duplicate role, so the natural fix was unavailable and users had to copy the column to itself | Three parts. (1) `pmx_roles()` now permits one column to be both `cmt` and `dvid`; every other collision stays an error. (2) `.require_endpoint_key()` refuses to generate when observations occupy more than one compartment and no `dvid` is declared, naming the fix in the message — refusing beats inferring, since `CMT` is not reliably the endpoint. (3) Independently, `.attendance_key()` numbers repeated `endpoint@time` cells and `.apply_attendance()` gives the n-th repeat the n-th nearest source row, so no collapse can delete an observation even where `dvid` *is* declared and a source holds duplicate records. `tests/testthat/test-avatar-dose-attendance.R` covers all three; `tests/testthat/test-roles-validation.R` pins the permitted overlap |
| `SIM-035` | AVATAR attendance sampling | With `nominal_time` and `occasion` declared, every generated observation row carried the *first* visit's metadata: on `pmx_simulated_fixture()` the nominal column collapsed from 14 distinct values to two (0 and 0.25) while `TIME` was correct, so `TIME` and `NTIME` disagreed on nearly every row, and the second dose interval was labelled occasion 1. Found 2026-07-30 while extending the README example to declare every role | `.apply_attendance()` cloned **one template row per endpoint** — `candidate[1L]`, the anchor's first observation of that endpoint — for every visit in the sampled pattern, then overwrote only `roles$time`. Every other row-varying column rode along from that one row. `TAD` hid the defect because `.recompute_tad()` rebuilds it afterwards, and `MDV` because `.derive_standard_mdv()` does; nothing rebuilds nominal time or occasion. The mechanism landed 2026-07-29 (`REV-026`, `min_pattern_share`) and is on by default, so any dataset declaring those roles was affected | The clone is now the anchor's observation of that endpoint **nearest in time** to the wanted visit, so row-varying columns arrive from the visit they belong to; where the source sits on its nominal grid the nominal column follows the sampled time exactly. `tests/testthat/test-avatar-dose-attendance.R` asserts `TIME == NTIME` row for row on output, that occasion 2 covers exactly the second dose interval, and that the nominal column keeps more than ten distinct values |
| `SIM-044` | AVATAR dose truncation | The `SIM-043` guarantee protected a patient who stopped dosing early by declining to build on them, which removed the regimen from the output entirely: nineteen patients on three doses, one on two and one on one came out as twenty-one on three. Raised 2026-08-03 by the owner, who asked for the case to be resolved rather than merely caught | Dose events are deliberately never redrawn, because a drawn regimen could be one no protocol permits. But there is one edit that escapes the objection: truncating a schedule at one of its OWN dose times yields a regimen the study actually gave someone. That was not being used | `.dose_truncation_plan()` treats a schedule that is a PREFIX of the cohort's full one as a depth, and re-truncates a depth held by one patient to one held by several or by nobody, walking **downward only** -- truncation removes dose rows and cannot add one, so a deeper target silently leaves the schedule where it was and puts the avatar back on the uniquely-held depth (found in testing). On a 32-patient fixture with four patients stopping at depths 4, 6, 7 and 9, output depths become 5, 8 and 10 with zero identifying schedules, where re-anchoring alone had put every avatar on the full schedule. It does not rescue every cohort and the arithmetic says which: depths `{1:1, 2:1, 3:19}` has no free depth between the two singletons, so that case still falls back to M2 and the report shows the regimen coverage lost |
| `SIM-043` | AVATAR dose-schedule guarantee | The visit-set guarantee (`SIM-042`) covered observation times and left dose times unguarded, so a patient whose set of dose times nobody shares handed it to every avatar built on them. Raised 2026-08-03 by the owner after the gap was reported: "please do fill that gap with regards to dose as well" | Attendance sampling deliberately never touches dose events -- resampling them risks emitting a regimen the protocol never permitted, which is a worse error than the one it fixes -- so the dose side had no mechanism at all. Measured after coarsening: `theo_md`, `warfarin` and `mavoglurant` have one shared dose schedule each and are unaffected, `wbcSim` has 3 of 45 alone, and `nimoData` has 12 of 12 | The regimen is not rewritten; a patient whose dose schedule nobody shares is simply not built upon, and the avatar is re-anchored. Where **every** patient is in that position, re-anchoring cannot help, so it is not attempted and the run alerts once instead of silently doing nothing. `identifying_dose_schedules` records the result and is 0 for every public dataset except `nimoData`, which honestly reports 12 of 12. Two false-check traps were found on the way and are worth not repeating: comparing recorded times directly always reports zero, because generation resamples the deviations coarsening removed; and snapping the finished table back onto the grid produces false alarms, because those deviations are comparable to the grid spacing on an irregular study (`wbcSim`). Both checks are therefore recorded per avatar, on what was actually applied. `tests/testthat/test-avatar-dose-attendance.R` recomputes both from the tables independently of the recorded settings |
| `SIM-042` | AVATAR visit-set guarantee | One avatar in 21 was emitted carrying a visit set exactly one real patient holds -- the property `min_pattern_share` exists to guarantee. Reported as GitHub issue #2. The immediate cause was a schedule group of ONE: groups are (stratum x endpoint set), and a single patient measuring a different set of endpoints from everybody else forms a group that can never hold a shared visit set, so there was no pool to draw from and nothing to substitute. The deeper cause is the pattern the owner named -- "not sure why things keep slipping through" -- which is that every mechanism here was verified on its own terms and nothing asked the whole question of the finished table | Three gaps, each individually reasonable. (1) The visit-set pool is grouped by stratum as well as endpoint set, so a small arm holds too few patients to share anything, even though donors already cross strata freely. (2) A group of one can never be rescued by any pooling. (3) The generation loop copied the anchor's own set as its last resort, without asking whether that set was one nobody else had | Three parts. (1) A second pool keyed on the **endpoint set alone**, borrowed from when a group has none of its own; the endpoint set is never crossed, since an avatar must not be handed a set naming an endpoint its anchor lacks. (2) The visit set is decided **before** the avatar is built, and an anchor that cannot be masked causes that avatar to be **re-anchored** rather than the patient to be dropped -- every source patient stays a donor and stays available to anchor others. An earlier attempt excluded such patients up front and removed 48% of a cohort to fix one avatar, which is why the decision is made per avatar. (3) An end-to-end check on the finished table, recorded as `identifying_visit_sets`, which must be 0 and is what the guarantee is now stated in terms of. Zero across `theo_md`, `warfarin`, `wbcSim`, `nimoData` and `mavoglurant`; `tests/testthat/test-avatar-dose-attendance.R` pins it on a lone-endpoint-set fixture and recomputes it from the tables independently of the recorded setting |
| `SIM-041` | Identifiability screen | `flag_identifiable_subjects()` flagged 10 of 21 avatars on follow-up time alone, including two whose follow-up agreed to three decimal places (109.0458 and 109.0438). Raised 2026-08-03 by the owner: "there are a bunch with similar follow-up times, especially at the long end. I wonder if the criteria are still too strict" | `.modified_z()` scales by the median absolute deviation, and the MAD is **zero** whenever more than half the values are identical. On trial data that is the ordinary case, not a degenerate one: most patients complete the protocol and stop at the same visit. The zero branch then assigned `sign(x - median) * Inf` to every other value, and `flag_identifiable_subjects()` flags infinite scores unconditionally, so every patient who stopped even slightly early was an outlier by construction | The standard Iglewicz-Hoaglin fallback: when the MAD is zero, scale by `1.253314 *` the **mean** absolute deviation from the median, which is zero only when every value is identical -- and then nothing departs from it to flag. On the reported cohort shape this goes from 9 of 20 flagged to 0, with the most extreme score at 2.61 against a threshold of 3.5, while a genuine extreme (a 4580-hour follow-up in a ~650-hour cohort) still scores 12.5. `mavoglurant` goes from 63 of 120 flagged to 17; `wbcSim` and `theo_md` are unchanged, since their spread is real and the MAD path still applies. `tests/testthat/test-flag-identifiable.R` pins the clustered case, the genuine extreme, and the no-variation case |
| `SIM-039` | AVATAR attendance placement | 86% of avatars on a real 21-patient study kept their **anchor's own visit set** -- one real patient's exact pattern of absences, copied -- which is the single thing `min_pattern_share` exists to prevent. The run reported it only as a fraction; nothing alerted. Raised 2026-08-03 by the owner asking why the number was so high when `nominal_time` was already declared | Two independent faults, both hit by ordinary dropout. (1) `.place_attendance()` retried a **deterministic** placement: `trailing` (which is every discontinuation) has exactly ONE arrangement for a given miss count, so when a real patient had dropped out at that depth the placement was correctly rejected as theirs and the other 23 tries re-proposed it unchanged. (2) Shapes are keyed on (kind, count), and under staggered discontinuation every count has exactly one holder, so no shape cleared the floor either and the group got **no pool at all** -- even though "these patients dropped out" is plainly shared by all of them | (1) Placements are enumerated rather than resampled, and the miss count walks outward (k, k+1, k-1, ...) when nothing at the wanted count is free, with the complete pattern as a last resort; `pattern_shifted_fraction` records when the count had to move. (2) A third, coarser abstraction -- the kind alone -- is reached only when the (kind, count) level clears nothing. On a 21-patient staggered-discontinuation fixture, avatars carrying a real patient's exact visit set go from 21 of 21 to 0 of 21. A saturated grid, where every arrangement is genuinely taken, still falls back, and `synpmx_avatar()` now alerts past 10%. `tests/testthat/test-avatar-dose-attendance.R` pins the rescue, the saturated case, and the no-pool case separately |
| `SIM-040` | AVATAR dose basis | A weight-based study came through with amounts copied verbatim, so every avatar's implied mg/kg was its anchor's rather than its own -- 23 distinct mg/kg levels out of a protocol with three -- and the copied milligrams still encoded one real patient's weight. Raised 2026-08-03: "It's actually missing that the dose was given on a weight-based basis" | `.detect_dose_basis()` must *infer* proportionality, and infers conservatively by design: the dose-to-covariate ratio has to collapse onto a handful of protocol levels within 2%. A study that escalates within patient by factors that are not identical across patients, and dispenses in vials, produces ratios that do not cluster. Refusing is the right call for an inference engine -- rewriting amounts on a study that is not proportional would be worse -- but there was no way for a caller who *knows* the study is weight-based to say so | `pmx_roles(dose_covariate = )` names the covariate outright and skips inference. The declared path holds each dose row's **own** ratio instead of snapping to shared levels, so intra-patient escalation survives exactly. It must also appear in `covariates`, since the amount is rebuilt from the avatar's *blended* value. Inference remains the default and now records why it declined, so "not detected" and "never attempted" are distinguishable in a run report |
| `SIM-038` | AVATAR time keys | Nearly every patient looked as though they held a visit set nobody else shared, and the "unique observation times" alert fired on cohorts already sitting on a shared grid. On `warfarin` the coarsened source reported 2 patients with a one-off observation time when the true count is 0, and a real study reported 15 of 17 visit sets discarded as too rare to reuse. Found 2026-08-03 by drawing `plot_pmx_schedule()` next to `skeleton_uniqueness()` and finding the picture disagreed with the table | `format(x, digits = 12, trim = TRUE)` fixes ONE layout for the whole vector it is given, so hour 12 keys as `"12"` for a patient sampled only on the hour and `"12.00000000000"` for a patient who also has a 1.2142857 sample. Both `.attendance_key()` (built per subject, then compared across subjects) and `skeleton_uniqueness()`'s per-time lookup did exactly that. One real visit split across several keys, so patterns had one holder each, `min_pattern_share` discarded them, and their avatars fell back to the anchor's own visit set — the outcome the mechanism exists to prevent | `.time_key()` formats element-wise with `%.12g`, which cannot depend on the rest of the vector, and both call sites use it. `tests/testthat/test-avatar-coarsen.R` pins that two patients with identical visit sets key identically and that a shared time is never counted as a one-off, in both cases with a fractional time present in one patient's vector and not the other's |
| `SIM-032` | AVATAR donor compatibility | Profiles with different numeric dose magnitudes could exchange trajectories because signatures retained only event sign | Event compatibility omitted exposure scale | Include rounded AMT/RATE magnitude in signatures; deterministic two-dose regression in `test-avatar.R` |

When an issue is fixed, keep its row. The registry is a permanent record of
the behavior that must continue to work.

### Related setup and documentation issues

Several problems encountered during the same review are important preflight
checks but are not simulation-fidelity metrics:

- The OpenDP production backend is an R package dependency. The evaluation
  runner should print `dp_backend_status()` and clearly distinguish an
  unavailable production backend from the explicitly nonprivate public-fixture
  backend.
- Vignette PDF/manual builds require a working `pdflatex`; TinyTeX installation
  and PATH discovery are environment checks rather than generator tests.
- Privacy concepts and the simulation algorithm were initially mixed in one
  method vignette. The package now requires separate privacy-introduction and
  simulation-method vignettes, plus the practical demo.
- The practical vignette must explain near the top that the implementation is
  a fixed-grid population generator, not an ODE, NLME likelihood, or spline.
- Evaluation should report missing optional packages and rendering tools as
  explicit preflight failures or skips, never as simulation successes.

### Existing coverage reused

The evaluation layer builds on, rather than silently replacing, these existing
checks:

- `tests/testthat/test-integration-nlmixr2data.R` covers inferred Theophylline
  dosing/sampling and single peaks, late Warfarin CP, WBC infusion/recovery,
  cohort size, schema, and timing non-copying.
- `tests/testthat/test-private-fitting.R` covers automatic grid bases,
  unoccupied-cell interpolation, and the separate occasion-presence and
  conditional-count release.
- `tests/testthat/test-generation.R` covers cohort-size defaults, timing-cell
  selection, schema restoration, repeated-dose versus study-time behavior,
  chronological coherence, reproducibility, subject-property/regimen
  coherence, and assigned-dose reconstruction.
- `tests/testthat/test-censoring.R` covers generated censoring conventions.
- `scripts/demo_nlmixr2data.R` currently performs cohort, point-count, and
  follow-up checks and constructs the visual comparison panels.

The shared evaluator described below closes the earlier gaps in metric reuse,
systematic multi-seed evaluation, machine-readable results, plot-semantic
assertions, and one report spanning all datasets and known gates. Future
coverage should extend this layer rather than create a competing metric path.

## Implemented evaluation

Use one shared metric layer with two callers.

### `tests/testthat/helper-simulation-evaluation.R`

This helper contains no test expectations and writes no files. It provides:

- a dataset registry returning source data, roles, endpoints, public bounds,
  public schema semantics, contribution limits, and budget allocation;
- row classifiers that distinguish dose/event, infusion-stop, and observed
  endpoint rows;
- regimen summaries by subject: dose count, interval, amount, rate, and
  infusion duration;
- sampling summaries by endpoint and occasion: activation, conditional count,
  total count, timing-cell coverage, first/last time, and late-follow-up
  coverage;
- trajectory summaries on study-time or TAD bins: bounded quantiles, peak or
  nadir time, directional peak count, decline/recovery indicators, and broad
  range;
- schema, ID, endpoint, timing-vector-copy, and PMX-validity checks; and
- plot-data checks for facet order, group identifiers, event exclusion, axis
  scale, and consistent dataset labels/colors.

Metrics should be ordinary data frames with stable column names so the same
code can feed `testthat` expectations and the longer report script.

### `tests/testthat/test-simulation-evaluation.R`

This is the fast deterministic regression gate run by `devtools::test()` and
`R CMD check`. It:

1. Runs package-owned fixtures unconditionally.
2. Runs `nlmixr2data` cases with `skip_if_not_installed("nlmixr2data")`.
3. Uses the guarded noiseless public-fixture backend and fixed generation
   seeds.
4. Fits without source-derived regimen or sampling overrides.
5. Asserts all hard invariants and the dataset-specific gates below.
6. Avoids writing plots or evaluation artifacts.

Some overlapping assertions remain in
`tests/testthat/test-integration-nlmixr2data.R` as focused regression tests.
Consolidate them only after equivalent coverage is demonstrated, and avoid
duplicating two independent implementations of the same metric.

### `scripts/evaluate_simulations.R`

This is the longer evaluation runner. A typical invocation is:

```sh
Rscript scripts/evaluate_simulations.R \
  --datasets=all \
  --seeds=101:200 \
  --backend=public \
  --output=output/simulation-evaluation
```

The script uses base-R argument parsing to avoid a new dependency. It:

1. Loads the registry from package/test helper code without duplicating dataset
   definitions.
2. Fits each public dataset once for an ordinary generation-seed sweep.
3. Generates every requested seed and computes the shared metrics.
4. Records hard failures immediately but finishes the run so all problems are
   visible together.
5. Produces source-above/synthetic-below study-time and scientific-clock plots,
   endpoint sampling panels, regimen summaries, and released-curve panels.
6. Writes only derived metrics and figures by default, not source or synthetic
   row-level datasets.
7. Exits nonzero when a hard gate fails.

Outputs under the ignored `output/` directory are:

```text
output/simulation-evaluation/
  run-manifest.txt
  metrics-by-seed.csv
  gate-results.csv
  failures.csv
  regimen-by-fit.csv
  sampling-by-fit.csv
  subject-properties-by-fit.csv
  summary.html
  figures/
    censoring-study-time.png
    theo_md-study-time.png
    theo_md-tad.png
    warfarin-study-time.png
    wbcSim-study-time.png
    nimoData-tad.png
    mavoglurant-tad.png
```

The manifest should record package version, Git commit when available, dirty
worktree status, R and dependency versions, dataset package versions, backend,
bounds, contribution limits, generation seeds, timestamp, and platform. It
must also state whether the run used the public-fixture or formal DP backend.

### Optional privacy-utility sweep

A separate mode may repeatedly fit the public/simulated datasets across
epsilon values and production-backend noise draws. Because fitting intentionally
does not accept a user-controlled privacy-noise seed, this mode should report
distributions rather than exact snapshots. It must never be pointed at a
confidential dataset merely for convenient tuning.

## Initial gates and tolerances

### Hard gates for every dataset

- `validate_pmx(..., strict = TRUE)` passes.
- Generated IDs are disjoint from source IDs.
- Column names, order, practical classes, declared factor levels, and endpoint
  set are restored.
- Every observation has finite DV, a declared endpoint, and zero or missing
  amount; every event has coherent event fields.
- Times are finite and bounded. They are nondecreasing within subject, or
  within subject and declared occasion when the source clock resets by OCC.
- No complete generated timing vector is identical to a source timing vector.
- The fitted model contains no raw source rows, identifiers, or unnoised
  aggregates.
- Named demo configurations contain no `dose_times`, `dose_interval`,
  `n_doses`, `dose_amount`, `dose_rate`, `infusion_duration`,
  `endpoint_grids`, or `endpoint_occasion_grids` derived from the source.
- With the public-fixture backend and omitted `n_subjects`, generated and
  source cohort sizes are equal.
- Every expected endpoint is represented in every dataset-level result.

### General utility gates

- Mean observations per subject and endpoint differ by no more than
  `max(1, 25% of the source mean)` in the deterministic public-fixture run.
- First and last endpoint observation times differ by no more than 20% of the
  declared public time span, unless a stricter dataset-specific gate applies.
- Event rows and observation rows are counted and reported separately.
- Multi-seed evaluation reports median, 5th, and 95th percentiles rather than
  relying on one favorable seed. Hard-invariant failures are never tolerated;
  soft utility gates should pass for at least 95% of the planned seeds.

These are initial engineering tolerances, not scientific equivalence margins.
Any change must be justified in this file and must not be made only to excuse a
new regression.

### `theo_md`

- Fitted generalized dose count is 7 and interval is 24 hours within 0.1 hour.
- Generated subjects retain all seven dose events.
- Released sampling probabilities identify occasions 1 and 7 as active,
  occasion 2 as sparse, and occasions 3--6 as inactive.
- Intensive occasions contain approximately 10--11 observations and the
  sparse occasion contains zero or one according to its fitted activation
  draw.
- Dose-relative observations stay strictly within their generating occasion.
- Every intensive generated profile has at most one directional peak after
  zero-length/flat differences are removed.
- Peak occurs after the first observation and is followed by an overall
  decline to the last observation.

### `warfarin`

- Both `cp` and `pca` are present with the lower-case source schema.
- Generated cohort and patients represented per endpoint equal the source
  cohort for the public fixture.
- CP observations per subject differ from source by no more than one in the
  deterministic run.
- At least one generated CP observation is later than 24 hours and the median
  subject-specific final CP time is at least 72 hours.
- Pca remains a study-time trajectory and is not restarted as a dose-relative
  excursion.

### `wbcSim`

- Positive infusion-start and negative infusion-stop rows both exist.
- Paired generated amount and rate fields are coherent and bounded.
- The exceptional 4580-hour source schedule is not reproduced.
- Every subject has a post-baseline value below baseline and a later recovery
  above its nadir when enough observations are generated.
- Late WBC follow-up remains within the general coverage tolerance.

### `nimoData`

- `subject_property_summary()` contains public DOS strata 50, 100, 200, and
  400 with ten inferred doses in every stratum.
- Every generated subject has one constant DOS value, ten positive dose events,
  and positive AMT equal to DOS.
- Positive infusion starts and negative generated stop rows both exist; fitted
  duration can be inferred from source AMT/RATE when no source stop row exists.
- The final occasion may extend beyond one weekly interval, while observations
  on occasions 1--9 stay before the following dose.
- WGT is absent from generated schema because it is explicitly excluded; BSA,
  AGE, and HGT remain constant within subject.

### `mavoglurant`

- The 120-subject public cohort is preserved by default.
- TIME ordering is valid within ID/OCC even though the source clock resets.
- Every generated subject has two positive dose occasions in the deterministic
  evaluation case.
- Generated DOSE is finite on every row, constant within ID/OCC, and equals
  positive event AMT.
- Numeric-coded SEX uses only its declared categorical levels.
- The current generalized-regimen case does not promise preservation of the
  source crossover-sequence distribution when no ACTARM/TRT/sequence property
  is present; that remains an explicit utility limitation, not a passed gate.

### Censoring fixture

- The source fixture contains uncensored, left-, right-, and interval-censored
  examples; generated states use only valid codes and include censoring.
- DV, CENS, and LIMIT combinations satisfy the declared convention.
- Censoring is applied to generated latent values, not copied source rows.

## Plot and vignette evaluation

Pixel snapshots are likely to be brittle across R, graphics-device, and font
versions. Automated tests should inspect plot semantics instead:

- plot data contain observations only;
- grouping connects a subject's chronological observations across study-time
  gaps;
- facet levels place Source above Synthetic and endpoints in columns;
- colors map consistently to Source and Synthetic;
- demo DV scales are linear; and
- every plot contains all expected subjects and endpoints.

The long evaluation script should still render PNGs for human review. Reviewers
should look for disappearing patients, lost late samples, unintended repeated
profiles, false peaks or nadirs, axis compression, connections through event
rows, and misleading summary lines. A visual concern should be converted into
an objective regression metric whenever possible.

Vignettes must be rendered from a clean temporary library containing the
current source build. Do not validate a stale HTML file or an older namespace
already loaded in an interactive R session.

## Where the tests should live

The core regression checks belong in `tests/testthat/`. Your instinct is right:
patient counts, endpoint coverage, inferred occasion sampling, late Warfarin
CP, infusion coherence, one-peak Theophylline profiles, schema restoration,
and timing-vector non-copying are package behavior. They should fail during
ordinary development and `R CMD check`, close to the code change that caused
the regression.

The complete evaluation simulation should not live only in the test suite.
Multi-seed runs, repeated privacy-noise fits, figure rendering, HTML reports,
and artifact writing are slower and partly judgment-based. Put that runner in
`scripts/evaluate_simulations.R` and keep its outputs under ignored `output/`.

The recommended division is therefore:

| Location | Responsibility |
|---|---|
| `tests/testthat/` | Fast deterministic invariants and minimal regressions for every issue that can be asserted objectively |
| `scripts/evaluate_simulations.R` | Multi-seed, multi-backend, distributional, and visual evaluation with persisted reports |
| `design/TEST_SIM.md` | Dataset registry, issue history, gates, tolerances, and rationale |

If only one location were allowed, choose `tests/testthat/`, because a test
that is not run automatically will eventually be forgotten. In practice the
hybrid is stronger: every issue first receives the smallest reliable automated
test, while the evaluation script answers broader questions that are too slow
or visual for package checks.

## Maintenance workflow

1. Reproduce a newly reported problem with a public or simulated dataset.
2. Add it to the issue registry before changing thresholds or implementation.
3. Add the smallest deterministic failing test.
4. Fix the implementation and run the focused test.
5. Run the complete `testthat` suite and `R CMD check`.
6. Run the multi-seed evaluation for affected datasets and inspect its plots.
7. Record intentional threshold or dataset changes in this document and NEWS.

No issue is considered closed solely because one plot or one seed looks good.
