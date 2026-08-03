# Checks of the synthetic data

**Status: specification for a vignette not yet written.** Everything below is
the plan for `vignettes/synthetic-data-checks.Rmd`, plus the reasoning behind
the taxonomy. Raised by the owner 2026-08-03, after a day of finding leaks one
at a time on a real study:

> There should be a vignette of all the checks we do of the synthetic data.
> [...] This is worth being really, really clear about because this could help
> anyone who's developing a synthetic data method.

That second sentence sets the scope. This is not only a user guide to
`synpmx`'s diagnostics — it is the *category list* a person building any
synthetic-data method should work through, with `synpmx` as the worked example.
The categories are what generalize; the function names are not.

Internal record, per `AGENTS.md`: nothing shipped may cite this path.

## Why this document exists

Every defect fixed on 2026-08-03 (`SIM-038` through `SIM-044`, plus the
identifiability screen) was found by *looking at output*, not by reasoning about
the algorithm. Several had been in place for weeks behind mechanisms that were
individually correct. The list below is an attempt to make that systematic
rather than lucky.

## Where the vignette belongs

`vignettes/`, not `vignettes/articles/`. It is method reference rather than
supporting evidence, and it should be rebuilt by `R CMD check` on every
behavioral change — a check that stops being true is exactly what needs to fail
loudly. This makes it the fifth shipped vignette, which is a real recurring
cost; that cost is the argument for keeping it tight and example-driven rather
than exhaustive.

Audience: someone who has generated a synthetic dataset and is deciding whether
to use it, and someone building their own generator who wants the question list.

## The taxonomy

Six categories. The owner proposed A/B-ish, C and D; B was split by *what*
identifies a patient, because the remedies differ completely, and E and F were
added.

### A. Is it a valid dataset at all?

The tier that is easy to skip and embarrassing to fail. Structure before
statistics: a table that is not a legal PMX dataset cannot be assessed for
anything else.

- `validate_pmx()` — schema, column classes, factor levels, event grammar,
  time monotonicity within subject, `CENS`/`DV` coherence, infusion start/stop
  pairing.
- Endpoint survival. `SIM-036` was an entire endpoint silently vanishing —
  108 rows in, 60 out, every `CMT = 3` row gone, no warning. Row counts alone
  would not have caught it; endpoint *sets* would.
- Patient and row counts against the source.

Worth saying in the vignette: `validate_pmx()` on the **source** is also a
check, and on `xgxr::case1_pkpd` it correctly refused a `cens` role because that
study's `CENS` is meaningful only for its PK endpoint (the PD effect is signed).
A validator that only ever passes is not evidence of anything.

### B. Is any single patient singled out?

The privacy tier, and the one worth splitting, because "identifying" is not one
property. A patient can be re-identified through any of these independently, and
fixing one does nothing for the others.

**B1 — timing and structure.** Who was observed when, and who was dosed when.
- `skeleton_uniqueness()` (`n_share_schedule`, `n_share_rarest_time`,
  `n_share_obs_count`, `n_share_dosing`, `nearest_set_diff`, per-endpoint
  breakdown), `plot_pmx_schedule()`.
- The two guarantees, measured on the output: `identifying_visit_sets` and
  `identifying_dose_schedules` in `pmx_settings`. Both must be 0.
- **Read the near-miss distance, not just the count.** Exact-set equality is a
  harsh test: on a real 21-patient study, 15 of 21 patients had a "unique"
  observation schedule and every one of them was a *single* missing sample away
  from somebody else. The count alone cannot distinguish that from a cohort of
  genuinely ad-hoc schedules.
- **Say which endpoint is responsible.** A schedule is only as shared as its
  least-shared part; on that same study the PK sampling accounted for all of it
  and the biomarker for none.

**B2 — extremes.** Who stands out from the cohort.
- `flag_identifiable_subjects()` on follow-up time, dose count, dose magnitude,
  peak DV; `remediate_identifiable_subjects()` to act on it.
- Cautionary note that belongs in the vignette: this screen was flagging every
  patient who stopped early, because its robust scale (the median absolute
  deviation) is **zero** whenever more than half a cohort shares one exact value
  — which on trial data is the ordinary case, not a degenerate one (`SIM-041`).
  A screen's own statistics need checking too.

**B3 — values and geometry.** Whether a synthetic patient's numbers sit too
close to a real patient's.
- `compare_pmx_proximity()`: nearest-neighbour adversarial accuracy against a
  split-half null built from the source itself.
- The null is wide at pharmacometric cohort sizes. "Inside the interval" means
  *nothing was detected*, never *nothing is there*.

**B4 — exact copies.** The crudest check and the one worth keeping.
- `SIM-014`: no generated time vector may be identical to a source subject's.
- Extend to DV vectors and to covariate rows. Cheap, and it is the check that
  catches a whole class of plumbing mistakes.

**B5 — rare covariate combinations. NOT IMPLEMENTED; a genuine gap.**
The only 85-year-old female on arm C is identifiable by her covariates alone,
whatever the schedule machinery did. This is ordinary k-anonymity on the
cross-tabulation of `covariates` and `subject_properties`, `synpmx` does not
check it, and it should. Blending makes it *less* likely than in the source —
covariates are new numbers nobody had — but categorical covariates are
resampled, not blended, so a rare combination can survive intact.

**B6 — deterministic proxies.** A column that is a function of a covariate
discloses that covariate exactly.
- The worked case: under mg/kg dosing an uncorrected `AMT` reveals the anchor's
  weight to the gram. `dose_covariate` (`SIM-040`) is the fix.
- The general rule is worth stating, because every dataset has candidates: BSA
  from height and weight, BMI, creatinine clearance, any derived exposure
  column carried through with `keep`.

### C. Is it still the same study?

The plausibility tier — the owner's "general properties of the dataset are
maintained". These are not privacy checks and failing them does not endanger
anyone; it produces data nobody can develop against.

- **Semantic ordering.** A trough sample must stay a trough. Concretely: the
  sign and rough magnitude of time-after-dose per observation, dose/observation
  ordering within subject, occasion assignment. This is the owner's own example
  and it is the check that would have caught the first two attempts at the
  dose-authoritative grid (see `design/TODO.md`), both of which put a dose
  *before* the sample that preceded it.
- **Regimen validity.** No invented regimen. Dose counts, intervals and
  amounts must be ones the protocol permitted — which is why `synpmx` never
  resamples dose times and only ever truncates a real schedule (`SIM-044`).
- **Endpoint coverage.** All endpoints present, comparable observations per
  patient per endpoint, comparable follow-up.
- **What coverage was lost.** The masking mechanisms buy safety with fidelity,
  and the cost is invisible unless reported: visit sets not reused, arms
  dropped below the donor floor, dose regimens not represented. A cohort of
  nineteen patients on three doses, one on two and one on one came out as
  twenty-one on three, silently, until `pmx_masking_report()` was made to say
  "1 of 3 regimens represented".

### D. Are the numbers still right?

The utility tier — the owner's "distributions of the covariates and of the
dependent variables".

- `compare_pmx_distributions()` per endpoint and per covariate.
- **Between-subject variability shrinks under blending**, necessarily: an
  avatar averaged over `k` donors retains roughly `sum(w^2)` of individual
  variance. `mean_effective_donors` and `cap_binding_fraction` are how much.
  This should be stated as an expected property, not discovered as a surprise.
- **Do the covariate–DV relationships survive?** The one a modeller actually
  cares about, and the least covered. A weight–exposure slope, a dose–response
  relationship, a treatment effect. `vignettes/articles/example-avatar-PKPD-
  covariate-treatment-effect.Rmd` is the existing evidence and should be cited
  rather than duplicated.
- Both failure directions matter. Too close is a privacy failure; too far is a
  utility failure. `compare_pmx_proximity()` is the one diagnostic that reports
  both from a single statistic, and the vignette should use it to make the
  point that most checks have two tails.

### E. Does it work in the workflow?

The purpose of the package is code development outside the environment holding
the real data, so the final check is whether it does that job.

- Does the pipeline that will consume the real data run unchanged on this?
- Does a control stream or model fit *execute* — not that estimates are
  correct, but that the plumbing works?
- Do joins, reshapes and ADaM-shaped derivations behave?

Cheap to state, easy to omit, and the one check that speaks to why anybody is
doing this.

### F. What these checks cannot tell you

A section, not an afterthought. The honest limits:

- Nothing here **bounds what an adversary learns**. These reduce the ways a
  real patient can be singled out; they do not limit how often that succeeds.
  That is what the differentially private modes are for.
- The guarantees are about **reproduction, not similarity**. An avatar whose
  visits sit *near* a real patient's is not covered by
  `identifying_visit_sets`.
- Most of these read the source and are therefore **restricted output**. The
  vignette must be explicit about which can leave the safe environment. Nearly
  all are marked `restricted_not_releasable`; `validate_pmx()` on the synthetic
  table and `flag_identifiable_subjects()` on the synthetic table are the ones
  that do not read the source.

## The section that generalizes: check the output, not the algorithm

If the vignette teaches one thing to somebody building their own method, this
is it, and 2026-08-03 supplies the evidence.

Every leak found that day sat behind a mechanism that was correct on its own
terms:

- `.place_attendance()` rejected identifying placements exactly as designed —
  and the caller then copied the anchor's own set when it returned nothing.
- The visit-set pool was built exactly as designed — and a schedule group with
  one member never got one, so there was nothing to draw from.
- `.apply_attendance()` returned the untouched skeleton when a drawn set could
  not be built on its anchor, so the caller believed a set had been applied
  while the avatar quietly kept its own.

Each was locally right and jointly wrong. The fix that generalizes is not any of
the individual patches: it is asking the question of the **finished table** —
*does any synthetic patient hold a visit set that fewer than `min_pattern_share`
real patients share?* — and recording the answer as a number that must be zero.

Two corollaries, both learned the hard way and both worth writing down:

1. **A check written from the mechanism's own vocabulary inherits its blind
   spots.** The first end-to-end check compared recorded times directly and
   reported "nothing identifying" on `nimoData`, where everything is — because
   generation resamples the deviations coarsening removed, so an avatar in the
   same grid cell never holds the same number. The second snapped back to the
   grid and produced a *false alarm* on `wbcSim`, because those deviations are
   comparable to grid spacing on an irregular study. Only the third — recording
   per avatar what was actually applied — was both exact and quiet.
2. **Passing on every available dataset is not sufficient evidence.** Attempt 2
   at the dose-authoritative grid measured *strictly better* on all five public
   datasets and still collapsed three real visits into one on a dense-PK
   fixture. Public data is a floor, not a proof.

## Inventory: what exists today

| Category | Function or field | Reads source? |
|---|---|---|
| A | `validate_pmx()` | no |
| A | endpoint set comparison | yes |
| B1 | `skeleton_uniqueness()`, `plot_pmx_schedule()` | yes |
| B1 | `identifying_visit_sets`, `identifying_dose_schedules` | (recorded at generation) |
| B2 | `flag_identifiable_subjects()`, `remediate_identifiable_subjects()` | no |
| B3 | `compare_pmx_proximity()` | yes |
| B4 | `SIM-014` gate — in tests only, no exported helper | yes |
| B5 | **nothing** | — |
| B6 | `dose_basis` / `dose_basis_note` in `pmx_masking_report()` | (recorded) |
| C | `pmx_masking_report()`, `compare_pmx()` | yes |
| C | semantic ordering — **nothing exported** | — |
| D | `compare_pmx_distributions()`, `mean_effective_donors`, `cap_binding_fraction` | yes |
| E | **nothing** — by nature the user's own pipeline | — |

## What to build alongside the vignette

Writing it will be the fastest way to find what is missing. From the inventory,
three gaps are already visible:

1. **B5, rare covariate combinations.** The clearest gap and probably the most
   valuable single addition. A cross-tabulation of `covariates` and
   `subject_properties` with a minimum cell count, source and synthetic side by
   side.
2. **C, semantic ordering.** An exported check that time-after-dose sign and
   dose/observation ordering are preserved per patient. It is the check the
   dose-grid work needed and did not have.
3. **B4 as an exported helper** rather than a test-only gate, so a user can run
   it on their own output.

Do not build these before the vignette. Write the vignette against what exists,
let the gaps be visible in the prose, and add them where the prose is
embarrassing to write.

## Suggested shape of the vignette

One worked dataset the whole way through — `xgxr::case1_pkpd` is the best
candidate, being the only public dataset shaped like a real study report (180
patients, declared nominal time, six arms, two endpoints, baseline weight) — and
one deliberately awkward one where checks *fail*, since a document in which
everything passes teaches nothing. `nlmixr2data::pheno_sd` is the honest failing
example: 59 real patients whose individualised neonatal dosing genuinely cannot
be masked, and the run says so.

Order the sections A through F. Lead each with the question in plain language,
then the check, then how to read a bad answer. Keep the "check the output, not
the algorithm" section near the end, where it reads as the lesson rather than
as a preamble.
