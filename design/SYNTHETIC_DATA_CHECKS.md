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

**B5 — rare categories and rare combinations. NOT IMPLEMENTED; a genuine gap,
and not the one it first looks like.** Measured 2026-08-03 rather than assumed,
after the owner asked why blending does not already handle it.

*The mechanism.* `.synthesize_covariates()` treats the two kinds of covariate
completely differently, and only one of them is blended:

- **Numeric** covariates are a weighted mean of the donors' values plus noise.
  The result is a new number nobody had.
- **Categorical** covariates are `sample()`d from the donors' values. That is
  not blending and there is no averaging available: a synthetic patient's
  category is always **some real patient's actual category**, copied.
- **`subject_properties`** (arm, dose group) are copied verbatim from the
  anchor — exact by design, since they are protocol facts.

*The first experiment was wrong about the risk, in an instructive way.* A
fixture with a singleton category (one patient of forty coded `RACE = "Other"`)
and an isolated numeric pair (two patients at 138 and 141 kg among 60--80)
produced **zero** leakage into 200 avatars: no avatar carried the singleton
category, and no avatar's weight came within 5 kg of the isolated pair. The
reason is worth understanding, because it inverts the intuition. Donors are the
*nearest* patients in profile space, and a patient who is unusual on their
covariates is nobody's nearest neighbour — so they are rarely selected as a
donor, and a patient is always excluded from their own donor set. **Being an
outlier is self-protecting under this design.**

*The real risk is the opposite shape, and it is the owner's example.* A patient
who is **typical in every way except one rare category** — an ultra-rare
mutation status, say — sits right in the middle of the profile space, is
selected as a donor constantly, and their category is copied out. The same
fixture with two otherwise-ordinary patients of forty carrying
`MUT = "TP53-R248W"` put that value on 2 of 200 avatars. Blending never touched
it, because there is nothing to average.

*Why that matters more than dataset-uniqueness.* The disclosure is not "this
value is unique in the dataset". It is that **the value may be rare in the
world**. If five people alive carry a mutation, a synthetic dataset containing
it discloses that someone with that mutation was in this study, which for a
named trial with public inclusion criteria can be nearly identifying on its own.
No amount of cohort size helps, and `min_pattern_share` — which protects visit
sets — has no analogue here.

*Combinations make it worse and unenumerable.* Any subset of covariates can be a
quasi-identifier: arm x sex x age band x mutation. With `d` covariates there are
`2^d` subsets and no principled way to know in advance which ones single someone
out in the presence of external data. One mitigation is already accidentally
present: each covariate is sampled **independently** in the loop, so `SEX` may
come from one donor and `RACE` from another, and a rare *joint* combination is
less likely to be reproduced whole than any of its parts. That is a fidelity
cost (real correlations between covariates are broken) doing double duty as a
weak privacy benefit, and it should be stated as such rather than claimed as a
mechanism.

*What to build.* A cross-tabulation of `covariates` and `subject_properties`
with a minimum cell count, source and synthetic side by side, plus a per-level
report of the rarest categories reaching the output. It cannot enumerate all
combinations, and the vignette should say so plainly rather than implying the
check is complete.

*Two reproducible fixtures* for whoever builds it, both under 40 lines and both
run on 2026-08-03:

```r
# (a) the self-protecting case: a singleton category on an OUTLIER patient.
#     Expect 0 of 200 avatars to carry it -- outliers are nobody's donor.
# (b) the leaking case: a rare category on TYPICAL patients.
#     Expect it to propagate at roughly its source frequency.
mk <- function(i) {
  t <- c(0, 1, 2, 4, 8, 24)
  data.frame(
    ID = i, TIME = c(0, t), NTIME = c(0, t),
    DV = c(NA, round(10 * exp(-0.1 * t) + rnorm(length(t), 0, .3), 3)),
    AMT = c(100, rep(0, length(t))), EVID = c(1L, rep(0L, length(t))),
    CMT = c(1L, rep(2L, length(t))),
    WT = round(runif(1, 60, 85), 1),                       # (b): all typical
    MUT = if (i %in% c(17, 23)) "TP53-R248W" else "wild-type",
    stringsAsFactors = FALSE)
}
raw <- do.call(rbind, lapply(1:40, mk))
roles <- pmx_roles(id = "ID", time = "TIME", dv = "DV", amt = "AMT",
                   evid = "EVID", cmt = "CMT", nominal_time = "NTIME",
                   covariates = c("WT", "MUT"))
syn <- synpmx_avatar(raw, roles, n_subjects = 200, seed = 9)
# 2 of 40 source patients carry it; 2 of 200 avatars did.
```

For (a), give the rare patient an outlying weight as well (138 and 141 kg among
60--80) and a singleton `RACE` level; nothing propagates.

*Does differential privacy cover this?* **Yes, and it is the strongest argument
for the DP modes.** This is precisely the difference between the two families.
AVATAR's protections are *enumerated* — visit sets, dose schedules, extremes,
proximity — so a risk nobody named is a risk nobody covered, and B5 is exactly
such a risk. An (epsilon, delta) guarantee is *unenumerated*: it bounds the
influence of any one individual on **any** function of the output, so it holds
for combinations nobody thought to check, including ones an attacker constructs
later with data that does not exist yet.

Concretely, `synpmx`'s DP path never copies a category. `pmx_covariate(levels =)`
releases a noisy **count vector** over the levels (L1 sensitivity 1, since one
subject occupies one level), normalises it to probabilities, and generation
draws from that distribution. A level held by two patients out of forty is
swamped by noise at any sensible epsilon, and `.support_threshold()` gates
weakly-supported cells out entirely.

Three honest caveats, all of which belong in the vignette:

1. **The level set must be public.** `pmx_covariate()` requires `levels` and a
   `source` citation for where that public knowledge came from — precisely
   because a domain derived from the data leaks the domain. For mutation status
   that is a real constraint: "the mutations observed in this trial" is itself
   disclosive, so the level set has to be pre-specified from outside.
2. **Delta matters here more than usual.** (epsilon, delta)-DP permits the
   guarantee to fail with probability delta, and for a category held by one or
   two people that is the failure that matters. Keep delta well below 1/n.
3. **DP bounds the mechanism, not the world.** If a trial's inclusion criteria
   are public and only a handful of people could qualify, membership may be
   inferable regardless. DP still holds — that is its point — but a large
   epsilon buys little absolute protection in that setting.

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
  ordering within subject, occasion assignment.

  **`TAD` is the carrier of this check**, and the vignette needs to be clear
  about a thing that surprises people: `tad` is an **output, not an input**.
  `synpmx_avatar()` recomputes it from the generated times and dose rows and
  overwrites whatever the source held. Declaring the role says which column to
  overwrite and to carry through; the source's values never generate anything.

  `validate_pmx()` is the one place they are read. It reports, as a non-fatal
  warning, where the declared column disagrees with time since the most recent
  dose row — and a disagreement is a real finding every time, because it means
  one of: the study measures TAD from the end of an infusion, or from a nominal
  dose time, or from an assigned occasion rather than the most recent dose; or
  the source column is wrong; or our derivation is wrong for that study. It
  cannot tell you which, and the message says so rather than guessing.

  Worth showing in the vignette because it is live in our own registry:
  `nlmixr2data::nimoData` disagrees on **143 of 321 observation rows (45%)**, by
  up to 311.9 hours, and the synthetic column follows the derivation rather than
  nimoData's convention. That is exactly the kind of silent semantic change this
  whole category exists to surface.

  Two limits to state: a sample taken before any dose is reported as TAD 0
  because `validate_pmx()` refuses a negative, not because it is genuinely zero
  hours after a dose it precedes; and where `addl`/`ii` are declared the
  derivation cannot see the doses they imply, so the check is skipped and says
  so.

  **Still missing from category C**: dose/observation ordering and occasion
  assignment have no check. TAD covers the first only indirectly. This is the owner's own example
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
| B5 | **nothing.** Mechanism pinned by `tests/testthat/test-avatar-relationships.R` | — |
| B6 | `dose_basis` / `dose_basis_note` in `pmx_masking_report()` | (recorded) |
| C | `pmx_masking_report()`, `compare_pmx()` | yes |
| C | `validate_pmx()` `tad_agreement` — TAD against the derivation | yes (source only) |
| C | dose/observation ordering, occasion assignment — **nothing** | — |
| D | `compare_pmx_distributions()`, `mean_effective_donors`, `cap_binding_fraction` | yes |
| E | **nothing** — by nature the user's own pipeline | — |

## What to build alongside the vignette

Writing it will be the fastest way to find what is missing. From the inventory,
three gaps are already visible:

1. **B5, rare categories and combinations.** The clearest gap and probably the
   most valuable single addition. A cross-tabulation of `covariates` and
   `subject_properties` with a minimum cell count, source and synthetic side by
   side, plus a per-level report of the rarest categories reaching the output.
   Two tests in `tests/testthat/test-avatar-relationships.R` already pin the
   mechanism it has to change: a categorical covariate is copied from a donor
   and never blended, and an outlying patient is self-protecting only as a side
   effect of donor selection.
2. **C, semantic ordering.** Partly done: `validate_pmx()` now reports where a
   declared `TAD` disagrees with the derivation. What is still missing is the
   ordering itself — that a dose never moves past a sample that preceded it, and
   that occasion assignment survives. That is the check the dose-grid work
   needed and did not have (`design/TODO.md`).
3. **B4 as an exported helper** rather than a test-only gate, so a user can run
   it on their own output.

Do not build these before the vignette. Write the vignette against what exists,
let the gaps be visible in the prose, and add them where the prose is
embarrassing to write.

## Suggested shape of the vignette

**List the datasets.** The vignette should carry a short table of the public
datasets it uses, one line each, saying what each one is for — a reader who
wants to reproduce a check needs to know which dataset shows it and why that
one. `design/TEST_SIM.md` holds the full inventory of what is available,
including the candidates not yet used and the ones deliberately skipped; the
vignette should show only what it actually uses and not duplicate the registry.

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
