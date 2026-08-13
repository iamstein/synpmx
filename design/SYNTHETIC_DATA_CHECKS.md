# Checks of the synthetic data

**Status: written 2026-08-04 as `vignettes/scorecard-synthetic-data-checks.Rmd`.** The
taxonomy below is what shipped, sections A--F, worked on `xgxr::case1_pkpd` with
`nlmixr2data::pheno_sd` as the failing case, as suggested at the end of this
document. What writing it exposed is recorded under "Findings from writing the
vignette" below; three claims in this specification were wrong and are corrected
in place. Everything else below is the reasoning behind the taxonomy and remains
the internal record. Raised by the owner 2026-08-03, after a day of finding leaks
one at a time on a real study:

> There should be a vignette of all the checks we do of the synthetic data.
> [...] This is worth being clear about because this could help
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
loudly. This makes it the fourth shipped vignette, which is a real recurring
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
- **`strata`** (arm, dose group) are copied verbatim from the
  anchor — exact by design, since they are protocol facts.

*The first experiment was wrong about the risk, in an instructive way.* A
fixture with a singleton category (one patient of forty coded `RACE = "Other"`)
and an isolated numeric pair (two patients at 138 and 141 kg among 60--80)
produced **zero** leakage into 200 avatars: no avatar carried the singleton
category, and no avatar's weight came within 5 kg of the isolated pair. The
reason is worth understanding, because it inverts the intuition. Donors are the
*nearest* patients in profile space, and a patient who is unusual on their
covariates is nobody's nearest neighbour — so they are rarely selected as a
donor, and a patient is always excluded from their own donor set.

**Corrected 2026-08-04 while writing the vignette.** That reading was right
about the numeric axis and wrong about the categorical one, because the fixture
varied two things at once: case (a) had **one** holder who was also an outlier,
and case (b) had **two** holders who were typical. Varying the holder count
alone, with everything else held fixed, gives the whole picture:

| holders of the level | of 40 | avatars carrying it, of 200 |
|---|---|---|
| 1 | 2.5% | 0 (0.0%) |
| 2 | 5.0% | 3 (1.5%) |
| 3 | 7.5% | 4 (2.0%) |
| 5 | 12.5% | 17 (8.5%) |
| 10 | 25.0% | 41 (20.5%) |
| 20 | 50.0% | 97 (48.5%) |

**The operative variable is the number of holders, not outlier status.** Making
the two holders numeric outliers as well changes nothing on the categorical
axis: 2 of 200 either way. The mechanism is that `.build_profiles()` one-hot
encodes each categorical level as its own feature, so a level's sole holder sits
alone on that axis and is nobody's nearest neighbour; a second holder gives each
of them a near neighbour, and the level propagates between them and out. By ten
holders the output tracks the source frequency closely.

So "being an outlier is self-protecting" holds for **numeric** covariates, which
are blended, and the numeric half of the original experiment reproduces exactly
(synthetic weight tops out at 116 kg against a source maximum of 141). It does
not hold for categories, where two holders is already a leak and nothing
enforces the protection at one.

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

*The checks, enumerated.* This is the actionable list, in the order a reader
should run them. None exists yet; each is a check the vignette should present
and the package should eventually provide.

1. **Level census, source against synthetic.** For every categorical covariate,
   tabulate the levels on both sides. Report any level that reaches the output,
   together with **how many source patients held it**. A level held by one or
   two real patients is the whole risk, and it is invisible in a marginal
   distribution comparison that only checks the proportions look similar.
2. **Minimum holders, as a threshold.** Flag any level appearing in the output
   that fewer than `min_pattern_share` source patients held. This is the direct
   analogue of what already protects visit sets, and the number should be the
   same one, for the same reason.
3. **Cross-tabulation with `strata`.** Arm x category, then arm x
   category x any other categorical. `strata` are copied verbatim
   from the anchor, so the arm is exact and any joint cell involving it is as
   rare as its rarest part. Report cells below a minimum count.
4. **Numeric covariates: is anyone in a cluster of their own?** Not a
   uniqueness question but a neighbourhood one — the nearest-neighbour distance
   in covariate space, source against synthetic. `compare_pmx_proximity()`
   answers this over the *whole* profile and would mask a covariate-only
   effect; a covariate-only version is worth having.
5. **Levels the user declares rare in the population.** The checks above can
   only see the dataset. Rarity *in the world* — the actual risk for a mutation
   — is knowledge the package does not have and cannot infer. It has to be
   declared, and nothing in `pmx_roles()` accepts such a declaration today.
   That is the gap worth closing first, because it is the only one where the
   right answer may be to **suppress or coarsen the level** rather than report
   it: collapse `TP53-R248W` to `TP53 mutation`, or to `mutation present`, or
   drop the column.

Checks 1--4 are mechanical and could ship as one function. Check 5 needs an API
decision — a `rare_levels` argument on `pmx_covariate()`, or a sensitivity flag
per covariate — and is the one that changes what the generator does rather than
only what it reports.

*What it still cannot do.* Enumerate all combinations. With `d` covariates
there are `2^d` subsets, checks 1--3 cover the ones involving `strata`
and single levels, and an attacker's chosen combination need not be any of them.
The vignette must say this plainly rather than implying the list is complete —
and it is the point at which the honest answer is the differentially private
mode, for the reasons under "Does differential privacy cover this?" below.

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
60--80) and a singleton `RACE` level; nothing propagates — but see the
correction above: what makes (a) safe is that the level has **one** holder, not
that the holder is an outlier. Set `n_holders` and vary only that; the vignette
ships this as a runnable chunk.

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

### D. How close are the parameter distributions to the original values

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


###  What these checks cannot tell you

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
| B2 | `flag_identifiable_subjects()` (scored within `strata`), `remediate_identifiable_subjects()` | no |
| B3 | `compare_pmx_proximity()` | yes |
| B4 | `SIM-014` gate — in tests only, no exported helper | yes |
| B5 | **nothing.** Mechanism pinned by `tests/testthat/test-avatar-relationships.R` | — |
| B6 | `dose_basis` / `dose_basis_note` in `pmx_masking_report()` | (recorded) |
| C | `pmx_masking_report()`, `compare_pmx()` | yes |
| C | `validate_pmx()` `tad_agreement` — TAD against the derivation | yes (source only) |
| C | dose/observation ordering, occasion assignment — **nothing** | — |
| A/C | dose-count fidelity — **nothing**; see finding 2 | — |
| D | `compare_pmx_distributions()`, `mean_effective_donors`, `cap_binding_fraction` | yes |

## Findings from writing the vignette (2026-08-04)

Writing it was the fastest way to find what is missing, as predicted. Four
things the specification did not have:

1. **`flag_identifiable_subjects()` was not stratified by `strata`, and
   misread dose-ranging designs. FIXED 2026-08-04.** On `xgxr::case1_pkpd` it flags **59 of
   180** avatars, 31 of them on "dose magnitude" and 24 on "number of doses".
   Those are not 59 outliers: the study is six arms from 3 mg to 300 mg, so dose
   magnitude is spread by protocol and a cohort-wide robust screen flags the ends
   of the range. This is the same class of error as the median-absolute-deviation
   failure under B2 — the screen's own statistics need checking — and the fix is
   the same idea as `min_pattern_share`: score within stratum. Now implemented:
   scores are computed inside each declared stratum, strata under five subjects
   fall back to the cohort, and `case1_pkpd` goes from 59 flagged to 1 while a
   patient given twice their arm's dose is still caught
   (`tests/testthat/test-strata-balance.R`).

2. **Dose-count fidelity is not reported anywhere.** On `nlmixr2data::pheno_sd`
   the generated cohort keeps its observations (155 rows to 149, 2.6 per patient
   to 2.5) and loses **half its dose rows** (589 to 284). The median real infant
   received twelve doses; the median avatar receives one. Nothing warns. The
   cause is legitimate — dose schedules there are nearly all unique, so
   truncating each to a shared depth collapses most of them to the shortest —
   and `pmx_masking_report()` does say how many regimens are represented, but no
   check states the consequence in the unit a user cares about. Since `SIM-048`
   the masking itself succeeds -- `identifying_dose_schedules` is 0 -- and the
   dosing cost is *worse*, 3 of 56 regimens and 1.2 doses a patient, because
   reaching 0 is exactly what truncates them. So the privacy tier now reports
   nothing at all on this dataset and A5 is the only place it shows. That is the
   vignette's argument for why this dataset should not be shipped, and the
   reason a row that must be 0 is never read on its own.

3. **B5's mechanism was misdiagnosed.** Corrected in place above: the operative
   variable is the number of patients holding a level, not whether they are
   outliers. One holder never propagates, two holders do. The measured curve is
   in the table under B5 and runs in the vignette.

4. **Arm balance was not preserved. FIXED 2026-08-04.** Anchors are sampled
   with replacement, so `case1_pkpd`'s exact 30-per-arm design came back between
   21 and 39 per arm across seeds. Cohort size and arm membership were both
   exact; only the balance moved, and a stratified analysis pipeline run against
   the output would notice. `synpmx_avatar(preserve_strata_balance = TRUE)` is
   now the default. **The privacy caveat is why there is a floor:** exactly
   reproducing a stratum's size discloses that size, so strata under three
   source patients are left stochastic. That is the same reasoning as
   `min_pattern_share`, applied to cell counts instead of visit sets.

Not a finding but worth recording: category A earns its place. Splitting the row
count by event type is what surfaced finding 2, and a total row count hides it
completely.

## What to build alongside the vignette

Writing it will be the fastest way to find what is missing. From the inventory,
three gaps are already visible:

1. **B5, rare categories and combinations.** The clearest gap and probably the
   most valuable single addition. Five checks are enumerated under B5; the
   first four are mechanical and could ship as one function, and the fifth --
   letting a user declare which levels are rare *in the population*, which the
   package cannot infer — needs an API decision and is the only one that
   changes what the generator does rather than what it reports. Two tests in
   `tests/testthat/test-avatar-relationships.R` already pin the mechanism any
   of this has to work against: a categorical covariate is copied from a donor
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

One worked dataset the whole way through — `xgxr::case1_pkpd` is the best
candidate, being the only public dataset shaped like a real study report (180
patients, declared nominal time, six arms, two endpoints, baseline weight) — and
one deliberately awkward one where checks *fail*, since a document in which
everything passes teaches nothing. `nlmixr2data::pheno_sd` is the honest failing
example: 59 real patients whose individualised neonatal dosing can only be masked
by truncating it away, so B1 passes and A5 is where the study is lost.

Order the sections A through F. Lead each with the question in plain language,
then the check, then how to read a bad answer. Keep the "check the output, not
the algorithm" section near the end, where it reads as the lesson rather than
as a preamble.

---

# Proposed additions, from the literature (2026-08-04)

Written after reading the literature review article — since split into
`vignettes/articles/synthetic-data-generation-review.Rmd` and
`vignettes/articles/synthetic-data-checking-review.Rmd` — against the
published evaluation frameworks it cites, plus the general synthetic-data
evaluation literature it does not. **Nothing below is implemented and nothing
below is agreed** — this is a proposal for the owner to cut down. Sources are
listed at the end of this section.

The organizing observation: the vignette's taxonomy is a *disclosure-route*
taxonomy, which is unusual and good — the literature almost always organizes by
fidelity / utility / privacy, and route-based is more actionable. What the
route-based framing loses is the standard machinery the literature has already
built for measuring each route, and four of the gaps below already have
published definitions and reference implementations that we could adopt rather
than invent.

## P0. The structural gap: there is no control group

**This is the single largest gap and it changes how every B-tier check is
read.**

Every privacy check in the vignette compares source against synthetic where the
whole source was available as donors. The literature has converged, from four
independent directions, on the same correction: you cannot distinguish *what an
attacker learns about this individual* from *what an attacker learns about the
population* without a control group the generator never saw.

- Anonymeter (Giomi et al., PoPETS 2023) makes the control dataset the center of
  the framework: risk is reported only when the attack succeeds better against
  training data than against control data.
- Yale et al.'s nearest-neighbour adversarial accuracy — which is the statistic
  `compare_pmx_proximity()` implements — is defined in its original form as the
  *difference* between accuracy computed against the training set and against a
  holdout. We compute only the training half, against a split-half null of the
  source.
- Raab et al. (2024) correct every disclosure measure by what is disclosive in
  the original data anyway (`DiO`), for exactly this reason.
- Ganev & De Cristofaro (IEEE S&P 2025, "the DCR delusion") show that
  similarity-based metrics without a holdout can be passed by datasets that
  leak, because those metrics are average-case and an attacker is worst-case.

**Proposal.** A documented holdout workflow, and eventually a `holdout =`
argument: partition the source into a donor set and a holdout set, generate from
donors only, then run the proximity and uniqueness checks *twice* — synthetic
against donors, synthetic against holdout. The privacy signal is the difference,
not either number.

This is also the honest fix for the vignette's own sentence *"inside the
interval means nothing was detected, never nothing is there."* With a holdout,
"inside the interval" acquires a meaning it does not currently have.

Costs to state plainly, because they are real: a holdout shrinks the donor pool,
which at n = 20 is expensive and at n = 12 is not possible. Recommend it as a
per-study validation exercise run once, not as a per-run default.

## P1. AVATAR's own published privacy metrics, which we do not compute

`synpmx_avatar()` is in the AVATAR family and the AVATAR paper (Guillaudeux et
al., npj Digital Medicine 2023) defines two record-level privacy metrics
specific to patient-centric generators. Both fit `synpmx` exactly, and both are
things the current single-number `adversarial_accuracy` cannot do: they are
**per source subject**, so they name *which* patients are exposed.

- **Local cloaking.** For each source subject, the number of avatars that are
  closer to it than its own avatar is. The paper reports medians of 11 and 24.
  The actionable statistic for us is not the median but the **count of subjects
  with local cloaking 0** — their own avatar is their nearest avatar.
- **Hidden rate.** The percentage of source subjects whose own avatar is *not*
  their nearest avatar. The paper reports 93–94%.
- **Nearest-neighbour distance ratio (NNDR)**, `d1/d2`, which the same paper
  reports alongside distance to closest record with a >= 0.8 rule of thumb.
  `compare_pmx_proximity()` already computes nearest-neighbour distances, so
  this is nearly free.

**Prerequisite, and it is the interesting part.** Local cloaking and hidden rate
require the anchor -> avatar correspondence. `synpmx` has it internally but does
not retain it: no `anchor_id` appears in `pmx_settings` or on the output.
Retaining it is a small change and a large hazard — **the anchor map is the most
disclosive artifact in the whole pipeline**, so it must be computed inside the
restricted environment, marked `restricted_not_releasable` more emphatically
than anything else in the package, and never shipped with the data. Worth
considering a design where the map is passed to the check function rather than
attached to the returned table, so it cannot travel by accident.

Because these are per-subject, they are also the natural input to
`remediate_identifiable_subjects()`, which today acts only on
`flag_identifiable_subjects()` output.

## P2. B5 has a published answer already: RepU and DiSCO

The vignette says of B5 *"nothing in the package checks this, and it is the most
valuable thing missing"*, then enumerates five checks from first principles.
Four of the five are a rediscovery of measures Raab, Nowok and Dibben formalized
in 2024 and shipped in `synthpop` >= 1.8.1 as `disclosure()` and
`multi.disclosure()`. Adopting their vocabulary would give us definitions to
cite, thresholds others have argued about, and a reference implementation to
test against.

- **RepU (replicated uniques)** — records unique in the original on a chosen key
  set that are *also* unique in the synthetic. This is B5 checks 1–2 and B4,
  generalized from visit sets to arbitrary quasi-identifier sets.
- **DiSCO (Disclosive in Synthetic, Correct in Original)** — the percentage of
  original records for which a key combination has a unique target value in the
  synthetic *and* that value matches the original. This is the attribute
  disclosure measure, and it is what "arm x sex x age band -> mutation status"
  risk actually is. It is B5 check 3, made into a number.
- **DiO baseline correction** — subtract what is disclosive in the original data
  anyway. A relationship that is deterministic in the population is not a leak,
  it is the science. **This is the same insight as the B2 stratification
  finding** — "a screen that ignores assignment reports the protocol back to you
  as a privacy finding" — and `DiSCO - DiO` is that insight as arithmetic. Worth
  saying so in the vignette, because it is the same idea arrived at twice from
  different directions.

The fit is good for a reason worth stating: at *subject baseline* level a PMX
dataset **is** rectangular, one row per patient, which is exactly the shape
`synthpop`'s machinery assumes. The event table is where their tools break (the
literature review says so); the covariate table is where they apply directly.

Two concrete sub-findings while checking this:

1. **`compare_pmx_distributions()$covariates_categorical` already is the level
   census.** It is computed per subject via `.subject_baseline_values()` and
   reports level counts and proportions for source and synthetic side by side.
   What is missing is not the table but the **flag** — no minimum-holder
   threshold, and the two datasets are stacked long rather than joined. The
   vignette's claim that nothing exists is too strong; the accurate claim is
   that nothing *flags*.
2. **`strata` are excluded from that table.** `compare_pmx_distributions()`
   loops over `roles$covariates` only, and `strata` are a separate role. So the
   one categorical axis that is **copied verbatim from the anchor** — the axis
   where a two-patient cell is reproduced exactly — is the axis missing from the
   census. That is why the vignette had to hand-write `level_census()` for
   `TRTACT`. Including `strata` in `compare_pmx_distributions()` is a small fix
   with a real finding behind it.
3. The vignette's `level_census()` chunk builds its level set from the
   **synthetic** side only (`levels_out`), so a source level that was *dropped*
   is invisible. That is a coverage failure as well as a privacy one, and the
   join should be a full outer join. Small, but it is a check that cannot see
   one of its two failure directions.

## P3. Linkability is missing entirely

Anonymeter operationalizes the Article 29 Working Party's three criteria for
factual anonymization: **singling out, linkability, inference**. The vignette's
B1–B6 cover singling out thoroughly and inference partly (B5, B6), and cover
linkability **not at all**.

Linkability is: an attacker holding *part* of a record from another source — the
covariates, say, obtained from a registry — links it to the right synthetic
record. For a record-based generator with an anchor map this is directly
testable and it is the attack the method is most exposed to.

**Proposal:** a mapping paragraph in B's preamble giving B1–B6 against the three
regulatory criteria, which (a) lets a reader arriving from a data-protection
officer navigate, and (b) makes the missing one visible as missing rather than
absent by omission.

## P4. Category D is thin, and PMX is where we should be ahead

Every pharmacometric benchmark in the literature — Destere (daptomycin popPK),
Woillard (pharmacogenetics), Jiang (PK/PD deep generative) — evaluates on the
*analysis*, not on the marginals, because the analysis is what a user of the
data does. Category D is currently two marginal-summary tables and a pointer to
an article. Four additions, in ascending cost:

1. **Non-compartmental analysis (NCA) summary, per subject, source against
   synthetic.** Cmax, Tmax, AUC(0-last), terminal half-life, and the
   between-subject CV of each. Model-free, trapezoid rule, roughly twenty lines,
   no new dependency. This is the pharmacometric idiom for "did the exposure
   survive", and it renders the variance-shrinkage discussion in the units a
   pharmacometrician thinks in rather than as a standard deviation on pooled DV.
   **Probably the highest value-per-line addition in this whole document.**
2. **Confidence interval overlap** (Karr et al. 2006) as the named metric for
   the parameter-recovery check. `example-avatar-PKPD-covariate-treatment-
   effect.Rmd` already does the work; D should name the statistic rather than
   gesturing at "the relationship survives". Destere's benchmark is exactly
   this: fit the same structural model to both, compare fixed effects, omega,
   sigma.
3. **Within-subject trajectory shape.** Marginals can match perfectly while
   every individual profile is scrambled — which is precisely the failure the
   literature review attributes to column-wise synthesizers, and which nothing
   in the vignette would catch. The time-series literature checks this with
   autocorrelation-function agreement and a discriminative score. Cheapest
   useful version for us: compare **within-subject** spread around each
   subject's own smooth against **between-subject** spread. Blending should
   shrink between-subject variance and leave within-subject residual structure
   alone; if it shrinks both, the trajectories are being flattened and the D
   tables cannot see it.
4. **pMSE / propensity-score utility** (Snoke & Raab; the measure Woillard uses
   for this purpose). Fit a CART classifier to distinguish source subjects from
   synthetic ones; pMSE at its null means indistinguishable. It is a general
   utility measure with a known null distribution — the same shape as
   `compare_pmx_proximity()` — and it complements the specific checks. Note it
   reads in both directions like proximity does: too distinguishable is a
   utility failure, indistinguishable *plus* adversarial accuracy near zero is
   memorization.

## P5. Coverage, not only closeness

Alaa et al. (ICML 2022) split fidelity into three axes, and we measure roughly
one of them:

- **alpha-precision** — are synthetic records typical of real ones. Roughly what
  `compare_pmx_distributions()` reports.
- **beta-recall** — are real records *covered* by synthetic ones. **Not
  measured.** This is the "did we lose the tails" direction, which today is
  visible only as a standard deviation shrinking in a table with no pass
  criterion attached.
- **authenticity** — the fraction of synthetic records that are not
  near-copies of a training record. This is the "too close" tail that
  `compare_pmx_proximity()` covers with one aggregate number; authenticity is
  per record, and per record is what remediation needs.

Cheap proxies worth having even without the full estimators: **range coverage**
per covariate and per endpoint (does the synthetic span the source's range), and
**category coverage** (levels present in the source but absent from the output).
Both are one line each and both have a stateable pass criterion.

## P6. Smaller PMX-specific checks that nothing covers

6a. **Dropout / follow-up curve.** A Kaplan-Meier of time to last observation,
source against synthetic. The vignette treats follow-up length only as a B2
*outlier axis* — whether one patient stands out — and never asks whether the
dropout **process** survived. This is close to indefensible given what the
package advertises: the literature review's own pitch is "the patient who
withdrew at week 12", and nothing checks that withdrawal behaves like the
source's. Category C.

6b. **Missingness and BLOQ pattern.** Proportion of below-limit-of-quantification
observations per endpoint per arm, source against synthetic, and the
missing-visit pattern. `cens` is already a declared role so this is nearly free.
`pmx_masking_report()` reports visit-set reuse but not whether the censoring
*rate* was preserved — and censoring rate is a first-class property of a PK
dataset that any assembly script has to handle.

6c. **Value plausibility bounds.** `validate_pmx()` checks `dv_finite`, not
plausibility. DV blending adds AR(1) residual noise after the transform, so a
generated value can land outside anything the source contained. Numeric
covariates are already floored above zero when all donors are positive
(`.synthesize_covariates()`), so the exposure is DV and any covariate whose
donor set is mixed-sign. An A-tier range check — synthetic min/max per endpoint
and covariate against the source's, or against declared `pmx_bounds` — has a
one-line pass criterion and belongs with the other structural checks.

6d. **Dose/observation ordering and occasion assignment.** Already on the
vignette's own gap list; noted here only so the C-tier list is complete in one
place.

## P7. One addition to section E

Section E says "nothing here bounds what an adversary learns", which is correct
but under-specified. The literature now supplies the precise failure mode: these
are **average-case** statistics and an attacker is **worst-case**, and datasets
that pass similarity-based metrics have been shown to leak under membership
inference (Ganev & De Cristofaro 2025). Naming that strengthens the DP argument
B5 already makes, and it costs two sentences.

---

# Proposed reorganization of `vignettes/scorecard-synthetic-data-checks.Rmd`

The short verdict: **the taxonomy is right and should not be restructured.**
A–F by disclosure route is better than the literature's fidelity/utility/privacy
split for someone deciding whether to ship a dataset. Seven changes, roughly in
order of value.

1. **Add a scorecard table near the top.** One row per check: category,
   function, reads source?, **pass criterion**, status (shipped / gap). Today
   the pass criteria are scattered through prose — "both must be 0", "zero, as
   required", "0.5 is the target", "must be zero" — and a reader deciding
   whether to ship has to reconstruct the list by reading the whole document.
   This is the highest-value change and it is purely additive.

2. **Fold the two gap lists into one.** B5's numbered 1–5 and the closing "What
   is still missing" table overlap. If the scorecard has a status column, the
   closing table becomes a filter of it rather than a second list to maintain.

3. **Move the B5 propagation experiment out, to `avatar-algorithm.Rmd`.** It is
   the largest code block in the vignette — a 50-line fixture plus a
   propagation loop — and it teaches the *mechanism of this generator*, not a
   check the reader should run on their own data. Keep in the checks vignette:
   the conclusion (one holder never propagates, two do), the census check, and
   the DP argument. This is the single biggest length reduction available and it
   sharpens the section from "a study of our own generator" into "here is the
   check, here is the gap".

4. **Make the release status a scorecard column, not a section-E paragraph.**
   It is per-check metadata and it is what determines whether a diagnostic can
   leave the environment — which is an operational question a reader asks per
   check, not once at the end.

5. **Name the fidelity / utility / privacy frame once**, in the preamble, and
   map A–D onto it in one sentence. Not a reorganization: a signpost, so a
   reader arriving from Anonymeter, `synthpop`, or SDMetrics can find the
   corresponding section. Costs a paragraph, buys the entire external audience
   the article is partly written for.

6. **Add the singling out / linkability / inference mapping to B's preamble**,
   per P3, so the regulatory vocabulary connects and the missing route is
   visible.

7. **Grow D.** It is currently the thinnest section relative to its importance —
   two tables and a pointer — and it is where P4's NCA comparison and
   confidence-interval overlap belong.

What should **not** change: sections E and F. F is the best writing in the
package and its position at the end is correct.

## Sources for this section

- Giomi M, Boenisch F, Wehmeyer C, Tasnádi B. *A Unified Framework for
  Quantifying Privacy Risk in Synthetic Data.* PoPETS 2023(2). (Anonymeter;
  singling out / linkability / inference; the control dataset.)
- Guillaudeux M, Rousseau O, Petot J, et al. *Patient-centric synthetic data
  generation.* npj Digital Medicine 2023;6. (Local cloaking, hidden rate, NNDR.)
- Raab GM, Nowok B, Dibben C. *Privacy risk from synthetic data: practical
  proposals.* arXiv:2409.04257, 2024, and *Practical privacy metrics for
  synthetic data*, arXiv:2406.16826. (RepU, DiSCO, DiO; `synthpop::disclosure()`.)
- Ganev G, De Cristofaro E. *The Inadequacy of Similarity-Based Privacy
  Metrics.* IEEE S&P 2025. ("The DCR delusion".)
- Alaa A, van Breugel B, Saveliev E, van der Schaar M. *How Faithful is your
  Synthetic Data?* ICML 2022. (alpha-precision, beta-recall, authenticity.)
- Yale A, Dash S, Dutta R, et al. *Generation and evaluation of privacy
  preserving synthetic health data.* Neurocomputing 2020. (Nearest-neighbour
  adversarial accuracy, train-versus-holdout form.)
- Snoke J, Raab GM, Nowok B, Dibben C, Slavkovic A. *General and specific
  utility measures for synthetic data.* JRSS-A 2018. (pMSE.)
- Karr AF, Kohnen CN, Oganian A, Reiter JP, Sanil AP. *A framework for
  evaluating the utility of data altered to protect confidentiality.* The
  American Statistician 2006. (Confidence interval overlap.)
- Destere A, Lombardi R, Labriffe M, et al. medRxiv 2026, and Woillard JB,
  Benoist C, et al. CPT:PSP 2025 — both already cited in the literature review.

---

## Owner's decisions on the proposal (2026-08-05)

Reviewed by the owner. What was accepted, and — more usefully — the scope
statement that came out of the review.

**The scope statement, which should govern everything above.** The use case is
that the data still *looks like* pharmacometric data and still protects privacy.
It is explicitly **not** required that distributions and processes are
maintained exactly, and it is **not** about scientific discovery. That decides
several of the proposals without further argument:

- **P4 specific-utility measures are out of scope by design.** Confidence
  interval overlap, parameter recovery, and non-compartmental exposure agreement
  (Cmax, Tmax, AUC) are not pursued. Not because they are unimportant, but
  because the package does not claim the property they measure. Stating that is
  more useful than reporting them badly.
- **P4's pMSE and P5's coverage measures are interesting, not binding.** Worth
  knowing about; not worth chasing a number for. Coverage is worth reporting
  only because it lets a shrinking spread be *attributed* — masking working
  versus the generator collapsing to the mean look identical in a standard
  deviation.
- **P6a, the dropout curve, is dropped.** The owner does not require that the
  dropout process survive exactly.

The reasoning is written up for readers in `vignettes/articles/literature-
review.Rmd`, under "why almost none of this is `synpmx`'s problem".

**Accepted and implemented 2026-08-05** (see `design/TODO.md` for the entry):

1. The scorecard, in both forms — static index at the top of the checks
   vignette, runnable version at the end of `demo.Rmd`.
2. B5's propagation experiment moved to `avatar-algorithm.Rmd` step 8.
3. Release status as a scorecard column, with `validate_pmx()` as rows A1/A2.

Plus, unrequested but required by the accuracy rule in `AGENTS.md`: the two
corrections under P2 (the census exists but does not flag; `strata` are missing
from it).

**Still open, unagreed**, in the order the owner is most likely to want them:
the holdout (P0), local cloaking and hidden rate (P1), `strata` in
`compare_pmx_distributions()`, B5 as one function in RepU/DiSCO vocabulary (P2),
and linkability (P3).

**The tutorial.** The owner asked for the checking literature to be written up
as a tutorial rather than a citation list, and for it to live in the literature
review rather than here. That is now
`vignettes/articles/synthetic-data-checking-review.Rmd` — originally the second
half of `literature-review.Rmd`, split into its own article on 2026-08-11.
Everything above stays here as the internal record.
