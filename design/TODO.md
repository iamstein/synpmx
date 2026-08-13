# TODO

Living task list. One line per actionable item, newest scope at the top.

How this relates to the other documents. `AGENTS.md` explains the three tiers
and which one new writing belongs in.

Internal design record (`design/`, cited by nothing shipped):

- `design/TODO.md` (this file) — **what to do next.** The working queue.
- `design/REVIEW_BACKLOG.md` — **why**, for defects and design findings. `REV-###`.
- `design/TEST_SIM.md` — **evidence**, for simulation defects and their gates. `SIM-###`.
- `design/METHOD_DISCUSSION.md` — **tradeoffs**, AVATAR blending vs formal DP,
  and why AVATAR is the trajectory-level analogue of synadam.
- `design/PRIVACY_EVAL.md` — **program**, the proposed empirical
  disclosure-risk evaluation. Phased, not yet started.
- `design/PROTOTYPE_SPEC.md` — **contract**, the specification being implemented.


## Open 2026-08-13: discrete endpoints come back continuous (`SIM-045`)

Found while rewriting `avatar-evaluation-public-data.Rmd` onto the demo's shape.
`xgxr::mad` carries binary, count and ordinal PD endpoints in the same numeric
`LIDV` column as PK, and all three come out continuous: the 0/1 endpoint takes
600 distinct values from -0.13 to 1.08. Every scorecard row passes, so nothing
catches it.

- [ ] Decide the fix: round generated values back onto the source level set per
      endpoint at emit, or refuse and say so. Class restoration is per column and
      cannot see an endpoint, so neither is free.
- [ ] Regression test on a fixture with one integer-level endpoint alongside a
      continuous one.
- [ ] The `mad` section of the evaluation vignette records it as a known
      limitation; update it when the fix lands.

## Done 2026-08-03: run-report readability, and `SIM-038` behind it

Raised by the owner reading `INTERNAL_STUDY`'s output: the alerts were too
long to read, the exposure table's column names were unintelligible, the
mechanism table's rows could not be interpreted, and there was no way to tell
whether weight-based dosing had been detected.

- [x] Alerts wrapped to a fixed width, split `ALERT` / `NOTE`, emitted once
      (signalled rather than raised, so the text stopped printing twice).
- [x] `skeleton_uniqueness()` columns renamed to one consistent question, given
      `coarsen_time`, and reprinted as a verdict plus two summary tables. It now
      states which side of coarsening it scored, which was the owner's question.
- [x] New `plot_pmx_schedule()`. A count cannot distinguish one-off sampling
      times from ordinary dropout; the picture can.
- [x] New `pmx_masking_report()`, replacing five hand-built copies of the same
      table. Says outright whether dose amounts were recomputed and, when not,
      why not.
- [x] `SIM-038`: `format(x, digits = 12)` is vector-dependent, so identical
      visit sets keyed differently and nearly every pattern looked unique.
      Found by drawing the schedule next to the table and seeing them disagree.

## Done 2026-08-03 (second pass): the 86% fallback, and declared dosing

From the owner reading the alerts on `INTERNAL_STUDY`.

- [x] `SIM-039`: 86% of avatars kept their anchor's own visit set. Deterministic
      `trailing` placements retried 24 times, and staggered discontinuation gave
      every (kind, count) shape a single holder so the group got no pool. Both
      fixed; 21 of 21 -> 0 of 21 on a staggered fixture.
- [x] `SIM-040`: `pmx_roles(dose_covariate = )`, so a weight-based study can be
      declared instead of inferred. Holds each dose row's own ratio, so
      intra-patient escalation survives.
- [x] Alerts printed twice in knitted reports (the condition carried the
      `warning` class, which knitr renders alongside the message).
- [x] Alerts no longer advise declaring `nominal_time` to someone who has.
- [x] Report wording: "floor" named, "stay in the cohort as donors" clarified.

Open on the dose side, from 2026-08-03:

- [x] **Dose truncation as a shape.** Done as `SIM-044`.
- [ ] **Dose times on their own grid — THREE ATTEMPTS, none shippable.** The
      goal is right and the diagnosis is certain: `nimoData` is 12 of 12
      identifying purely because dose times are recorded actuals. Doses are
      weekly but 165.70 / 167.24 / 168.05 are three cells, and 86 of 97 dose
      times after coarsening are held by one patient. A dose-only grid merges
      them perfectly — 104 distinct times collapse to 10 weekly cells — and
      every attempt below fixes `nimoData` to 0. What none of them survives is
      what happens to the observations around a dose that has just moved.

      1. *Independent observation grid.* An observation at 167.16 snaps to
         167.10 while its dose at 167.20 snaps to 166.73, so the dose lands
         before the sample that preceded it. Fires on all 12 subjects; guarded
         per subject it reverts all 12 and is a no-op everywhere.
      2. *Clamp observations between neighbouring dose cells.* Measured
         strictly better on all five public sets (`nimoData` dose 12 to 0,
         `wbcSim` observation-unique 17 to 13, nothing worse) — but it drags
         every POST-dose sample forward onto the dose when the dose moves far,
         collapsing three real visits into one. The existing spurious-collision
         repair then reverts those rows, leaving a mixed state that can itself
         reorder. Passes the public numbers, fails a dense-PK fixture.
      3. *Carry samples with their dose, preserving time-after-dose.* The
         conceptually right one, and what the owner described. Fixes `nimoData`
         and produces identical schedules on a dense-PK fixture. Regresses
         `theo_md` from 0 to 12 unique observation schedules, because the
         elapsed-time axis is unstable at interval boundaries: a trough drawn a
         minute before a dose is 23.9 h after the previous one, the same trough
         drawn a minute after is 0.02 h after this one. Referring each sample to
         the NEAREST dose rather than the preceding one was tried and does not
         fix it.

      Attempt 3 is the one to finish. The remaining problem is entirely the
      boundary: it needs the elapsed-time grid derived per dose interval, or a
      reference rule that cannot straddle. Whatever is tried must be measured on
      all five public sets plus `case1_pkpd` before it ships — attempt 2 shows
      that passing the public numbers is not sufficient evidence.


      What would work is making **doses authoritative**: derive the dose grid
      first, then clamp each observation into the interval between its
      neighbouring snapped dose times, so a pre-dose sample ties with its dose
      rather than crossing it. That keeps dose cells shared, which is the whole
      point. It touches coarsening for every dataset, so it needs its own
      before/after on all five public sets plus the SIM-014 timing gates before
      it can be trusted. Not attempted yet.

Still open from this pass:

- [x] Documentation for `dose_covariate`: `README.md` (and its pinned test), the
      demo's theophylline section, the algorithm vignette's Step 1 and M5, and
      the four-methods vignette. `theo_md` is the public test case.
- [ ] Run the four `scripts_private/` templates against the real studies again
      now that `SIM-038`, `SIM-039` and `SIM-040` are fixed. Check the
      "avatars keeping their anchor's own visit set" row is near 0, and set
      `dose_covariate` on the weight-based studies.
- [ ] Consolidate the four templates into one shared body plus a per-study
      config chunk. Agreed in principle; the owner wants one template in good
      shape first.
- [x] `compare_pmx_distributions()` gained a `knit_print()` method and the four
      local `kable_distributions()` copies are gone.
- [x] Template audit: `N_SUBJECTS` was set in three templates and never passed
      to `synpmx_avatar()`; a relative `OUT_DIR` put source-derived CSV at the
      repository root, outside the `scripts_private/` ignore, when the chunks
      are run interactively (now covered by `/output_*/`); the ignore allowlist
      still named `try_avatar.R`, which no longer exists.

## Why is there a `tad` role at all?

Raised by the owner 2026-08-03: "Why is there a TAD role? Can't that just be
rederived? I guess it could perform a check if it was derived correctly?"

Investigated rather than answered from memory, and the answer is worse than the
question assumes: **the role means three different things in three places.**

- **AVATAR ignores the declared values entirely.** `.recompute_tad()` overwrites
  the column from generated `TIME` and the generated dose rows. Declaring `tad`
  does exactly two things: names the column to overwrite, and keeps it in the
  output. The source's TAD values are never read.
- **`validate_pmx()` barely checks it.** Finite and non-negative on observation
  rows. It does *not* check that TAD agrees with `TIME` minus the most recent
  dose time — which is the only thing that makes it a TAD.
- **The DP path trusts it.** `.subject_clock()` in `representation.R` derives
  TAD from time and dose times and then *overrides* the derived value with the
  declared one wherever it is finite.

So the two engines disagree with each other about whether the column is input or
output, and nothing reconciles them.

**Measured on `nimoData`, which is in our own registry: the declared TAD
disagrees with what `.recompute_tad()` derives on 143 of 321 observation rows
(45%), with a maximum disagreement of 311.94 hours.** Example rows:

| ID | TIME | declared TAD | derived TAD |
|---|---|---|---|
| 3 | 166.27 | 0.00 | 0.57 |
| 3 | 332.35 | 166.08 | 166.65 |
| 3 | 333.10 | 0.00 | 0.70 |

The cause is not diagnosed. The offsets look systematic rather than random, so
plausible explanations are that nimoData measures TAD from the end of an
infusion rather than its start, or from a nominal dose time rather than the
recorded one. Whatever it is, generating from nimoData today silently replaces
one TAD convention with another on nearly half its rows, and no check notices.

Two further problems found while looking:

- **`.recompute_tad()` is wrong for any dataset using `ADDL`/`II`.** It reads
  explicit dose *rows* via `findInterval`, and `synpmx` accepts `addl`/`ii`
  without expanding them. TAD would then be measured from the last written
  dose rather than the last actual one. No dataset under test exercises this;
  `nlmixr2data::nmtest` would (see the inventory in `design/TEST_SIM.md`).
- **Pre-dose samples get TAD 0.** The code comments that zero is "the only
  value `validate_pmx()` accepts and the only one that is not a fiction". A
  baseline sample taken before any dose has *no* time-after-dose, and 0 is a
  fiction; `validate_pmx()` rejecting negatives is what forces it. Worth
  deciding rather than inheriting.

Decide what the role is *for*, then make all three places agree:

- [x] **Option A — check-only. DONE 2026-08-03.** `validate_pmx()` gains a
      `tad_agreement` check reporting, non-fatally, where the declared column
      disagrees with time since the most recent dose row; `pmx_roles()`
      documents that `tad` is an output rather than an input; and one
      `.derived_tad()` now serves both the check and `.recompute_tad()`, so the
      two cannot drift. Fixing this also made warnings on a *valid* report
      visible — `print()` returned on the valid branch before reaching them, so
      any non-fatal finding on a well-formed dataset had been unreachable.
Left open, and deliberately not done for now:

- [ ] **Option B — authoritative.** The declared column is the truth and the
      generator carries it through the same transformation as `TIME`. Harder,
      and it makes the DP path right and the AVATAR path wrong today.
- [ ] **Option C — drop the role.** Derive TAD always, name the column with
      `keep` if you want it carried. Simplest, and it loses the check.

The `ADDL`/`II` case is handled by refusal: the agreement check is skipped and
says why, since the derivation cannot see doses that were never written as rows.
`.recompute_tad()` itself is still wrong in that case and nothing under test
exercises it — `nmtest` would. The pre-dose convention is now stated in the
`tad` roxygen rather than only in a code comment.

## Next big piece: a vignette of the checks on synthetic data

Scoped 2026-08-03 with the owner. **The specification is
`design/SYNTHETIC_DATA_CHECKS.md`** — six categories (valid dataset, who is
singled out and by what, is it still the same study, are the numbers right,
does it work in the workflow, and what these checks cannot tell you), the
inventory of what already exists against each, the three gaps that inventory
makes visible, and the section that generalizes beyond this package: check the
finished table, not the mechanism.

- [x] Write `vignettes/synthetic-data-checks.Rmd` from that specification.
      Done 2026-08-04.
- [ ] Then, and only then, close the gaps the prose makes obvious — rare
      covariate combinations (B5, nothing exists), semantic ordering such as a
      trough staying a trough (C, nothing exported), and the `SIM-014` exact-copy
      gate as a user-runnable helper (B4, test-only today).

### Done 2026-08-05: the scorecard, and the checking literature

Owner reviewed the literature-derived proposal at the end of
`design/SYNTHETIC_DATA_CHECKS.md` and accepted three of its organizational
changes, plus a tutorial. Scope decision recorded with it: **specific-utility
measures are out of scope by design** — the package does not claim that a
modeller reaches the same scientific conclusion, so confidence interval overlap,
parameter recovery and exposure-metric agreement are deliberately not pursued.
Distributions and processes need not be maintained exactly; the output must be
the same *kind* of object and must protect privacy.

- [x] Scorecard. Static index table at the top of
      `vignettes/synthetic-data-checks.Rmd` (question, what to run, what it
      reads, pass criterion, gaps marked), and a **runnable** version computed at
      the end of `vignettes/demo.Rmd` on the real run. Release status is a column
      rather than a section-E paragraph, and `validate_pmx()` is rows A1/A2.
- [x] B5's propagation experiment moved to `vignettes/avatar-algorithm.Rmd`
      step 8, where it belongs: it explains the generator's mechanism rather than
      a check a user runs. The checks vignette keeps the conclusion, the census,
      and the differential privacy argument.
- [x] `vignettes/articles/literature-review.Rmd` restructured into two halves —
      how synthetic data is made, and **how it is checked** — with the second
      half written as a tutorial: the patient-versus-population confound,
      adversarial accuracy and why it needs a holdout, local cloaking and hidden
      rate, keys/targets/RepU/DiSCO/DiO, the WP29 criteria and linkability,
      the similarity-metric critique, alpha-precision/beta-recall/authenticity,
      pMSE versus confidence interval overlap, and a scope statement for why
      most utility measurement is not this package's problem. **Superseded
      2026-08-11:** the two halves became two articles, see below.

### Done 2026-08-11: the literature review split in two

- [x] `vignettes/articles/literature-review.Rmd` split at the halfway point into
      `synthetic-data-generation-review.Rmd` (the four generation families,
      differential privacy as a guarantee rather than a family, where `synpmx`
      fits) and `synthetic-data-checking-review.Rmd` (the checking tutorial).
      Each carries only the references it cites; the few citations both need are
      repeated rather than shared. Both appear under **Background** in the
      navbar, `articles/literature-review.html` redirects to the generation
      half, and `vignettes/synthetic-data-checks.Rmd` now links to each by name.

Two accuracy corrections found while doing it, both now fixed in the checks
vignette:

- [x] B5's claim that "nothing in the package checks this" was too strong.
      `compare_pmx_distributions()$covariates_categorical` *is* a per-patient
      level census of source against synthetic. Nothing **flags**, which is the
      accurate claim.
- [x] That census loops over `roles$covariates` only, so **`strata` are absent
      from it** — the one categorical axis copied from the anchor verbatim.
      Worth fixing in the function, not only in the prose: **open.**

Open, from the same proposal, in the owner's stated priority order (nothing
below is agreed work yet):

- [ ] **Holdout support.** The largest gap. Every source-reading privacy check
      compares the synthetic data against patients the generator was allowed to
      use, which confounds population capture with memorization. A `holdout =`
      argument plus running the proximity checks twice would fix it. Expensive at
      small n; recommend it as a per-study validation exercise, not a default.
- [ ] **Local cloaking and hidden rate.** Per-patient, published for exactly this
      generator family, and the natural input to
      `remediate_identifiable_subjects()`. Prerequisite is retaining the
      anchor -> avatar map, which is itself the most disclosive artifact in the
      pipeline and must never leave the trusted environment.
- [ ] **`strata` in `compare_pmx_distributions()`**, per the correction above.
- [ ] **B5 as one function** over the subject-level table, in `synthpop`'s
      RepU/DiSCO vocabulary rather than invented terms.
- [ ] **Linkability**, which nothing measures.

## Owver's next steps: 2026-08-03

- Try out on real data and apply checks 
- If happy with real data, then Go though the demo file and apply to all 5 examples
- If happy with all above, then go through full avatar algorithm page
- Identify a few places to spot check the code. Think about if/how to document.  Or think about other validation mechanisms
- Email authors of recent paper AVATAR to share

- In parallel - Meet with an internal colleague (mid august)
- Make poster (due Sept 28 for online submission)

## The owner's five next steps (2026-07-29)

Stated by the owner, in their order. Everything below this section is detail
that one of these five owns; if a task further down does not serve one of them,
it is not next.

1. [ ] **Land the dose-fingerprint masking and prove it on real data.** The
       mechanisms are built (anchor screen, `coarsen_time`, `min_pattern_share`,
       dose recomputed from the avatar's own covariate); what is unproven is all
       of it on a real study, and what is still copied is the *sequence* of dose
       levels an anchor climbed. Detail in "Open on the two new mechanisms"
       below; measurement is `skeleton_uniqueness()` and
       `compare_pmx_proximity()`.
2. [ ] **Schedule the internal review meeting.** Presentation following their
       template; a second internal group possibly after.
3. [ ] **Simulation study to evaluate privacy.** Scoped 2026-07-29 into
       "Next: empirical disclosure-risk evaluation" below and
       `design/PRIVACY_EVAL.md`.
4. [ ] **Trinity Metrics blog post: how to try `synpmx`.** New, nothing written
       yet. Source belongs in `communications/` (`.Rbuildignore`d, and see
       `communications/README.md` — public or fixture data only, internal review
       before it goes out). Largest reusable input is `README.md`'s runnable
       example; the post is the shortest path from "installed" to "ran it on my
       own data", not another method explanation.
5. [ ] **Make the ACoP poster.** `communications/2026-acop-poster-notes.md`.
       Still gated on the abstract's acceptance status and the deadline.

## Next Steps: adoption, external validation, and the ACoP poster

The package, its documentation, and the live site are all in place. The work now
is to get the method used by others, understood well enough to defend, and
presented. The ACoP (American Conference on Pharmacometrics) poster is the
culminating deliverable, and it paces the rest.

Try it on real data — no mode has run on a real study yet, and that is the test
that matters: role declaration against a real schema, event grammar the template
sampler has not seen, and whether the output is actually useful for workflow
development. Keep all of it in `scripts_private/`.

- [ ] Run AVATAR and the calibrated path on the internal INTERNAL_STUDY data.
- [ ] Ask Alex to try it on his own dataset.
- [ ] Consider asking Anwesha or Bambang to try it as well.
- [ ] Turn tester friction into `REV-###` entries in `REVIEW_BACKLOG.md` as it
      surfaces, so external use feeds back into the package.

The AVATAR extreme-structure work (`REV-025`/`REV-026`) is settled for now and
should be judged on real data before any more is built. Current design (see
`design/METHOD_DISCUSSION.md` §6a; guiding principle is good enough, not
perfect — output must not look extreme, and simplicity is valued):

- **Default, on:** `synpmx_avatar(screen = TRUE)` never anchors on a subject
  whose follow-up or dose count exceeds **twice the cohort's 90th percentile**,
  so no avatar looks structurally extreme; clean datasets are byte-identical to
  `screen = FALSE`. This went z-score → 2× median → 2× p90; p90 is the one that
  keeps ordinary spread. One toggle, no other knobs.
- **Manual layer:** `flag_identifiable_subjects()` (four axes) and
  `remediate_identifiable_subjects()` (truncate long / drop others / replace).
- **Parked:** skeleton sampling (`REV-026`) — the "perfect" cure, deliberately
  not built.

- **Default, on:** `synpmx_avatar(coarsen_time = TRUE)` collapses source times
  onto a shared visit grid before generation and resamples pooled deviations
  back afterwards, so no avatar carries one real subject's exact visit schedule
  (`SIM-033`). Exact when a `nominal_time` role is declared; inferred and
  best-effort otherwise, with a loud alert when it cannot collapse a subject
  (`SIM-034`). `time_jitter` is *not* an alternative -- its Voronoi clamp holds
  every time within half a gap of the source value at any magnitude.
- **Default, on:** dose recomputed from the avatar's own blended covariate where
  dosing is covariate-proportional within a `strata` stratum
  (`REV-027`). Fixes a coherence defect as well as the disclosure.
- **Default, on:** `min_pattern_share = 2` draws each avatar's attended-visit
  pattern from ones at least two subjects share (`REV-026`, partly closed), so no
  synthetic patient carries a schedule unique to a real one. Rare patterns are
  discarded rather than approximated, and the loss is reported per run.
- **Measurement:** `skeleton_uniqueness()` on the source, and
  `scripts/measure_skeleton_uniqueness.R` for the before/after table.

Open on the two new mechanisms:
- [ ] Confirm `min_pattern_share = 2` on real data. Attendance patterns on the
      public sets are distributed 18/2/1/1/1 -- one common pattern and a tail of
      singletons, no middle -- so the floor bites hard: on `warfarin` 2 keeps 2 of
      14 patterns and 3 keeps 1. If INTERNAL_STUDY and the oncology study have a
      populated middle, a higher floor may be affordable; if they look like
      `warfarin`, 2 is the only usable value and option D in the design
      discussion (decompose the pattern into miss-count plus placement rather
      than copying it whole) becomes the way to keep interruptions at all.
- **Done 2026-07-29.** The shape fallback closed most of the loss: across the
  public sets, patterns discarded fell from 12/17/12/54 to 2/4/2/1. What is left
  to confirm on real data is whether the *resolution* cost matters -- which
  specific visits were missed is no longer preserved, only how many and of what
  kind.
- **Done 2026-07-29.** `compare_pmx_proximity()` gives donor blending the
  measurement it lacked. Still open: whether the null is usefully narrow at
  INTERNAL_STUDY's cohort size, or whether it is so wide that only a blatant leak is
  detectable there too.
- [ ] `nimoData` has no pattern shared by even two subjects, so sampling has no
      pool. Constructing a nominal time (rounding the roughly-weekly infusions to
      their protocol week) would fix both this and its inferred-grid failure --
      worth demonstrating, since real studies will hit the same shape.
- [ ] The dose *amount* is now handled at any number of levels, intra-patient
      escalation included. What is still copied is the *sequence of levels* an
      anchor climbed; where escalation is outcome-adaptive that sequence encodes
      the subject's own response. The oncology study will show whether it
      matters.

Open, only if real data shows a need:
- [ ] Revisit the screen cut (2× the 90th percentile) if it over- or
      under-trims on a real study.
- [ ] Declare `nominal_time` in the private templates' roles blocks once the
      real schemas are confirmed to carry it. The placeholder is in
      `scripts_private/try_avatar*.qmd`; without it the inferred grid is
      best-effort and, on the public datasets, collapses nothing at all for
      warfarin and mavoglurant.
- [ ] Record the skeleton-exposure numbers per study in
      `scripts_private/README.md`'s study inventory, including the oncology
      repeated-dosing/intra-patient-escalation study, which no public dataset
      and neither INTERNAL_STUDY part covers.
- [ ] Consider randomized dropping of *observations* (not doses) if
      `n_obs_alone` stays high after coarsening on real data. Dropping doses was
      considered and rejected: it is protocol-invalid under fixed dose
      escalation, and the regime cannot be inferred reliably when intra-patient
      escalation makes dose amount vary within subject either way. Observation
      drops are protocol-neutral and cover a larger share of the residual.
- [ ] Expose the screen multiplier as an argument if per-dataset tuning is
      wanted (kept hardcoded now for simplicity).

Understand it well enough to defend it.

- [ ] Read the original AVATAR paper (`references/Guillaudeux23.pdf`). Feeds the
      novelty positioning owed under "Verification owed" below.
- [ ] Write my own mental map of how AVATAR works, checked against
      `vignettes/avatar-algorithm.Rmd`, to confirm my understanding.

Documentation the owner asked for (2026-07-25). Both are writing tasks that
depend on inputs not yet in hand, so they are queued, not started:

- [ ] A simplified, made-basic-and-clear explanation of how AVATAR works, from
      the 2023 paper (`references/Guillaudeux23.pdf`; now readable — poppler
      installed 2026-07-25). Likely lives in the method vignette or a new short
      article. Note the paper's k = 20 headline / k = 4 floor when writing the
      k-choice justification.
- [ ] A "what synpmx adds to AVATAR" section: the BLOQ/censoring handling, the
      full event-table (dosing + observation grammar) rather than a flat matrix,
      and — once designed and built — the minimum-donor pooling and rare-event
      pooling from `REV-025`/`REV-026`. Do not document the pooling as a feature
      until it exists.

Present it.

- [ ] Follow up with the internal colleague: a presentation following their
      template. Then potentially with a second internal group.
- [ ] Trinity Metrics blog post — instructions for trying the package. Draft in
      `communications/`, aimed at someone who has not installed it yet: install,
      declare roles on their own dataset, run `synpmx_avatar()`, read the
      warnings. Point at the pkgdown site for the method rather than restating
      it. Public or fixture data only, and internal review before posting.
- [ ] Confirm the ACoP abstract status and the poster deadline. The abstract is
      drafted (`communications/2026-acop-abstract.md`) — is it submitted and
      accepted, and by when is the poster due? This gates everything above.
- [ ] Make the ACoP poster (`communications/2026-acop-poster-notes.md`). Public
      or fixture data only, and internal review before submission — see
      `communications/README.md`; any internal-study figure stays in
      `scripts_private/` and never on the poster.

Deferred.

- [ ] Decide how date/datetime columns should be handled. Low priority; dates
      are rarely analysis-relevant. Today `time` must be numeric elapsed time,
      and a raw `RFSTDTC`-style datetime column is either converted by the user
      beforehand or dropped as undeclared. Not needed for INTERNAL_STUDY.

## Done: documentation reorganization (2026-07-23)

Decided and executed. The reasoning, the audience analysis, and the rationale
for each call are in `design/DOCUMENTATION_SCOPE.md`; delete that file once
this section is stale. `AGENTS.md` now records the resulting three-tier rule.

- [x] Adopt pkgdown. `_pkgdown.yml` with a grouped reference index over all 31
      exports, plus a GitHub Actions workflow deploying to `gh-pages`.
- [x] `README.Rmd` → `README.md` as the entry point: pitch, one runnable
      example, the four-mode table, and the documentation map. 247 lines of
      API reference and limitations came out.
- [x] Vignette set cut from five to three: `synpmx-method` (all four modes,
      high level), `synpmx-demo`, `synpmx-privacy`. `synpmx-intro` and
      `synpmx-epsilon-exploration` were merged away.
- [x] Deep AVATAR mathematics moved out of the method vignette into
      `vignettes/avatar-algorithm.Rmd`.
- [x] Five design documents moved to `vignettes/articles/`: feasibility,
      privacy background, privacy argument, model and data elicitation.
- [x] Every citation into `design/` removed from shipped documentation and
      roxygen comments; roxygen now links to the website.
- [x] Deleted `design/METHODS_VIGNETTE_SPEC.md` (stale) and `scripts/README.md`
      (stale). `NEWS.md` reduced to a stub.

### Earlier scope: one entry point, four named modes

- [x] An entry point covering all four modes applied to `theo_md`, a properties
      table, and a table mapping environments (trusted / partner / published)
      to acceptable modes. Now the four-mode tour in `synpmx-method`.
- [x] Privacy vignette: explain what AVATAR *is* and what DP *is* before
      comparing them, with the formal `(epsilon, delta)` definition and the
      kind-not-degree table.
- [x] Method vignette: all options at the top, then the default AVATAR
      algorithm, then the model-based alternatives at the end.
- [x] Demo vignette: state the four modes up front and run the model-based path
      on theophylline, as `scripts/demo_nlmixr2data.R` does.
- [x] `README.md`: table of contents explaining how the documentation set is
      organized and which document answers which question. Superseded by the
      full `README.Rmd` rewrite above.
- [x] House style: spell out every acronym on first use in a document
      (`AGENTS.md`). The word "mock" is gone: prose says "synthetic data", the
      `compare_pmx()` argument and outputs are `synthetic`, and generated
      character/factor IDs are `syn_001` rather than `mock_001`.
- [x] `design/DOCUMENTATION_SCOPE.md` — inventory of all 23 documents with
      guessed audiences, since rewritten as the decision record.

## Version 4 — return to AVATAR blending as the primary method

Scope decision (2026-07-22): after comparing to Novartis's `synadam` (which
resamples each column marginally from the data with no formal guarantee), AVATAR
is the trajectory-level analogue of the same governance-based approach, and it is
the right default when the synthetic data reaches no one the source data could
not. The DP (v2) and structural
(v3) engines are kept as superseded alternatives for when a formal guarantee is
required, not removed.

- [x] `design/METHOD_DISCUSSION.md`: the AVATAR vs DP tradeoff essay, the
      trajectory-is-a-fingerprint asymmetry, and the synadam parity argument.
- [x] `PROTOTYPE_SPEC.md`: Version 4 section (history) plus a section 0 banner
      making AVATAR the default and the trust-boundary the decision rule.
- [x] Restore the AVATAR engine as `synpmx_avatar()` (renamed from the v1
      `mock_pmx`), `synthesis.R`, `profiles.R`, plus the ported `utils.R`
      helpers. Exported and working; no name collisions with v2/v3.
- [x] AVATAR tests, including all five nlmixr2data datasets. `test-avatar.R`
- [x] Rebuilt the demo vignette around AVATAR, keeping all five nlmixr2data
      datasets, and the method vignette from the Version 1 "How synpmx
      Works".
- [x] Method vignette explains the (epsilon, delta) vs AVATAR distinction and
      the synadam parity argument.
- [x] Slimmed the two DP vignettes to short "formal-privacy alternative" asides.
- [x] `./build.sh` clean, all four vignettes knit, 395 tests pass.

Version 4 is complete.

Hardening completed after review: `REV-018` and `REV-019` are fixed with
regression tests in `test-structural-v3.R` and `test-avatar.R`. The remaining
open findings below apply to the superseded formal-DP v2 path.

## Superseded: v3 low-dimensional structural generator (kept as an alternative)

The scope decision (2026-07-22): prefer small trials over pooled corpora. That
requires releasing a handful of parameters against public structural priors
instead of a dense grid. See `vignettes/articles/feasibility.Rmd` section 8 and
`design/PROTOTYPE_SPEC.md` "Version 3 scope".

### Core: model in, correction out

- [x] `pmx_structural_model()` — public structural model with built-in analytic
      1-cmt IV/oral/infusion PK and direct-effect or indirect-response PD.
      Built-ins need no compiler; `rxode2` is accepted and validated but not yet
      wired through. `R/structural.R`
- [x] `pmx_trial_design()` — dose levels, cohort sizes, protocol sampling
      schedule, dosing interval, infusion duration, visit windows.
- [x] Multiplicative correction release rather than absolute parameters.
      `synpmx_calibrated()`, `d = 2-3`.
- [x] `pmx_prior(range, source)` / `pmx_priors()` with mandatory provenance.
- [x] Per-subject NCA estimator. The correction is the ratio of predicted to
      observed AUC on the subject's own grid, which avoids needing `F` or an
      extrapolation to infinity.
- [x] Sampling schedule from the protocol. The timing release group is gone.
- [x] `at_prior_boundary` diagnostic, warned on and shown by `print()`.
- [x] Prior mode is first-class: `pmx_generate(model, design)` takes no data and
      no epsilon.
- [ ] Support ladder levels 3-4 (intra-patient escalation, titration):
      occasion-varying dose applied to *generated* subjects from the public
      rule. Never replay a source subject's sequence.
- [ ] `REV-017` Record realized trial-design quantities in `proof_assumptions`
      with a source field.
- [x] Two-compartment PK (`2cmt_iv`, `2cmt_oral`), analytic and verified
      against `Dose/CL`.
- [ ] Transit-absorption option.
- [x] Covariates: **out of scope**, stated in `PROTOTYPE_SPEC.md` section 6.
      Users join their own; covariate-handling code is not exercised by this
      package's output.
- [ ] Schema flexibility. Output is a fixed column set; real datasets are not.

### Guardrails

- [ ] Enforce public model selection ergonomically: the model is built and
      passed in before `synpmx_empirical()` sees data, so there is no easy route
      to fit, look, and revise.
- [ ] Guard against the prior range being set from the data without budget.
      Mean +/- 2 SD is the right target for a private stage 1 to *discover*;
      computing it directly voids the guarantee undetectably.
- [ ] Document that confidential data must never be sent to an external service,
      including an LLM. Outside the accounting entirely.

### Optional, if priors turn out too wide

- [ ] Two-stage range-finding: coarse log-scale histogram at `epsilon_0` (~20%
      of budget, L1 sensitivity 1 regardless of bin count), then clip and
      release means at `epsilon_1`. Measured worth with an absolute prior:
      2.44-fold -> 1.52-fold at N = 20, epsilon 1. Likely unnecessary once the
      release is a correction factor, since that prior is already tight.
- [ ] Small budget for realized messiness rates (dropout, missed dose, BLQ) if
      public assumptions prove too crude.

### Verification owed

- [ ] Verify the `rxode2` templates in `vignettes/articles/model-elicitation.Rmd` compile and
      produce sensible profiles. Authored without a working C compiler and never
      executed.
- [ ] Measure the correction-factor parameterization. The ~1.37-fold estimate at
      N = 20, epsilon 1 is arithmetic from the error law, not a measurement.
- [x] Measure the PD correction. Exact without residual error; biased low by
      about a third with 15% residual on a small deviation. Documented in
      `vignettes/articles/feasibility.Rmd`; PD is experimental.
- [x] Improve the PD estimator. Solved by changing the endpoint rather than the
      statistic: simple time-course shapes take a level correction (ratio of
      means), which is unbiased under residual error. Exposure-driven PD keeps
      the signed-area statistic and stays experimental.
- [x] Retire the exposure-driven PD shapes. Removed entirely; PD is now a
      simple exposure-independent time course.
- [ ] Literature check before claiming novelty: DP + non-compartmental analysis,
      DP + popPK, DP synthetic data under informative structural priors.

## Next: empirical disclosure-risk evaluation (proposed 2026-07-29)

Raised by the owner. Measure whether an attacker can actually recover anything
from released synthetic data — membership inference, singling out, linkability,
attribute inference — rather than only arguing from the trust boundary. The
program, its phasing, what already exists, and why a cohort of 30 limits what
any of it can conclude are in `design/PRIVACY_EVAL.md`. Not started, and it
competes with the ACoP poster and the real-data validation above for the same
weeks.

- [ ] Decide what it is for — internal defensibility, the poster, or a paper.
      That sets the depth, and nothing else here starts until it is answered.
- [ ] P1: read the attack literature and Guillaudeux 2023's own metrics, and
      settle whether the attacks run on the profile vector or get restated for
      event tables. Shares the paper read already queued below.
- [ ] P2: train/holdout membership inference over the public datasets, with the
      verbatim-copy positive control and `synpmx_prior()` as negative control.
      The cheap slice that says whether the rest is worth building.

## Next: covariate relationships and covariance

Raised by the owner 2026-07-28, after the pyrazinamide example
(`scripts/example_prior_pyrazinamide.qmd`) had to omit covariates entirely.
**Status measured, not assumed** --- the numbers below are from a 80-subject
source with a known allometric WT effect on CL, 15 synthetic replicates.

| Mode | Covariate -> PK/PD | Covariate <-> covariate |
|---|---|---|
| `synpmx_avatar()` | **Retained, mildly diluted.** Measured over 30 runs in the worked article: allometric exponent 0.75 -> 0.68 (92%), dose effect 104%, treatment effect 99%, exposure-response slope 89%. All within ~10% of source. | Retained by the same mechanism. |
| `synpmx_prior()` | **None.** | **None.** |
| `synpmx_calibrated()` | **None.** | **None.** |
| `synpmx_empirical()` | **None.** | **None.** |

**Correction, 2026-07-28.** An earlier version of this entry recorded the
covariate effect as *amplified* to ~143%. That was an artifact of a confounded
measurement: log AUC was regressed on log WT without dividing out the dose, so a
chance weight/arm imbalance leaked into the slope. Dose-normalised, the source
slope is -0.747 -- the allometric exponent 0.75 itself -- and AVATAR returns
0.683, a mild dilution. There is no amplification. Every relationship measured
moves the same way: slightly toward the null.

Why AVATAR keeps it: `.synthesize_covariates()` and `.synthesize_trajectories()`
are handed *the same* `donors$indices` and `donors$weights`
(`R/synthesis.R`), so a subject whose concentrations came mostly from donor A
also gets mostly donor A's weight. The relationship survives because it is never
taken apart. The correlation attenuates --- blending plus `subject_noise_sd`
adds noise to both sides --- but the *slope*, which is what a covariate model
estimates, comes through close to intact. That is a stronger result than
expected and worth not breaking.

Why the model-based modes lose it: `.generate_structural()` draws parameters at
line 118 and simulates the whole profile from them; covariates are drawn at line
184 by `.draw_covariate_table()` and `merge()`d onto the finished table. The
concentration exists before the weight does. And `.draw_covariate_table()` loops
over covariates one at a time, so there is no joint structure between them
either --- two correlated covariates come out independent.

- [ ] **Decide whether the model-based modes should carry covariate effects at
      all.** Owner's position (2026-07-28) is that they probably should not: the
      space of PMX models is unbounded, every study wants something
      idiosyncratic, and a user is better served writing the model they want in
      `rxode2`/`nlmixr2` --- with an LLM and a skill file as the assistant ---
      than by this package growing a model language. `vignettes/synpmx-4-methods.Rmd`
      now says so outright ("The built-in models are illustrative, and
      deliberately so"). If that holds, this item closes as *documented, not
      built*, and the work is to keep the documentation honest.
- [ ] **If it is built anyway, the contained version** is a `covariate_effects`
      argument on `pmx_structural_model()`, drawing covariates *before*
      parameters and scaling each subject's parameters by their covariate terms.
      Needs the functional form pinned down first: whether a categorical effect
      of -0.40 means `x (1 - 0.40)` or `x exp(-0.40)`, and whether a continuous
      effect is a power on a centred covariate or a linear term. Would benefit
      `synpmx_calibrated()` for free, since it shares the generator.
- [ ] **Covariate--covariate covariance is a separate and cheaper problem.** Even
      without parameter effects, two covariates that are correlated in reality
      (weight and height, age and renal function) come out independent from the
      model-based modes, and from `pmx_covariates_auto()`, which resamples each
      column marginally in the `synadam` style. A joint draw would fix that
      without touching the structural model at all.
- [x] **Protect what AVATAR already does.** Done 2026-07-28:
      `test-avatar-relationships.R` pins the covariate slope, the arm dose
      ratio, arm/dose coherence for a `keep` column, and structurally that
      `.select_donors()` is called once per subject and feeds both synthesisers.
      Bounds are deliberately wide -- these catch destruction, not drift.
      Worked evidence in
      `vignettes/articles/example-avatar-PKPD-covariate-treatment-effect.Rmd`.
- [ ] **Revisit the per-endpoint donor weighting idea in that light.** Reweighting
      donors per endpoint was floated to stop a densely-sampled endpoint
      dominating donor selection (see the multi-endpoint section of the AVATAR
      article). It would also break the single shared donor draw that every
      relationship above depends on. The two goals are in direct tension; decide
      deliberately rather than discovering it in a test failure.

## Next: correctness and privacy findings

- [ ] `REV-003` Data-dependent `stop()` on confidential rows before any noise is
      applied. An unaccounted output channel in the core DP claim.
- [ ] `REV-016` Pooling breaks the stated adjacency when one person appears more
      than once (rollover, extension, crossover, re-enrollment). Group privacy
      silently degrades the guarantee to `k * epsilon`.
- [ ] `REV-004` `delta` is validated and requires justification but is never
      spent. Either restrict to `delta = 0` or wire it to a real mechanism.
- [ ] `REV-002` Pre-flight feasibility check, so an infeasible `(N, epsilon, d)`
      is refused before budget is spent rather than discovered in a plot.
- [ ] `REV-030` `addl`/`ii` are accepted, carried, and never read, so a source
      encoding repeated dosing as `ADDL` presents to AVATAR as a single-dose
      study. Refuse the roles with a message telling the user to expand first,
      or expand and re-collapse around synthesis. Found 2026-07-30.
- [ ] `REV-020` `pmx_structural_model(rx = )` is stored but never used; it now
      warns. Either wire it through `rxode2::rxSolve()` with a regression test
      against the analytic solution, or reject it outright. (Was stranded in the
      completed documentation-reorganization list.)
- [x] `REV-025` **Owner-flagged, 2026-07-25; closed 2026-07-27.** AVATAR
      borrows the nearest donors across dose/schedule groups to blend `k` (=5)
      real patients into every avatar, with a loud red alert when the source has
      fewer than `k + 1` subjects (`.select_donors()`, `.loud_warn()`).
      Discovery: weight-based dosing made nearly every subject its own
      singleton, so this was the common case, not an edge case. Outlier detector
      done too: `flag_identifiable_subjects()` screens follow-up time, dose
      count, dose magnitude, and DV. Closed out 2026-07-27 with the route
      barrier (IV/bolus/oral never blended; short route arms dropped from the
      anchor pool) and `max_donor_weight` 0.80 -> 0.50 enforced on every donor.
      0.50 was chosen on how often the cap fires (68%, effective donors 2.90):
      0.30 fires on 99% and *is* the weighting scheme rather than a guardrail,
      and it cost too much between-subject variability; 0.80 fires on 15% and
      still lets one real patient be most of a synthetic one.
      Tests in `test-avatar-pooling.R`, `test-avatar-weights.R`; algorithm
      written out in `vignettes/avatar-algorithm.Rmd` Steps 6-7.
      Still open: whether the headline floor should exceed 5 -- now the weaker
      of the two constraints, since the cap binds first. Judge on INTERNAL_STUDY.
- [ ] `REV-026` **Consider, not committed (2026-07-25 decision).** "Skeleton
      sampling." Today each avatar copies one real anchor's *event skeleton* —
      the number of doses, their sizes, and the observation times — verbatim,
      and only the DV values are blended. So a source subject with a unique
      structure (the long-followed `wbcSim` subject, a one-off dose) yields an
      avatar with that same unique structure, which is identifying even though
      its values are blended. Skeleton sampling would instead *draw* each
      avatar's schedule (follow-up length, dose count, observation times) from
      the cohort distribution, so no avatar carries any one real subject's
      unique structure and there is nothing to flag. It is the thorough,
      preventive cure.
      **Why it is paused:** `flag_identifiable_subjects()` +
      `remediate_identifiable_subjects()` (with replacement) already remove
      structural outliers after generation and refill the cohort, which covers
      the same privacy goal from the other direction. The owner chose the
      detect-and-remediate path for now; skeleton sampling stays a possible
      future direction, to weigh against the diversity it would add versus its
      cost. Design notes in `design/METHOD_DISCUSSION.md` §6a; would build on
      the `REV-025` relaxed pooling.

## Then: utility headroom in the existing dense-grid path

Keep this path for pooled corpora; it is not superseded by v3.

- [ ] `REV-005` Derive sensitivity from the declared contribution limits instead
      of matrix width.
- [ ] `REV-006` Gaussian/zCDP accounting. **Note:** only worthwhile at high
      dimension. At `d` around 6 pure-DP Laplace is the better mechanism, so
      this belongs to the dense-grid path, not to v3.
- [ ] `REV-007` Low-dimensional trajectory representation. Largely subsumed by
      the v3 work above.
- [ ] `REV-008` Name and document the post-processing rescue constants; add a
      diagnostic for how much output shape comes from the release versus the
      fallbacks.

## Hygiene

- [ ] `REV-009` `.public_input_manifest()` ignores all six arguments and returns
      a constant. Derive it or delete it.
- [ ] `REV-010` Replace the nine-name payload denylist with an allowlist.
- [ ] `REV-011` Add `print.pmx_private_validation` and `print.pmx_backend_tests`.
- [ ] `REV-012` `run_dp_backend_tests()` advertises a `privacy_map` test it does
      not actually assert.
- [ ] `REV-013` `enable_features("contrib")` mutates global OpenDP state on every
      backend resolve.
- [ ] `REV-015` Convert load-bearing vignette prose into assertions or tests.

## Done

- [x] `REV-001` / `SIM-020` Scale-aware support threshold in decoding. `44db89f`
- [x] `REV-014` Land the working tree in reviewable commits.
- [x] `./build.sh` for `R CMD check` and clean-library vignette rendering. `778848b`
- [x] `vignettes/articles/feasibility.Rmd` scoping assessment. `d5b0e30`
- [x] `REV-023` Session-level `synpmx_enable_dp_engines()` gate on the DP
      engines' unaudited status, so it is enforced rather than only documented.
- [x] `REV-024` Trust boundary reframed as organizational rather than
      geographic, so taking synthetic data to a local machine for code
      development is a supported use rather than a boundary crossing.
- [x] Publish the pkgdown site. Live and verified 2026-07-24: landing page,
      articles, and reference all render, and `AGENTS.html` / `CLAUDE.html` are
      pruned. The stringfish/RcppParallel ABI break that failed every build is
      fixed by keeping `rxode2` out of the docs job. `f1326fe`
- [x] Delete `.github/workflows/r.yml`, which could never pass (R 3.6.3 leg
      against `DESCRIPTION` R >= 4.1.0); `./build.sh` checks locally. `2521f8e`
- [x] One README. Deleted `README.Rmd`; `README.md` is the only entry point,
      its example pinned by `tests/testthat/test-readme.R`. `e09011d` `215cf50`
