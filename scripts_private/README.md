# Private company-data work

Scratch space for trying the calibrated generator on a real modeling dataset,
inside the safe computing environment.

## This folder does not go to git

The root `.gitignore` ignores everything here **except** `README.md` and
`try_dp_calibrated.R`. Real datasets, generated tables, fitted models, and every
figure stay local. Before adding any new tracked file, confirm it contains no
patient data and allow-list it explicitly.

- Put the real dataset in `scripts_private/data/` (ignored).
- Output lands in `scripts_private/output/` (ignored).
- Never set `public_source = TRUE` for confidential data.
- The source-vs-synthetic comparison is a **restricted** diagnostic. It is derived
  from the real data and must stay in this environment.

## The three inputs, and their privacy status

The workflow separates a trial into three inputs that differ in how sensitive
they are. Keeping them separate is what lets most of the work happen outside the
privacy budget.

**Config** — structural model, priors, role mapping, sampling schedule
- *Privacy status:* public, from protocol + preclinical prediction
- *Who can produce it:* a pharmacometrician, or an assisting AI reading the
  protocol and the dataset headers

**Regimen skeleton** — one row per subject: dose(s), cohort, timing
- *Privacy status:* design if prespecified; outcome if adaptive
- *Who can produce it:* see below

**Observed DV** — the actual concentrations / responses
- *Privacy status:* confidential
- *Who can produce it:* only `synpmx_calibrated()` ever reads it, and only to
  extract a correction factor

The config can be drafted from public documents alone, so an AI that sees only
the protocol and the column headers (not the data) can produce it. That is the
intended path; see `vignettes/articles/model-elicitation.Rmd` and
`vignettes/articles/data-elicitation.Rmd`, which are written to be worked by a person or an
agent.

## The regimen skeleton

The most faithful way to describe a real trial's dosing is not a parametric
`dose_levels` / `cohort_sizes` summary — real dosing is messier than that. It is
a **one-row-per-subject (or few-rows-per-subject) table of the dose regimen**:
who got what dose, when, on which occasion.

Whether that table is a public input depends on one question:

- **Prespecified design** — fixed cohorts, protocol-defined titration. The
  regimen is the protocol applied to each subject, and can be used directly as
  a public template. Using the real one reproduces the realized cohort sizes and
  dose distribution exactly; that is a realized-design disclosure, treated as
  public by assertion and recorded (see `REV-017` and
  `vignettes/articles/data-elicitation.Rmd` section 3).
- **Outcome-adaptive** — the dose a subject received depends on their own
  tolerability or response (adaptive escalation, response titration). Then the
  per-subject dose *sequence* encodes that subject's outcome and must **not** be
  copied row-for-row. Generate the regimen from the public rule instead, or work
  from the planned design.

This step is being shaped; see the note at the top of `try_dp_calibrated.R` and the
project TODO. For a first pass the template uses the parametric design, which is
correct for fixed cohorts.

## Templates, one decision

Pick by the trust boundary (see the `synpmx-privacy` vignette for the decision
rule):

- **`try_avatar.qmd`** — start here. AVATAR blending (`synpmx_avatar()`) for
  output that reaches no one the source data could not — which includes taking
  it to your own workstation under the same access controls. Simpler, more
  faithful, no formal privacy guarantee. Fill in the data path and the column roles. Runs
  chunk-by-chunk like a script, or renders a source-vs-synthetic report. Prefer
  a plain `.R`? `knitr::purl("try_avatar.qmd")` writes one out — there is a
  single template and the script is generated from it.
- **`try_dp_calibrated.R`** — the differentially private structural path, for
  output that may REACH SOMEONE THE SOURCE DATA COULD NOT (a partner, a vendor,
  a publication) and needs a formal (epsilon, delta) guarantee. Needs a structural model, priors, and a trial design. This is the
  secondary, provided-as-is path (see the package README's maintenance status).

`try_avatar_pit565a1.qmd` and `try_avatar_pit565b1.qmd` are filled-in AVATAR
examples on a real study, kept as a reference for what a completed roles block
looks like.

## The study inventory

Each real study exercises a different part of the generator, and no public
dataset covers some of them at all. This table is the working record of what is
available and what has actually been checked on it. **Keep it current** — it is
the private counterpart to the dataset registry in `design/TEST_SIM.md`, and the
only place the coverage is written down.

Name studies by their internal identifier only. No patient-level facts, no
enrollment figures, no results belong in this file: it is the one tracked file in
this folder, so anything written here leaves the environment.

### PIT565 A1
- **Template** — `try_avatar_pit565a1.qmd`
- **Design shape** — Phase 1 dose escalation, with subsequent weekly dosing.
  Dose interruptions exist.
- **Why it is here** — first real schema
- **Checked** — roles, validation

### PIT565 B1
- **Template** — `try_avatar_pit565b1.qmd`
- **Design shape** — Phase 1, dose escalation
- **Why it is here** — fixed 3 doses per subject with intrapatient dose
  escalation
- **Checked** — roles, validation

### Oncology study *(to add)*
- **Template** — not yet written
- **Design shape** — repeated dosing with intra-patient escalation
- **Why it is here** — the case fixed-dose studies do not cover: dose amount
  varies *within* subject, so the event signature stays unique after coarsening
  and a dropped-dose mechanism would be protocol-invalid
- **Checked** — not yet run

### What to record per study

Run `scripts/measure_skeleton_uniqueness.R`'s `measure()` against the study and
record the structural properties below. `pmx_masking_report(synthetic, raw,
roles)` reports the same quantities for one run with an explanation beside each
number, and `plot_pmx_schedule(raw, roles)` shows the schedule the counts are
describing. They are what decide whether the default
mechanisms do anything, and they are cheap to compute:

Design, which decides whether the mechanisms can do anything at all:

- **Subjects** — every mechanism's strength scales with it
- **`nominal_time` declared?** — an exact snap to the protocol grid versus the
  best-effort inferred one. The single most important field here
- **Times nominal or actual?** — decides whether the schedule is exposed at all
- **Dose regime** — fixed, escalating, or intra-patient escalating
- **Doses per subject** — fixed or variable

Exposure, in the same terms the demo vignette's table uses:

- **`alone` before coarsening** — subjects whose visit schedule no other subject
  shares. Under actual recorded times, expect nearly all of them
- **`alone` after coarsening** — the number you would consider dropping
- **`unique_moment`** — of those, the ones observed at a time nobody else was.
  This is the grid's job, and declaring `nominal_time` is what fixes it
- **`unique_pattern`** — of those, the ones whose every time is shared and only
  whose combination of attended visits is unique. Dropout; no grid touches it
- **`unique_dosing`** — unique dose structure or amount. Weight-based dosing
  keeps this high regardless of what the grid does
- **`dose_basis`** — whether dose amounts were recomputed from a covariate. The
  **Dose** section of `pmx_masking_report()` answers this in words and, when
  the answer is no, states which covariates were tried and what failed. "No" on
  a study you believe is weight-based usually means the dose-to-covariate ratio
  did not collapse onto a small number of protocol levels; check whether a dose
  group needs declaring in `subject_properties`
- **Avatars keeping their anchor's own visit set** — the complement of
  `pattern_sampled_fraction`. Above 0 means some stratum had no visit set
  shared by `min_pattern_share` subjects, and the run will have alerted. Those
  avatars carry one real patient's absences, so read it against the discarded
  count directly above it in the report

Also record the existing utility metrics on the same row —
`compare_pmx_distributions()`, `flag_identifiable_subjects()` counts,
`cap_binding_fraction`, `mean_effective_donors` — so a privacy change and its
utility cost can be read off one line rather than reconstructed.

## Running the AVATAR template

1. Edit the two configuration chunks: the data path, and the column roles.
2. Run. Undeclared columns are dropped (the run says which); a failed validation
   names the role and column to fix.
3. Inspect the restricted comparison figure before treating anything as usable,
   or exporting it.

## Running the DP script

1. Fill in the `CONFIG` block: the structural model, priors, and trial design.
   Nothing in it should come from looking at the data — see the cardinal rule in
   the model-elicitation article.
2. Run with `DRY_RUN <- TRUE` first. This exercises the whole pipeline on public
   simulated data and never touches your dataset.
3. Read the pre-flight verdict; if `f >= 1` the release will not beat the prior
   and you should not spend budget.
4. Set `DRY_RUN <- FALSE`, point `DATA_PATH` at the real dataset, and run.
5. Inspect the output and the restricted comparison before exporting anything.
