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

| Input | What it is | Privacy status | Who can produce it |
|---|---|---|---|
| **Config** | structural model, priors, role mapping, sampling schedule | Public. From protocol + preclinical prediction | A pharmacometrician, or an assisting AI reading the protocol and the dataset headers |
| **Regimen skeleton** | one row per subject: dose(s), cohort, timing | Design if prespecified; outcome if adaptive | See below |
| **Observed DV** | the actual concentrations / responses | Confidential | Only `synpmx_calibrated()` ever reads it, and only to extract a correction factor |

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
  data that may CROSS A TRUST BOUNDARY and needs a formal (epsilon, delta)
  guarantee. Needs a structural model, priors, and a trial design. This is the
  secondary, provided-as-is path (see the package README's maintenance status).

`try_avatar_pit565.qmd` is a filled-in AVATAR example on a real study, kept as a
reference for what a completed roles block looks like.

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
