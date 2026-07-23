# Documentation scope

Working document for the documentation-reorganization decision in
`design/TODO.md` item 1. Created 2026-07-23.

**Status: an inventory with guessed audiences, not a decision.** The audience
column is my inference from each document's content and tone, not something
anyone has stated. Correct it before we use it to decide anything — a wrong
guess here produces a wrong merge later.

The package currently ships **23 documents totaling roughly 6,800 lines**, of
which about 2,300 lines are user-facing (vignettes and `README.md`) and 4,000
are internal design records. That ratio is the first thing worth discussing.

---

## 1. User-facing: vignettes

Shipped with the package; visible to anyone who installs it.

| Document | Lines | Guessed audience | Job it does today | Concern |
|---|---:|---|---|---|
| `pmxSynthData-intro.Rmd` | 352 | A pharmacometrician who has just heard of the package and does not yet know what "differentially private" buys them | Big picture; all four modes run on `theo_md`; properties table; environment→mode table | New (2026-07-23). Overlaps `demo` on theophylline and `privacy-intro` on the decision rule — deliberately, as an entry point, but the boundary needs agreeing |
| `pmxSynthData-demo.Rmd` | 587 | A pharmacometrician who has chosen a mode and wants to run it on their own study | Worked AVATAR workflow over five public datasets with structural checks, plus one model-based example | Half its length is plotting helper functions in the setup chunk that no reader needs to see |
| `pmxSynthData-simulation-method.Rmd` | 1119 | Someone who must defend, debug, or review the generator — a methods-minded colleague or an internal reviewer | Step-by-step AVATAR algorithm with the mathematics, edge cases, and a worked example; alternatives at the end | **The longest document in the package.** Is this a vignette or a specification? It reads like both |
| `pmxSynthData-privacy-intro.Rmd` | 202 | Someone deciding whether a release is allowed to leave the environment — possibly not a modeler at all | What AVATAR is, what DP is, the formal definition, and the trust-boundary rule | Its audience may be the one least likely to open an R vignette |
| `pmxSynthData-epsilon-exploration.Rmd` | 74 | Someone who has already committed to a DP mode and is choosing epsilon | `f = d/(epsilon N)`, the measured frontier table, `pmx_preflight()` | Shortest by far. Is it a vignette or a section of `privacy-intro`? |

**Open questions.** Is five the right number? Does the method vignette want to
become a specification appendix instead? Does `epsilon-exploration` survive as
its own document, or fold into `privacy-intro` now that `intro` carries the
decision rule?

## 2. User-facing: repository front matter

| Document | Lines | Guessed audience | Job it does today | Concern |
|---|---:|---|---|---|
| `README.md` | 247 | Someone browsing the GitHub repository who has not installed anything | Currently doing four jobs at once: pitch, privacy contract, API reference, limitations list, documentation map | **Needs the rewrite.** The API-reference and limitations material is deep internal detail on a page whose first job is to explain what this is and why to care |
| `NEWS.md` | 144 | An existing user upgrading | Version history by feature | Reasonable as-is. Written for users who do not yet exist |
| `scripts/README.md` | 15 | A developer running the public demonstrations | Describes `demo_nlmixr2data.R` and `evaluate_simulations.R` | Refers to `test-avatar.R`, which is not in `scripts/`; stale |
| `scripts_private/README.md` | 87 | Whoever runs this inside the safe environment — realistically you | The trust-boundary rules for the private scripts folder and what may be committed | Longer than the public scripts README, and the most operationally important page in the repository |
| `references/README.md` | 9 | A contributor adding a paper | Where to put references and what not to commit | Fine |

## 3. Contributor-facing: repository conventions

| Document | Lines | Guessed audience | Job it does today | Concern |
|---|---:|---|---|---|
| `AGENTS.md` | 100 | An AI coding agent, and secondarily a human contributor | Repository conventions: where code goes, documentation-synchronization rules, the acronym rule, testing discipline | Doubles as the de facto contributor guide. There is no `CONTRIBUTING.md` |
| `CLAUDE.md` | 3 | Claude Code specifically | Points at `AGENTS.md` | Correct and minimal |

## 4. Internal: design record

Not shipped with the package. This is where most of the writing lives.

| Document | Lines | Guessed audience | Job it does today | Concern |
|---|---:|---|---|---|
| `PROTOTYPE_SPEC.md` | 1091 | Future you, and any reviewer asking "what is this supposed to do?" | The implementation contract, with four versions of scope history newest-first | Longest document in the repository. Carries three superseded designs alongside the current one |
| `FEASIBILITY.md` | 665 | Someone asking "will this work at N = 20?" | The measured evidence behind every scope decision; DP utility by cohort size | The strongest document in the set. Arguably deserves to be user-facing, at least in part |
| `TEST_SIM.md` | 415 | Whoever is fixing a simulation defect | The living evaluation specification: dataset registry, `SIM-###` issues, metrics, gates | Working document; healthy |
| `REVIEW_BACKLOG.md` | 300 | Whoever is fixing a mechanism or API defect | `REV-###` findings from code review, with status | Working document; healthy |
| `MODEL_ELICITATION.md` | 300 | Someone preparing the public inputs for a DP fit | The interview producing a public structural model and priors without touching data | Only relevant to modes 2–4. Would a user ever find it? |
| `DATA_ELICITATION.md` | 291 | The same person, for trial design | The complexity ladder and which parts of a protocol are genuinely public | Same question |
| `PRIVACY_BACKGROUND.md` | 250 | Someone who wants the arithmetic to make sense | `d`, `f`, sensitivity, the error law, worked examples | Tutorial in tone. **This reads like a vignette that ended up in `design/`** |
| `TODO.md` | 244 | You and me, this week | The working queue, and the index to every other design document | Currently the entry point to `design/`. Works |
| `PRIVACY_ARGUMENT.md` | 180 | A privacy reviewer who will not read R code | The formal mechanism-level argument for the v2 engine | Explicitly labeled as not independently reviewed. Scoped to the superseded v2 path |
| `METHOD_DISCUSSION.md` | 169 | Someone asking "why is AVATAR the default?" | The AVATAR-vs-DP essay and the `synadam` parity argument | Its content is now partly duplicated across three vignettes |
| `METHODS_VIGNETTE_SPEC.md` | 56 | Whoever maintains the vignettes | Vignette maintenance notes | **Stale.** Describes "Version 2 has four documents"; we now have five and a different structure. Candidate for deletion |
| `DOCUMENTATION_SCOPE.md` | this | You and me | This inventory | Delete once the reorganization is decided and executed |

---

## 5. What I would want to decide, in order

1. **Who is the actual audience?** Everything above is a guess. The plausible
   set is: (a) you, (b) pharmacometrician colleagues inside the company,
   (c) an internal privacy or governance reviewer, (d) the open-source R
   pharmacometrics community. These want very different documents, and the
   current set reads as though it is serving (a) and (c) while being formatted
   for (d).
2. **Is `design/` shipped, linked, or private?** Right now it is in the public
   repository but not in the package, and the vignettes cite it by path. If a
   reader cannot follow those citations, they are noise; if they can, several
   design documents are effectively user-facing and should be written that way.
3. **What is the entry point?** `README.md` and `intro` currently compete.
4. **What gets deleted?** My candidates: `METHODS_VIGNETTE_SPEC.md` (stale),
   the superseded version history in `PROTOTYPE_SPEC.md` (move to `NEWS.md` or
   drop), and the duplicated AVATAR-vs-DP argument in at least one of the three
   places it now appears.
5. **What moves from `design/` to `vignettes/`?** `PRIVACY_BACKGROUND.md` and
   parts of `FEASIBILITY.md` are the strongest candidates: both are written to
   teach, and both answer questions a user will actually ask.

## 6. Constraints worth remembering

- The package renames to `synpmx` (`TODO.md` item 2), which changes all five
  vignette filenames and their `\%\VignetteIndexEntry` values. **Decide the
  vignette set before renaming**, so the filenames are only churned once.
- `R CMD check` rebuilds every vignette, so each one is a maintenance cost paid
  on every behavioral change, not just a document.
- `AGENTS.md` requires that vignette prose be verified against the code rather
  than preserved for diff-minimization. Any document we keep, we are agreeing
  to keep true.
