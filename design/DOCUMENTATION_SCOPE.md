# Documentation scope

Working document for the documentation-reorganization decision in
`design/TODO.md` item 2. Created 2026-07-23.

**Status: an inventory with guessed audiences, not a decision.** The audience
line under each document is my inference from its content and tone, not
something anyone has stated. Correct it before we use it to decide anything —
a wrong guess here produces a wrong merge later.

The package carries 23 documents. Far more of that writing is internal design
record than user-facing explanation, and that imbalance is the first thing worth
discussing.

---

## 1. User-facing: vignettes

Shipped with the package; visible to anyone who installs it.

### `synpmx-intro.Rmd`

- **Audience:** a pharmacometrician who has just heard of the package and does
  not yet know what "differentially private" buys them.
- **Job:** the big picture; all four modes run on `theo_md`; the properties
  table; the environment-to-mode table.
- **Concern:** new (2026-07-23). It overlaps `demo` on theophylline and
  `privacy-intro` on the decision rule. That is deliberate for an entry point,
  but the boundary needs agreeing.

### `synpmx-demo.Rmd`

- **Audience:** a pharmacometrician who has chosen a mode and wants to run it on
  their own study.
- **Job:** the worked AVATAR workflow over five public datasets with structural
  checks, plus one model-based example.
- **Concern:** a large fraction of it is plotting helper functions in the setup
  chunk that no reader needs to see.

### `synpmx-simulation-method.Rmd`

- **Audience:** someone who must defend, debug, or review the generator — a
  methods-minded colleague or an internal reviewer.
- **Job:** the step-by-step AVATAR algorithm with the mathematics, edge cases,
  and a worked example; the alternatives at the end.
- **Concern:** the longest document in the package. Is it a vignette or a
  specification? It currently reads like both.

### `synpmx-privacy-intro.Rmd`

- **Audience:** someone deciding whether a release is allowed to leave the
  environment — possibly not a modeler at all.
- **Job:** what AVATAR is, what differential privacy is, the formal definition,
  and the trust-boundary rule.
- **Concern:** its audience may be the one least likely to open an R vignette.

### `synpmx-epsilon-exploration.Rmd`

- **Audience:** someone who has already committed to a private mode and is
  choosing epsilon.
- **Job:** `f = d / (epsilon N)`, the measured frontier table, `pmx_preflight()`.
- **Concern:** the shortest by a wide margin. Is it a vignette, or a section of
  `privacy-intro`?

**Open questions.** Is five the right number? Does the method vignette want to
become a specification appendix instead? Does `epsilon-exploration` survive as
its own document now that `intro` carries the decision rule?

---

## 2. User-facing: repository front matter

### `README.md`

- **Audience:** someone browsing the GitHub repository who has not installed
  anything.
- **Job:** currently four jobs at once — the pitch, the privacy contract, an API
  reference, a limitations list, and the documentation map.
- **Concern:** **this is the one that needs the rewrite.** The API-reference and
  limitations material is deep internal detail on a page whose first job is to
  say what this is and why to care.

### `NEWS.md`

- **Audience:** an existing user upgrading.
- **Job:** version history by feature.
- **Concern:** reasonable as-is, though written for users who do not yet exist.

### `scripts/README.md`

- **Audience:** a developer running the public demonstrations.
- **Job:** describes `demo_nlmixr2data.R` and `evaluate_simulations.R`.
- **Concern:** stale. It refers to `test-avatar.R`, which is not in `scripts/`.

### `scripts_private/README.md`

- **Audience:** whoever runs this inside the safe environment — realistically
  you.
- **Job:** the trust-boundary rules for the private scripts folder and what may
  be committed.
- **Concern:** longer than the public scripts README, and probably the most
  operationally important page in the repository.

### `references/README.md`

- **Audience:** a contributor adding a paper.
- **Job:** where to put references and what not to commit.
- **Concern:** none. It is fine.

---

## 3. Contributor-facing: repository conventions

### `AGENTS.md`

- **Audience:** an AI coding agent, and secondarily a human contributor.
- **Job:** repository conventions — where code goes, documentation-
  synchronization rules, the acronym rule, testing discipline.
- **Concern:** it doubles as the de facto contributor guide. There is no
  `CONTRIBUTING.md`.

### `CLAUDE.md`

- **Audience:** Claude Code specifically.
- **Job:** points at `AGENTS.md`.
- **Concern:** none. Correct and minimal.

---

## 4. Internal: design record

Not shipped with the package. This is where most of the writing lives.

### `PROTOTYPE_SPEC.md`

- **Audience:** future you, and any reviewer asking "what is this supposed to
  do?"
- **Job:** the implementation contract, with four versions of scope history
  newest-first.
- **Concern:** the longest document in the repository, and it carries three
  superseded designs alongside the current one.

### `FEASIBILITY.md`

- **Audience:** someone asking "will this work at N = 20?"
- **Job:** the measured evidence behind every scope decision; private-mode
  utility by cohort size.
- **Concern:** the strongest document in the set. It arguably deserves to be
  user-facing, at least in part.

### `TEST_SIM.md`

- **Audience:** whoever is fixing a simulation defect.
- **Job:** the living evaluation specification — dataset registry, `SIM-###`
  issues, metrics, gates.
- **Concern:** none. A healthy working document.

### `REVIEW_BACKLOG.md`

- **Audience:** whoever is fixing a mechanism or API defect.
- **Job:** `REV-###` findings from code review, with status.
- **Concern:** none. A healthy working document.

### `MODEL_ELICITATION.md`

- **Audience:** someone preparing the public inputs for a private fit.
- **Job:** the interview producing a public structural model and priors without
  touching data.
- **Concern:** relevant only to modes 2–4. Would a user ever find it?

### `DATA_ELICITATION.md`

- **Audience:** the same person, for trial design.
- **Job:** the complexity ladder, and which parts of a protocol are genuinely
  public.
- **Concern:** same question as above.

### `PRIVACY_BACKGROUND.md`

- **Audience:** someone who wants the arithmetic to make sense.
- **Job:** `d`, `f`, sensitivity, the error law, worked examples.
- **Concern:** tutorial in tone. **This reads like a vignette that ended up in
  `design/`.**

### `TODO.md`

- **Audience:** you and me, this week.
- **Job:** the working queue, and the index to every other design document.
- **Concern:** it is currently the entry point to `design/`, and that works.

### `PRIVACY_ARGUMENT.md`

- **Audience:** a privacy reviewer who will not read R code.
- **Job:** the formal mechanism-level argument for the v2 engine.
- **Concern:** explicitly labeled as not independently reviewed, and scoped to
  the superseded v2 path.

### `METHOD_DISCUSSION.md`

- **Audience:** someone asking "why is AVATAR the default?"
- **Job:** the AVATAR-versus-DP essay and the `synadam` parity argument.
- **Concern:** its content is now partly duplicated across three vignettes.

### `METHODS_VIGNETTE_SPEC.md`

- **Audience:** whoever maintains the vignettes.
- **Job:** vignette maintenance notes.
- **Concern:** **stale.** It describes "Version 2 has four documents"; there are
  now five, in a different structure. A deletion candidate.

### `DOCUMENTATION_SCOPE.md`

- **Audience:** you and me.
- **Job:** this inventory.
- **Concern:** delete it once the reorganization is decided and executed.

---

## 5. What I would want to decide, in order

1. **Who is the actual audience?** Everything above is a guess. The plausible
   set is: (a) you, (b) pharmacometrician colleagues inside the company, (c) an
   internal privacy or governance reviewer, (d) the open-source R
   pharmacometrics community. These want very different documents, and the
   current set reads as though it is serving (a) and (c) while being formatted
   for (d).
2. **Is `design/` shipped, linked, or private?** Right now it is in the public
   repository but not in the package, and the vignettes cite it by path. If a
   reader cannot follow those citations, they are noise; if they can, several
   design documents are effectively user-facing and should be written that way.
3. **What is the entry point?** `README.md` and `intro` currently compete.
4. **What gets deleted?** My candidates: `METHODS_VIGNETTE_SPEC.md` (stale), the
   superseded version history in `PROTOTYPE_SPEC.md` (move to `NEWS.md` or
   drop), and the duplicated AVATAR-versus-DP argument in at least one of the
   three places it now appears.
5. **What moves from `design/` to `vignettes/`?** `PRIVACY_BACKGROUND.md` and
   parts of `FEASIBILITY.md` are the strongest candidates: both are written to
   teach, and both answer questions a user will actually ask.

---

## 6. Constraints worth remembering

- `R CMD check` rebuilds every vignette, so each one is a maintenance cost paid
  on every behavioral change, not just a document.
- `AGENTS.md` requires that vignette prose be verified against the code rather
  than preserved to minimize a diff. Any document we keep, we are agreeing to
  keep true.
- The package rename to `synpmx` is done, so the vignette filenames are settled
  for now; changing the vignette *set* is the remaining churn.
