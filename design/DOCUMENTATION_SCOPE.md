# Documentation scope — decision record

Decisions for the documentation reorganization in `design/TODO.md` item 2.
Inventory created 2026-07-23; decided 2026-07-23 after a scoping conversation.

**Status: decided. Execution in progress.** Delete this document once the
reorganization is complete and `design/TODO.md` item 2 is ticked.

---

## 1. The problem this solves

The package carried 23 documents: five vignettes, a `README.md` doing five jobs
at once, and twelve internal design documents. Far more of the writing was
internal design record than user-facing explanation, the README and the
introduction vignette competed to be the entry point, and the vignettes cited
`design/` by path in 16 places — paths a reader who installed the package
cannot follow.

Three costs drove the decision:

- **`R CMD check` rebuilds every vignette.** Each shipped vignette is a
  maintenance cost paid on every behavioral change, not just a document. The
  vignette set was 2334 lines, 1119 of them in one document.
- **`AGENTS.md` requires vignette prose be verified against the code.** Anything
  kept, we are agreeing to keep true.
- **Citations into `design/` are dead ends** for anyone who is not standing in
  a clone of the repository.

## 2. Audience

Three audiences, in the order they arrive:

- **You**, and pharmacometrician colleagues inside the company — served by
  `design/` plus the website.
- **A privacy or governance reviewer** — will not read R code and is not
  expected to go through the material in detail. Served by a link to the
  website. This is the first arrangement that actually reaches them; an `.Rmd`
  vignette never did.
- **The open-source R pharmacometrics community**, ultimately — served by the
  README plus three vignettes.

## 3. The mechanisms that made this possible

Two standard R facilities do most of the work, both of which sidestep the
`R CMD check` bottleneck:

- **`README.md`, written directly.** Tables and figures on the GitHub landing
  page, and GitHub renders LaTeX math in markdown. Superseded 2026-07-24: the
  `README.Rmd` that used to knit this file is gone, so that there is one README
  and no question which to edit. Its example is now pinned by
  `tests/testthat/test-readme.R` instead of by re-execution.
- **`vignettes/articles/`.** Excluded from the build via `.Rbuildignore`, so
  `R CMD check` never touches these and they are not shipped in the tarball —
  but pkgdown renders them into the website as Articles. Full `.Rmd`, zero check
  cost.

pkgdown was adopted. The one thing that usually blocks it — an exported
function with no documented topic fails the reference index — does not apply
here: 31 exports against 33 man pages. Note that pkgdown *executes* article
code, so a broken article fails the site build; that is a separate signal from
`R CMD check`, which is the decoupling that buys the cheap article space.

## 4. The decided structure

### Entry point

**`README.md`.** The pitch, installation, one minimal
example, a short tour naming the four modes, and the documentation map. High
level throughout. The API reference and the limitations list come out — they
were deep internal detail on a page whose first job is to say what this is and
why to care.

### Vignettes — shipped, rebuilt by `R CMD check`, kept small

| Vignette | Job |
|---|---|
| `synpmx-method.Rmd` | All four modes at a high level: what each does and when to use it. Renamed from `synpmx-simulation-method.Rmd`. |
| `synpmx-demo.Rmd` | The worked workflow over the public datasets. Plot helpers hidden behind code folding. |
| `synpmx-privacy.Rmd` | What differential privacy guarantees, the trust-boundary decision rule, and choosing epsilon. Absorbs `synpmx-epsilon-exploration.Rmd`. |

### Articles — pkgdown only, no check cost

| Article | Source |
|---|---|
| `avatar-mathematics.Rmd` | The deep AVATAR algorithm, mathematics, and edge cases, lifted out of the old method vignette. |
| `privacy-background.Rmd` | `design/PRIVACY_BACKGROUND.md` — `d`, `f`, sensitivity, the error law, worked examples. |
| `feasibility.Rmd` | `design/FEASIBILITY.md` — measured utility by cohort size. |
| `privacy-argument.Rmd` | `design/PRIVACY_ARGUMENT.md` — the formal mechanism-level argument, for a reviewer. |
| `model-elicitation.Rmd` | `design/MODEL_ELICITATION.md` |
| `data-elicitation.Rmd` | `design/DATA_ELICITATION.md` |

### `design/` — internal record, cited by nothing shipped

`TODO.md`, `REVIEW_BACKLOG.md`, `TEST_SIM.md`, `PROTOTYPE_SPEC.md`,
`METHOD_DISCUSSION.md`.

## 5. Rationale for the harder calls

**Why the method vignette covers all four modes rather than only AVATAR.**
The first instinct was to narrow it to AVATAR alone. Covering all four instead
makes one document the canonical method reference rather than scattering
prior-only, calibration, and empirical across the privacy and demo vignettes.
The size problem this would otherwise create is solved by moving the deep
mathematics to `avatar-mathematics.Rmd`: the vignette stays high level and
short, and stops reading as half-vignette, half-specification.

This also absorbs part of the old "alternatives" section. Once the vignette
compares all four modes, *why AVATAR is the default* is in scope there. What
stays internal in `METHOD_DISCUSSION.md` is the `synadam` parity argument.

**Why teaching material is not uniformly moved into vignettes.** There is a
real line:

- *Teaching required to use the package correctly* → vignette. You cannot
  choose an epsilon without understanding `f = d / (epsilon N)`. That earns its
  place and its rebuild cost.
- *Evidence for the package's claims* → article. `FEASIBILITY.md` is the
  strongest document in the set, but it is justification, not instruction. A
  user does not read it to do a task; a reviewer reads it to believe the claims.

**Why the elicitation documents moved out of `design/`.** Their stated audience
is "someone preparing the public inputs for a private fit" — a user doing a
task, which is the definition of an article, not a design record. The open
question in the inventory was "would a user ever find it?" Where they sat, no.

**Why the introduction vignette goes away.** Its job was the four-mode tour on
`theo_md`. That tour now lives in `synpmx-method.Rmd`, and the README carries
the pitch and one runnable example. Keeping both meant README and intro
duplicating each other, which was the original complaint.

**Why `PRIVACY_BACKGROUND.md` became an article rather than folding into the
privacy vignette.** It is 250 lines of arithmetic and worked examples — deep
math by the same standard that moved the AVATAR mathematics out. The privacy
vignette carries the decision rule and enough arithmetic to choose an epsilon,
and links to the article for the full derivation.

## 6. Deletions

- `design/METHODS_VIGNETTE_SPEC.md` — stale. Describes "Version 2 has four
  documents"; there were five, in a different structure.
- `scripts/README.md` — stale. Refers to `test-avatar.R`, which is not in
  `scripts/`.
- `vignettes/synpmx-intro.Rmd` — superseded by the README and
  `synpmx-method.Rmd`.
- `vignettes/synpmx-epsilon-exploration.Rmd` — merged into the privacy vignette.
- `NEWS.md` — emptied to a stub rather than removed. The package has no users
  yet, so the version history is written for people who do not exist; the file
  is kept because recreating it at first release costs more than keeping it.

**`scripts_private/README.md` is kept**, against the initial call to delete it.
It holds the trust-boundary rules for what may be committed out of the safe
environment — the thing standing between this repository and committed
patient-level data. The inventory itself called it "probably the most
operationally important page in the repository." `.Rbuildignore` already
excludes `scripts_private`, so it ships to no one. If the goal is that it not be
public, that is an argument about the git remote, not about deleting the rules.

## 7. Execution order

Relocation precedes citation removal, or the 16 `design/` links become dead
ends mid-flight.

1. Set up pkgdown: `_pkgdown.yml`, GitHub Action, `.Rbuildignore` entry for
   `vignettes/articles`.
2. Move to articles: feasibility, privacy argument, privacy background, both
   elicitation documents.
3. Split the deep AVATAR mathematics out of the method vignette into
   `avatar-mathematics.Rmd`.
4. Rename `synpmx-simulation-method.Rmd` → `synpmx-method.Rmd` and restructure
   to the four modes at a high level.
5. Merge `synpmx-epsilon-exploration.Rmd` into `synpmx-privacy.Rmd`.
6. Write `README.Rmd`, knit `README.md`, delete `synpmx-intro.Rmd`.
7. Repoint every citation at an article; verify nothing shipped points into
   `design/`.
8. Apply the deletions in section 6.
9. `./build.sh`.
