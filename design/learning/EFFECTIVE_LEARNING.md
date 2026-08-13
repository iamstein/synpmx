# Learning plan for synpmx

Applying the A-B-C-D-E-F-G
checklist from <https://sites.google.com/site/andrewsteinphd/effective-learning>
to the work in this repository.

Evidence comes from the git history over 2026-07-17 to 2026-08-11,
`design/_TODO_owner.md`, `design/_THEORY_OF_MIND.md`, and the file sizes of what
has been produced.

## Aspiration

**The stated aspiration and the growth plan are the same object, and that is the
central problem with the current setup.**

1. Develop experience directing AI agents in technical work that would be time consuming do alone.  Learn to tell when the output is wrong.
2. Learn statistical disclosure control well enough to make and defend a privacy
   claim in front of pharmacometricians and those who enforce privacy boundaries.

## Where you are, letter by letter

| | Status |
|---|---|
| A Aspiration | Gain practice with AI coding tools and become conversant in privacy literature |
| B Broad overview | Two review articles |
| C Consult experts | Claude only; no human synthetic data researchers so far |
| D Deconstruct | AVATAR algorithm; Basic privacy concepts and metrics |
| E Routine | Minimum of one learning block and one coding block per week |
| F Feedback | Will come from internal and external presentations |
| G Growth plan | Write `synpmx` package |

## Claude assessment of what has gone well

Three things are working and should not be disturbed.

**The context system is a real AI-direction skill.** `AGENTS.md`, the three
documentation tiers, and the `SIM-###` / `REV-###` registries are how you get
consistent output from an agent across 26 days without repeating yourself. Most
people do not build this. It is the concrete evidence that skill 1 of the
aspiration is developing.

**Your questions are at the mechanism level.** O1 in `_THEORY_OF_MIND.md`
records that a method name, an author, a year, and a one-line gloss are not
accepted as an explanation. Asking "what is computed, from what, and what does a
good value look like" is the correct instinct and it is the reason the pMSE and
adversarial-accuracy passages ended up usable.

**You located the prior art before the poster.** O3 records the reaction to
finding that the hand-derived B5 checks were most of `synthpop`'s RepU and
DiSCO. Finding that in August rather than at ACoP in a question from the floor
is worth more than any section of prose in the repository.

## Problem 1: the package has outrun the reader

The artifact is 27,000 lines across `R/`, `tests/`, and `vignettes/`. Nobody has
read all of it.

The failure mode of learning by building with an agent is that the build does
not slow down when comprehension does. The repository looks healthier every
week either way, so the gap between what the package does and what its author
can explain widens without a signal.

In this repository it shows up in two places. `_THEORY_OF_MIND.md` records that
the second half of the literature review and most of `scorecard-synthetic-data-checks.Rmd`
have not been read by the owner. And "identify a few places to spot check the
code" has not started, so the 11,217 lines in `R/` have had no line-by-line
human read.

**No new vignette or article prose is drafted until the queue in
`_TODO_owner.md` is empty.** Bug fixes, tests, and code changes continue.
Generating explanation does not.

## Problem 2: every feedback loop points at the code

`R CMD check`, the test suite, `TEST_SIM.md`, and `REVIEW_BACKLOG.md` are a
serious feedback on the coding. They measure package progress, but not learning.

`_THEORY_OF_MIND.md` is feedback for Claude, but not feedback for author's learning. 

Because there is no feedback mechanism yet on the author's understanding, the question of how efficiently is the author learning cannot be answered from inside the current setup. You are asking whether the approach is working, and there is no
measurement that could tell you.

The three instruments below are the smallest set that fixes this.

**Instrument 1: explain-back before re-reading.** For each document in the review
queue, before opening it, write five sentences from memory on what it claims.
Then read it. The delta is the measurement. Keep the deltas in the log at the
bottom of this file.

**Instrument 2: `design/GLOSSARY.md`, written by you and not by Claude.** One
line each, in your own words, for every term that arrived from outside
pharmacometrics: pMSE, linkability, WP29, authenticity, local cloaking,
adversarial accuracy, control set, epsilon, sensitivity, RepU, DiSCO. O11
already established that a mangled term is a reliable signal the concept has not
landed, so the terms you cannot write are the syllabus. Writing them is retrieval
practice, which is the one approach from your own page that appears nowhere in
this project.

**Instrument 3: predict-then-read on three functions.** Pick
`compare_pmx_proximity()`, `skeleton_uniqueness()`, and one masking mechanism
from `synpmx_avatar()`. Write what you expect the code to do from the vignette
alone, then read the source. Every disagreement is either a defect or a gap in
your understanding, and both are worth finding. This is the only exercise that
tests both halves of the aspiration at once, and it is already on your own list.

## Problem 3: no teaching, no retrieval

Your page lists five learning approaches. Reading, Organizing, Teaching, Doing, Recalling.

Reading is heavy. Organizing is heavy, arguably too heavy: `design/` holds 12
files, four of which are documentation about documentation
(`WRITING_STYLE.md`, `_THEORY_OF_MIND.md`, `DOCUMENTATION_SCOPE.md`,
`BUILD_DOCUMENTATION.md`) and serve Claude's output quality rather than your
understanding. Doing is partial, since the directing is yours and the doing is
Claude's. Teaching is at zero. Recalling is at zero.

Teaching and recalling are where retention comes from, and you already have two
scheduled teaching events that are being treated as deliverables downstream of
the package rather than as the instruments they are:

- **The internal review meeting, mid-August.** Four days out.
- **The ACoP poster, submission 2026-09-28.** Seven weeks out.

Reframe both. The poster is not a report on the package. It is the Feynman test
for the whole project, with a fixed date and an audience that can ask questions.
Working backwards from it orders everything else.

## D. Deconstruction

**Stack 1, statistical disclosure control.** Deep theory needed in two places
only:

1. *The AVATAR mechanism itself.* Blending, `k` donors, local cloaking, and the
   six masking mechanisms. 

   **Status, 2026-08-11: solid at the mechanism level.** Stated unaided as snap
   to a grid, PCA to reduce dimensionality, weighted blending of patients with
   randomness, unsnap, plus identification checks that are not complete. That
   matches the code: `stats::prcomp` at `R/profiles.R:223`, components retained
   by cumulative variance at `:229`, donors taken nearest-first by Euclidean
   distance in the retained coordinates at `R/synthesis.R:1607`, weights drawn
   at `R/profiles.R:306`. Two details were missing and both matter for the
   privacy claim: it is **one PCA over all endpoints together**, not one per
   endpoint, and there is a **per-donor weight cap of 0.50** with water-filling
   redistribution, which is what stops an avatar from being one real patient
   with noise on top. Remaining work here is practice and delivery, not study.

2. *Disclosure risk measurement.* Adversarial accuracy, uniqueness, RepU,
   DiSCO, the training-versus-control split. This is what an internal review
   could question. Learn it to the point where you can compute it by hand on ten
   rows.

   **Status, 2026-08-11: this is the gap, and it is the whole gap.** "Have I
   checked enough" resists study because it is not a fact to be looked up. It is
   a decision that stays open until two things are fixed: who receives the data,
   and what counts as a failure. Section E of `scorecard-synthetic-data-checks.Rmd`
   already states the second half — a leak cannot be detected from the synthetic
   data alone — which means output checks cannot close the question no matter
   how many are added. The train-versus-control holdout is what closes it,
   because without a baseline for how close a synthetic patient lands to a real
   one by chance, every proximity number is uninterpretable. Build the holdout
   before adding another check.

Recognition level is enough for the rest: formal differential privacy (know
what epsilon buys and that this package's implementation is illustrative), CART
and synthpop's generation side, utility measures such as pMSE. Being able to say
what each one is for, and that synpmx does not do it, covers every use you have.

**Stack 2, directing AI.** Ranked by what currently limits you:

1. *Reviewing and verifying agent output.* The binding constraint. Everything
   else in the repository is gated on it, and it has had the least practice.
2. *Knowing when the agent is confidently wrong.* Related but distinct, and only
   trainable by finding cases. Instrument 3 is the drill.
3. *Context engineering.* Already competent. No further investment needed.
4. *Prompting.* Not the constraint. Ignore.

The 80/20: two topics in stack 1, one skill in stack 2. Everything else in this
repository is background you can read once and look up later.

## E. Routine

The plan is a few hours a week. The pattern is long single-day pushes on
Jul 22, Jul 23 and Aug 3, then five quiet days from Aug 6 to Aug 10, then work
again on Aug 11.

Bursts are fine for building and bad for retention, and they have a specific
side effect here: in a burst session, building beats reading, because building
produces visible output and reading does not. That is the mechanism that leaves
the review queue where it is. It will not fix itself by intending to read more.

Proposed shape, two blocks a week, ninety minutes each:

- **Block 1, comprehension.** Explain-back, then read, then glossary entries.
  No agent-directed generation in this block.
- **Block 2, building.** Whatever `design/TODO.md` says.

If only one block happens in a given week, it is block 1. That ordering is the
whole point of the change.

## Sequence to 2026-09-28

**This week, before the internal review.** Write the one-page poster outline from memory.
Read `Guillaudeux23.pdf` and `Destere26.pdf`. Drain the first two items in
`_TODO_owner.md`: `scorecard-synthetic-data-checks.Rmd` and `demo.Rmd`. Start the
glossary. Go into that meeting able to state the privacy claim and its limits
without reading from a document.

**Ask for one specific thing at that meeting:** to try to break the claim that an avatar cannot be traced to a real patient. Do not ask whether the work is good. A
question that can only be answered yes is not feedback.

**Aug 18 to Aug 31.** Another iteration of the checks on the internal datasets,
and the three predict-then-read code spot checks. This is the block where the
two stacks are exercised together.

Decide the stopping condition for that iteration before starting it. Each pass
over real data has so far taught something and changed the package, which is
productive and has no natural end. The train-versus-control holdout is the end
condition: it is the first result that can come out either way, so it is the
first one that can say the checking is finished rather than merely continuing.

**Sept 1 to Sept 14.** Fill out the August outline into the poster, from the
glossary and from memory, before showing it to Claude. Claude critiques; Claude
does not draft. Anything you cannot write unaided is a gap the poster will
expose in front of an audience, and finding it in this window is the point of
the exercise.

**Sept 15 to Sept 28.** Iterate and submit.

The remaining items on `_TODO_owner.md`, the eight public-dataset examples and
the full avatar-algorithm page, fit in whichever comprehension blocks are free.
The literature review stays off the critical path, where you already put it.

## Should Claude keep writing tutorials

**Yes. Keep the loop exactly as it is for material you do not know yet.**

It is working. The evidence is that the questions moved from "what is this
called" to "what is computed and what does a good value look like," which is the
transition that matters. Claude drafts, the owner reads and asks, Claude
revises. That is a fast B and a reasonable D, and nothing on the reading list
replaces it. Drafting a tutorial on a method you have not learned is not
possible, and asking for one is not a failure of independence.

The rule applies to one artifact only, and an earlier version of this file
stated it too broadly. **A tutorial is input. The poster is output.** The poster
is what an audience will judge, so the poster is the one thing to draft unaided
and hand to Claude for critique rather than for drafting. Vignettes, articles,
and `design/` material carry the owner's name too and are not covered by this;
Claude drafts those.

The one real limit on the tutorial loop is that reading is passive, so it
produces recognition rather than recall. That is what Instruments 1 and 2 (explaining back and writing a glossary) are for, and they cost minutes. The loop itself stays.

## Teaching something you have not learned yet

**Write the poster outline this week, not in September.** One page, claims and
the evidence for each, from memory. The blanks are the reading list, generated by
you rather than assigned. The meeting then becomes the rehearsal rather than a
status update.

[Author note: poster outline is in progress in communication/2026-acop-poster-notes.md]

## Reading list

Ordered. Four items, and the first two are the ones that matter.

1. **`references/Guillaudeux23.pdf`** and **`references/Destere26.pdf`**. Both
   are already in the repository. You are shipping this method under your name.
   Reading these in full is not optional and nothing else on the list outranks them.
2. **The disclosure-risk part of `synthetic-data-checking-review.Rmd` only.**
   Adversarial accuracy, the training-versus-control split, RepU and DiSCO. Skip
   the utility and distribution material;
3. **One synthpop paper on RepU and DiSCO**, because your B5 checks are that and
   naming the prior art is worth more on a poster than the checks themselves.
4. **`R/` for three functions**, by predict-then-read. This is reading too, and
   it is the only item that converts directing the build into knowing the build.

Off the list until after 2026-09-28: the differential-privacy articles,
`feasibility.Rmd`, both elicitation articles, and the generation review beyond a
skim. None of them appear on the poster.

## Directing a build is not the same as knowing it

You are right that you are not building the package, and it is worth being exact
about what you are doing instead. You are writing the specification, making the
scope calls, and deciding when output is wrong. Those are the skills in stack 2
and they are developing.

They produce knowledge of *what* the package does. They do not produce knowledge
of *how*. Predict-then-read is
the only activity on this list that converts one into the other, which is why it
sits above everything except the two papers.

## What Claude cannot do for you

You have named Claude as mentor and asked who else it would be. An honest
account of the limits.

Claude is a good B, a good D, and a usable C for locating prior art. O3 is the
strongest case for that last one.

Claude is a bad F, for two structural reasons. Claude wrote the explanations, so
testing your understanding with the same source that produced the material is
not an independent check. And Claude will tend to agree with you, which means
any self-assessment run through this channel is biased toward passing.

Independent feedback has to come from outside. Two candidates:

- **The internal review, mid-August.** The nearest thing to a mentor you have. Use the meeting
  as an adversarial review of the privacy claim, not a status update.
- **The ACoP audience, September.** An audience that can ask questions is the
  only F that cannot be gamed.

## Log

One entry per comprehension block. Keep it short; the value is in the second
line.

### YYYY-MM-DD

- **Could not explain unaided:**
- **Closed by:**
