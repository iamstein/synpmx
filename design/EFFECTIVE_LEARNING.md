# Learning plan for synpmx

Written 2026-08-11 by Claude, at the owner's request, applying the A-B-C-D-E-F-G
checklist from <https://sites.google.com/site/andrewsteinphd/effective-learning>
to the work in this repository.

Evidence comes from the git history (238 commits, 2026-07-17 to 2026-08-11),
`design/_TODO_owner.md`, `design/_THEORY_OF_MIND.md`, and the file sizes of what
has been produced.

## Aspiration

**The stated aspiration and the growth plan are the same object, and that is the
central problem with the current setup.**

The stated goal is "write the package." Writing the package is G. It is the
long-term application that cements the learning. It is not A.

The actual aspiration, from the surrounding context, is two things that have not
been written down separately:

1. Become able to direct AI agents to do technical work you could not do alone,
   and to tell when the output is wrong.
2. Learn statistical disclosure control well enough to make and defend a privacy
   claim in front of pharmacometricians.

Neither of those is measured by the package existing. When A and G collapse into
one object, progress on the artifact reads as progress on the learner, and the
two can diverge without any signal. They have already diverged, and the size of
the gap is measurable. See "The package has outrun the reader" below.

Restating A separately is the highest-value correction in this document. Every
other recommendation follows from it.

## Where you are, letter by letter

| | Status | Evidence |
|---|---|---|
| A Aspiration | Misidentified | Stated as the package, which is G |
| B Broad overview | Done | Two review articles, 1,263 lines |
| C Consult experts | Partial | Claude only; no human has reviewed the method |
| D Deconstruct | Not done | See the deconstruction below |
| E Routine | Stated but not kept | 45 commits on Jul 23, zero from Aug 6 to Aug 10 |
| F Feedback | Absent for you, strong for the code | `R CMD check`, `SIM-###`, `REV-###` measure the package; nothing measures what you know |
| G Growth plan | Running ahead of the rest | 11,217 lines of R, 6,207 of tests, 10,050 of vignettes |

## What has gone well

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

Of 238 commits, 44 carry your message style and 194 carry Claude's. The
artifact is 27,000 lines across `R/`, `tests/`, and `vignettes/`.

Against that, `_TODO_owner.md` is a six-item review queue in which nothing is
ticked, and `_THEORY_OF_MIND.md` records that the second half of the literature
review and most of `synthetic-data-checks.Rmd` have not been read by you. The
item "identify a few places to spot check the code" has not started, so the
11,217 lines in `R/` have had no human read at all.

This is the specific failure mode of learning by building with an agent. The
build does not slow down when comprehension does, so the gap widens silently and
the repository looks healthier every week. Growth is only real to the extent it
is absorbed, and right now the package knows more than you do.

The correction is a stop rule, not more effort:

**No new vignette or article prose is drafted until the queue in
`_TODO_owner.md` is empty.** Bug fixes, tests, and code changes continue.
Generating explanation does not.

## Problem 2: every feedback loop points at the code

`R CMD check`, the test suite, `TEST_SIM.md`, and `REVIEW_BACKLOG.md` are a
serious F. They measure the package.

`_THEORY_OF_MIND.md` looks like an F for you but runs the other way. It measures
how Claude should write so that you understand it. It tunes the explainer, not
the learner. There is no instrument anywhere in the repository that answers
"what does Andy now know that he did not know in July."

That asymmetry is why the efficiency question cannot be answered from inside the
current setup. You are asking whether the approach is working, and there is no
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

Your page lists five approaches. Reading, Organizing, Teaching, Doing, Recalling.

Reading is heavy. Organizing is heavy, arguably too heavy: `design/` holds 12
files, four of which are documentation about documentation
(`WRITING_STYLE.md`, `_THEORY_OF_MIND.md`, `DOCUMENTATION_SCOPE.md`,
`BUILD_DOCUMENTATION.md`) and serve Claude's output quality rather than your
understanding. Doing is partial, since the directing is yours and the doing is
Claude's. Teaching is at zero. Recalling is at zero.

Teaching and recalling are where retention comes from, and you already have two
scheduled teaching events that are being treated as deliverables downstream of
the package rather than as the instruments they are:

- **David, mid-August.** Four days out.
- **The ACoP poster, submission 2026-09-28.** Seven weeks out.

Reframe both. The poster is not a report on the package. It is the Feynman test
for the whole project, with a fixed date and an audience that can ask questions.
Working backwards from it orders everything else.

## Deconstruction

You have not done D, and the tutorials are not D. A tutorial is B. Deconstruction
is naming the sub-skills and deciding which need theory and which need only
recognition.

There are two stacks, and treating them as one is part of why the work feels
chaotic.

**Stack 1, statistical disclosure control.** Deep theory needed in two places
only:

1. *Disclosure risk measurement.* Adversarial accuracy, uniqueness, RepU,
   DiSCO, the training-versus-control split. This is what the poster claims and
   what David will attack. Learn it to the point where you can compute it by
   hand on ten rows.
2. *The AVATAR mechanism itself.* Blending, `k` donors, local cloaking, and the
   six masking mechanisms. You are shipping it under your name.

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

## Routine

The plan is a few hours a week. The history is 36 commits on Jul 22, 45 on
Jul 23, 29 on Aug 3, then five days of nothing, then 6 on Aug 11.

Bursts are fine for building and bad for retention, and they have a specific
side effect here: in a burst session, building always beats reading, because
building produces visible output and reading does not. That is the mechanism
that produced the unstarted review queue. It will not fix itself by intending to
read more.

Proposed shape, two blocks a week, ninety minutes each:

- **Block 1, comprehension.** Explain-back, then read, then glossary entries.
  No agent-directed generation in this block.
- **Block 2, building.** Whatever `design/TODO.md` says.

If only one block happens in a given week, it is block 1. That ordering is the
whole point of the change.

## Sequence to 2026-09-28

**This week, before David.** Drain the first two items in `_TODO_owner.md`:
`synthetic-data-checks.Rmd` and `demo.Rmd`. Start the glossary. Go into that
meeting able to state the privacy claim and its limits without reading from a
document.

**Ask David for one specific thing:** to try to break the claim that an avatar
cannot be traced to a real patient. Do not ask whether the work is good. A
question that can only be answered yes is not feedback.

**Aug 18 to Aug 31.** The real-data run (pit, eci) with the checks applied, and
the three predict-then-read code spot checks. This is the block where the two
stacks are exercised together.

**Sept 1 to Sept 14.** Write the poster yourself, from the glossary and from
memory, before showing it to Claude. Claude critiques; Claude does not draft.
Anything you cannot write unaided is a gap the poster will expose in front of an
audience, and finding it in this window is the point of the exercise.

**Sept 15 to Sept 28.** Iterate and submit.

The remaining items on `_TODO_owner.md`, the eight public-dataset examples and
the full avatar-algorithm page, fit in whichever comprehension blocks are free.
The literature review stays off the critical path, where you already put it.

## Should Claude keep writing tutorials

**Keep the method, demote it, and reverse who drafts.**

It is working. The evidence is that your questions moved from "what is this
called" to "what is computed and what does a good value look like," which is the
transition that matters. Reading generated explanation is an efficient B and a
reasonable D.

Two limits. It is passive, so it produces recognition rather than recall, which
is why you can read a document and still not be able to state its claim four days
later. And it scales badly against an agent, which is how 10,050 lines of
explanation accumulated faster than they could be read.

So: Claude drafts internal `design/` material and code. For anything with your
name on it, you draft and Claude critiques. That inverts the current
relationship, and it converts the same hours from reading into teaching.

## What Claude cannot do for you

You have named Claude as mentor and asked who else it would be. An honest
account of the limits.

Claude is a good B, a good D, and a usable C for locating prior art. O3 is the
strongest case for that last one.

Claude is a bad F, for two structural reasons. Claude wrote the explanations, so
testing your understanding with the same source that produced the material is
not an independent check. And Claude will tend to agree with you, which means
any self-assessment run through this channel is biased toward passing.

Independent feedback has to come from outside. Three candidates:

- **David, mid-August.** The nearest thing to a mentor you have. Use the meeting
  as an adversarial review of the privacy claim, not a status update.
- **The ACoP audience, September.** An audience that can ask questions is the
  only F that cannot be gamed.
- **The methods authors.** Guillaudeux and colleagues on AVATAR, Destere and
  colleagues on the pop-PK benchmark, and the synthpop group (Nowok, Raab,
  Dibben). A package that implements someone's method and measures it on real
  clinical data is a legitimate reason to write to them. This is the C step your
  own page describes, and it has not been attempted.

## Log

One row per comprehension block. Keep it short; the value is in the middle
column.

| Date | Could not explain unaided | Closed by |
|---|---|---|
| | | |
