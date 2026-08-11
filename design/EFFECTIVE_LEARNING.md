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

**This week, before David.** Write the one-page poster outline from memory.
Read `Guillaudeux23.pdf` and `Destere26.pdf`. Drain the first two items in
`_TODO_owner.md`: `synthetic-data-checks.Rmd` and `demo.Rmd`. Start the
glossary. Go into that meeting able to state the privacy claim and its limits
without reading from a document.

**Ask David for one specific thing:** to try to break the claim that an avatar
cannot be traced to a real patient. Do not ask whether the work is good. A
question that can only be answered yes is not feedback.

**Aug 18 to Aug 31.** The real-data run (pit, eci) with the checks applied, and
the three predict-then-read code spot checks. This is the block where the two
stacks are exercised together.

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

## Teaching something you have not learned yet

The objection: you cannot explain the tool until you know it, you cannot know it
until it is built, and it is not built yet. So the teaching step keeps sliding
into the future and the reading has no focus.

The dependency runs the other way. The poster's claim decides what you have to
know, and the claim is already written. `communications/2026-acop-abstract.md`
asserts four things:

1. Structure can be separated from sensitive patient-level values.
2. What survives is enough to develop ingestion, exploratory analysis, NLME
   specification, and Shiny simulation code.
3. Code developed on synthetic data transfers into the secure environment and
   runs on the real data.
4. Confidentiality is preserved well enough for the synthetic data to leave.

That is the syllabus. It is not 27,000 lines. Claims 1 to 3 are about the
workflow and you already know them from doing the work; claim 4 is the only one
that requires study, and it is the only one an audience will attack. The
question from the floor will be "how do you know no patient is identifiable,"
and everything you need to answer it is disclosure-risk measurement.

You also already have something to explain. Read those four claims and see which
you could state, unprompted, to David. That test takes ten minutes and it
replaces the feeling of not being ready with a list.

**Write the poster outline this week, not in September.** One page, claims and
the evidence for each, from memory. The blanks are the reading list, generated by
you rather than assigned. David then becomes the rehearsal rather than a status
update.

## Reading list

Ordered. Four items, and the first two are the ones that matter.

1. **`references/Guillaudeux23.pdf`** and **`references/Destere26.pdf`**. Both
   are already in the repository. You are shipping this method under your name
   and have written to one of the authors. Reading these is not optional and
   nothing else on the list outranks them.
2. **The disclosure-risk part of `synthetic-data-checking-review.Rmd` only.**
   Adversarial accuracy, the training-versus-control split, RepU and DiSCO. Skip
   the utility and distribution material; claim 4 does not depend on it.
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
of *how*, and the poster's claim 4 needs the second kind. Predict-then-read is
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

Independent feedback has to come from outside. Three candidates:

- **David, mid-August.** The nearest thing to a mentor you have. Use the meeting
  as an adversarial review of the privacy claim, not a status update.
- **The ACoP audience, September.** An audience that can ask questions is the
  only F that cannot be gamed.
- **The methods authors.** Destere was written to; no response as of
  2026-08-11. Guillaudeux and colleagues on AVATAR and the synthpop group
  (Nowok, Raab, Dibben) have not been. A package that implements someone's
  method and measures it on real clinical data is a legitimate reason to write,
  and a poster with results in it is a second reason to follow up.

## Log

One row per comprehension block. Keep it short; the value is in the middle
column.

| Date | Could not explain unaided | Closed by |
|---|---|---|
| | | |
