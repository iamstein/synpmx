## Develop Theory for How Andy Learns and Tailor documents for him

I often ask Claude to create good explanatory documents: Readme, literature reviews, etc., as in the vignettes for this project.  

Usually, the first few drafts, I do not fully understand, I have a lot of questions.  There may be technical material I don't know.  Or just definitions.  And it requires a lot of back and forth with Claude to get something that answers my question.

I'd like the AI here to start developing a theory of mind for me around how to provide good explanations in vignettes and articles that will be accessible by humans.  To do that, start by using this document.  Keep a structured set of observations based on the questions I ask and the documents I create.  

Start by looking at the existing vignettes and the literature review article.  

In the coming days, I will be focusing in particular on literature-review.md and synthetic-data-checking.md, so the questions I ask and changes in August-September 2026 may be of particular value.

**Companion file: `design/WRITING_STYLE.md`.** This file records how Andy reads
and what he asks for. That one turns the same evidence into rules for writing —
headings, transitions, pronouns, punctuation, and the AI tics to keep out of a
draft. Update both when an edit teaches something new.

## Structured Observations.

*First pass written 2026-08-05 by Claude, after the session that added the
checking-literature tutorial to `literature-review.Rmd` and the scorecard to
`synthetic-data-checks.Rmd` / `demo.Rmd`.*

*(`literature-review.Rmd` was split on 2026-08-11 into
`synthetic-data-generation-review.Rmd` and `synthetic-data-checking-review.Rmd`;
references to it below are to the document as it stood when they were written.)*

**Standing caveat on everything below.** Andy has **not yet reviewed** the
second half of `literature-review.Rmd` or most of
`synthetic-data-checks.Rmd` himself. Observations that cite those documents are
evidence about *what he asked for*, not evidence that the result worked. The
strongest evidence in this file is instead his **own edits** to Claude prose,
because those are decisions rather than requests. Re-read this file after he has
reviewed those two documents and mark what turned out to be wrong.

### Where the evidence comes from

Three sources, in decreasing order of reliability.

1. **His own commits editing Claude's prose.** `f8c22f1 "simplify readme"`
   (−123/+31 lines) and the README rewrite in `0edffdc`. These are revealed
   preference and outrank anything he says about his preferences.
2. **The verbatim questions in the 2026-08-04/05 session.** A dense list of
   "I don't know X" / "I'm not sure I get Y" / "I don't care about Z".
3. **Rules he has written into `AGENTS.md`** — especially the acronym rule and
   the README-ownership rule. Each of those is a scar from earlier friction.

---

### O1 — He asks for the mechanism, not the citation

**Evidence.** "I don't know the measures by Raab/Nowok/Dibben in 2024 that
shipped in synthpop. I don't know what Linability is. Or the WP criterion."
"I'm not sure what pMSE is." "What is authenticity?" In every one of those
cases the draft had already given the name, the authors, the year, and a
one-sentence gloss. That was not enough.

**What to do.** A named method is not explained until the reader could either
compute it or recognize a bad value. The working test: *what is calculated, from
what inputs, and what does a good number look like?* The pMSE passage that
survived is four numbered mechanical steps; the version that failed was
"the standard measure is pMSE (Snoke et al. 2018)".

### O2 — Implementing a method is not understanding it

**Evidence.** "I haven't really learned about adversarial accuracy yet (even
though you implemented it)." `compare_pmx_proximity()` had been in the package,
documented and tested, for weeks.

**What to do.** Never assume the package's own statistics are understood because
the package computes them. Explain in-house methods with the same care as
external ones. The parenthetical also reads as slightly self-deprecating, so
explain without any framing that implies he should have known already.

### O3 — Telling him he has reinvented something is welcome, not a criticism

**Evidence.** "I see I'm inventing something rather than learning from what's
been done before." This was his reaction to being shown that his hand-derived B5
checks were most of `synthpop`'s RepU/DiSCO. The tone is relief, and the
follow-up was to ask for *more* prior art, not less.

**What to do.** When a home-grown idea has a published name, say so early and
plainly: "the thing you already do is called X." Treat locating prior art as a
deliverable in its own right, not as a politeness. This is probably the single
highest-value thing Claude can do for him on a literature document.

### O4 — "I don't care about X" and "I don't get X" arrive in one sentence and need opposite responses

**Evidence.** One message contained "I honestly don't care that much about the
distributions" (a **scope** decision), "I'm not sure what pMSE is" (an
**explanation** request), and "Regarding coverage, I'm not sure I get what
you're talking about. On did we lose the tails, is the idea we want to?"
(**both** — explain the concept, then rule on whether it matters).

**What to do.** Separate the two explicitly. A scope statement should change what
the document *claims*; an explanation request should change what the document
*teaches*. Conflating them produces the worst outcome — cutting a section he
didn't understand but did want, or explaining at length something he had already
ruled out.

### O5 — He wants scope stated negatively, and stated loudly

**Evidence.** "Maybe I want to be explicit about this is not about scientific
discovery." His own README rewrite replaced a four-row table of hedged verdicts
with four bolded lines carrying ✅ / ❌ markers, two of which are flat "No"s.

**What to do.** Every explanatory document should carry an explicit
**out-of-scope** statement, and it should be blunt rather than qualified. The
"why almost none of this is `synpmx`'s problem" subsection exists because of
this observation; check whether it survives his review.

### O6 — He converts hedged verdicts into binary ones

[Andrew Stein comment: I do not agree with this one.  ACtually, on reading this, I think I made a mistake.  I like having the ❌ but believe I should change **No** back to **Only with care** and I do that now.]

**Evidence.** Claude wrote: "**Only with care.** This needs a formal guarantee;
see the privacy modes below, which are illustrative rather than audited." Andy
replaced it with: "**❌ Send data past a trust boundary (No).** ... this package
should not be used."

**What to do.** Lead with the verdict, then the nuance. Never lead with the
nuance. The nuance is not deleted in his version — "the formal privacy-protecting
methods are illustrative, but not audited" is still there — it just comes second
and no longer softens the headline. This is the clearest stylistic finding in the
file, because it is a direct rewrite of the same content.

### O7 — He deletes exhaustive detail from entry-point documents

**Evidence.** `f8c22f1` cut a fully annotated 13-argument `synpmx_avatar()` call
down to four arguments, and removed the multi-line hanging comments explaining
`dvid` and `dose_covariate`. He also deleted the conceptual coinage "The
declaration is also the **manifest of what survives**" in favour of the literal
"The function drops every column that is not described."

**What to do.** At the entry point, show the minimum that works. Detail belongs
one document deeper. And prefer literal description to coined framing — the
metaphors Claude reaches for to make something memorable are what he deletes
first. (Consistent with the existing memory note "Thin templates and docs".)

### O8 — He is an expert in one half of every document and a novice in the other

**Evidence.** He needed no explanation of PK, PD, BLOQ, trough samples, mg/kg
dosing, arms, or dropout. He needed full definitions of linkability, WP29, pMSE,
authenticity, local cloaking, and adversarial accuracy. The line falls exactly at
the pharmacometrics / statistical-disclosure-control boundary.

**What to do.** Make that boundary conscious rather than accidental. In any
document that crosses it, the pharmacometric side can be terse to the point of
shorthand and the privacy/statistics side must be taught from zero. The
`AGENTS.md` acronym rule is a blunt instrument aimed at this problem; the sharper
version is **expand every acronym, and fully explain every named method that
comes from outside pharmacometrics.**

This is the *expertise reversal effect* (see the reading below): the scaffolding
that helps him in the unfamiliar half actively annoys him in the familiar half.
Uniform explanation depth is wrong in both directions at once.

### O9 — "Where do I start?" is a literal request

**Evidence.** "I don't know these different methods you described. I'm not sure
where to start."

**What to do.** Any survey should end with a ranked entry path: read these N
things, in this order, and one clause on why each. The four-item list closing the
literature review exists for this reason. A flat alphabetical reference list does
not answer the question he actually asked.

### O10 — He optimizes for fewer review rounds, which fights O7

**Evidence.** "I'd first like to make them more complete so that I iterate less."
He will accept a longer first draft if it reduces the number of passes.

**Resolution of the conflict.** The tolerance for length is a function of the
document *tier*, not of his mood. Internal `design/` documents and proposals
should be exhaustive — that is where "iterate less" applies. Shipped entry points
(`README.md`, and the top of each vignette) must be thin — that is where O7
applies. Getting this backwards is what produced the README he had to cut in
half.

### O11 — A mangled term is a reliable signal that the concept has not landed

**Evidence.** "local cleaking", "Linability", "the WP criterion",
"synpmx_validate" (the function is `validate_pmx()`).

**What to do.** Do not silently correct and move on. A term he has half-absorbed
is one he met once and has not yet used. Also note the `synpmx_validate` slip: he
holds the *concept* ("is the dataset even validated") and reconstructs the name
from it. So documents should lead with what a function does and let the name
follow, not the other way round.

### O12 — He reviews by running things and by looking at output

**Evidence.** `_TODO_owner.md`: "Try out on real data and apply checks (pit,
eci), see if I can follow all steps." The render-and-open-in-browser loop he set
up for `avatar-algorithm.Rmd`. And the entire thesis of the checks vignette —
every defect was found by looking at output, not by reasoning about the
algorithm.

**What to do.** Prefer a document that *computes* its claims on a real dataset
over one that asserts them. The scorecard he accepted immediately has exactly
this shape: a static table of claims in one document, and a runnable version that
fills it in on a real run in another. That pairing is worth reusing.

### O13 — When more than about four comparable things are enumerated, he wants a table

**Evidence.** "I like your idea of the scorecard table near the top" — accepted
without modification, unlike almost everything else in that message. He then
immediately extended it ("maybe one part of the scorecard is just
synpmx_validate"), which is the response of someone who has understood a
structure well enough to add to it.

**What to do.** Use a table when the items answer the same question, and make
every row answer it in the same form. Prose comparison of five things is a
reliable way to lose him.

### O14 — Unresolved: the long narrative sections drew no questions

**Evidence.** He asked nothing at all about sections E ("what these checks cannot
tell you") and F ("check the output, not the algorithm") of the checks vignette —
the two most essay-like, least tabular parts of the corpus.

**Two readings, and this file cannot distinguish them.** Either they worked, or
he skipped them. This matters a lot for how much narrative future documents
should carry. **Ask him directly.**

### O15 — One document, one question. A document with "two halves" is two documents

**Evidence.** 2026-08-11, unprompted: *"I think the literature review should be
split into two files"* — generation into one article, checking into another,
both under **Background** in the navbar. The article had explicitly announced
itself as having two halves since it was written, and it survived one review
pass in that form before he cut it. In the same working-tree pass he had already
renamed its two top headings from "How synthetic data is made" / "How synthetic
data is checked" to the plainer **"Algorithms for Generating Synthetic Data"** /
**"Checking Synthetic Data"**, and deleted a hedging clause about the claim
being "about the applications cited here rather than about the whole
machine-learning literature."

**What to do.** Treat "this article has two halves" as a defect report Claude
wrote about its own draft. When an outline needs that sentence, propose two
documents instead — navigation is cheap, and a reader arriving for the checking
tutorial should not scroll through the generation survey. Related to O7 (thin
entry points) and O13 (tables over enumeration): the shared instinct is that
**each artifact should do one thing, at the smallest size that does it**. Note
also the heading edits: he prefers a heading that names the subject flatly over
one phrased as a narrative beat, which is O7's "literal over coined" applied to
structure.

---

### Open questions to put to him

Things that cannot be inferred from the evidence available, roughly in order of
how much they would change future drafts.

1. **Did you read sections E and F of the checks vignette, and did they land?**
   (Resolves O14 — that is, whether argumentative prose earns its length.)
2. **Is the mathematics helpful or is it noise?** The literature review now
   carries the adversarial-accuracy formula in display math. It could as easily
   be three sentences of English.
3. **Inline definitions or a glossary?** Right now every term is defined inline
   at first use, which is why that section is long.
4. **What is the length ceiling for a shipped vignette?** `avatar-algorithm.Rmd`
   is ~2,100 lines. Is that already past the point of usefulness, or is it fine
   because it is reference material nobody reads end to end?
5. **Does the ✅/❌ marker style from your README rewrite generalize?** If you
   want it in the vignettes too, that is a cheap and consistent change.

---

### What the outside literature says about this

Browsed 2026-08-05. Three findings map onto the observations above closely
enough to be worth naming.

**Mixing documentation modes is the most common cause of confusing docs.** The
[Diátaxis framework](https://diataxis.fr/) splits documentation into four kinds
serving four different needs — *tutorial* (learning-oriented), *how-to*
(goal-oriented), *reference* (information-oriented), and *explanation*
(understanding-oriented) — and its central claim is that a page which tries to
be more than one of them serves none of its readers.

This diagnoses O9 precisely. The literature review was written in **explanation**
mode: a survey that assumes you already have the vocabulary and want the
landscape. Andy needed **tutorial** mode: start from zero, one idea at a time, in
an order chosen for learning rather than for taxonomy. That is why "make it a
tutorial" was the fix, and it suggests a general rule — *check which of the four
modes a document is in before writing it, and say so in its first paragraph.*
The corpus already mostly does this: `demo.Rmd` is a tutorial,
`avatar-algorithm.Rmd` is reference, `synthetic-data-checks.Rmd` is a how-to
guide wrapped around an explanation, and the literature review is now explanation
plus tutorial in two labelled halves.

**The [expertise reversal effect](https://en.wikipedia.org/wiki/Expertise_reversal_effect)**
— instructional support that measurably helps a low-knowledge reader measurably
*hurts* a high-knowledge one, because for the expert it is redundant material
competing for the same working memory. This is O8 with a name and an evidence
base, and it is the strongest argument against uniform explanation depth.

**The [curse of knowledge](https://earthly.dev/blog/curse-of-knowledge/)** in
technical writing is attributed to *fluency misattribution* — the writer finds
the material easy to retrieve and misreads that ease as the reader's. Worth
naming because the failure mode here is not carelessness: the drafts that failed
were carefully written, by a writer for whom "pMSE" retrieves a full definition.

**On the style side**, the commonly catalogued LLM writing pathologies —
[hedging and deferential qualifiers](https://passo.uno/whats-wrong-ai-generated-docs/),
verbosity, and over-elaboration — match Andy's edits almost exactly. What he
deleted from the README was hedging (O6) and elaboration (O7). This is not a
personal idiosyncrasy; it is the general complaint, and it should be treated as a
default to correct rather than a preference to accommodate.

---

## Proposed changes to this document

Suggestions, for Andy to accept or reject. Nothing here has been done.

1. **Split raw material from conclusions.** This file currently mixes them. Add a
   `## Question log` section holding his questions **verbatim and dated**, and
   keep `## Structured Observations` for distilled findings that cite the log.
   Verbatim wording matters — "I'm not sure I get what you're talking about" and
   "I don't care about that" look similar in a summary and mean opposite things
   (O4).

2. **Give every observation a fixed schema.** *Observation / Evidence / What to
   do*, as above. An observation without evidence is a guess, and an observation
   without an action does not change any document.

3. **Add a disagreement marker.** Right now this file is written by Claude about
   Andy, with no way for Andy to push back inside it. Suggest a convention —
   `**[AS: no, actually...]**` inline — so a wrong model gets corrected in place
   rather than re-derived every session. **This is the most valuable change on
   the list**, because most of the file is currently inference from indirect
   evidence.

4. **Add a document contract table.** One row per shipped document: audience,
   Diátaxis mode, target length, and how thin it must be. That turns the theory
   into something checkable before writing rather than a description of what went
   wrong afterwards. It also resolves the O7/O10 conflict mechanically instead of
   by judgement each time.

5. **Sharpen the `AGENTS.md` hook.** It currently reads "Consult the
   design/_THEORY_OF_MIND.md and update." Suggest naming *when*: before drafting
   or substantially revising any vignette, article, or README, and after any
   session in which he asks more than a couple of clarifying questions. As
   written it is easy to satisfy trivially.

6. **Record what was tried and rejected.** If a way of explaining something fails
   review, that is more informative than the version that passed, and it is
   currently lost. Two lines per rejection is enough.

7. **Revisit after the reviews.** The two documents this file has the most to say
   about are the two he has not read yet. Re-run this analysis after he works
   through the two literature-review articles and `synthetic-data-checks.Rmd`; the questions
   he asks then will be better evidence than anything above.