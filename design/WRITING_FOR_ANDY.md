# Writing for Andy

Who the reader is, and how to write a document he will not have to send back.

**When to read this.** When asked to draft a new document, or to substantially
rewrite an existing one — a vignette, an article, `README.md`, or a long
summary. Not for code, test, or registry work, and not as a session preamble.
Part 3 opens with the short list of prose habits that applies to *everything*,
including conversation, commit messages and comments; the rest of the file is
calibrated for one reader.

**Who updates it.** Andy, by asking for a review of recent history. Do not
append to it during a session. If something here looks wrong, or contradicted by
what just happened, say so in the conversation and let him decide.

**Status.** Merged 2026-08-14 from `WRITING_STYLE.md` and `_THEORY_OF_MIND.md`.
Observations first written 2026-08-05; rules accumulated through 2026-08-14. The
merge dropped `WRITING_STYLE.md`'s tics list and its search step, and both were
restored to Part 3 on 2026-08-14 at Andy's request.
Sections are numbered in two independent series — `O###` for observations about
the reader, plain numbers for the writing rules — and the gaps in the `O###`
series are deliberate, because those identifiers are cited elsewhere. Do not
renumber to close them.

---

## The original charter, in his words

> *Written 2026-08-05. Kept verbatim; the document names have since changed.*
>
> I often ask Claude to create good explanatory documents: Readme, literature
> reviews, etc., as in the vignettes for this project.
>
> Usually, the first few drafts, I do not fully understand, I have a lot of
> questions. There may be technical material I don't know. Or just definitions.
> And it requires a lot of back and forth with Claude to get something that
> answers my question.
>
> I'd like the AI here to start developing a theory of mind for me around how to
> provide good explanations in vignettes and articles that will be accessible by
> humans. To do that, start by using this document. Keep a structured set of
> observations based on the questions I ask and the documents I create.

---

# Part 1 — The reader

The rules in Part 2 read as arbitrary until you know where they came from. This
part is that.

## Where the evidence comes from

Three sources, in decreasing order of reliability.

1. **His own commits editing Claude's prose.** `f8c22f1 "simplify readme"`
   (−123/+31 lines) and the README rewrite in `0edffdc`. These are revealed
   preference and outrank anything he says about his preferences.
2. **The verbatim questions in the 2026-08-04/05 session.** A dense list of
   "I don't know X" / "I'm not sure I get Y" / "I don't care about Z".
3. **Rules he has written into `AGENTS.md`** — especially the acronym rule and
   the README-ownership rule. Each of those is a scar from earlier friction.

**Standing caveat.** Andy had **not reviewed** the checking-literature article or
most of `scorecard-synthetic-data-checks.Rmd` when the observations below were
written. Anything citing those documents is evidence about *what he asked for*,
not evidence that the result worked. His own edits are the strong evidence
because they are decisions rather than requests. Mark what turns out to be wrong.

## O1 — He asks for the mechanism, not the citation

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

## O2 — Implementing a method is not understanding it

**Evidence.** "I haven't really learned about adversarial accuracy yet (even
though you implemented it)." `compare_pmx_proximity()` had been in the package,
documented and tested, for weeks.

**What to do.** Never assume the package's own statistics are understood because
the package computes them. Explain in-house methods with the same care as
external ones. The parenthetical also reads as slightly self-deprecating, so
explain without any framing that implies he should have known already.

## O3 — Telling him he has reinvented something is welcome, not a criticism

**Evidence.** "I see I'm inventing something rather than learning from what's
been done before." This was his reaction to being shown that his hand-derived B5
checks were most of `synthpop`'s RepU/DiSCO. The tone is relief, and the
follow-up was to ask for *more* prior art, not less.

**What to do.** When a home-grown idea has a published name, say so early and
plainly: "the thing you already do is called X." Treat locating prior art as a
deliverable in its own right, not as a politeness. This is probably the single
highest-value thing Claude can do for him on a literature document.

## O4 — "I don't care about X" and "I don't get X" arrive in one sentence and need opposite responses

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

## O8 — He is an expert in one half of every document and a novice in the other

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

This is the *expertise reversal effect* (Part 4): the scaffolding that helps him
in the unfamiliar half actively annoys him in the familiar half. Uniform
explanation depth is wrong in both directions at once.

## O9 — "Where do I start?" is a literal request

**Evidence.** "I don't know these different methods you described. I'm not sure
where to start."

**What to do.** Any survey should end with a ranked entry path: read these N
things, in this order, and one clause on why each. The four-item list closing the
literature review exists for this reason. A flat alphabetical reference list does
not answer the question he actually asked.

## O10 — He optimizes for fewer review rounds, which fights the thinness rule

**Evidence.** "I'd first like to make them more complete so that I iterate less."
He will accept a longer first draft if it reduces the number of passes.

**Resolution of the conflict.** The tolerance for length is a function of the
document *tier*, not of his mood. Internal `design/` documents and proposals
should be exhaustive — that is where "iterate less" applies. Shipped entry points
(`README.md`, and the top of each vignette) must be thin — that is where rule 12
applies. Getting this backwards is what produced the README he had to cut in
half.

## O11 — A mangled term is a reliable signal that the concept has not landed

**Evidence.** "local cleaking", "Linability", "the WP criterion",
"synpmx_validate" (the function is `validate_pmx()`).

**What to do.** Do not silently correct and move on. A term he has half-absorbed
is one he met once and has not yet used. Also note the `synpmx_validate` slip: he
holds the *concept* ("is the dataset even validated") and reconstructs the name
from it. So documents should lead with what a function does and let the name
follow, not the other way round.

## O12 — He reviews by running things and by looking at output

**Evidence.** `_TODO_owner.md`: "Try out on real data and apply checks (pit,
eci), see if I can follow all steps." The render-and-open-in-browser loop he set
up for `avatar-algorithm.Rmd`. And the entire thesis of the checks vignette —
every defect was found by looking at output, not by reasoning about the
algorithm.

**What to do.** Prefer a document that *computes* its claims on a real dataset
over one that asserts them. The scorecard he accepted immediately has exactly
this shape: a static table of claims in one document, and a runnable version that
fills it in on a real run in another. That pairing is worth reusing.

## O14 — Unresolved: the long narrative sections drew no questions

**Evidence.** He asked nothing at all about sections E ("what these checks cannot
tell you") and F ("check the output, not the algorithm") of the checks vignette —
the two most essay-like, least tabular parts of the corpus.

**Two readings, and this file cannot distinguish them.** Either they worked, or
he skipped them. This matters a lot for how much narrative future documents
should carry. **Ask him directly.**

## O16 — He asks for the design tradeoff to be argued, not just implemented

**Evidence.** 2026-08-13, on the discrete-endpoint defect: *"Implement a check
and a fix for when LIDV is binary or ordinal or categorical. And actually, give
a thought of whether it's reasonable for this to be determined from the data or
whether it should be specified in the pmxroles somehow."* The instruction to
build came first and was unambiguous; the second sentence reopened the interface
question that the first sentence's phrasing had already implied one answer to.
He wanted the alternative weighed before the code existed, not after.

Note also what he did not ask for: he named the symptom in his own terms ("when
LIDV is binary or ordinal or categorical") rather than the ID of the issue that
had just been filed, and he asked for **a check and a fix**, both, in one
sentence. Consistent with O12 (he reviews by running things): a fix he cannot
see fire is half a delivery.

**What to do.** When a task has an inference-versus-declaration fork — or any
comparable interface choice — state the fork and the answer in the reply, with
the reason, in a few sentences. Do not present it as an open question to be
resolved before starting, and do not bury it in a code comment either. The
answer that fit here was *both*: infer by default where the data answers the
question outright, and provide the declaration as the override, which is the
same shape as `dose_covariate` one step further along. Cite the existing
precedent in the package when there is one; it is the strongest argument
available and it keeps the API consistent.

---

# Part 2 — The rules

Every rule below is derived from an edit he actually made. The cited commits and
diffs are the evidence, and when a rule and the evidence disagree, the evidence
wins.

## 1. Do not describe the document inside the document

Cut on 2026-08-11 from `synthetic-data-checking-review.Rmd`:

- "is separate from the literature on generating it, has its own vocabulary,
  and is easy to reinvent badly if you have not read it"
- "written as a tutorial rather than a survey"
- "and says where `synpmx` fits among them"
- "Current as of August 2026"

What survived: "This article is a tutorial on these methods and at each step it
states which measures `synpmx` implements, which it does not, and why." One
sentence, and it describes content rather than the author's intentions. A
reader who is reading the document does not need to be sold on it.

## 2. A heading names its subject

| Deleted | Replaced with |
|---|---|
| The fix, and it is the same fix everywhere | Training vs Control Set |
| What it costs in pharmacometrics, honestly | Applications of control set to pharmacometrics datasets |
| Why this is hard: population facts versus patient facts | Introduction |
| How synthetic data is made | Algorithms for Generating Synthetic Data |

No stance words (`honestly`, `worth having`), no narrative beats (`The fix`,
`Now the part that needs the holdout`), no colon-plus-restatement. A reader
scanning the table of contents should be able to find a topic, not follow a
plot.

## 3. Cut transitional stage directions

Deleted: "Here is why the separation is hard." / "Now make it a measurement
problem." / "If that feels familiar, it should:" / "Here is the part that closes
the loop." These announce a move instead of making it. Delete the sentence and
start the paragraph at its content.

## 4. Repeat the noun; do not lean on a pronoun

He *adds* words for this, while deleting words everywhere else.

- "**You cannot tell from the synthetic data alone.**" became "**You cannot tell
  from the synthetic data alone if a leak of individual information occurred.**"
- "...only person who could have qualified for it, it is." became "...only person
  who could have qualified for it, then it is a leak."

A sentence that depends on the previous sentence for its subject fails when the
reader arrives from the table of contents, and this reader does arrive from the
table of contents.

## 5. Do not grade the evidence

[Author Note: Collect more data to confirm if this is truly a signal, I'm not sure.  Might just be an idiosyncratic preference here]

Deleted: "which is the strongest signal in this whole area that it is not
optional" / "the raw value of either one on its own is close to meaningless" /
"This is a real reason the technique is rare in this field and not a reason it
is wrong." Report what the four lines of work do. He will decide how strong that
is.

## 6. State a recommendation as an instruction

"**Hold some patients out.** Split the real cohort in two: a *training* set the
generator is allowed to use, and a *control* set it never sees." became "**To
assess a synthetic data generating algorithm, split real cohort into a training
and control (holdout) set.**" The bolded clause carries the action and the
definitions ride along inside it. See rule 10: the verdict comes first, the
nuance second.

## 7. Two items is a sentence, three is a list, five is a table

He collapsed a two-bullet list ("the generator captured the population, which is
success" / "the generator memorized individual patients, which is failure") back
into one clause. Bullets start at three items, or where each item is long enough
to be hard to hold.

Past about four comparable items, use a table instead. "I like your idea of the
scorecard table near the top" was accepted without modification, unlike almost
everything else in that message, and he immediately extended it ("maybe one part
of the scorecard is just synpmx_validate") — the response of someone who has
understood a structure well enough to add to it. Use a table when the items
answer the same question, and make every row answer it in the same form. Prose
comparison of five things is a reliable way to lose him.

## 8. Punctuation

- Em dashes are for a genuine parenthesis, not for a rhetorical pause. Most of
  his edits replace one with a period or a comma. Two per page, not two per
  paragraph.
- He types `--` and two spaces after a period. Leave both alone, and do not
  reflow his paragraphs to 80 columns when editing next to them.

## 9. Numbers survive; adjectives do not

"a *median local cloaking of 11*", "253 patients", "`k` = 5 donors" all survive
review untouched. "considerably", "genuinely", "wildly", "close to meaningless",
"the elegant part" are what gets cut. Prefer the measurement to the
characterization of the measurement.

## 10. A verdict leads, and the verdict may be a hedge

The four use-case lines in `README.md` each open with a marker and a verdict in
parentheses: **✅ Develop code**, **✅ Teaching tool (Yes)**, **⚠️ Send data past
a trust boundary (Use Caution)**, **❌ Answer scientific questions (No)**. The
nuance follows within the same line — "the formal privacy-protecting methods
provided with this package are illustrative, but not audited" is still there,
after the verdict rather than softening it.

**An earlier version of this rule said he converts hedged verdicts into binary
ones. He rejected it in writing and reverted the example.** The trust-boundary
row was Claude's "Only with care", became his "No", and is now "⚠️ Use Caution".
What the evidence supports is the *marker* first so the line is scannable. Do not
flatten a real "it depends" into a "No" to satisfy this rule.

## 11. State the out-of-scope, and state it bluntly

"Maybe I want to be explicit about this is not about scientific discovery." His
own `README.md` rewrite replaced a four-row table of hedged verdicts with four
bolded lines carrying ✅ / ⚠️ / ❌ markers. Every explanatory document should
carry an explicit out-of-scope statement, unqualified. The "why almost none of
this is `synpmx`'s problem" subsection in
`vignettes/articles/synthetic-data-checking-review.Rmd` exists because of this.

## 12. Thin at the entry point; literal over coined

`f8c22f1` cut a fully annotated 13-argument `synpmx_avatar()` call down to four
arguments and removed the multi-line hanging comments explaining `dvid` and
`dose_covariate`. The same commit deleted the coinage "The declaration is also
the **manifest of what survives**" in favour of "The function drops every column
that is not described."

Show the minimum that works; detail belongs one document deeper. And prefer
literal description to coined framing — the metaphors reached for to make
something memorable are what he deletes first.

In tension with O10, his preference for fewer review rounds: thinning a document
can cost a round when the cut detail was the answer to his next question. The
tier decides it — internal `design/` documents exhaustive, shipped entry points
thin.

## 13. One document, one question

2026-08-11, unprompted: *"I think the literature review should be split into two
files"* — generation into one article, checking into another. The article had
announced itself as having two halves since it was written and survived one
review pass in that form before he cut it.

Treat "this article has two halves" as a defect report a draft wrote about
itself. When an outline needs that sentence, propose two documents instead;
navigation is cheap, and a reader arriving for the checking tutorial should not
scroll through the generation survey. The shared instinct behind this rule,
rule 7 and rule 12 is that **each artifact should do one thing, at the smallest
size that does it**.

---

# Part 3 — Before handing over a draft

## The machine-prose tics

Widely reported markers of machine-written prose, several of which appear in
this repo's own drafts. Unlike the rest of this file, the list applies to
*everything* — conversation, commit messages and code comments as much as
documents.

- **Contrastive antithesis.** "It is not X, it is Y." "This is not about X, it is
  about Y." "not a bug, a feature."
- **The rule of three** applied to everything: three adjectives, three-clause
  sentences, three-item lists where two items exist.
- **Significance announcements.** "crucial", "pivotal", "key insight",
  "fundamental", "underscores", "highlights", "showcases", "testament to",
  "the important thing is".
- **Hedge-then-assert.** "It is worth noting that", "It is important to
  understand", "arguably", "in many ways".
- **Vocabulary.** delve, leverage, robust, seamless, landscape, realm, tapestry,
  navigate the complexities, deep dive, at its core, in essence, that said.
- **The closing summary** that repeats what the section just said. "In
  conclusion", "Ultimately", "The takeaway is". If the section needs a summary,
  it is too long.
- **Bold as emphasis spray.** Bold marks a term being defined or a verdict.
  Three bolded phrases in one paragraph mark none of them.
- **Symmetric sentence pairs.** "X does A. Y does B." repeated for rhythm rather
  than for content.
- **Second-person coaching.** "Let us break this down", "Think of it like",
  "Here is the thing".
- **Invented framing.** A coined metaphor ("the manifest of what survives") in
  place of a literal description. Deleted in `f8c22f1`. Rule 12 is the same
  finding stated as a rule.
- **Sycophancy.** No "great question", no praising the request, no announcing
  that something is a strong idea before doing it.

## The checklist

1. Every heading names a subject (rule 2).
2. No sentence describes the document, its structure, or its own difficulty
   (rule 1).
3. Each paragraph's first sentence stands alone with the topic explicit (rule 4).
4. Every acronym expanded at first use, and every method from outside
   pharmacometrics explained from zero (O8).
5. Verdicts lead; nuance follows, and a hedge is a legitimate verdict (rule 10).
6. Out-of-scope stated explicitly and bluntly (rule 11).
7. Numbers where a claim is measurable; no adjective standing in for one
   (rule 9).
8. Nothing says the document has two halves (rule 13). More than four comparable
   items are in a table (rule 7).
9. Read the em dashes. Keep the parentheticals, cut the pauses (rule 8).
10. Every named method answers *what is calculated, from what inputs, and what
    does a good number look like* (O1).
11. A survey ends with a ranked entry path (O9).
12. Claims are computed on a real dataset where they can be (O12).
13. The tics above are *searched for*, not read for: crucial, key, worth,
    honestly, delve, leverage, robust, seamless, "not just", "it is worth
    noting", "here is why", "the point is".

---

# Part 4 — What the outside literature says

Browsed 2026-08-05. Three findings map onto the observations closely enough to
be worth naming.

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
`avatar-algorithm.Rmd` is reference, `scorecard-synthetic-data-checks.Rmd` is a
how-to guide wrapped around an explanation.

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
verbosity, and over-elaboration — match his edits almost exactly. What he deleted
from the README was hedging and elaboration (rules 10 and 12, though note he
later restored one hedge and rule 10 records why). This is not a personal
idiosyncrasy; it is the general complaint, and it should be treated as a default
to correct rather than a preference to accommodate.

---

# Part 5 — Open questions to put to him

Things that cannot be inferred from the evidence available, roughly in order of
how much they would change future drafts.

1. **Did you read sections E and F of the checks vignette, and did they land?**
   (Resolves O14 — that is, whether argumentative prose earns its length.)
2. **Is the mathematics helpful or is it noise?** The literature review carries
   the adversarial-accuracy formula in display math. It could as easily be three
   sentences of English.
3. **Inline definitions or a glossary?** Right now every term is defined inline
   at first use, which is why that section is long.
4. **What is the length ceiling for a shipped vignette?** `avatar-algorithm.Rmd`
   is ~2,100 lines. Is that already past the point of usefulness, or is it fine
   because it is reference material nobody reads end to end?
5. **Does the ✅/⚠️/❌ marker style from your README rewrite generalize?** If you
   want it in the vignettes too, that is a cheap and consistent change.

---

# Part 6 — Proposed changes to this document

Suggestions, for Andy to accept or reject. Nothing here has been done.

1. **Split raw material from conclusions.** Part 1 currently mixes them. Add a
   question log holding his questions **verbatim and dated**, and keep the
   observations for distilled findings that cite the log. Verbatim wording
   matters — "I'm not sure I get what you're talking about" and "I don't care
   about that" look similar in a summary and mean opposite things (O4).

2. **Give every observation a fixed schema.** *Observation / Evidence / What to
   do*, as above. An observation without evidence is a guess, and an observation
   without an action does not change any document.

3. **Adopt an explicit disagreement marker.** Most of Part 1 is Claude's
   inference about Andy, and it needs a way for him to push back in place rather
   than have a wrong model re-derived each time. This has happened twice
   informally — the author note in rule 5 and the rejection recorded in rule 10 —
   and both are the most useful lines in their sections. A convention such as
   `**[AS: no, actually...]**` would make it routine.

4. **Add a document contract table.** One row per shipped document: audience,
   Diátaxis mode, target length, and how thin it must be. That turns the theory
   into something checkable before writing rather than a description of what went
   wrong afterwards. It also resolves the O10 conflict with rule 12 mechanically
   instead of by judgement each time.

5. **Record what was tried and rejected.** If a way of explaining something fails
   review, that is more informative than the version that passed, and it is
   currently lost. Two lines per rejection is enough.

6. **Revisit after the reviews.** The documents this file has the most to say
   about are ones he had not read when it was written. Re-run the analysis after
   he works through the two literature-review articles and
   `scorecard-synthetic-data-checks.Rmd`; the questions he asks then will be
   better evidence than anything above.
