# Writing style for synpmx documents

Rules for prose in `README.md`, vignettes, and articles. Every rule below is
derived from an edit Andy actually made to Claude's draft, or from a widely
reported failure mode of AI-written prose. Cited commits and diffs are the
evidence; when a rule and the evidence disagree, the evidence wins.

Related: `design/_THEORY_OF_MIND.md` records *how he reads*. This file records
*how to write*. O7, O6 and O15 there are the same findings stated as
observations rather than as rules.

## Rules from his edits

### 1. Do not describe the document inside the document

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

### 2. A heading names its subject

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

### 3. Cut transitional stage directions

Deleted: "Here is why the separation is hard." / "Now make it a measurement
problem." / "If that feels familiar, it should:" / "Here is the part that closes
the loop." These announce a move instead of making it. Delete the sentence and
start the paragraph at its content.

### 4. Repeat the noun; do not lean on a pronoun

He *adds* words for this, while deleting words everywhere else.

- "**You cannot tell from the synthetic data alone.**" became "**You cannot tell
  from the synthetic data alone if a leak of individual information occurred.**"
- "...only person who could have qualified for it, it is." became "...only person
  who could have qualified for it, then it is a leak."

A sentence that depends on the previous sentence for its subject fails when the
reader arrives from the table of contents, and this reader does arrive from the
table of contents.

### 5. Do not grade the evidence

[Author Note: Collect more data to confirm if this is truly a signal, I'm not sure.  Might just be an idiosyncratic preference here]

Deleted: "which is the strongest signal in this whole area that it is not
optional" / "the raw value of either one on its own is close to meaningless" /
"This is a real reason the technique is rare in this field and not a reason it
is wrong." Report what the four lines of work do. He will decide how strong that
is.

### 6. State a recommendation as an instruction

"**Hold some patients out.** Split the real cohort in two: a *training* set the
generator is allowed to use, and a *control* set it never sees." became "**To
assess a synthetic data generating algorithm, split real cohort into a training
and control (holdout) set.**" The bolded clause carries the action and the
definitions ride along inside it. Compare O6: verdict first, nuance second.

### 7. Two items is a sentence, not a list

He collapsed a two-bullet list ("the generator captured the population, which is
success" / "the generator memorized individual patients, which is failure") back
into one clause. Bullets start at three items, or where each item is long enough
to be hard to hold. (O13: past four comparable items, use a table.)

### 8. Punctuation

- Em dashes are for a genuine parenthesis, not for a rhetorical pause. Most of
  his edits replace one with a period or a comma. Two per page, not two per
  paragraph.
- He types `--` and two spaces after a period. Leave both alone, and do not
  reflow his paragraphs to 80 columns when editing next to them.

### 9. Numbers survive; adjectives do not

"a *median local cloaking of 11*", "253 patients", "`k` = 5 donors" all survive
review untouched. "considerably", "genuinely", "wildly", "close to meaningless",
"the elegant part" are what gets cut. Prefer the measurement to the
characterization of the measurement.

## General AI tics to avoid

Widely reported markers of machine-written prose. Several appear in this repo's
own drafts.

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
  place of a literal description. Deleted in `f8c22f1`.
- **Sycophancy.** No "great question", no praising the request, no announcing
  that something is a strong idea before doing it.

## Checklist before handing over a draft

1. Every heading names a subject.
2. No sentence describes the document, its structure, or its own difficulty.
3. Each paragraph's first sentence stands alone with the topic explicit.
4. Every acronym expanded at first use in that document (`AGENTS.md`).
5. Verdicts lead; nuance follows.
6. Out-of-scope stated explicitly and bluntly (O5).
7. Numbers where a claim is measurable; no adjective standing in for one.
8. Search the draft for: crucial, key, worth, honestly, delve, leverage, robust,
   seamless, "not just", "it is worth noting", "here is why", "the point is".
9. Read the em dashes. Keep the parentheticals, cut the pauses.
