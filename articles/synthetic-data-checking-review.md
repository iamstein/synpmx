# Checking Synthetic Data: A Tutorial on the Published Methods

### Introduction

The literature on checking synthetic data is spread across three
disciplines: statistics, machine learning, and data protection law. This
article is a tutorial on these methods and at each step it states which
measures `synpmx` implements, which it does not, and why.

Its companion, [the generation
review](https://iamstein.github.io/synpmx/articles/synthetic-data-generation-review.html),
surveys the four families of generation method in use.

A synthetic dataset is *supposed* to tell you about the population but
it must not tell you about a patient and each privacy check is an
attempt to mask patient identifying information.

Masking patient IDs is not necessarily sufficient for protecting
privacy. For example, suppose your trial found that patients over 80 kg
cleared the drug slowly, and the synthetic data reproduces that – that
is a fact about the population and ok to be reproduced. Now suppose
exactly one patient in your trial weighed 138 kg, and the synthetic data
contains a 137 kg patient. Is that a leak of individual information? It
might be.

**You cannot tell from the synthetic data alone if a leak of individual
information occurred.** If people of 138 kg are common in the world, it
is not. If your trial’s 138 kg patient is nearly the only person who
could have qualified for it, then it is a leak. Looking at the synthetic
data alone cannot answer the question of whether there was a leak of
private information. Any statistic computed from the source and the
synthetic dataset together (a distance, distribution, or comparison) can
measure both properties of the population and properties of individual
patients.

### Training vs Control Set

**To assess a synthetic data generating algorithm, split real cohort
into a training and control (holdout) set.** Then ask your similarity
question twice — synthetic against training, synthetic against control.

Whatever the generator learned about the **population** shows up in both
comparisons, because the control patients are drawn from the same
population. Whatever it learned about **specific training patients**
shows up only in the first. The **difference between the two is the
leak**. This process is exactly the train/test split from machine
learning, run for the opposite purpose. In modelling, accuracy on
training data confounds a “learned pattern” with “memorized, overfit
examples”, and a test set separates them. Privacy measurement is the
other side of this problem — here memorization/over-fitting is the thing
you are trying to *detect* rather than to avoid. Four independent lines
of work use this framework.

- **Anonymeter** (Giomi and colleagues) builds the control dataset into
  the framework and reports risk only where an attack succeeds better
  against training data than against control data.
- **Adversarial accuracy** (Yale and colleagues) is *defined* as a
  training-versus-holdout difference; see below.
- **synthpop’s disclosure measures** (Raab and colleagues) subtract what
  is disclosive in the original data anyway, which is the same
  correction applied to a computation rather than to a dataset.
- **Differential privacy** is this idea taken to its limit: rather than
  measuring the difference one patient makes, it *bounds* it in advance,
  for every patient and every possible statistic.

#### Applications of control set to pharmacometrics datasets

For Phase 2-3 studies, removing 20% of 200 patients in a two arm trial
is reasonable, but removing 20% of patients when there are only 12-24 in
an early dataset, or when cohort sizes are only 3 patients could remove
an entire cohort from the study. Thus a synthetic data generator can be
characterized using larger studies, but for smaller studies, further
caution is needed to ensure manage the individual-patient-level
information making it into the synthetic dataset.

### Adversarial Accuracy

[`compare_pmx_proximity()`](https://iamstein.github.io/synpmx/reference/compare_pmx_proximity.md)
implements this statistic.

Put every real subject and every synthetic subject into one space, where
each subject is a point. For `synpmx` that space is built from
covariates and trajectory features. For each point, ask one question:
*is my nearest neighbour in my own dataset, or in the other one?*

Write $`d_{TS}(i)`$ for the distance from real subject $`i`$ to the
nearest synthetic subject, and $`d_{TT}(i)`$ for the distance from real
subject $`i`$ to the nearest **other real** subject; define
$`d_{ST}(j)`$ and $`d_{SS}(j)`$ the same way from the synthetic side.
Adversarial accuracy is the average, over both sides, of how often the
*other* dataset is farther away than your own:

``` math
\mathrm{AA}=\frac12\left[\frac1n\sum_i \mathbf{1}\{d_{TS}(i)>d_{TT}(i)\}
+\frac1n\sum_j \mathbf{1}\{d_{ST}(j)>d_{SS}(j)\}\right].
```

It is called *adversarial* because it is the success rate of the
simplest possible adversary: one who is handed a record and guesses
which dataset it came from by looking at what is nearest to it.

Both ends of the scale are failures:

- **AA near 1.** Every subject’s nearest neighbour is in its own
  dataset. The two clouds are separated and the synthetic data is
  nothing like the real data. This is a *utility* failure.
- **AA near 0.** Every subject’s nearest neighbour is in the *other*
  dataset. The two clouds are interleaved so tightly that each real
  patient has a synthetic partner closer than any of their real peers.
  This is the signature of **memorization**, a privacy failure.
- **AA near 0.5.** A synthetic subject is no more like a real subject
  than two real subjects are like each other. This is the target, and it
  is a target rather than a maximum.

One number therefore reports two failure directions.

#### Adversarial accuracy with a holdout

In Yale and colleagues’ original formulation, adversarial accuracy is
computed **twice** — once against the training data the generator saw
($`\mathrm{AA}_{\text{train}}`$) and once against a holdout it did not
($`\mathrm{AA}_{\text{test}}`$) — and the reported privacy loss is the
difference:

``` math
\text{privacy loss}=\mathrm{AA}_{\text{test}}-\mathrm{AA}_{\text{train}}.
```

A generator that captured the population produces synthetic subjects
that sit equally close to training and holdout patients, so the two
accuracies agree and the difference is zero. A generator that memorized
its training patients produces synthetic subjects that sit *unusually*
close to the training patients specifically, which drags
$`\mathrm{AA}_{\text{train}}`$ toward 0 while leaving
$`\mathrm{AA}_{\text{test}}`$ where it was. A positive difference is
leakage, in units you can compare across datasets.

#### What synpmx computes

[`compare_pmx_proximity()`](https://iamstein.github.io/synpmx/reference/compare_pmx_proximity.md)
computes $`\mathrm{AA}_{\text{train}}`$ only, and calibrates it against
a null interval built by recomputing the same statistic on random halves
of the source cohort. At 20 patients in a 30-dimensional profile space,
nearest-neighbour statistics are unstable and 0.5 is *not* the right
expectation, so a raw adversarial accuracy is uninterpretable without
that null. The null interval states what this statistic does on data of
this size and shape when nothing is wrong.

The null interval corrects for cohort size and dimension. It does not
correct for the population/patient confound, because both halves of the
source were seen by the generator. A holdout is what would add that
correction, and its absence is the largest gap in the package’s privacy
measurement today.

### Local Cloaking and Hidden Rate

Guillaudeux and colleagues built two privacy metrics for a
**patient-centric** generator, one where each synthetic record is built
from one identified original record.
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
is such a generator: every avatar has an anchor. A global generator such
as a GAN cannot compute these metrics, because no record-to-record
correspondence exists to compute them over.

Both metrics answer one question:

> For this particular real patient, how well hidden is their own avatar?

**Local cloaking.** Take real patient $`i`$ and their avatar $`a(i)`$.
Count how many *other* avatars lie closer to $`i`$ than $`a(i)`$ does.
That count is patient $`i`$’s local cloaking: the number of other
avatars an attacker would pass through before reaching this patient’s
own.

- **Local cloaking 0** is the worst case: your own avatar is the avatar
  nearest to you, so an attacker who has your real record and guesses
  “the closest avatar is derived from me” is right.
- Higher is better. The published values are a **median local cloaking
  of 11** on an AIDS clinical trial dataset and **24** on a breast
  cancer dataset. Local cloaking cannot exceed the number of other
  avatars, so on a 12-patient study the ceiling is 11 and those
  published values are not reachable as targets. The scale of this
  metric is set by cohort size.

**Hidden rate.** The percentage of real patients whose own avatar is
*not* the avatar nearest to them, which is the percentage with local
cloaking of at least 1. It aggregates the same measurement, and it is
the failure probability of the nearest-neighbour attack. The published
values are **93% and 94%**.

#### Granularity: a cohort number against a patient list

The two metrics answer a different question from adversarial accuracy,
at a different granularity.

- Adversarial accuracy is a **cohort-level** statistic: *does this
  dataset as a whole look memorized?* It gives you one number and no
  names.
- Local cloaking is a **per-patient** statistic: *is this patient
  exposed?* It gives you a list.

For remediation, a list is what you need.
[`remediate_identifiable_subjects()`](https://iamstein.github.io/synpmx/reference/remediate_identifiable_subjects.md)
acts on individual subjects, and today it can only act on what
[`flag_identifiable_subjects()`](https://iamstein.github.io/synpmx/reference/flag_identifiable_subjects.md)
found — structural outliers. A patient with local cloaking 0 is exposed
for a completely different reason and is invisible to that screen.

Hidden rate is a **membership inference** measure: the error rate of an
adversary trying to link a known individual to the release. Membership
inference is the strongest empirical privacy test in the general
literature.

#### Limits of both metrics

**The pairing map is the most disclosive artifact in the whole
pipeline.** Local cloaking and hidden rate require knowing which avatar
came from which patient. That map is a re-identification key: anyone
holding it does not need to attack anything. It has to be computed
inside the trusted environment, used, and destroyed. It must never be
attached to the released data or written next to it. `synpmx` does not
retain the anchor of each avatar, which is why these metrics are not
implemented. Implementing them means creating that map deliberately and
handling it more carefully than anything else in the package.

**Both measure one specific attack.** Hidden rate is the failure rate of
an attacker who holds a real record and the whole avatar set and guesses
nearest-neighbour. A smarter attacker exists. Both metrics are concrete
and neither is a bound, which is the caveat on every measure in this
article except differential privacy.

The same paper reports two generic distance measures:

- **Distance to closest record (DCR).** For each synthetic record, the
  distance to the nearest real one. Larger is safer; zero is a copy.
- **Nearest-neighbour distance ratio (NNDR).** The ratio of the nearest
  distance to the second-nearest. Near 1 means the record sits in a
  crowd with no single conspicuous partner. Near 0 means one real record
  is much closer than everything else, which is the shape of a copy. The
  paper treats 0.8 and above as satisfactory.

Both are cheap to compute. The section on similarity metrics below gives
the reason not to rely on either.

### Disclosure Risk in Official Statistics: RepU and DiSCO

The most developed disclosure-risk framework for synthetic data comes
from official statistics rather than from machine learning, and it ships
as working code. Raab, Nowok and Dibben formalized it in 2024 and
implemented it in `synthpop` (version 1.8.1 and later) as `disclosure()`
and `multi.disclosure()`.

#### Keys and targets

- A **key**, called a *quasi-identifier* in the older literature, is a
  set of variables an attacker plausibly **already knows** about the
  person they are looking for. Age band, sex, treatment arm, site, date
  of first visit. Not identifying on their own; identifying in
  combination.
- A **target** is the variable they want to **learn**: a genotype, a
  diagnosis, a lab value.

Disclosure is then a checkable event: *the key lets them pin down the
target.* The choice of key is a judgement call by the data holder rather
than something a package can infer, which the paper states explicitly
and which is the framework’s main practical weakness. New external data
sources create new keys, so an assessment has a shelf life.

#### Identity disclosure: RepU

*Can an attacker work out which record is this person?*

**RepU — replicated uniques** — is the percentage of original records
that are unique on the key set in the **original** data *and* also
unique on that same key in the **synthetic** data.

The two conditions do different jobs. Being unique in the original is
the patient’s pre-existing exposure, which the release did not cause.
Being *replicated* as a unique in the synthetic data is what the release
**added**: the attacker who knows the key can now find a single record,
and its other values are asserted about that person. `synthpop` will
optionally remove replicated uniques from the output, which is a
disclosure control rather than a diagnostic.

#### Attribute disclosure: DiSCO

*Can an attacker learn a value about this person, without necessarily
picking out their record?*

Attribute disclosure does not require identifying anybody. If every
synthetic patient in a given key cell has the same target value, then an
attacker who knows a real person is in that cell learns their value
without ever locating their record.

**DiSCO — Disclosive in Synthetic, Correct in Original** — is the
percentage of original records whose key combination maps to a
**single** target value in the synthetic data, and where that value is
the one the real person **actually has**. Two conditions again: the
synthetic data has to make a confident assertion, and the assertion has
to be right.

#### DiO, the control computation

**DiO — Disclosive in Original** — is the same computation run on the
original data. If the key already determines the target in the real
data, then the attacker did not need your synthetic dataset: they could
have derived the relationship from published literature, from the
protocol, or from clinical knowledge. That is a **population fact**. It
is the science.

So the quantity that measures what the *release* added is not DiSCO but
**DiSCO − DiO**.

This is the patient-versus-population separation again. Anonymeter uses
a control **dataset**, Yale uses a holdout, and Raab uses a control
**computation** on the original data. One correction, three
implementations.

`synpmx` arrived at the same correction independently: the
identifiability screen scores each patient **within their treatment
arm** rather than against the whole cohort, because on a six-arm
dose-ranging study a cohort-wide screen flags thirty patients for having
received the dose their arm was assigned. That is a protocol fact
reported back as a privacy finding, which is DiO under another name.
`DiSCO − DiO` states the same correction as arithmetic.

#### Fit to pharmacometric data

The framework fits the baseline table and does not fit the event table.

At the level of **subject baselines** — one row per patient, covariates
and strata — a pharmacometric dataset *is* rectangular, which is the
shape `synthpop`’s machinery assumes. The covariate table can go
straight into `disclosure()` with keys chosen from arm, sex, age band,
site and so on. This is the natural home for the rare-category question.

At the level of the **event table** it does not apply, for the reasons
the [generation
review](https://iamstein.github.io/synpmx/articles/synthetic-data-generation-review.html)
gives: many rows per subject, ordered in time, with a grammar the tools
do not model. That is where the timing and structural checks have to do
the work instead, and where there is essentially no external literature
to borrow from.

### The Article 29 Working Party Criteria

#### WP29 and Opinion 05/2014

The **Article 29 Data Protection Working Party** — universally “WP29” —
was the advisory body of the European national data protection
authorities under the 1995 Data Protection Directive. It was replaced by
the European Data Protection Board in 2018 when the General Data
Protection Regulation (GDPR) took effect, but its **Opinion 05/2014 on
Anonymisation Techniques** is still the reference document practitioners
and regulators reach for, and it is short.

Its contribution is a three-part test. A dataset is only effectively
anonymized if an attacker can do **none** of the following:

1.  **Singling out** — isolate some or all records that identify an
    individual within the dataset. *“The 63-year-old man at site 4 on
    the 300 mg arm”* picks out one row, whether or not you know his
    name.
2.  **Linkability** — link at least two records concerning the same
    individual, either within one dataset or across two datasets.
3.  **Inference** — deduce, with significant probability, the value of
    an attribute from the values of other attributes.

A data protection officer applies these three criteria, and they
partition the space differently from an engineer’s intuition.

#### Linkability

Linkability is the ability to connect two records that belong to the
same person, and **it does not require identifying anybody.**
Establishing that record 47 in your synthetic dataset and record 12 in
some other release describe the same person is a failure on its own,
because such links accumulate: each one narrows the field, and enough of
them identify.

The concrete attack for synthetic data is that an adversary holds a
*partial* record from elsewhere — a covariate row from a registry, a
table in a publication, or their own knowledge as the person’s treating
clinician — and asks which synthetic record corresponds to it.
Anonymeter implements this attack: split the attributes into two sets,
give the attacker one set, and measure whether they recover the right
partner.

**A record-based generator is more exposed to linkability than any other
family**, because a link between an avatar and its anchor exists *by
construction*. That link is what local cloaking measures. The two
literatures describe the same risk from opposite ends. The `synpmx`
checks vignette covers singling out thoroughly, inference partially, and
linkability not at all.

#### Anonymeter

Giomi and colleagues turned the three criteria into a working evaluation
suite, one attack per criterion, with the control dataset built into
each of them. Risk is reported only where the attack performs better
against the training data than against the control, which stops the tool
from reporting the population back as a leak.

### Similarity Metrics Do Not Bound Attack Success

Ganev and De Cristofaro (2025) attack the whole family of
similarity-based metrics: distance to closest record, nearest-neighbour
ratios, and their relatives. Synthetic datasets that pass these metrics
comfortably were still attacked successfully, with training records
reconstructed.

The reason is structural rather than a matter of tuning:

> These are **average-case** statistics over a whole dataset. An
> attacker targets **one** record, and only needs to succeed once.

A synthetic dataset can sit, on average, a respectable distance from its
source and still contain one record that is a near-copy of one patient.
The mean does not see it. A better threshold does not fix this, because
the statistic answers a different question from the one being asked.

Four conclusions follow:

- **Distance metrics are a smoke alarm, not a safety certificate.** They
  catch blatant copying: hand
  [`compare_pmx_proximity()`](https://iamstein.github.io/synpmx/reference/compare_pmx_proximity.md)
  a verbatim copy of the source and it goes to zero and objects, which
  the package’s own regression test requires of it.
- **Prefer per-record measures over cohort averages** wherever both
  exist. This is the argument that makes local cloaking more actionable
  than adversarial accuracy.
- **The empirical gold standard is a membership inference attack**,
  which determines whether a specific individual was in the training
  data. It is expensive, it needs a holdout, and distance metrics are a
  cheap proxy for it.
- **The only claim that survives all of this is a formal one**, which is
  the argument for the differentially private modes.

### Fidelity, Diversity, and Authenticity

Alaa and colleagues (2022) split “is the synthetic data good?” into
three questions. Collapsing them into one hides which of the three
broke:

- **α-precision — fidelity.** Are synthetic records *typical* of real
  ones? The failure it catches is invention: patients the generator made
  up who could not exist.
- **β-recall — diversity, or coverage.** Do the synthetic records
  *cover* the real ones? The failure it catches is called **mode drop**:
  the generator produces only the comfortable middle of the
  distribution, and the unusual patients have no counterpart at all.
- **Authenticity — generalization.** Are the synthetic records *new*, or
  are they near-copies of training records? This is the privacy axis,
  and it is computed **per record**.

#### Losing the tails

Losing the tails is mostly intended.

Coverage and privacy are in direct opposition **at the tails**, because
in a clinical dataset the tail *is* the identifiable patient. The
patient with the longest follow-up, the highest dose, the extreme weight
— these are the ones a masking mechanism is designed to remove or blend
away. β-recall falling at the extremes is the masking **working**.

Coverage is therefore not a target to hit. It is useful for **telling
two causes apart**, and the two causes have opposite implications:

- The tail is missing because the outlier was deliberately screened out
  or blended away. That is a designed cost, and `synpmx` reports it
  directly as anchors screened out, patterns dropped, and regimens not
  represented.
- The tail is missing because the generator collapsed toward the mean
  and produces nothing but average patients. That is a defect, and for a
  package whose entire purpose is producing the *awkward* patients that
  break software, it is close to a total failure.

**In a standard deviation these two causes look identical.** Measure
coverage so that a shrinking spread can be attributed, not so that it
can be optimized. Report it, do not chase it.

#### Authenticity

Authenticity is the “too close” direction, per record. It asks, of each
synthetic record individually, whether the record is a copy of some
training record’s neighbourhood rather than a new draw.

**Authenticity is the machine-learning literature’s version of local
cloaking.** Different field, different derivation, same question: *is
this particular output too close to a particular input?* Three
traditions, three names, one measurement.

### Utility Measures

#### General utility: pMSE

**General** — or *broad* — utility asks whether the synthetic dataset
resembles the real one at all, without reference to any particular
analysis. The standard measure is **pMSE**, the propensity score mean
squared error (Snoke, Nowok, Raab, Dibben and Slavkovic, 2018).

The mechanism, which appears in the literature under four different
names:

1.  Stack the real and the synthetic records into one table, and label
    each row with which dataset it came from.
2.  Fit a classifier — usually a CART tree — to predict that label from
    the data.
3.  Look at the predicted probabilities, the *propensity scores*. If the
    two datasets are indistinguishable, every record’s predicted
    probability sits near the mixing proportion (0.5 for equal sizes),
    because the classifier cannot do better than guessing.
4.  pMSE is the mean squared deviation of those probabilities from that
    proportion. **Low is good**, and it has a known null value to
    compare against.

The same idea is called the *discriminative score* in the time-series
literature and the *classifier two-sample test* in statistics. It has
the same shape as adversarial accuracy: a “can you tell them apart?”
statistic with a calibrated null.

#### Specific utility: confidence interval overlap

**Specific** — or *narrow* — utility asks whether **the analysis you
actually intended to run** gives the same answer on the synthetic data
as on the real data. The standard measure is **confidence interval
overlap** (Karr and colleagues, 2006): run the analysis on both
datasets, and measure how much the two confidence intervals for the same
quantity overlap.

In pharmacometrics that means fitting the same population model to both
and comparing the fixed effects, the between-subject variance terms, and
their intervals. Destere’s daptomycin benchmark is this measure, and it
is the question a modeller asks first.

#### Broad and specific utility are weakly correlated

Comparing marginal distributions asks whether the columns look right
**one at a time**. Comparing the analysis asks whether the
**conclusion** is the same. Every column can match while the
relationships between them, which is what a model estimates, are
destroyed.

Woillard and colleagues measured this: broad utility and specific
utility were only **weakly correlated** in their benchmark, and data
that looked right column by column did not support the same predictions.
**Marginal agreement is weak evidence**, and a document reporting only
distribution comparisons is reporting the easy half.

#### Utility measures out of scope for synpmx

`synpmx` exists so that **software** can be developed against realistic
pharmacometric data. Under that use case, “would a modeller reach the
same scientific conclusion from this data?” is not a requirement that
has been weakened — it is **not a requirement at all**, because nobody
should be reaching scientific conclusions from this data in the first
place, and the package says so in every document it ships.
Between-subject variance is knowingly shrunk by blending. A dose amount
is rebuilt from each avatar’s own covariate only where the study is
weight- or body-surface-area-based and that was declared with
`dose_covariate` or inferred. These are documented properties, not
defects awaiting a fix.

The output does have to be **the same kind of object** as a real study
dataset: a legal event grammar, dosing histories that a protocol could
have produced, several endpoints on one visit grid, patients who missed
visits and withdrew early, and a schema an assembly script can be
pointed at.

That reorders the list above:

- **Structural and validity checks matter most**, and they are the part
  this literature barely addresses — because tabular and
  machine-learning tools do not have event tables, so nobody outside
  pharmacometrics has needed to check one. There is little to borrow
  here and it is where the package’s own work has to be original.
- **Privacy measures matter in full**, unchanged. The source data is
  real regardless of what the output is used for, and none of the
  reasoning in the preceding sections is affected by the use case.
- **General utility measures such as pMSE are not binding.** A poor pMSE
  would prompt a look at *why*, not a change to hit a number.
- **Specific utility measures are out of scope by design**: confidence
  interval overlap, parameter recovery, non-compartmental exposure
  agreement. The package does not claim the property they measure.
  Stating that is more useful than reporting them badly and letting a
  reader infer a claim that was never made.

The one binding requirement in the utility column appears nowhere in
this literature, and no function can compute it: **does the pipeline
that will consume the real study run unchanged against this?** Only the
person holding that pipeline can answer it.

### Can Enumerated Checks Be Sufficient?

`synpmx` protects the output with a list of enumerated conditions rather
than with a bound. No visit pattern held by fewer than
`min_pattern_share` = 2 source patients reaches an avatar. A stratum
holding fewer than three source patients is deliberately not
size-balanced, because reproducing its size exactly would disclose it: a
cell holding one real patient would otherwise receive exactly one avatar
on every seed. Structurally extreme subjects never anchor. Values are
blended across `k` = 5 donors with no donor exceeding `max_donor_weight`
= 0.50.

#### These checks have published names

Each of those conditions is a known disclosure control, which is worth
knowing because the published versions come with published attacks.

| The check in `synpmx` | Its name in the literature |
|----|----|
| No visit pattern held by fewer than 2 patients is released | **k-anonymity** with k = 2, on the visit pattern as quasi-identifier |
| Strata under three patients are not size-balanced | **Small-cell suppression** |
| Structurally extreme subjects never anchor | **Special-uniques suppression** |
| No donor exceeds half of any avatar | **Bounded per-record influence** |

#### The published answer is that enumeration does not close

k-anonymity (Sweeney, 2002) was defeated by the homogeneity and
background-knowledge attacks, which motivated l-diversity
(Machanavajjhala and colleagues, 2007). l-diversity was defeated by
skewness and similarity attacks, which motivated t-closeness (Li and
colleagues, 2007). Ganta and colleagues (2008) then showed that two
releases each satisfying these conditions can be composed to identify
people neither release identifies alone.

Every step in that sequence added a condition in response to an attack
the previous enumeration did not cover. Differential privacy was
formulated in the same period as the alternative to continuing it.

Three reasons the enumeration does not close:

- **Combinations outrun enumeration.** Checks are written per field or
  per small group of fields. Identification comes from combinations,
  which grow exponentially in the number of fields.
- **Auxiliary information is outside your control and changes over
  time.** Raab and colleagues make this point about the choice of key: a
  new registry, a new publication, or a treating clinician’s own
  knowledge creates a key that did not exist when the assessment was
  written.
- **The attacker needs one success.** Enumerated checks and average-case
  statistics both describe the dataset. The attack targets one record.

#### When enumerated checks are sufficient anyway

They are sufficient when three conditions hold together, and none of the
three is a property of the data.

1.  **The quasi-identifiers are closed and known.** You can list what an
    attacker plausibly knows, and the list is stable over the life of
    the release.
2.  **The release is one-shot.** No second release from the same source,
    because composition is what defeats condition 1.
3.  **The recipient is inside a trust boundary.** Someone is
    contractually obliged not to attack the data and not to pass it on.

Under those three conditions the checks are a coverage argument that a
human evaluates, which is how most clinical data sharing is in fact
authorized. What the checks do not produce is anonymous data, which is a
property that would survive leaving the boundary. This is why
`README.md` marks sending output past a trust boundary as ⚠️ rather than
✅.

#### Why the blending cannot be made formal

Differential privacy needs three things: per-record influence bounded by
a known sensitivity, randomization calibrated to that sensitivity, and
accounting across everything released.

[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
has the first, in `max_donor_weight` = 0.50, and it randomizes the donor
weights. The randomization is not calibrated to the sensitivity, and
nothing accounts across the release, so no epsilon follows. The bound
itself is also loose: one donor may move an avatar by half its value,
and noise calibrated to a sensitivity that large would destroy the data.
Tightening the cap toward 1/`k` and adding calibrated noise would
approach a mechanism with worse utility than
[`synpmx_empirical()`](https://iamstein.github.io/synpmx/reference/synpmx_empirical.md),
which already releases differentially private summaries.

The recommendation is to keep the two paths separate rather than to
formalize the blending, which is what the package already does. State
each check as an invariant together with the assumption it needs, so
that the argument is auditable:

> No released visit pattern is held by fewer than two source patients,
> **assuming** attendance is not observable from outside the trial.

The invariant is verifiable and the assumption is not. Naming both is
what distinguishes an argument a data protection officer can evaluate
from a claim.

### What To Do At Each Cohort Size

Cohort size decides which checks are interpretable, not only how
expensive they are. Sizes below are source patients, and the mechanism
thresholds are the
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
defaults.

| Source cohort | What breaks | Which checks mean anything | Recommended |
|----|----|----|----|
| **Fewer than 6 on a route** | No legal donor set exists, since `k` + 1 = 6 subjects are needed. `on_donor_shortfall = "drop"` omits those subjects and that arm is not represented | None | [`synpmx_prior()`](https://iamstein.github.io/synpmx/reference/synpmx_prior.md) only. Every other mode reduces to noising one patient |
| **6-24** (Phase 1 SAD/MAD) | Holdout impossible. `k` = 5 is reached only by borrowing donors across dose groups, and the run warns that measurements are blended across doses. Strata under three patients are not size-balanced | Adversarial accuracy against its null interval. Structural and validity checks in full. Local cloaking ceiling is n − 1, so the published values of 11 and 24 are unreachable | [`synpmx_prior()`](https://iamstein.github.io/synpmx/reference/synpmx_prior.md) for anything leaving the trust boundary. [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md) inside it, with the masking report read in full |
| **25-60** | Holdout would remove a whole cohort of a dose-ranging design. Differential privacy measured as Tier A only at N = 60 | As above, with same-schedule donors more often reaching `k` without crossing dose groups | [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md) inside a trust boundary. No formal guarantee is available at this size except epsilon = 0 |
| **60-200** | Differential privacy still not usable; N = 100 measured as unusable | A one-off holdout characterization becomes affordable at the top of this range: 20% of 200 is 40 patients | [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md), characterized once with a holdout, then shipped without one |
| **200-1000** | Differential privacy marginal at N = 500 and usable at N = 1000 only at a weak epsilon | Everything above, plus `synthpop::disclosure()` on the baseline table, since key cells now hold several patients | [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md) with measured privacy, or a differentially private mode if the release must cross a boundary |
| **Above 1000** | Nothing specific to size | The published literature applies as written. N = 10000 is the first size where differential privacy works at epsilon 1 | Any mode, chosen by where the data is going |

The differential privacy verdicts are measured rather than extrapolated;
the frontier is in the [feasibility
article](https://iamstein.github.io/synpmx/articles/feasibility.html).

Two consequences of this table are worth stating separately.

**Below about 60 patients, no measurement in this article distinguishes
a good generator from a bad one.** The holdout is unaffordable, the null
intervals are wide, and the per-patient metrics have ceilings set by
cohort size. Characterize a generator on a large public dataset, then
run it on the small study. The argument for a small study is inherited
from that characterization plus the mechanism counts on the run, not
established on the study itself.

**A Phase 1 dose-escalation design is the hardest case in this table**,
and not because of its size alone. Cohorts of 3 to 6 mean that donors
are borrowed across dose levels, so an avatar keeps its anchor’s dose
while its measurements blend patients dosed at other levels, which
flattens the exposure-dose relationship. A single-patient cohort is
worse: its dose level identifies that patient by construction, and the
package refuses to represent a regimen only one patient received. For
those studies the useful question is not which check to run but whether
[`synpmx_prior()`](https://iamstein.github.io/synpmx/reference/synpmx_prior.md)
from a public model is enough, since it reads no data and needs no check
at all.

### Reading Order

Four papers, in this order:

1.  **WP29 Opinion 05/2014**, section 2.2. The three criteria, in about
    three pages, and the vocabulary every other document assumes.
2.  **Giomi and colleagues (Anonymeter).** The control-group idea and
    one attack per criterion. Short, implemented, and framed the way
    regulators frame it.
3.  **Raab, Nowok and Dibben (2024).** Keys and targets, RepU, DiSCO,
    and the DiO correction. Then run `synthpop::disclosure()` on
    something.
4.  **Guillaudeux and colleagues.** Local cloaking and hidden rate, the
    patient-centric measures that apply directly to this package.

Then **Ganev and De Cristofaro** for the critique of items 2 through 4.

#### Summary of measures

| Measure | Question it asks | Granularity | Needs a holdout? |
|----|----|----|----|
| Adversarial accuracy | Is this cohort memorized, or too dissimilar? | Cohort | Yes, to mean anything |
| Local cloaking | Is *this patient’s* record too close to *their own* output? | Per patient | No |
| Hidden rate | What fraction of patients are exposed that way? | Cohort | No |
| DCR / NNDR | Is any synthetic record suspiciously close to a real one? | Per record | Helps |
| RepU | Does the release make a unique record findable again? | Per record | No |
| DiSCO − DiO | Does the release let someone infer a value they could not before? | Per record | No (DiO is the control) |
| Singling out / linkability / inference | Can the three regulatory attacks succeed? | Per attack | Yes |
| Membership inference | Was this individual in the training data? | Per patient | Yes |
| pMSE | Can a classifier tell real from synthetic? | Cohort | No |
| Confidence interval overlap | Does the intended analysis reach the same answer? | Per estimate | No |
| α-precision / β-recall | Are outputs typical / is the real data covered? | Cohort | No |
| Authenticity | Is this output a near-copy of an input? | Per record | No |
| Differential privacy (ε) | How much can *any* one patient change *any* release? | Guarantee | Not applicable |

The last row is the only one that is a bound rather than a measurement.
Why differential privacy is a guarantee rather than a method, and what
it costs at pharmacometric cohort sizes, is in the [generation
review](https://iamstein.github.io/synpmx/articles/synthetic-data-generation-review.html).

## Where to go next

- [Checks of the synthetic
  data](https://iamstein.github.io/synpmx/articles/synthetic-data-checks.html)
  — `synpmx`’s own worked checks, and the gaps between them and the list
  above.
- [Generating synthetic pharmacometric
  data](https://iamstein.github.io/synpmx/articles/synthetic-data-generation-review.html)
  — the companion review: the four families of generation method, and
  where `synpmx` sits among them.
- [Privacy](https://iamstein.github.io/synpmx/articles/synpmx-privacy.html)
  — when the enumerated protections are not enough, and what a formal
  guarantee buys instead.

## References

- Article 29 Data Protection Working Party. *Opinion 05/2014 on
  Anonymisation Techniques.* WP216, adopted 10 April 2014. (The singling
  out / linkability / inference criteria; section 2.2.)

- Sweeney L. *k-anonymity: a model for protecting privacy.*
  International Journal of Uncertainty, Fuzziness and Knowledge-Based
  Systems. 2002;10(5):557-570.

- Machanavajjhala A, Kifer D, Gehrke J, Venkitasubramaniam M.
  *l-diversity: privacy beyond k-anonymity.* ACM Transactions on
  Knowledge Discovery from Data. 2007;1(1):3.

- Li N, Li T, Venkatasubramanian S. *t-closeness: privacy beyond
  k-anonymity and l-diversity.* IEEE International Conference on Data
  Engineering (ICDE). 2007:106-115.

- Ganta SR, Kasiviswanathan SP, Smith A. *Composition attacks and
  auxiliary information in data privacy.* ACM SIGKDD. 2008:265-273. (Two
  releases that each satisfy a syntactic condition can be composed.)

- Giomi M, Boenisch F, Wehmeyer C, Tasnádi B. *A Unified Framework for
  Quantifying Privacy Risk in Synthetic Data.* Proceedings on Privacy
  Enhancing Technologies. 2023;2023(2):312-328. doi:
  [10.56553/popets-2023-0055](https://doi.org/10.56553/popets-2023-0055).
  (Anonymeter, and the control dataset.)

- Raab GM, Nowok B, Dibben C. *Privacy risk from synthetic data:
  practical proposals.* arXiv:2409.04257, 2024, and *Practical privacy
  metrics for synthetic data*, arXiv:2406.16826, 2024. (Keys and
  targets, RepU, DiSCO, DiO; implemented as `synthpop::disclosure()`.)

- Yale A, Dash S, Dutta R, Guyon I, Pavao A, Bennett KP. *Generation and
  evaluation of privacy preserving synthetic health data.*
  Neurocomputing. 2020;416:244-255. doi:
  [10.1016/j.neucom.2019.12.136](https://doi.org/10.1016/j.neucom.2019.12.136).
  (Nearest-neighbour adversarial accuracy, in its train-versus-holdout
  form.)

- Guillaudeux M, Rousseau O, Petot J, et al. *Patient-centric synthetic
  data generation, no reason to risk re-identification in biomedical
  data analysis.* npj Digital Medicine. 2023;6. doi:
  [10.1038/s41746-023-00771-5](https://doi.org/10.1038/s41746-023-00771-5).
  (Local cloaking, hidden rate, DCR, NNDR.)

- Ganev G, De Cristofaro E. *The Inadequacy of Similarity-Based Privacy
  Metrics: Privacy Attacks Against “Truly Anonymous” Synthetic
  Datasets.* IEEE Symposium on Security and Privacy, 2025.

- Alaa A, van Breugel B, Saveliev ES, van der Schaar M. *How Faithful is
  your Synthetic Data? Sample-level Metrics for Evaluating and Auditing
  Generative Models.* Proceedings of the 39th International Conference
  on Machine Learning (ICML), PMLR 162:290-306, 2022. (α-precision,
  β-recall, authenticity.)

- Snoke J, Nowok B, Raab GM, Dibben C, Slavkovic A. *General and
  specific utility measures for synthetic data.* Journal of the Royal
  Statistical Society Series A. 2018;181(3):663-688. doi:
  [10.1111/rssa.12358](https://doi.org/10.1111/rssa.12358). (pMSE.)

- Karr AF, Kohnen CN, Oganian A, Reiter JP, Sanil AP. *A framework for
  evaluating the utility of data altered to protect confidentiality.*
  The American Statistician. 2006;60(3):224-232. doi:
  [10.1198/000313006X124640](https://doi.org/10.1198/000313006X124640).
  (Confidence interval overlap.)

- Woillard JB, Benoist C, et al. *To be or not to be, when synthetic
  data meet clinical pharmacology: A focused study on pharmacogenetics.*
  CPT Pharmacometrics Syst Pharmacol. 2025. doi:
  [10.1002/psp4.13240](https://doi.org/10.1002/psp4.13240). (Broad
  versus specific utility, weakly correlated.)

- Destere A, Lombardi R, Labriffe M, et al. *Can synthetic data overcome
  the privacy and fidelity bottleneck in Pharmacometrics? A comparative
  benchmark using a daptomycin population pharmacokinetic model.*
  medRxiv preprint, posted June 2, 2026. doi:
  [10.64898/2026.05.30.26354512](https://doi.org/10.64898/2026.05.30.26354512).
  (The daptomycin confidence-interval-overlap benchmark.)
