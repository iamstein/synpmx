# Synthetic Pharmacometric Data: How It Is Made, How It Is Checked, and How synpmx Fits

Synthetic datasets in pharmacometrics (PMX) are useful for open-source
software development, AI-assisted software development, teaching, and
external collaborations. Several different approaches exist for
generating them, and they optimize for genuinely different things —
which means the question is never “which method is best” but “best at
what.”

This article has two halves. The first surveys the four generation
families in use, explains why differential privacy sits across all of
them rather than inside one, and says where `synpmx` fits. The second is
a tutorial on the other literature — the one about **checking**
synthetic data once it exists, which is a separate field with its own
vocabulary, and which is easy to reinvent badly if you have not read it.
Current as of August 2026.

## How synthetic data is made

### Mechanistic simulation

Trial simulation has long been the standard approach within
pharmacometrics, using tools like `rxode2` or `mrgsolve`. A structural
pharmacokinetic (PK) or pharmacodynamic (PD) model, a set of parameters
with between-subject variance, and a trial design produce as many
virtual subjects as you like.

**Strengths**

- No patient record is read, *provided the model is public* (see below)
- Mechanistically interpretable: every feature of the output traces to a
  model assumption you chose
- Unlimited dataset size, and any design you can describe

**Limitations**

- Simulations often lack the irregularities of operational clinical
  datasets: unexpected dosing histories, protocol deviations, and
  inconsistent observation schedules. These mechanisms *can* be built
  in, but doing so makes the simulation considerably more complex, and
  accurate models of the deviation processes are themselves hard to
  develop
- The output is only as good as the model, and a wrong structural
  assumption produces confidently wrong data

The privacy claim deserves care, because it is often overstated.
Simulation is perfectly private only when the model it simulates from is
public. A model *fitted to the trial you are trying to protect* carries
information about those patients in its parameter estimates, and
publishing simulations from it is a release of that information, not an
alternative to one. This distinction is why `synpmx` separates a mode
that reads no data at all from modes that read data under a budget.

#### Copulas: realistic covariates to simulate from

A structural model needs covariates to simulate on — age, body weight,
renal function, liver function, laboratory values, disease
characteristics — and drawing each one independently from its own
marginal distribution produces patients who do not exist: 45 kg adults
with the creatinine clearance of a 90 kg one.

A **copula** solves this by separating a joint distribution into two
parts that can be handled independently: the marginal distribution of
each variable on its own, and the *dependence structure* linking them.
Fit the dependence structure on real data, keep the marginals, and you
can then draw whole covariate vectors in which the relationships between
variables are realistic even though no vector belongs to a real patient.
Zwep and colleagues developed this for virtual patient simulation, and
Guo and colleagues for realistic virtual adult populations.

The important thing about copula methods, for the purposes of this
survey, is what they produce: **one row per virtual patient — a
covariate vector, not a longitudinal profile.** They do not compete with
mechanistic simulation, they supply its input. A published copula also
travels well: it can be shared without the underlying data, which is
what makes the combination a genuinely strong option for a fully
synthetic study.

### Deep generative models

Recent advances in machine learning have introduced deep generative
models for synthetic PK/PD data generation. Methods explored include
generative adversarial networks (GANs), TimeGAN for sequences,
variational autoencoders (VAEs), diffusion models, and probabilistic
autoregressive networks. Jiang and colleagues benchmark several of these
on PK/PD data; Gadgil and colleagues apply diffusion models to virtual
populations and pharmacometric simulation.

These models aim to reproduce the statistical properties of longitudinal
PK datasets while maintaining realistic relationships among patients and
observations.

**Strengths**

- Flexible nonlinear modeling with no structural assumption imposed in
  advance
- Learn dependence structures nobody specified, including ones nobody
  noticed
- Actively developed, with a large methods literature to draw on

**Limitations**

- Fidelity is an aim rather than a result. Where it has been measured in
  pharmacometrics the answer is mixed: Jiang and colleagues find
  performance varies substantially by method and scenario, and Woillard
  and colleagues found TVAE clearly behind CT-GAN and AVATAR on the same
  data
- Data-hungry in a field whose cohorts are often 12 to 200 subjects
- In the pharmacometric applications published so far, the input is
  usually a simulated dataset or a simplified longitudinal concentration
  profile rather than an operational trial dataset. This is the gap this
  package is aimed at, and it is worth being precise that the claim is
  about the applications cited here rather than about the whole
  machine-learning literature

### Sequential conditional models: synthpop and tabular tools

Outside pharmacometrics, the most widely used synthetic-data approach is
neither mechanistic nor deep. The `synthpop` R package (Nowok, Raab and
Dibben) synthesizes a table **column by column**: choose an order for
the variables, model each one conditional on the variables already
synthesized — usually with a classification and regression tree (CART) —
and draw from that fitted conditional. Repeat until every column is
synthetic.

It is the default tool in official statistics and social science, and
deserves to be better known in this field. The same niche is occupied by
the tabular deep-learning tools — CT-GAN and TVAE, from the Synthetic
Data Vault project — which Woillard and colleagues benchmark against
AVATAR.

**Strengths**

- No distributional assumptions, and mixed variable types handled
  natively
- Fast, mature, well documented, and equipped with its own utility and
  disclosure diagnostics
- Sequential conditioning preserves relationships between columns by
  construction

**Limitations**

- It is built for **rectangular, one-row-per-unit data**. A
  pharmacometric event table is not that: it has many rows per subject,
  ordered in time, with a grammar in which `EVID`, `AMT`, `CMT`, and
  `DV` mean different things on different rows
- Synthesizing such a table column by column treats each row as an
  independent unit. That destroys the within-subject trajectory and
  readily produces illegal event sequences — an observation before the
  first dose, a dosing row carrying a measurement, a subject whose times
  do not increase

You can flatten a study to one row per subject and synthesize that
successfully — which is exactly what Woillard and colleagues do, on 253
patients with a single measurement each. What you cannot do that way is
keep the longitudinal endpoint, which for pharmacometrics is the data.

### Record-based blending

Another class of algorithms generates synthetic records directly from
existing ones, by blending or sampling among neighboring individuals
within the original dataset. The original method here is AVATAR, due to
Guillaudeux and colleagues, in which each synthetic record is built from
the local neighborhood of real records. “AVATAR” is a method name rather
than an initialism. Destere and colleagues benchmark a modified AVATAR
against differentially private alternatives on a population PK model,
and Woillard and colleagues include a simplified AVATAR in the
pharmacogenetics comparison above.

**Strengths**

- High realism, because the output is made of real data rather than of a
  model of it
- Empirical distributions and complex, unmodeled relationships are
  preserved without anyone having to notice them first
- Works at small cohort sizes, where deep methods struggle

**Limitations**

- No formal privacy guarantee. Protection is a set of mechanisms whose
  coverage must be argued and measured, not a bound that can be stated
- Real datasets contain operational characteristics that identify people
  on their own — actual dosing times, actual amounts, observation
  schedules, protocol deviations — and blending measurement *values*
  does nothing about any of them. These require separate masking
- The output is assembled from real trajectories, so it inherits the
  source data’s handling obligations wherever it goes

### Differential privacy is a guarantee, not a family

Differential privacy (DP) is frequently listed alongside the approaches
above, which is a category error worth undoing. DP is not a way of
generating data. It is a property a *release mechanism* can have: a
bound, epsilon, on how much the presence or absence of any single
individual can change the distribution of what is released. Any of the
families above either has such a bound or does not.

Three things follow that matter in practice.

- **It is the only claim that survives a determined recipient.** Every
  other protection in this article is an argument about what an attacker
  is likely to manage. DP is a statement about what is possible.
- **It composes, and it is spent.** Two releases from the same data cost
  the sum of their epsilons. A budget is a finite resource, which is why
  it is allocated deliberately across the quantities being released.
- **It is not free.** The guarantee is purchased with noise, and at
  pharmacometric cohort sizes the price is steep. A study of 12 subjects
  cannot hide one subject cheaply.

The practical consequence is that DP and record-based realism sit at
opposite ends of one axis, and choosing between them is a question about
where the data is going rather than about which method is better.
Destere and colleagues measure exactly this tradeoff. `synpmx` treats it
as the primary design decision; the reasoning is in the [privacy
article](https://iamstein.github.io/synpmx/articles/synpmx-privacy.html),
with the formal argument and the cohort-size feasibility analysis
alongside it.

### Where synpmx fits

`synpmx` generates **event-based pharmacometric datasets** — dosing and
measurement tables with the schema, event grammar, and rough behavior of
a real study — for software development, package testing, reproducible
examples, educational materials, AI-assisted programming, and workflow
validation.

The emphasis is on the structural characteristics that software
encounters when processing real trial data. Concretely, that means a
synthetic dataset should contain the patient who missed the week-4
visit, had a dose interrupted and restarted, and withdrew at week 12 —
because the assembly script, the diagnostic plot, and the model-run
plumbing all have to survive that patient, and a clean simulation never
contains them.

#### It is in the record-based family, and inherits its problems

The honest placement is the one that costs the most to defend: the
default method,
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md),
is record-based blending, and every limitation listed in that section
applies to it. What distinguishes the package is not escaping those
limitations but treating them as the specification.

- Values are protected by blending across at least `k` = 5 donors with
  no donor exceeding `max_donor_weight` = 0.50 of any avatar
- Observation times are coarsened onto the declared protocol grid, which
  is what turns a cohort of individually identifiable sampling schedules
  into a shared one
- A pattern of attended visits held by fewer than two real patients is
  refused rather than copied, so no avatar wears one person’s history of
  absences
- Structurally extreme subjects are screened out of the anchor pool, and
  donors are never blended across routes of administration
- Every one of these mechanisms reports what it actually did on the run,
  in counts, so the argument is auditable rather than asserted

The mechanisms are specified in the [AVATAR algorithm
article](https://iamstein.github.io/synpmx/articles/avatar-algorithm.html),
measured across eight public datasets in the
[evaluation](https://iamstein.github.io/synpmx/articles/avatar-evaluation-public-data.html),
and the questions to ask of any generated dataset — from any method —
are in [the checks
article](https://iamstein.github.io/synpmx/articles/synthetic-data-checks.html).

#### The four modes across the families

The package offers four generation modes, and they do not all sit in the
same family. Choosing among them is a privacy decision, not a technical
preference.

| Mode | Family | Formal guarantee |
|----|----|----|
| [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md) | Record-based blending | None; mechanisms, measured |
| [`synpmx_prior()`](https://iamstein.github.io/synpmx/reference/synpmx_prior.md) | Mechanistic simulation from a public model | epsilon = 0 — no data is read |
| [`synpmx_calibrated()`](https://iamstein.github.io/synpmx/reference/synpmx_calibrated.md) | Mechanistic, with magnitude corrected under a budget | (epsilon, delta) DP |
| [`synpmx_empirical()`](https://iamstein.github.io/synpmx/reference/synpmx_empirical.md) | Statistical, rebuilt from DP summaries | (epsilon, delta) DP |

#### What is genuinely different

The families answer different questions, and the third column is where
this package lives.

| Family | Primary goal | Typical methods |
|----|----|----|
| Mechanistic | Simulate the biology | `mrgsolve`, `rxode2`, copulas for covariates |
| Statistical | Match the observed distributions | GANs, diffusion models, `synthpop`, CT-GAN |
| Structural | Preserve how the dataset behaves | AVATAR-style blending on event tables |

Mechanistic methods are indispensable for simulation studies.
Statistical methods excel at reproducing distributions and supporting
data sharing. What neither is built to do is keep the event table intact
— the dosing history, the irregular visit, the multiple endpoints on one
grid, the deviation — because for their purposes those are nuisance
rather than signal. For software engineering and methodological
development they are the entire point.

`synpmx` is not intended to replace any of these. It fills the gap
between mechanistic simulation and statistical synthetic data
generation, and it is explicit about what that costs: blending shrinks
variance, dose amounts are not always rescaled to each avatar’s own
covariates, and no output of the default mode is anonymous data. It is
built for developing software against realistic data, not for estimating
parameters from it.

## How synthetic data is checked

Everything above is about *making* the data. This half is about
*checking* it, and it is a separate literature — one that is scattered
across three fields that barely cite each other: official statistics,
machine learning, and data protection law. Each has its own name for the
same idea, which is the main reason the area is hard to enter.

It is worth reading before inventing anything, and this section is
written as a tutorial rather than a survey for that reason. The last two
parts say which of these measures matter for `synpmx`’s use case and
which are deliberately out of scope.

### The one idea underneath all of it: the patient versus the population

A synthetic dataset is *supposed* to tell you about the population. That
is what it is for. It must not tell you about a patient. Every privacy
check in this literature is an attempt to separate those two, and nearly
every naive check fails because it cannot separate them at all.

Here is why the separation is hard. Suppose your trial found that
patients over 80 kg cleared the drug slowly, and the synthetic data
reproduces that. Good — that is a population fact, and reproducing it is
the point. Now suppose exactly one patient in your trial weighed 138 kg,
and the synthetic data contains a 137 kg patient. Is that a leak?

**You cannot tell from the synthetic data alone.** If people of 138 kg
are common in the world, it is not. If your trial’s 138 kg patient is
nearly the only person who could have qualified for it, it is. The
synthetic dataset looks identical in both cases.

Now make it a measurement problem. Any statistic computed from the
source and the synthetic dataset together — any distance, any
distribution comparison, any “how similar are they” number — confounds
two completely different things:

- the generator **captured the population**, which is success, and
- the generator **memorized individual patients**, which is failure.

Both make the synthetic data resemble the source. A similarity statistic
cannot tell you which one you are looking at, because both push it the
same direction.

#### The fix, and it is the same fix everywhere

**Hold some patients out.** Split the real cohort in two: a *training*
set the generator is allowed to use, and a *control* (or holdout) set it
never sees. Then ask your similarity question twice — synthetic against
training, synthetic against control.

Whatever the generator learned about the **population** shows up in both
comparisons, because the control patients are drawn from the same
population. Whatever it learned about **specific training patients**
shows up only in the first. The **difference between the two is the
leak**, and the raw value of either one on its own is close to
meaningless.

If that feels familiar, it should: it is exactly the train/test split
from machine learning, run for the opposite purpose. In modelling,
accuracy on training data confounds “learned the pattern” with
“memorized the examples”, and a test set separates them. Privacy
measurement is the same problem with the sign flipped — here
memorization is the thing you are trying to *detect* rather than to
avoid, but the confound and the remedy are identical.

Four independent lines of work arrive at this correction, which is the
strongest signal in this whole area that it is not optional:

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

#### What it costs in pharmacometrics, honestly

A holdout is expensive at our cohort sizes. Removing 20% of 200 patients
is routine; removing 20% of 20 is not, and at 12 subjects it is not
possible at all — the donor pool is already the binding constraint. This
is a real reason the technique is rare in this field and not a reason it
is wrong.

The practical compromise is to treat it as a **per-study validation
exercise** run once, deliberately, to characterize a generator on a
dataset — not as something every production run does. You lose donors
for the run that measures privacy, and you keep them for the run you
ship.

### Adversarial accuracy, and why the holdout is what makes it mean something

This is the statistic
[`compare_pmx_proximity()`](https://iamstein.github.io/synpmx/reference/compare_pmx_proximity.md)
implements, so it is worth understanding properly.

**The setup.** Put every real subject and every synthetic subject into
one space, where each subject is a point — for `synpmx`, the profile
space built from covariates and trajectory features. Now, for each
point, ask a single question: *is my nearest neighbour in my own
dataset, or in the other one?*

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

**How to read it.** The scale is the useful part, and both ends are
failures:

- **AA near 1** — every subject’s nearest neighbour is in its own
  dataset. The two clouds are separated. The synthetic data is nothing
  like the real data: a *utility* failure.
- **AA near 0** — every subject’s nearest neighbour is in the *other*
  dataset. The two clouds are interleaved so tightly that each real
  patient has a synthetic partner closer than any of their real peers.
  That is the signature of **memorization**: a privacy failure.
- **AA near 0.5** — a synthetic subject is no more like a real subject
  than two real subjects are like each other. This is the target, and it
  is a target rather than a maximum.

This two-tailed reading is why the statistic is worth having: most
diagnostics have picked a side without saying so, and this one reports
both failure directions from one number.

**Now the part that needs the holdout.** In Yale and colleagues’
original formulation, adversarial accuracy is computed **twice** — once
against the training data the generator saw
($`\mathrm{AA}_{\text{train}}`$) and once against a holdout it did not
($`\mathrm{AA}_{\text{test}}`$) — and the reported privacy loss is the
difference:

``` math
\text{privacy loss}=\mathrm{AA}_{\text{test}}-\mathrm{AA}_{\text{train}}.
```

The logic is the confound from the previous section. A generator that
captured the population produces synthetic subjects that sit equally
close to training and holdout patients, so the two accuracies agree and
the difference is zero. A generator that memorized its training patients
produces synthetic subjects that sit *unusually* close to the training
patients specifically, which drags $`\mathrm{AA}_{\text{train}}`$ toward
0 while leaving $`\mathrm{AA}_{\text{test}}`$ where it was. A positive
difference is leakage, in units you can compare across datasets.

**What `synpmx` currently does, and what it does not.**
[`compare_pmx_proximity()`](https://iamstein.github.io/synpmx/reference/compare_pmx_proximity.md)
computes $`\mathrm{AA}_{\text{train}}`$ only, and calibrates it against
a null interval built by recomputing the same statistic on random halves
of the source cohort. That null does real work: at 20 patients in a
30-dimensional profile space, nearest-neighbour statistics are wildly
unstable and 0.5 is *not* the right expectation, so a raw value is
uninterpretable without it. The null interval tells you what this
statistic does on data of this size and shape when nothing is wrong.

What the null interval cannot do is tell you how much of the resemblance
was legitimate. It corrects for cohort size and dimension. It does not
correct for the population/patient confound, because both halves of the
source were seen by the generator. That is what a holdout would add, and
it is the single largest gap in the package’s privacy measurement today.

### What the AVATAR authors measure: local cloaking and hidden rate

Guillaudeux and colleagues built two privacy metrics that only make
sense for a **patient-centric** generator — one where each synthetic
record is built from one identified original record.
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
is exactly that: every avatar has an anchor. A global generator like a
GAN cannot compute these at all, because there is no correspondence to
compute them over.

The question they ask is different from adversarial accuracy’s, and
better suited to acting on:

> For this particular real patient, how well hidden is their own avatar?

**Local cloaking.** Take real patient $`i`$ and their avatar $`a(i)`$.
Count how many *other* avatars lie closer to $`i`$ than $`a(i)`$ does.
That count is patient $`i`$’s local cloaking. The picture is a crowd:
your own avatar is somewhere in it, and local cloaking is how many other
avatars an attacker would have to walk past to reach yours.

- **Local cloaking 0** is the worst case: your own avatar is the avatar
  nearest to you, so an attacker who has your real record and guesses
  “the closest avatar is derived from me” is right.
- Higher is better. The published values are a **median local cloaking
  of 11** on an AIDS clinical trial dataset and **24** on a breast
  cancer dataset.

**Hidden rate.** The percentage of real patients whose own avatar is
*not* the avatar nearest to them — that is, the percentage with local
cloaking of at least 1. It is the aggregate of the same measurement, and
it is exactly the failure probability of that simple attack. The
published values are **93% and 94%**.

#### Why these are worth having even though we already have adversarial accuracy

They answer a different question at a different granularity, and the
granularity is the point.

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

Hidden rate is also, in flavour, a **membership inference** measure: it
is the error rate of an adversary trying to establish a link between a
known individual and the release. That connects it to the strongest
empirical privacy tests in the general literature.

#### Two honest caveats

**The pairing map is the most disclosive artifact in the whole
pipeline.** Local cloaking and hidden rate require knowing which avatar
came from which patient. That map is a re-identification key: anyone
holding it does not need to attack anything. It has to be computed
inside the trusted environment, used, and destroyed — never attached to
the released data, never written next to it. `synpmx` does not currently
retain the anchor of each avatar at all, which is why these metrics are
not implemented; adding them means creating that map deliberately and
handling it more carefully than anything else in the package.

**They measure one specific attack.** Hidden rate is the failure rate of
an attacker who has a real record, has the whole avatar set, and guesses
nearest-neighbour. A smarter attacker exists. These are useful,
concrete, and not a bound — the same caveat that applies to every
measure in this section except differential privacy.

The same paper also reports two generic distance measures worth knowing
by name, since they appear everywhere:

- **Distance to closest record (DCR)** — for each synthetic record, the
  distance to the nearest real one. Larger is safer; zero is a copy.
- **Nearest-neighbour distance ratio (NNDR)** — the ratio of the nearest
  distance to the second-nearest. Near 1 means the record sits in a
  crowd with no single conspicuous partner; near 0 means it has one real
  record much closer than everything else, which is the shape of a copy.
  The paper treats 0.8 and above as satisfactory.

Both are cheap. Both are also, on their own, weaker than they look — see
the critique two sections down.

### The official-statistics tradition: keys, targets, RepU and DiSCO

The most developed disclosure-risk framework for synthetic data comes
from official statistics rather than from machine learning, and it ships
as working code: Raab, Nowok and Dibben formalized it in 2024 and
implemented it in `synthpop` (version 1.8.1 and later) as `disclosure()`
and `multi.disclosure()`.

#### The vocabulary: keys and targets

Everything in this framework is organized around a distinction worth
adopting regardless of which tools you use.

- A **key** — the older literature says *quasi-identifier* — is a set of
  variables an attacker plausibly **already knows** about the person
  they are looking for. Age band, sex, treatment arm, site, date of
  first visit. Not identifying on their own; identifying in combination.
- A **target** is the variable they want to **learn**: a genotype, a
  diagnosis, a lab value.

Disclosure is then a concrete, checkable event: *the key lets them pin
down the target.* And the choice of key is a judgement call by the data
holder, not something a package can infer — which the paper is explicit
about, and which is the framework’s main practical weakness. New
external data sources create new keys, so an assessment has a shelf
life.

#### Identity disclosure: RepU

*Can an attacker work out which record is this person?*

**RepU — replicated uniques** — is the percentage of original records
that are unique on the key set in the **original** data *and* also
unique on that same key in the **synthetic** data.

The two conditions do different jobs, and this is the elegant part.
Being unique in the original is the patient’s pre-existing exposure — it
is not something the release caused. Being *replicated* as a unique in
the synthetic data is what the release **added**: the attacker who knows
the key can now find a single record, and its other values are being
asserted about that person. `synthpop` will optionally remove replicated
uniques from the output, which is a disclosure control rather than a
diagnostic.

#### Attribute disclosure: DiSCO

*Can an attacker learn a value about this person, without necessarily
picking out their record?*

This is the more important one, and it is subtler because it does not
require identifying anybody. If every synthetic patient in a given key
cell has the same target value, then an attacker who knows a real person
is in that cell learns their value without ever locating their record.

**DiSCO — Disclosive in Synthetic, Correct in Original** — is the
percentage of original records whose key combination maps to a
**single** target value in the synthetic data, and where that value is
the one the real person **actually has**. Two conditions again: the
synthetic data has to make a confident assertion, and the assertion has
to be right.

#### The correction that matters most: DiO

Here is the part that closes the loop back to the first section.

**DiO — Disclosive in Original** — is the same computation run on the
original data. If the key already determines the target in the real
data, then the attacker did not need your synthetic dataset: they could
have derived the relationship from published literature, from the
protocol, or from clinical knowledge. That is a **population fact**. It
is the science.

So the quantity that measures what the *release* added is not DiSCO but
**DiSCO − DiO**.

This is the patient-versus-population separation again, in a third
guise. Anonymeter does it with a control **dataset**; Yale does it with
a holdout; Raab does it with a control **computation** on the original
data. Same idea, three implementations, and once you see it you see it
everywhere.

It also has a direct analogue already in `synpmx`, arrived at
independently: the identifiability screen scores each patient **within
their treatment arm** rather than against the whole cohort, because on a
six-arm dose-ranging study a cohort-wide screen flags thirty patients
for having received the dose their arm was assigned. That is a protocol
fact being reported back as a privacy finding — DiO, in a different
costume. `DiSCO − DiO` is that same insight expressed as arithmetic.

#### How well it fits pharmacometric data

Better than you would expect, in one place, and not at all in another.

At the level of **subject baselines** — one row per patient, covariates
and strata — a pharmacometric dataset *is* rectangular, which is exactly
the shape `synthpop`’s machinery assumes. The covariate table can go
straight into `disclosure()` with keys chosen from arm, sex, age band,
site and so on. This is the natural home for the rare-category question.

At the level of the **event table** it does not apply, for the reasons
the first half of this article gives: many rows per subject, ordered in
time, with a grammar the tools do not model. That is where the timing
and structural checks have to do the work instead, and where there is
essentially no external literature to borrow from.

### The data-protection tradition: the Article 29 Working Party criteria

#### What the Article 29 Working Party is

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

The reason to know these three is that they are the criteria a data
protection officer will actually apply, and they cut the space
differently — and more completely — than an engineer’s intuition does.

#### Linkability, since it is the least intuitive

Linkability is the ability to connect two records that belong to the
same person. It is worth stating plainly that **it does not require
identifying anybody.** Establishing that record 47 in your synthetic
dataset and record 12 in some other release describe the same human
being is a failure on its own, because such links accumulate: each one
narrows the field, and enough of them identify.

The concrete attack for synthetic data: an adversary holds a *partial*
record from elsewhere — a covariate row from a registry, a table in a
publication, or simply their own knowledge as the person’s treating
clinician — and asks which synthetic record corresponds to it.
Anonymeter implements this directly: split the attributes into two sets,
give the attacker one set, and see whether they can recover the right
partner.

**A record-based generator is more exposed to this than any other
family**, because a link between an avatar and its anchor exists *by
construction*. That is not a hypothetical: it is precisely the link that
local cloaking measures. The two literatures are describing the same
risk from opposite ends, and the `synpmx` checks vignette currently
covers singling out thoroughly, inference partially, and linkability not
at all.

#### Anonymeter

Giomi and colleagues turned the three criteria into a working evaluation
suite — one concrete attack per criterion — with the control dataset
from the first section built into every one of them. Risk is reported
only where the attack performs better against the training data than
against the control, which is what stops the tool from reporting the
population back to you as a leak.

If you read one thing about attack-based privacy evaluation, this is the
one to read: it is short, it is implemented, and the framing is the one
that regulators already use.

### The critique everybody should read before trusting a distance

Ganev and De Cristofaro (2025) attack the whole family of
similarity-based metrics — distance to closest record, nearest-neighbour
ratios, and their relatives. Their result is that synthetic datasets
which pass these metrics comfortably can still be attacked successfully,
with training records reconstructed.

The reason is structural rather than a matter of tuning, and it is easy
to state:

> These are **average-case** statistics over a whole dataset. An
> attacker targets **one** record, and only needs to succeed once.

A synthetic dataset can be, on average, a perfectly respectable distance
from its source and still contain one record that is a near-copy of one
patient. The mean does not see it. Nothing about a better threshold
fixes this, because the statistic is answering a different question from
the one being asked.

The practical conclusions are worth adopting even if you never read the
paper:

- **Distance metrics are a smoke alarm, not a safety certificate.** They
  reliably catch blatant copying — hand
  [`compare_pmx_proximity()`](https://iamstein.github.io/synpmx/reference/compare_pmx_proximity.md)
  a verbatim copy of the source and it goes to zero and objects, which
  is what the package’s own regression test requires of it. That is a
  real and useful job.
- **Prefer per-record measures over cohort averages** wherever both
  exist. This is the same argument that makes local cloaking more
  actionable than adversarial accuracy.
- **The empirical gold standard is a membership inference attack**:
  determine whether a specific individual was in the training data. It
  is expensive, it needs a holdout, and it is the thing distance metrics
  are a cheap proxy for.
- **The only claim that survives all of this is a formal one**, which is
  the argument for the differentially private modes.

### The machine-learning tradition: fidelity, diversity, and authenticity

The generative-modelling literature has one framing worth importing,
from Alaa and colleagues (2022). Their observation is that “is the
synthetic data good?” is three questions that get collapsed into one,
and collapsing them hides which thing broke:

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

#### “Did we lose the tails” — and do we want to?

Mostly **yes, we want to**, and the reason is worth being precise about.

Coverage and privacy are in direct opposition **at the tails**, because
in a clinical dataset the tail *is* the identifiable patient. The
patient with the longest follow-up, the highest dose, the extreme weight
— these are the ones a masking mechanism is designed to remove or blend
away. β-recall falling at the extremes is the masking **working**.

So coverage is not a target to hit. What it is good for is **telling two
causes apart**, and they are causes with opposite implications:

- The tail is missing because the outlier was deliberately screened out
  or blended away. That is a designed cost, it is on purpose, and
  `synpmx` already reports it directly — anchors screened out, patterns
  dropped, regimens not represented.
- The tail is missing because the generator collapsed toward the mean
  and produces nothing but average patients. That is a defect, and for a
  package whose entire purpose is producing the *awkward* patients that
  break software, it is close to a total failure.

**In a standard deviation these two look identical.** That is the
argument for measuring coverage at all: not to optimize it, but so that
a shrinking spread can be attributed. Report it, do not chase it.

#### Authenticity

Authenticity is the “too close” direction, per record: it asks, of each
synthetic record individually, whether it is essentially a copy of some
training record’s neighbourhood rather than a new draw.

Note where it lands: **authenticity is the machine-learning literature’s
version of local cloaking.** Different field, different derivation, same
question — *is this particular output too close to a particular input?*
— and the same advantage over an aggregate distance. Three traditions,
three names, one measurement.

### Utility measures, and which of them are not your problem

This is where the phrase “evaluates on the analysis, not the marginals”
comes from, so here is the distinction it refers to.

#### General utility: pMSE

**General** — or *broad* — utility asks whether the synthetic dataset
resembles the real one at all, without reference to any particular
analysis. The standard measure is **pMSE**, the propensity score mean
squared error (Snoke, Nowok, Raab, Dibben and Slavkovic, 2018).

The mechanism is simple and worth knowing because the same trick appears
under four different names:

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

You will meet the same idea as the *discriminative score* in the
time-series literature and as the *classifier two-sample test* in
statistics. It has the same logical shape as adversarial accuracy — a
“can you tell them apart?” statistic with a calibrated null — which is a
good sign that this shape is the right one for the question.

#### Specific utility: confidence interval overlap

**Specific** — or *narrow* — utility asks whether **the analysis you
actually intended to run** gives the same answer on the synthetic data
as on the real data. The standard measure is **confidence interval
overlap** (Karr and colleagues, 2006): run the analysis on both
datasets, and measure how much the two confidence intervals for the same
quantity overlap.

In pharmacometrics that means fitting the same population model to both
and comparing the fixed effects, the between-subject variance terms, and
their intervals. Destere’s daptomycin benchmark is precisely this, and
it is the question a modeller asks first.

#### Why the distinction matters: they come apart

This is what “evaluates on the analysis, not the marginals” means.
Comparing marginal distributions asks whether the columns look right
**one at a time**. Comparing the analysis asks whether the
**conclusion** is the same. It is entirely possible for every column to
match beautifully while the relationships between them — which is what
any model estimates — are destroyed.

This is measured, not speculative: Woillard and colleagues found that
broad utility and specific utility were only **weakly correlated** in
their benchmark. Data that looked right column by column did not support
the same predictions. The lesson to carry away is the negative one:
**marginal agreement is weak evidence**, and a document that reports
only distribution comparisons is reporting the easy half.

#### And why almost none of this is `synpmx`’s problem

Now the scope statement, because it changes which of these measures are
worth implementing.

`synpmx` exists so that **software** can be developed against realistic
pharmacometric data. Under that use case, “would a modeller reach the
same scientific conclusion from this data?” is not a requirement that
has been weakened — it is **not a requirement at all**, because nobody
should be reaching scientific conclusions from this data in the first
place, and the package says so in every document it ships.
Between-subject variance is knowingly shrunk by blending. Dose amounts
are not always rescaled. These are documented properties, not defects
awaiting a fix.

What the output does have to be is **the same kind of object** as a real
study dataset: a legal event grammar, dosing histories that a protocol
could have produced, several endpoints on one visit grid, patients who
missed visits and withdrew early, and a schema an assembly script can be
pointed at.

That reorders the whole list above:

- **Structural and validity checks matter most**, and they are the part
  this literature barely addresses — because tabular and
  machine-learning tools do not have event tables, so nobody outside
  pharmacometrics has needed to check one. There is little to borrow
  here and it is where the package’s own work has to be original.
- **Privacy measures matter in full**, unchanged. The source data is
  real regardless of what the output is used for, and none of the
  reasoning in the first six sections is affected by the use case.
- **General utility measures such as pMSE are interesting but not
  binding.** Worth knowing; not worth chasing. A poor pMSE would prompt
  a look at *why*, not a change to hit a number.
- **Specific utility measures are out of scope by design** — confidence
  interval overlap, parameter recovery, non-compartmental exposure
  agreement. Not because they are unimportant, but because the package
  does not claim the property they measure. Saying that plainly is more
  useful than reporting them badly and letting a reader infer a claim
  that was never made.

The one thing in the utility column that *is* binding does not appear in
this literature at all, and no function can compute it: **does the
pipeline that will consume the real study run unchanged against this?**
Only the person holding that pipeline can answer it.

### Where to start

If you read four things, in this order:

1.  **WP29 Opinion 05/2014**, section 2.2 — the three criteria, in about
    three pages. It gives you the vocabulary everyone else assumes.
2.  **Giomi and colleagues (Anonymeter)** — the control-group idea, and
    one concrete attack per criterion.
3.  **Raab, Nowok and Dibben (2024)** — keys and targets, RepU, DiSCO,
    and the DiO correction. Then run `synthpop::disclosure()` on
    something.
4.  **Guillaudeux and colleagues** — local cloaking and hidden rate, the
    patient-centric measures that apply directly to this package.

Then **Ganev and De Cristofaro** for the critique of everything in items
2 through 4, which will stop you over-trusting any of it.

#### The vocabulary, in one table

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

The last row is the only one that is a bound rather than a measurement,
which is the distinction the differential privacy section above is
about.

`synpmx`’s own worked checks, and the gaps between them and this list,
are in [the checks
article](https://iamstein.github.io/synpmx/articles/synthetic-data-checks.html).

## References

#### Generating synthetic data

- Jiang Y, García-Durán A, Losada IB, Girard P, Terranova N. *Generative
  models for synthetic data generation: application to
  pharmacokinetic/pharmacodynamic data.* J Pharmacokinet
  Pharmacodyn. 2024. doi:
  [10.1007/s10928-024-09935-6](https://doi.org/10.1007/s10928-024-09935-6).

- Gadgil PK, Poojari SM, Ramanathan M. *Diffusion models for virtual
  populations and pharmacometric simulations.* J Pharmacokinet
  Pharmacodyn. 2026;53(5):45.

- Zwep LB, Guo T, Nagler T, Knibbe CAJ, Meulman JJ, van Hasselt JGC.
  *Virtual Patient Simulation Using Copula Modeling.* Clin Pharmacol
  Ther. 2024;115(4):795-804. doi:
  [10.1002/cpt.3099](https://doi.org/10.1002/cpt.3099).

- Guo T, et al. *Generation of realistic virtual adult populations using
  a model-based copula approach.* J Pharmacokinet Pharmacodyn. 2024.
  doi:
  [10.1007/s10928-024-09929-4](https://doi.org/10.1007/s10928-024-09929-4).

- Nowok B, Raab GM, Dibben C. *synthpop: Bespoke creation of synthetic
  data in R.* Journal of Statistical Software. 2016;74(11):1-26. doi:
  [10.18637/jss.v074.i11](https://doi.org/10.18637/jss.v074.i11).

- Xu L, Skoularidou M, Cuesta-Infante A, Veeramachaneni K. *Modeling
  tabular data using conditional GAN.* Advances in Neural Information
  Processing Systems 32 (NeurIPS 2019). 2019:7335-7345. (CT-GAN and
  TVAE, from the Synthetic Data Vault project.)

- Guillaudeux M, Rousseau O, Petot J, et al. *Patient-centric synthetic
  data generation, no reason to risk re-identification in biomedical
  data analysis.* npj Digital Medicine. 2023;6. doi:
  [10.1038/s41746-023-00771-5](https://doi.org/10.1038/s41746-023-00771-5).

- Destere A, Lombardi R, Labriffe M, et al. *Can synthetic data overcome
  the privacy and fidelity bottleneck in Pharmacometrics? A comparative
  benchmark using a daptomycin population pharmacokinetic model.*
  medRxiv preprint, posted June 2, 2026. doi:
  [10.64898/2026.05.30.26354512](https://doi.org/10.64898/2026.05.30.26354512).

- Woillard JB, Benoist C, et al. *To be or not to be, when synthetic
  data meet clinical pharmacology: A focused study on pharmacogenetics.*
  CPT Pharmacometrics Syst Pharmacol. 2025. doi:
  [10.1002/psp4.13240](https://doi.org/10.1002/psp4.13240).

- Dwork C, Roth A. *The Algorithmic Foundations of Differential
  Privacy.* Foundations and Trends in Theoretical Computer Science.
  2014;9(3-4). doi:
  [10.1561/0400000042](https://doi.org/10.1561/0400000042).

#### Checking synthetic data

- Article 29 Data Protection Working Party. *Opinion 05/2014 on
  Anonymisation Techniques.* WP216, adopted 10 April 2014. (The singling
  out / linkability / inference criteria; section 2.2.)

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
