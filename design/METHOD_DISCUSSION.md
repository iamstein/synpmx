# Method discussion: AVATAR blending vs formal differential privacy

Why `synpmx` uses AVATAR-style blending as its primary method, when it
also contains a fully differentially private engine. Written after building
both, measuring both, and comparing to Novartis's `synadam`.

The short version: for **synthetic data that reaches no one the source data
could not**, a resampling method is the right tool, and AVATAR is the
trajectory-level version of exactly what `synadam` already does column by
column. Formal differential privacy is the right tool only when the output
**does reach someone new**, and it is retained here for that case. The boundary
is organizational — who may see the data and under what obligations — not
geographic, so moving output onto a workstation covered by the same controls
stays on the AVATAR side of it.

---

## 1. The two methods in one paragraph each

**AVATAR blending** (the `synpmx_avatar()` engine). For each generated subject,
sample a compatible source subject's event skeleton as a template, then fill its
covariates and endpoint trajectories with a distance-weighted blend of a handful
of similar source subjects, plus subject-level and within-trajectory noise. The
output looks like real data because it is built from real data. There is **no
formal privacy guarantee**.

**Differential privacy** (the `synpmx_calibrated()` and `synpmx_empirical()`
engines). Compute a small number of aggregate statistics, add mathematically
calibrated noise so that no single subject can move any released number by more
than a bounded amount, and generate from those noised aggregates against public
structural priors. The output carries a provable `(epsilon, delta)` guarantee,
at the cost of utility that degrades sharply as the cohort shrinks.

---

## 2. What synadam does, and why it matters here

`synadam` is Novartis's own ADaM synthetic-data package. Its `simulate_vec`
generates each column independently:

- **continuous:** `runif(n, min = observed_min, max = observed_max)` — a uniform
  draw over the observed range;
- **categorical:** `sample(unique_values, replace = TRUE)` — resample the
  observed categories;
- **flags:** sampled in the observed proportions.

There is no differential privacy, no epsilon. The method preserves each column's
marginal support and nothing else, and relies on the surrounding governance —
ADaM data is already analysis-level de-identified, and it stays inside a
controlled environment — rather than on a mathematical guarantee.

This is a completely standard and defensible position. Most real clinical-data
sharing runs on governance and heuristic de-identification (HIPAA Safe Harbor,
Expert Determination), not on differential privacy. `synadam` shipping this way
is evidence that the governance-based model is accepted practice.

**The key observation:** `synadam` resamples each *column* from the data.
AVATAR resamples each *subject trajectory* from the data. They are the same idea
at two different granularities. If `synadam`'s privacy model is acceptable for
its use, AVATAR's is acceptable for the same use — with one important caveat in
section 4.

---

## 3. Why a resampling method is the right default for this package

The package's stated purpose (`design/PROTOTYPE_SPEC.md` section 1) is **synthetic
data for model-workflow exploration**: exercising cleaning, joins, reshaping,
plotting, control-file plumbing, and repeated-dose or longitudinal analysis
code. The accuracy bar is "vaguely right, structurally exact." It is explicitly
*not* parameter estimation, inference, or scientific conclusions.

For that purpose, AVATAR is simply better than the DP engine on every axis the
user cares about:

- **It works at any cohort size.** The DP engine's utility collapses below a few
  hundred subjects (`vignettes/articles/feasibility.Rmd`), which is most early-phase work.
  AVATAR produces plausible data from twelve subjects.
- **It preserves joint structure for free.** Real covariate correlations, real
  trajectory shapes, real timing patterns come through because whole real
  subjects come through (blended). The DP engine has to assert all of that from
  public structural models.
- **It needs no elicitation.** No public priors, no structural model, no
  protocol declaration. Point it at the data and go.
- **It spends no budget**, because there is no budget.

The DP engine, by contrast, spends most of its effort defending against an
adversary who wants to re-identify a patient from the output. If no such
adversary can reach the output — because it never leaves the obligations the
source data already carries — that effort buys nothing.

---

## 4. The asymmetry that makes AVATAR riskier than synadam

Honesty requires stating where the `synadam` analogy strains.

A resampled **covariate value** is weakly identifying. A weight of 72 kg is
shared by thousands of people; releasing it reveals almost nothing about any
individual. This is why `synadam`'s per-column resampling is low-risk, and why
`synpmx`'s bootstrap covariates (`pmx_covariates_auto()`) are a reasonable
default.

A resampled **subject trajectory** is strongly identifying. A full
concentration-time profile — with its particular sampling times, missed visits,
noise pattern, and idiosyncratic shape — is close to a fingerprint. Blending
several donors and adding noise mitigates this, but the mitigation is not free
and not formal: push the noise high enough to defeat a nearest-neighbor linkage
attack and you have destroyed the same signal the DP engine would have
destroyed, without the accounting to prove it. This is the exact failure mode
that made the Version 1 design uncomfortable in the first place
(`vignettes/articles/feasibility.Rmd` section 2).

So AVATAR sits a notch higher on the risk ladder than `synadam`'s column
resampling, for the same governance model. The conclusion is not that AVATAR is
unsafe; it is that **AVATAR's safety depends more heavily on the governance
context**, and the more identifying the trajectory, the more that dependence
matters.

---

## 5. The decision rule

The choice between the two engines is not about which is "more private" in the
abstract. It is about a single question:

> **Does the generated data reach anyone the source data could not?**

- **No** — the same organization, the same access controls and confidentiality
  obligations, wherever the file physically sits → **AVATAR.** It is more
  useful, works at any N, and its lack of a formal guarantee costs nothing
  because there is no adversary to guarantee against. This is `synadam`'s
  situation, and the reasoning is the same.
- **Yes** — shared with a partner, a vendor, a system outside those
  obligations, or published → **the DP engine.** A formal guarantee is the only thing that
  survives a determined adversary, and here there might be one.

Differential privacy is expensive precisely because it defends against someone
who *wants* to break it. Buying that defense when no one can touch the output is
paying for a threat you do not have. Refusing to buy it when the output will be
handed to strangers is negligence. The engines exist side by side so the user
can match the tool to the boundary.

---

## 6. What the package keeps, and why

Both engines remain, deliberately.

- **`synpmx_avatar()`** — AVATAR, the primary and default method. The right
  answer when the synthetic data reaches no one the source data could not, which
  is the common case.
- **`synpmx_calibrated()`** — the structural-correction DP engine. The right
  answer when a formal guarantee is needed and the cohort is small; it asserts
  shape from a public model and privately calibrates only the magnitude.
- **`synpmx_empirical()`** — the dense-grid DP engine. Retained for large pooled
  corpora where its cost is affordable.

Keeping all three is not indecision. It reflects that "synthetic clinical data"
is not one problem with one right method — it is a family of problems separated
by the trust boundary and the cohort size, and the honest package offers the
tool that fits each, clearly labeled with what it does and does not guarantee.

---

## 6a. Minimum donor pooling and event-skeleton sampling (design owed)

Raised by the package owner on 2026-07-25, tracked as `REV-025` (privacy) and
`REV-026` (coherence). This section scopes the problem and the options; nothing
here is implemented yet.

### The defect

AVATAR groups subjects by an **exact** event signature (`R/profiles.R`:
`EVID` + `CMT` + `DVID` + dose/rate magnitude, the count and gaps of dose
starts, and the observed endpoint set), and only blends within a group
(`.select_donors`, `R/synthesis.R:363`). Exact matching buys a real guarantee —
a synthetic subject's trajectory is always structurally coherent with its event
skeleton — but it fragments the cohort:

- A **singleton** group (one subject at a dose, or one unusual schedule such as
  the single long-followed subject in `wbcSim`) falls back to that subject as
  the sole donor, with only subject/residual noise (defaults 0.15/0.05) as
  perturbation. The synthetic subject is a **noised near-copy of one real
  person** — the fingerprint risk the privacy sections warn about, produced
  silently with only a warning.
- A **pair** group blends exactly two real subjects (`k` collapses to 1).

Sparse groups are usually the extreme dose arms and the off-protocol subjects —
exactly the individuals most identifiable. `compare_pmx_distributions()` cannot
catch this because the risk is per-subject, not distributional.

### Owner decisions (2026-07-25)

The owner settled several open questions; recorded here so implementation does
not relitigate them.

- **Floor = 5, and it is a hard floor: `stop()` with an error, not a warning.**
  Loudness is the point — see "Why it was silent" below.
- **No dose-rescaling.** The goal is structural realism for code development, not
  statistical fidelity, so donors pulled from a neighbouring dose are blended
  **raw**, even where that is not physiologically sensible. Optionally a future
  refinement could dose-normalise only DVIDs flagged as PK, but the default
  stays simple: nearest dose, blend as-is.
- **Observations need not be dose-supported.** The earlier draft claimed a
  synthetic observation cannot sit where no dosing supports it; the owner
  rejects that — PD and baseline observations routinely exist without a nearby
  dose. So schedule sampling is freer than first scoped.
- **Dropping unique subjects is acceptable.** When a subject cannot reach the
  floor even after pooling (e.g. the single long-followed `wbcSim` subject),
  dropping it from the source before synthesis is a fine outcome, not only
  erroring.

### On the floor value: what the AVATAR paper actually supports

Guillaudeux 2023 (`references/Guillaudeux23.pdf`) does **not** justify 5. Its
headline setting is **k = 20**; **k = 4 is the lowest value it tested** (range
4–750 for the AIDS dataset, 4–150 for WBCD). Privacy scales with k: *"lower k
values indicated … lower local cloaking; higher k values indicated more
protected individuals."* So 5 sits at the least-private end of the paper's
range. It is defensible here only because the paper's k is nearest-neighbours
across the **whole** dataset (hundreds–thousands of cross-sectional records),
whereas our k is donors **within an event-signature group** in a PMX cohort of
12–60 — you cannot require 20 donors when a dose arm has 6 patients. So 5 is
"about the largest floor feasible at PMX cohort sizes", not a privacy-optimal
choice. The paper names *"dynamic adaptation of k depending on … density"* as
future work, which is exactly the cross-cohort pooling below. Note also the
then-current `max_donor_weight = 0.80` cap: even with 5 donors, one real subject
could still be 80% of a synthetic one, so the floor alone does not bound
individual contribution — the weight cap does, and 0.80 was loose. *(The cap is
now 0.50 and enforced on every donor; see A2 below. This paragraph is left as
the 2026-07-25 record of why the cap was identified as the binding constraint.)*

### Why it was silent, and how to make it loud

The small-group fallback did emit a real `warning()` (`R/synthesis.R:586`), but
as a deferred `warn=0` warning printed unobtrusively after the call — easy to
miss, and trivially removed by `suppressWarnings()`. The shipped demo does
exactly that (`avatar-evaluation-public-data.Rmd`, wbcSim), so the website never shows it. The
hard-floor `stop()` is the fix: an error cannot be suppressed by
`suppressWarnings` and halts execution. Groups that still meet the floor but hit
lesser fallbacks (e.g. `k` reduced from 5 to 5-available) stay warnings.

### A. Reaching the donor floor (REV-025)

Replace signature *equality* with a signature *distance* and expand the donor
set outward — nearest dose first, then nearest schedule, then nearest endpoint
set — until 5 donors are gathered, blending raw (no rescaling, per above). If 5
cannot be reached even after pooling, the subject is either **dropped** from the
source or the run **errors**, by an explicit option; both are acceptable, drop
is the friendlier default for a dataset with a few oddballs like `wbcSim`.

### B. How many events a subject should have (REV-026)

A synthetic subject currently inherits its anchor's event skeleton verbatim
(`R/synthesis.R:531`), so the number and timing of observations equals one real
person's — the reason a unique-schedule subject reappears intact. The
count/timing should instead be **sampled** from the (relaxed) pool: draw a
skeleton from the donor group independently of the value donors, or draw the
observation/visit count from the cohort distribution. Since observations need
not be dose-supported (owner decision), this is less constrained than first
scoped — the main remaining constraint is that dosing events themselves stay
coherent with the regimen. Connects to protocol structure in
`vignettes/articles/data-elicitation.Rmd`.

### A2. Route barrier and donor weight cap — DONE (2026-07-27)

Owner decisions of 2026-07-27, closing the "relax exact compatibility" half of
`REV-025`:

- **Route of administration is an absolute barrier.** IV, bolus, and oral are
  never blended. The distinction the owner drew: a dose-size or schedule
  difference makes a donor a *worse match*, which the fallback may accept, but a
  route difference is a *different experiment* whose blend is a trajectory no
  protocol could produce. Route key (`.route_key()`, `R/profiles.R`) is the set
  of `(EVID, CMT, RATE != 0)` triples on dosing rows — a set, so dose count does
  not enter; NONMEM `RATE < 0` counts as an infusion.
- **No structural distance metric after all.** The earlier scoping proposed
  replacing signature *equality* with a signature *distance*. The owner rejected
  the added machinery: keep exact-signature-first, then nearest in profile
  space, and write the algorithm down explicitly instead. So structure enters as
  a two-stage ordering, not a score, and there is no new tuning knob. The
  explicit statement lives in `articles/avatar-algorithm.Rmd` Step 6.
- **A route arm below the floor is handled by `on_donor_shortfall`**, loudly in
  every branch. With no legal donor left to borrow there is no good answer, only
  a choice between omitting the arm and reproducing it, and only the caller
  knows which matters more. `"drop"` (default) omits those anchors --- the
  owner-approved outcome for subjects that cannot reach the floor. `"noise"`
  keeps them on whatever same-route donors exist plus noise; the owner asked for
  this escape hatch (2026-07-27) on the condition that the alert name it
  explicitly and mark it not recommended, which the `"drop"` and `"error"`
  messages both do. `"error"` refuses and names both alternatives. Not gated by
  `screen`: that is cosmetic, this is a privacy floor.
- **`max_donor_weight` 0.80 → 0.50**, exposed as an argument. The owner's point:
  `k` bounds how many patients are blended, but the cap is what bounds how much
  of any *one* patient lands in an avatar, so the cap — not `k` — is the real
  anonymization parameter. 0.80 let one donor be four-fifths of an avatar.
  Landed at 0.30 first, then revised to **0.50 on 2026-07-27** after the owner
  pushed back that 0.30 felt unnecessary and over-complicated. Both instincts
  were checked against measurement rather than argued: see "Choosing the cap"
  below.
- **The cap now applies to every donor**, by water-filling (`.cap_weights()`).
  Capping only `which.max()` was adequate at 0.80; at 0.30 the runner-up
  routinely exceeds the cap after the leader's excess is redistributed, so the
  documented maximum was violated by the donor the redistribution created. A cap
  below `1/K` relaxes to `1/K` (uniform), so small sources degrade rather than
  error.
- **The `2^(-R)` rank attenuation stays.** It is what stops a capped blend
  collapsing to a flat cohort average, which was the owner's stated worry about
  losing individual variability.

**Measured cost, and a correction.** The independence formula (a blend retains
`sum(w^2)` of individual variance) predicted between-subject SD falling from
~81% of source at cap 0.80 to ~49% at 0.30. That is wrong in practice, and the
error is worth recording so it is not re-derived: donors are *nearest
neighbours* and therefore strongly correlated, so averaging them destroys far
less variance than independence implies, and `subject_noise_sd` restores more.
Measured on `theo_md` (between-subject SD of log AUC, 20 seeds, source 0.273):

| cap | effective donors `1/sum(w^2)` | BSV retained |
|---|---|---|
| 0.80 | 2.50 | 72% |
| 0.50 | 2.93 | 75% |
| 0.30 | 3.91 | 78% |
| 0.25 | 4.43 | 74% |
| 0.20 | 5.00 | 68% |

So the tightening is close to free: effective donors rise 2.5 → 3.9 while BSV is
flat within noise. One dataset, one summary, 12 subjects — enough to justify the
default, not enough to call general. Re-measure on INTERNAL_STUDY. `mean_effective_donors`
and `min_effective_donors` are now recorded in `pmx_settings` so the question can
be asked of any run.

**Choosing the cap — settled 2026-07-27 at 0.50, on evidence.** Two questions
were asked and both were answered by simulating the raw weight formula rather
than by intuition.

*Is a cap needed at all?* Yes, and this is the number that settles it. Uncapped
at k = 5, the largest donor share has median **0.58** (IQR 0.47–0.72), exceeds
0.8 in 14% of avatars, and the effective donor count `1/sum(w^2)` is **2.37**.
So without a cap the floor is largely decorative: nominally five patients are
blended, effectively about two and a half. The cap is the only thing closing
that gap.

*What value?* The decisive diagnostic is **how often the cap fires**, because
that says what role it plays:

| cap | binds | effective donors |
|---|---|---|
| 0.30 | 99% | 3.92 |
| 0.40 | 89% | 3.28 |
| 0.50 | 68% | 2.90 |
| 0.60 | 47% | 2.64 |
| 0.80 | 15% | 2.40 |
| none | 0% | 2.37 |

Both ends are the wrong kind of parameter. A cap firing on 99% of subjects is
not a guardrail — it *is* the weighting scheme, and the inverse-distance term
beneath it stops mattering; that is what the owner was reacting to. A cap at
0.80 fires on 15%, trimming only the tail while still allowing one patient to be
four-fifths of an avatar. **0.50 is the chosen default**: it fires on about two
thirds of subjects, so it genuinely constrains without replacing the distance
weighting, and it states as one checkable sentence — *no single real patient is
more than half of any synthetic patient*. Honest cost, recorded so it is not
rediscovered: effective donors fall to 2.90, so "blends 5 patients" is really
"about 3".

**A simplification that fell out of 0.50.** Multi-pass water-filling is only
needed when redistribution lifts the runner-up over the ceiling, i.e. when
`(1-c)*w2/(1-w1) > c`. At c = 0.5 that requires `w1 + w2 > 1`, which is
impossible, so **at any cap >= 0.5 the loop provably runs at most once**
(simulation: at 0.5 the cap is untouched 33% and pinned once 67%, never twice;
at 0.3 two passes 53%, three passes 11%). The general routine is kept because
the cap is a user argument and below 0.5 a single pass is genuinely wrong, but
at the shipped default the simple thing and the correct thing coincide. This is
worth remembering before anyone "simplifies" `.cap_weights()` back to one pass.

**Both diagnostics are now recorded** in `pmx_settings`:
`cap_binding_fraction` and `mean_effective_donors`, so the question can be asked
of real data rather than of a simulation. On `theo_md` the default binds on 67%
of subjects with 2.82 effective donors, matching the simulation closely.

**Dose--exposure caveat surfaced, not changed.** The owner asked whether stage-2
selection is "just averaging DVs without caring about dosing". Near enough: AMT
is not a profile feature, so the fallback distance compares covariates and the
DV trajectory only. Dose enters *indirectly* --- a higher dose raises
concentrations, and concentrations are profile features, so different-dose
subjects land further apart --- which makes the ranking prefer similar doses
without being told to. But nothing rescales (the standing "no dose-rescaling"
decision), so an avatar carries its anchor's AMT with concentrations possibly
blended across doses, and the dose--exposure relationship is not guaranteed.
Fine for code development, wrong for parameter estimation. Now stated outright
in `avatar-algorithm.Rmd` Step 6 rather than left implicit. If this becomes a
problem in practice, the lever is dose-normalising PK-flagged DVIDs, which §6a
already scoped and deferred.

**Still open:** whether the floor `k` should exceed 5, unchanged by this work.
The cap is now the tighter of the two constraints, which weakens the argument
for raising `k`.

### D. Default anchor screen — DONE (2026-07-25), the "good enough" guard

The owner's guiding principle is *good enough, not perfect; output must not look
extreme; simplicity is valuable* (the canonical bad case is the lone 4580-hour
`wbcSim` avatar). The simplest expression of that is preventive and on by
default: `synpmx_avatar(screen = TRUE)` does not use a source subject whose
follow-up or dose count is more than twice the cohort's **90th percentile** as
an anchor, so no avatar inherits an extreme skeleton. Two iterations got here.
A MAD z-score (cutoff 3.5) over-excluded: on a tight core with a heavy tail
(wbcSim follow-up clusters near 480 h) the tiny MAD flags ordinary high-end
subjects, cutting the synthetic max to 576 and excluding 14 of 45. A median
multiple was better but still catches ordinary spread (2× the median can be
normal). Anchoring the cut on 2× the 90th percentile is the rule that only fires
well beyond the high end of normal: on wbcSim it excludes just the 1730/4580 h
subjects and keeps the ordinary tail (~1130 h and below).

Only follow-up length and dose count are screened — the axes that make a
skeleton look wrong. Dose magnitude is left alone (weight-based dosing makes
ordinary subjects look like dose outliers, and screening it broke `theo_md`),
and DV is blended rather than copied. Screening uses no randomness, so a source
with no extreme subject yields byte-identical output to `screen = FALSE`; only
datasets that actually contain an extreme structure change. This makes
not-extreme the default with no extra step, and leaves
`flag_identifiable_subjects()` / `remediate_identifiable_subjects()` as the
fuller, tunable manual layer.

Given this, skeleton sampling (B) is firmly a "consider, not do": the default
screen already removes the extreme structure simply, which is the bar the owner
set.

### C. A post-generation outlier detector (owner request) — DONE (2026-07-25)

`flag_identifiable_subjects(data, roles)` (`R/compare.R`) — the per-subject
counterpart to `compare_pmx_distributions()`, which checks distributions, not
individuals. Per the owner's steer, it screens the four **structural** axes that
make a subject easy to single out, each with a robust median/MAD modified
z-score (Iglewicz–Hoaglin, cutoff 3.5): **follow-up time** (last observation
time — the lone long-followed `wbcSim` subject), **number of doses**, **dose
magnitude**, and **DV value** (peak). A subject is flagged as an outlier on any
axis, with the offending axes named, riskiest first. This is the right tool for
the currently-live `REV-026` risk: because each avatar copies one anchor's event
skeleton, a structurally unique source subject yields a structurally unique
avatar even after the measurement blending. Run it on the synthetic output (or
the source) and drop or regenerate the flagged subjects. Tests in
`test-flag-identifiable.R`.

A near-copy-of-one-real-subject distance metric was prototyped and set aside:
the borrowing fix (A) now prevents near-copies by construction, so the live
residual risk is structural uniqueness, which the axis screen targets directly.

`remediate_identifiable_subjects()` (`R/compare.R`) acts on the flags with the
owner's policy (2026-07-25): a subject flagged **only** for a long follow-up is
*truncated* to the cohort's longest ordinary follow-up (the one axis a
value-level edit can fix), and a subject flagged for any other reason is
*dropped*, since an extreme-DV subject is elevated throughout and a rare dose
cannot be trimmed without breaking the regimen. An unusually *short* follow-up
is also dropped (nothing to truncate). When `source` is supplied, dropped
subjects are regenerated and refilled so the cohort keeps its size. `time` and
`other` options expose the policy. This is a stop-gap; skeleton sampling (B) is
the cure.

**Limitation found on `wbcSim` (2026-07-25).** Outlier screening is *relative*:
on a genuinely heavy-tailed axis (wbcSim follow-up runs from a few hours to
4580) the modified-z flags a fraction of subjects, removing the extreme tail
re-centres the distribution, and one pass does not fully converge. Remediation
still removes the individually extreme subjects (the 4580-hour avatar), which is
the identifiability that matters, but it cannot turn a heterogeneous cohort into
a homogeneous one, nor should it. This is direct evidence that skeleton sampling
(B) -- drawing follow-up length from the cohort distribution -- is the cleaner
fix for structural heterogeneity, where detect-and-remediate is a blunt tool.

### What was discovered on first implementation

The problem is far more pervasive than "sparse dose arms". Because the signature
carries the dose to 8 significant figures, any **weight-based (mg/kg) dosing**
gives nearly every subject a *distinct* dose and therefore its own singleton
group: `theo_md` has 11 distinct doses across 12 subjects. So the old code was
emitting a near-verbatim noised copy of essentially the **entire cohort** of any
individualized-dosing study, not just a few extreme-arm subjects. This makes the
borrowing fix load-bearing for the common case, not an edge-case guard.

### Phasing and status

1. **Cross-dose borrowing + loud alert — DONE (2026-07-25).** `.select_donors()`
   prefers same-signature donors, then borrows the nearest subjects from other
   groups to reach `k`; only a source smaller than `k + 1` cannot reach the
   floor, and that raises `.loud_warn()` — a red immediate `message()` (survives
   `suppressWarnings`) plus a warning condition. The owner chose the loud
   warning over a hard `stop()` for now. Tests in `test-avatar-pooling.R`.
2. **Post-generation outlier detector (C) — DONE (2026-07-25).**
   `flag_identifiable_subjects()` screens follow-up time, dose count, dose
   magnitude, and DV with a robust modified z-score.
3. **Event-skeleton sampling (B) — TODO.** Decouple schedule from value donors.
   The outlier detector (C) is the interim guard for the structural-uniqueness
   risk this would remove.
   Note the time-borrowing half already exists: `.interpolate_trajectory()`
   remaps a donor's trajectory onto out-of-range target times by proportional
   position, so a donor observed only early still contributes at later times.

### Still open

- Whether the *headline* floor should stay 5 given the paper leans higher — 5 is
  the feasible minimum, not a privacy recommendation; the owner may want to
  revisit once cross-dose pooling makes larger floors attainable.
- Tightening `max_donor_weight` below 0.80.

---

## 7. Reading guide

- `design/PROTOTYPE_SPEC.md` — the specification, with Version 4 (AVATAR) at the
  top and the DP versions retained below as alternatives.
- `vignettes/articles/feasibility.Rmd` — the measurements behind "DP utility collapses at
  small N" and the Version 1 re-identification analysis.
- `vignettes/articles/privacy-background.Rmd` — how the DP arithmetic (`d`, `f`, epsilon)
  works, for when the DP engine is the right choice.
