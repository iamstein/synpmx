# Empirical disclosure-risk evaluation (program)

**Status: proposed 2026-07-29 by the owner. Not started.**

This is a research program, not a queue item. `design/TODO.md` carries a pointer
and whichever slice is live; the reasoning, the phasing, and the honest costs
live here.

## 1. Why, and where it fits

`vignettes/articles/feasibility.Rmd` section 5 separates three things that get
called privacy:

1. **Formal guarantee** — a worst-case bound on one person's influence
   (differential privacy, DP). Provable, needs N.
2. **Empirical risk assessment** — measured attack success. No guarantee, but
   evidence.
3. **Governance** — trust boundary, contract, access control.

The package makes claim (1) for the DP engines and rests on (3) for AVATAR,
which is the default. There is no measured attack-success evidence for the
default mode. That is the gap this program fills.

The trust-boundary decision rule (`METHOD_DISCUSSION.md` section 5) is a
governance rule. It is defensible, and as stated it is also unfalsifiable —
nothing measures what an attacker could actually recover from a released
dataset. An empirical evaluation does not replace the rule; it says whether the
rule is doing any work.

## 2. What already exists — do not rebuild it

- **`compare_pmx_proximity()`** (`REV-029`). Nearest-neighbour adversarial
  accuracy between source and synthetic subjects in a shared principal
  component analysis (PCA) profile space. This is already a
  singling-out/linkability proxy, and it is attack number one of the framework
  rather than a precursor to it.
- **The split-half null.** Two halves of the *source* cohort through the
  identical statistic, so every small-sample artefact is present in null and
  observed alike and cancels. Every attack added here should get its null the
  same way. This is the part that makes a number at N = 30 mean anything.
- **The positive-control habit.** `test-proximity.R` leads with a verbatim copy
  that the metric must reject, on the reasoning that a privacy metric never seen
  to fire is an untested branch. The owner's "intentionally leaky datasets"
  point is already house rule for one metric; the work is to generalize it.
- **`skeleton_uniqueness()`.** Exact structural uniqueness — `n_share_rarest_time`,
  attendance patterns, dose counts. Not an attack, but the ground truth an
  attack can be checked against: where uniqueness is known to be total, an
  attack that scores at chance is a broken attack.
- **`flag_identifiable_subjects()`.** Four-axis outlier screen on the source.
- **`synpmx_prior()` is a genuine negative control.** It never reads the study,
  so its disclosure risk is zero by construction. Any attack scoring above
  chance against prior-mode output is measuring an artefact of the attack, not a
  leak. Most synthetic-data evaluations have no such control available.

## 3. The honest cost, stated before the plan

Pharmacometric cohorts are 12–60 subjects. `REV-029` already found what that
does: the split-half null is wide enough that only a blatant leak is detectable,
and the print method says so out loud. A membership-inference attack that splits
N = 30 into 15 training and 15 holdout subjects has almost no power. Two
consequences, and they are not caveats to append at the end:

- Every result is **"nothing detected"**, never "nothing there". That phrasing
  belongs in each print method and each conclusion, not only in the discussion.
- Power comes from replication across splits, seeds, and datasets rather than
  from any single run. Plan the simulation study around that, and expect the
  compute.

Decide whether the answer justifies the build before building all of it. Phases
1 and 2 are the cheap slice that answers exactly this.

## 4. Phases

### P1 — Literature, and one design decision that follows from it

Read and summarise: membership inference; singling out; linkability; attribute
inference; Anonymeter (the closest off-the-shelf framework, built around the
first three); and AVATAR's own metrics — hidden rate and local cloaking — from
Guillaudeux 2023 (`references/Guillaudeux23.pdf`, already owed a read under
`TODO.md`).

The decision this forces: **every framework named assumes one flat record per
subject.** Pharmacometric data is a variable-length event table. So either each
subject is first summarised into a fixed-length vector — which is what
`compare_pmx_proximity()` already does via the profile PCA — or the attack
definitions have to be restated for longitudinal event data. That choice is the
actual research contribution here, and it should be made deliberately rather
than inherited from whichever library is convenient.

Output: a section appended to this file. Not a new document.

### P2 — Train/holdout membership inference

Split subjects into training and holdout, generate from the training set only,
and ask whether an attacker can tell which subjects were used. Repeat across
random splits and generation seeds.

Reuse the profile space and the null machinery in `R/compare.R`. Deliverable is
one function plus a script over the public datasets. Nothing is believed until
both controls run: the positive controls (exact copies; lightly perturbed
copies) must be caught, and the negative control (`synpmx_prior()`) must score
at chance.

#### Added 2026-08-12, from writing the checking review

**The production run and the characterization run are different activities, and
only the second one splits.** A shipped dataset is generated from the whole
cohort: every patient should be represented and the donor pool is the binding
constraint. The split exists to produce evidence about the *algorithm*, once,
and its output is a statement about the generator rather than a dataset anyone
uses. Scoping P2 as a per-study step would be a mistake; it is a
characterization run, and it belongs on public data large enough for the
statistic to resolve something.

**The API is decided by a projection constraint, not by taste.**
`compare_pmx_proximity()` builds its PCA profile space from the *pair* it is
handed: `combined <- rbind(source, synthetic)` then `.build_profiles(combined,
...)`. Calling it twice — once with the training set, once with the control set
— therefore fits two different projections, and part of the resulting difference
would be a change of coordinate system rather than memorization. $T$, $C$ and
$S$ must be projected **together, in one call**. That rules out a naive public
`compute_adversarial_accuracy(a, b)` taking two arbitrary sets, because it hands
the trap to the caller. The safe shapes are:

- `compare_pmx_proximity(source, synthetic, roles, holdout = <ids or data>)`,
  projecting all three together, returning `aa_train`, `aa_control` and
  `privacy_loss`. Backward compatible; `holdout = NULL` behaves as today. The
  comparison size becomes `min(|T|, |C|, |S|)`.
- A `scripts/` driver for the repeated stratified split, per `AGENTS.md`:
  multi-seed and report-producing work lives there, sharing the one metric
  implementation above. Cost is `replicates` full generation runs.

**Measured noise floors, so section 3's warning has numbers.** Null interval
width from `compare_pmx_proximity()`, and the memorized fraction needed to move
$\mathrm{AA}$ clear of it, using $\delta \approx p/2$:

| Dataset | Subjects | Compared | Null width | Fraction needed |
|---|---|---|---|---|
| `theo_md` | 12 | 6 | 0.648 | 130% — nothing is detectable |
| `warfarin` | 32 | 16 | 0.385 | 77% |
| `case1_pkpd` | 180 | 90 | 0.141 | 28% |

The floor falls as $1/\sqrt{n}$ while one patient's contribution falls as
$1/n$, so the gap widens with cohort size: the procedure detects **systematic**
memorization at every size tested and never approaches single-patient
resolution. That is the argument for P3 and for per-record measures, and it is
also the honest answer to "should we run this per study" — no.

**A smaller defect found on the way.** The null is averaged over `replicates`
= 50, but the observed statistic is a **single** subsample:
`observed <- .adversarial_accuracy(take(fake, size), take(real, size))`. At
`case1_pkpd` sizes this is immaterial; at 12 subjects the observed value is one
draw of 6-vs-6 compared against an interval smoothed 50 times. Averaging the
observed over the same replicates, or reporting its spread, costs nothing and
makes the comparison symmetric. Worth doing before P2 rather than after, since
P2 differences inherit it.

### P3 — Which part of a subject leaks

The same attack restricted to one surface at a time: covariates only; dosing and
event history only; observations only; the combined fingerprint.

The hypothesis has evidence behind it already — the **event history should
dominate**. `skeleton_uniqueness()` reported nearly every subject unique on dose
before `REV-027`; `SIM-033` found every avatar reproducing one real subject's
exact observation-time vector; `REV-027` found a copied `AMT` disclosing a real
subject's weight exactly under milligram-per-kilogram dosing. If the
decomposition confirms it, that is the finding worth publishing, and it is one
the pharmacometrics literature has not stated.

### P4 — Compare what is being compared

Two axes, and the second is the more useful one:

- **Across modes:** AVATAR, `synpmx_prior()` (negative control), and the DP path
  (`synpmx_empirical()`) at a few epsilons.
- **Across mechanisms — the ablation.** AVATAR with `coarsen_time`,
  `min_pattern_share`, the dose-basis recomputation, the anchor screen, and
  `max_donor_weight` disabled one at a time. Each of those costs utility today
  and none has been shown to buy measurable privacy. The ablation is what values
  them.

### P5 — Privacy against utility

Pair each privacy number with the utility measures that exist:
`compare_pmx_distributions()`, and the covariate-slope and treatment-effect
evidence in `test-avatar-relationships.R` and
`vignettes/articles/example-avatar-PKPD-covariate-treatment-effect.Rmd`. Only
worth doing once P2–P4 produce a privacy number that actually moves.

## 5. How this is tracked

**No new identifier namespace.** Phases live in this file. When a phase produces
buildable work it becomes a `REV-###` in `REVIEW_BACKLOG.md`; when it produces a
regression gate it becomes a `SIM-###` in `TEST_SIM.md`. `design/TODO.md` holds
the pointer and the live slice only. Two registries are enough.

## 6. Open questions for the owner

- **What is this for?** Internal defensibility, the ACoP poster, or a paper. P1
  and P2 alone cover defensibility; the full program is a paper, and it competes
  directly with the poster deadline and the real-data validation currently
  pacing `TODO.md`.
- **"Differentially private synthesis, if available."** It is available —
  `synpmx_empirical()` — but gated behind `synpmx_enable_dp_engines()` as
  unaudited. Usable as a comparator, with that caveat attached to any number it
  produces.
