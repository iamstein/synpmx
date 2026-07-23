# Describing the trial: a data-structure elicitation guide

The companion to `design/MODEL_ELICITATION.md`. That guide asks **what the drug
does**; this one asks **what the dataset looks like**.

Both are needed, and they fail differently. A wrong PK prior gives you data with
the wrong numbers. A wrong design description gives you data with the wrong
*shape* — missing occasion columns, doses that never change when they should,
a `TIME` axis that resets when the real one does not. Against the accuracy bar
in `design/PROTOTYPE_SPEC.md` section 1, the second is the more serious failure:
structure must be exact, numbers need only be close.

Section 3 also settles which parts of a trial design are public. The short
answer: trial-level design is public and uncontroversial; an *individual's*
escalation or titration path is that individual's outcome and must never be
copied.

---

## 1. The complexity ladder

Find the lowest level that describes your trial. Each level adds structure the
generator must reproduce, and most add a privacy consideration.

### Level 0 — Single dose, single group, one endpoint

Every subject receives the same dose once. One PK endpoint.

- **Dataset signature:** one `EVID = 1` row per subject; `AMT` identical across
  subjects; no occasion column needed; `TIME` since dose.
- **Generator needs:** dose amount, sampling schedule, one endpoint.
- **Privacy:** everything about the design is in the protocol. Fully public.

### Level 1 — Dose escalation between cohorts (SAD)

Subjects are assigned to dose cohorts; each subject receives one dose level.

- **Dataset signature:** `AMT` varies between subjects, constant within subject.
  Usually a cohort or dose-group column.
- **Generator needs:** the dose levels, cohort sizes, and the dose-group
  assignment as a **subject property** (see `pmx_roles(subject_properties=)`).
  The structural model then gives dose-proportional exposure for free.
- **Privacy:** dose levels and cohort assignment are design facts and are
  treated as public. See section 3.

### Level 2 — Repeated dosing (MAD)

Multiple doses per subject, with accumulation, and sampling concentrated on
particular days.

- **Dataset signature:** multiple `EVID = 1` rows per subject, or `ADDL`/`II`.
  An occasion column. Rich profiles on day 1 and at steady state, troughs
  between.
- **Generator needs:** the dosing interval, number of doses, and **which
  occasions are sampled** — these are different questions, and conflating them
  was defect `SIM-003`. A Q24H regimen does not mean Q24H sampling.
- **Privacy:** all from the protocol. Public.

### Level 3 — Intra-patient dose escalation

A subject's dose changes over time, by protocol-defined rules.

- **Dataset signature:** `AMT` varies *within* subject across occasions. Often
  an assigned-dose column that changes at a defined visit.
- **Generator needs:** dose is no longer a subject property — it is
  occasion-varying. The `assigned_dose` role and occasion-conditioned regimen
  handle this, but the escalation *schedule* must be declared.
- **Privacy:** **an individual's escalation path is that individual's outcome**,
  not purely a design variable. Apply the rule to generated subjects; never
  replay a source subject's sequence. See section 3.

### Level 4 — Titration to effect, response-driven dosing

Dose is adjusted based on the individual's own measured response or toxicity.

- **Dataset signature:** like level 3, but dose changes correlate with the
  subject's own endpoint values.
- **Generator needs:** the titration rule as a public algorithm, applied to the
  *generated* response. Never replay a source dose sequence.
- **Privacy:** **the most sensitive structure in this ladder.** A subject's dose
  trajectory is a near-direct encoding of their response trajectory. Reproducing
  a real dose sequence discloses that subject's outcomes even if every DV is
  synthetic. Generate the sequence from the rule; never copy it.

### Level 5 — Crossover and sequence designs

Each subject receives multiple treatments in a randomized sequence.

- **Dataset signature:** period and sequence columns; `TIME` often resets within
  period or occasion; washout gaps; the same subject contributing several
  profiles.
- **Generator needs:** the sequence set and randomization ratio as public
  design; period/occasion clock handling. `mavoglurant` exercises the
  reset-clock case.
- **Privacy:** one person still contributes one record, so the adjacency holds —
  **provided the crossover is not pooled with a parallel study in which the same
  people also appear** (`REV-016`).

### Level 6 — Adaptive designs

Cohort sizes, dose levels, allocation ratios, or stopping are determined by
accumulating data.

- **Dataset signature:** irregular cohort sizes; dose levels that are not on the
  prespecified grid; early stopping; expansion cohorts.
- **Generator needs:** treat the *realized* design as an input, but see below.
- **Privacy:** the realized design is formally a function of the data, but is
  treated as public and is usually already disclosed. Recorded as an assumption
  in section 3, not as a risk to remediate.

---

## 2. Which levels are supported

| Level | Prior (no data) | Calibrated (v3 correction release) |
|---|---|---|
| 0-2 | Yes | Yes |
| 3 | Yes, with a declared escalation schedule | Yes, if the schedule is declared public |
| 4 | Yes, with a public titration rule | Possible, but the rule must be public and applied to generated responses |
| 5 | Yes | Yes |
| 6 | Yes | Yes, with the realized design recorded as a public input |

---

## 3. Which parts of the design are public

Throughout `design/PROTOTYPE_SPEC.md` the protocol is treated as a public input
costing no budget. That is right, with one narrow qualification that applies to
per-subject quantities rather than to trial-level ones.

### Trial-level realized design is public in practice

Realized dose levels, cohort sizes, and stopping points are **not** a meaningful
disclosure risk, for three independent reasons:

1. **They are usually already published.** Registry postings, conference
   presentations, and publications routinely disclose which dose levels were
   reached and how many subjects were enrolled, typically well before any
   dataset is shared. Differential privacy has no work to do protecting a fact
   that is already in the public domain.
2. **The inference is many-to-one and mostly confounded.** A trial stopping at a
   given dose could reflect a dose-limiting toxicity, but equally: target
   exposure or efficacy achieved, a prespecified maximum reached, enrollment
   difficulty, a data cutoff with later cohorts still ongoing, portfolio
   reprioritization, or a program halted for reasons unrelated to this study.
   Reading a stopping point as evidence about any individual is not sound.
3. **It is a property of the study, not of a person.** Escalation decisions are
   made by a committee weighing aggregate safety, exposure, and business
   context. Differential privacy protects the marginal contribution of one
   individual's record; a committee's decision about a program is not that.

Treat trial-level design as public, record its source in the provenance table,
and move on. Generating from the *planned* grid rather than the realized one is
a reasonable simplification if it is easier — synthetic data showing doses that were
never given is harmless for workflow testing — but it is a convenience, not a
privacy requirement.

### Per-subject escalation paths are individual outcomes

The qualification is here, and it is genuinely different in kind.

At levels 3 and 4, **an individual's dose trajectory is that individual's
outcome.** "This subject escalated at cycle 3, then reduced at cycle 5" is a
statement about one person's tolerability, with none of the many-to-one
confounding that makes the trial-level fact uninformative. In a titration design
the dose sequence is close to a direct encoding of the response sequence — it is
the variable the protocol adjusts *in order to* track the individual.

So the rule is not about the protocol. It is about what gets copied:

> **Apply the escalation or titration rule to generated subjects. Never replay a
> source subject's realized dose sequence.**

A dataset with entirely synthetic DV values but real per-subject dose
trajectories still discloses who tolerated treatment and who did not. This costs
nothing to get right — the rule is public, and applying it to generated
responses produces a more useful dataset anyway, since it stays coherent with
the generated exposures.

### Recorded as an assumption, not a defect

Strictly, realized design quantities are functions of the source dataset, so
treating them as fixed public inputs is an assumption rather than a derivation.
It belongs in the proof assumptions alongside the other public-input assertions,
where a reviewer can see it. It is not a leak to remediate. Tracked as
`REV-017`.

---

## 4. Dataset structure questions

Independent of design complexity, the generator needs the conventions your
datasets actually use. All of this is organizational convention, not patient
data, and costs nothing.

### Identifiers and grouping
1. **ID scheme.** Integer, character, site-subject composite?
2. **Cohort, group, or arm columns**, and their levels.
3. **Does any subject appear more than once** across the studies you intend to
   pool — rollover, extension, re-enrollment? This breaks the privacy adjacency
   (`REV-016`) and must be answered before pooling.

### Time
4. **Time units**, and whether `TIME` is continuous from first dose or resets.
5. **Does `TIME` reset within occasion or period?** Level 5 designs often do.
6. **Nominal time column?** Separate `NTIME`/`TAD`, or derived?
7. **Occasion column**, and exactly what an occasion means: a dose, a visit, a
   period, a cycle?

### Dosing
8. **`ADDL`/`II` or explicit dose rows?**
9. **Infusions:** explicit stop rows with negative `AMT`/`RATE`, or duration
   implied by `RATE`?
10. **Assigned-dose column** distinct from `AMT`? Does it vary by occasion
    (level 3+)?
11. **Dose interruptions, reductions, or holds** — how are they represented?

### Endpoints
12. **`DVID`/`CMT` coding** for each endpoint.
13. **`MDV` convention.**
14. **Censoring:** `CENS`/`LIMIT` Monolix-style, or a `BLQ` flag, or dropped?

### Covariates
15. **Baseline versus time-varying**, and which is which.
16. **How are missing covariates represented** — `NA`, `-99`, carried forward?

---

## 5. Output: a design declaration

The result of this interview, alongside the model declaration from
`design/MODEL_ELICITATION.md`:

```r
design <- pmx_trial_design(
  complexity      = "level_3_intrapatient_escalation",
  dose_levels     = c(10, 30, 100, 300),
  dose_levels_src = "protocol v4 section 5.1 (prespecified grid)",
  realized_levels = c(10, 30, 100),
  realized_src    = "public: ClinicalTrials.gov NCT01234567 results posting",
  escalation_rule = "escalate one level at cycle 3 if no grade 2+ event",
  cohort_sizes    = c(6, 6, 6),
  occasion_means  = "one cycle",
  time_resets     = FALSE
)
```

Every `_src` field records provenance. A realized quantity without a public
source must move to option 2 or option 3 of section 3.

### Provenance table addition

Extend the table in `design/MODEL_ELICITATION.md` with the design rows:

| Input | Value | Source | Data-independent? |
|---|---|---|---|
| Planned dose grid | 10/30/100/300 mg | Protocol v4 §5.1 | Yes |
| Realized levels reached | 10/30/100 mg | NCT01234567 results posting | Yes — publicly posted |
| Escalation rule | Cycle 3, no grade 2+ | Protocol v4 §5.3 | Yes |
| Realized per-subject escalation | — | **Not used**; generated from the rule | n/a |

---

## 6. Worked example: oncology dose escalation with intra-patient escalation

| Question | Answer |
|---|---|
| Level | **3.** Between-cohort escalation plus intra-patient escalation at cycle 3 |
| Planned grid | 10, 30, 100, 300 mg, 3+3 |
| Realized | Stopped at 100 mg; posted on ClinicalTrials.gov → public |
| Cohort sizes | 3, 3, 6 — publicly posted |
| Escalation rule | Escalate one level at cycle 3 if no grade 2+ event |
| Occasion | One cycle (21 days) |
| Time | Continuous from first dose; does not reset |
| Dosing | Explicit rows; `assigned_dose` column varies by occasion |
| Endpoints | PK (`cp`, dose-relative), platelet count (`plt`, study-time, IDR) |

**Generation approach.** Assign generated subjects to dose groups in the
realized cohort proportions. Apply the escalation rule to each *generated*
subject at cycle 3 — deterministically, or with a public probability, but never
copied from a source subject. Exposure follows the structural model, so a
subject escalating from 30 to 100 mg shows a step up in concentration
automatically.

**What is released privately.** Only the exposure correction factor and the
platelet-response correction. The entire design structure above is public.

**What must not happen.** Taking each real subject's escalation timing and
replaying it. Unlike the trial-level stopping point, an individual's escalation
timing is that individual's tolerability outcome, and copying it would disclose
which patients tolerated treatment — the primary endpoint of the study.
