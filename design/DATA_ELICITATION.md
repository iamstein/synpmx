# Describing the trial: a data-structure elicitation guide

The companion to `design/MODEL_ELICITATION.md`. That guide asks **what the drug
does**; this one asks **what the dataset looks like**.

Both are needed, and they fail differently. A wrong PK prior gives you data with
the wrong numbers. A wrong design description gives you data with the wrong
*shape* — missing occasion columns, doses that never change when they should,
a `TIME` axis that resets when the real one does not. Against the accuracy bar
in `design/PROTOTYPE_SPEC.md` section 1, the second is the more serious failure:
structure must be exact, numbers need only be close.

This guide also carries a privacy finding that does not appear anywhere else:
**in escalation and adaptive designs, part of the "protocol" is not public.**
See section 3.

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
- **Privacy:** the *planned* dose levels are public. **Which levels were
  actually reached is not** — see section 3.

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
- **Privacy:** **whether a given patient escalated depends on that patient's own
  tolerability.** The escalation path is an outcome, not purely a design
  variable. See section 3.

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
- **Privacy:** **the realized design is a function of the data.** This is the
  deepest issue on the ladder and section 3 is about it.

---

## 2. Which levels are supported

| Level | Prior (no data) | Calibrated (v3 correction release) |
|---|---|---|
| 0-2 | Yes | Yes |
| 3 | Yes, with a declared escalation schedule | Yes, if the schedule is declared public |
| 4 | Yes, with a public titration rule | Possible, but the rule must be public and applied to generated responses |
| 5 | Yes | Yes |
| 6 | Yes, if the realized design is asserted public | Requires an explicit decision — section 3 |

---

## 3. The privacy finding: not all of the protocol is public

Throughout `design/PROTOTYPE_SPEC.md` the protocol is treated as a public input
costing no budget. **For levels 0, 2, and 5 that is straightforwardly true.**
For levels 1, 3, 4, and 6 it is only partly true, and the distinction is easy to
miss because it all arrives in the same document.

The rule:

> **Prespecified design is public. Realized design that depended on the data is
> not.**

Concretely:

| Quantity | Status |
|---|---|
| The planned dose-escalation grid (10, 30, 100, 300 mg) | Public. In the protocol before enrollment |
| **Which of those levels were actually reached** | **Data-dependent.** Stopping at 100 mg implies a dose-limiting toxicity at 300 |
| Planned cohort size | Public |
| **Realized cohort size after a 3+3 expansion** | **Data-dependent.** Six subjects instead of three means a DLT occurred |
| The titration rule ("escalate if no grade 2 event") | Public. It is an algorithm |
| **Which subjects escalated, and when** | **Data-dependent.** It is their outcome |
| Planned visit schedule | Public |
| **Realized visit times, dropout, discontinuation** | **Data-dependent** |

This matters because escalation stopping points are informative about small
numbers of people. "The trial stopped at 100 mg" can be equivalent to "at least
one of the six subjects at 300 mg had a serious event" — a statement about a
handful of individuals, released with no accounting at all.

### What to do about it

Three acceptable options, in preference order:

1. **Declare the realized design public by assertion.** Often legitimate:
   escalation outcomes are frequently disclosed at conferences, in press
   releases, or on ClinicalTrials.gov before any dataset is shared. If the
   realized dose levels are already in the public domain, say so, record the
   source in the provenance table, and proceed. **This is usually the right
   answer and it costs nothing.**
2. **Generate from the planned design instead of the realized one.** Use the
   full prespecified escalation grid and planned cohort sizes. The mock data
   will show doses that were never given, which is harmless for workflow
   testing and discloses nothing.
3. **Budget it.** Treat realized cohort sizes as private counts and release them
   through the mechanism. Correct, but it consumes budget that would be better
   spent on exposure magnitude, and for a handful of cohorts the noise will
   dominate anyway.

Option 2 is underrated. The generator does not need the trial's actual
escalation history to produce a dataset that exercises an escalation-aware
analysis pipeline.

### What is never acceptable

Replaying a real subject's dose-escalation or titration sequence. At level 4 in
particular, the dose trajectory encodes the response trajectory. A dataset with
fully synthetic DV values but real dose sequences can still disclose who
responded and who did not.

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
| Realized levels reached | 10/30/100 mg | NCT01234567 results posting | Yes — publicly disclosed |
| Escalation rule | Cycle 3, no grade 2+ | Protocol v4 §5.3 | Yes |
| Realized per-subject escalation | — | **Not used**; generated from the rule | n/a |

---

## 6. Worked example: oncology dose escalation with intra-patient escalation

| Question | Answer |
|---|---|
| Level | **3.** Between-cohort escalation plus intra-patient escalation at cycle 3 |
| Planned grid | 10, 30, 100, 300 mg, 3+3 |
| Realized | Stopped at 100 mg; disclosed on ClinicalTrials.gov → assert public |
| Cohort sizes | 3, 3, 6 — the 6 implies a DLT, but is publicly posted |
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
replaying it. That would disclose which patients tolerated treatment, which is
the primary endpoint of the study.
