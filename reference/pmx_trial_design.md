# Declare a public trial design

Every field is a design fact from the protocol and consumes no privacy
budget. See the data elicitation guide at
<https://iamstein.github.io/synpmx/articles/data-elicitation.html> for
which parts of a realized design need a provenance note.

## Usage

``` r
pmx_trial_design(
  dose_levels = NULL,
  cohort_sizes = NULL,
  sampling,
  n_doses = 1L,
  dose_interval = 24,
  dose_escalation = NULL,
  dose_times = NULL,
  duration = 0,
  visit_window = 0.05,
  source
)
```

## Arguments

- dose_levels:

  Dose amounts, one per cohort. Omit when using `dose_escalation`.

- cohort_sizes:

  Planned subjects per cohort, recycled over `dose_levels`. Defaults to
  equal cohorts.

- sampling:

  Nominal sampling times after each dose, from the protocol.

- n_doses:

  Number of doses per subject in a parallel design.

- dose_interval:

  Time between doses when doses are equally spaced.

- dose_escalation:

  Per-occasion dose amounts for a within-subject escalation, for example
  `c(10, 30, 100)`. Applied to every subject. For a trial that also
  escalates between cohorts, supply a list of sequences, one per cohort,
  for example `list(c(1, 2, 4), c(2, 4, 8), c(4, 8, 16))`, with
  `cohort_sizes` giving the subjects in each. Every sequence must have
  the same length, since the dosing schedule is shared.

- dose_times:

  Explicit dose times, for example `c(0, 7, 14)`. Defaults to equally
  spaced times at `dose_interval`.

- duration:

  Infusion duration; zero for bolus or oral.

- visit_window:

  Fractional jitter applied to nominal times.

- source:

  Required provenance string.

## Value

A `pmx_trial_design`.

## Details

Two dosing patterns are supported. A parallel design gives each cohort
one dose level (`dose_levels`) repeated `n_doses` times. A
within-subject escalation gives every subject the same increasing
sequence of doses (`dose_escalation`), one per occasion; this is
prespecified design, not an outcome, when the escalation follows a fixed
protocol schedule.
