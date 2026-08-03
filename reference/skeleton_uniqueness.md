# Score how many patients share each patient's event skeleton

A source-side screen for patients that are **alone**, the complement to
[`flag_identifiable_subjects()`](https://iamstein.github.io/synpmx/reference/flag_identifiable_subjects.md),
which finds patients that are **extreme**. A patient can be perfectly
ordinary on every distribution and still hold the only copy of their
visit schedule, and
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
copies the anchor's event skeleton verbatim, so such a patient hands an
identifying schedule to every avatar anchored on it.

## Usage

``` r
skeleton_uniqueness(data, roles, coarsen_time = FALSE)
```

## Arguments

- data:

  A PMX dataset – normally the source.

- roles:

  Explicit roles from
  [`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md).

- coarsen_time:

  Score the coarsened visit grid
  [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
  would build (`TRUE`) or the recorded times as given (`FALSE`, the
  default).

## Value

A `pmx_skeleton_uniqueness` data frame, most-exposed first, one row per
patient: `subject_id`, `n_obs`, `n_doses`, `n_share_dosing`,
`n_share_schedule`, `n_share_rarest_time`, `n_share_obs_count` (each
counting the patients sharing that property, including this one),
`unique_schedule` (`TRUE` when `n_share_schedule == 1`), and
`why_unique`. Attributes `n_unique_schedule`, `n_unique_dose_signature`,
`n_unique_obs_count`, `n_unshared_time`, `min_class`, and `coarsened`
summarize the cohort; `summary_table` and `sharing_table` hold the two
tables [`print()`](https://rdrr.io/r/base/print.html) shows.

## Details

Three questions are asked of every patient, and the answer to each is a
count of how many patients share that property, the patient included. A
count of 1 means "nobody else":

- **`n_share_schedule`** – who else was observed at exactly this list of
  times? This is the fingerprint, because
  [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
  copies the anchor's event skeleton verbatim. Under nominal visit times
  the count is large, since the schedule is protocol-driven; under
  actual recorded times it is near-universally 1. Coarsening collapses
  the second case into the first.

- **`n_share_obs_count`** – who else has this many observations?
  Coarsening cannot change a count, so this is what survives it:
  dropout, early discontinuation, and missed visits.

- **`n_share_dosing`** – who else has this dose structure and these dose
  amounts? This is the full
  [`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md)
  event signature and it does *not* include observation times; it is the
  key donor compatibility uses. Weight-based dosing or per-patient
  titration makes it unique regardless of schedule, and coarsening does
  not change that either.

`n_share_rarest_time` splits the schedule count by cause, which matters
because the two causes have opposite remedies. A patient whose schedule
is unique **and** whose rarest single time was shared with nobody
(`n_share_rarest_time == 1`) was sampled at a one-off moment: a time
grid is meant to absorb that, and declaring `nominal_time` is the fix. A
patient whose schedule is unique while every individual time is shared
(`n_share_rarest_time >= 2`) is a dropout or missed-visit pattern, and
no grid at any resolution touches it. `why_unique` states which.

## Before or after coarsening

By default this scores the times **exactly as they appear in `data`**.
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
snaps the source onto a shared visit grid first (`coarsen_time = TRUE`,
its default) and the numbers it records in `pmx_settings` are therefore
post-coarsening. Pass `coarsen_time = TRUE` here to score the same grid
the generator would build, and run it both ways to see how much of the
exposure coarsening actually removed. The printed header always says
which of the two you are looking at.

Run it on the **source**, before generating. It is a heuristic screen,
not a privacy guarantee, and is marked `"restricted_not_releasable"`.

## See also

[`plot_pmx_schedule()`](https://iamstein.github.io/synpmx/reference/plot_pmx_schedule.md)
for the same information as a picture,
[`flag_identifiable_subjects()`](https://iamstein.github.io/synpmx/reference/flag_identifiable_subjects.md),
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md).

## Examples

``` r
data <- pmx_simulated_fixture(30)
roles <- pmx_roles(
  id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
  cmt = "CMT", dvid = "DVID", covariates = "WT"
)
skeleton_uniqueness(data, roles)
#> Restricted PMX schedule-uniqueness screen
#> Scored on the recorded times AS GIVEN, before any coarsening.
#> `synpmx_avatar()` coarsens first by default, so run this again with
#> `coarsen_time = TRUE` to see what the grid removes.
#> 
#> Every patient shares their observation schedule with somebody. Nothing to
#> do.
#> 
#>                    Patients whose ...  n % of cohort
#>  Observation schedule nobody else has  0           0
#>        ... a one-off observation time  0           0
#>        ... the set of visits attended  0           0
#>     Observation count nobody else has  0           0
#>                Dosing nobody else has 30         100
#> 
#> How crowded is each schedule (1 = nobody else has it):
#>  Patients sharing that schedule Patients % of cohort
#>                              30       30         100
#> 
#> One row per patient is in the returned data frame; `plot_pmx_schedule()`
#> draws the same cohort. Source-derived; not releasable unless separately
#> public or privately budgeted.
skeleton_uniqueness(data, roles, coarsen_time = TRUE)
#> Restricted PMX schedule-uniqueness screen
#> Scored AFTER coarsening, on the shared visit grid `synpmx_avatar()` builds.
#> These are the numbers a run reports.
#> 
#> Every patient shares their observation schedule with somebody. Nothing to
#> do.
#> 
#>                    Patients whose ...  n % of cohort
#>  Observation schedule nobody else has  0           0
#>        ... a one-off observation time  0           0
#>        ... the set of visits attended  0           0
#>     Observation count nobody else has  0           0
#>                Dosing nobody else has 30         100
#> 
#> How crowded is each schedule (1 = nobody else has it):
#>  Patients sharing that schedule Patients % of cohort
#>                              30       30         100
#> 
#> One row per patient is in the returned data frame; `plot_pmx_schedule()`
#> draws the same cohort. Source-derived; not releasable unless separately
#> public or privately budgeted.
```
