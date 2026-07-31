# Score how many subjects share each subject's event skeleton

A source-side screen for subjects that are **alone**, the complement to
[`flag_identifiable_subjects()`](https://iamstein.github.io/synpmx/reference/flag_identifiable_subjects.md),
which finds subjects that are **extreme**. A patient can be perfectly
ordinary on every distribution and still hold the only copy of their
visit schedule, and
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
copies the anchor's event skeleton verbatim, so such a subject hands an
identifying schedule to every avatar anchored on it.

## Usage

``` r
skeleton_uniqueness(data, roles)
```

## Arguments

- data:

  A PMX dataset – normally the source.

- roles:

  Explicit roles from
  [`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md).

## Value

A `pmx_skeleton_uniqueness` data frame, most-exposed first, one row per
subject: `subject_id`, `n_obs`, `n_doses`, `signature_class`,
`obs_time_class`, `n_obs_class` (each counting the subjects sharing that
key, including itself), and `alone` (`TRUE` when `obs_time_class == 1`).
Attributes `n_unique_schedule`, `n_alone_signature`, `n_alone_n_obs`,
and `min_class` summarize the cohort.

## Details

Three equivalence classes are scored per subject:

- **`obs_time`** – subjects sharing the exact observation time vector.
  This is the fingerprint, because
  [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
  copies the anchor's event skeleton verbatim. Under nominal visit times
  the class is large, since the schedule is protocol-driven; under
  actual recorded times it is near-universally of size one.
  `coarsen_time = TRUE` collapses the second case into the first, and
  `alone` reports this class.

- **`n_obs`** – subjects sharing the observation count. Coarsening
  cannot change a count, so this is what survives it: dropout, early
  discontinuation, and missed visits. That residual is the outlier
  screen's job rather than the grid's.

- **`signature`** – subjects sharing the full
  [`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md)
  event signature: dose structure, dose amounts, and endpoint set. Note
  this does *not* include observation times – it is the key donor
  compatibility uses. Weight-based dosing or per-subject titration makes
  it unique regardless of schedule, and coarsening does not change that
  either.

Run it on the **source**, before generating, to decide whether
coarsening is needed and what is left over once it is applied. It is a
heuristic screen, not a privacy guarantee, and is marked
`"restricted_not_releasable"`.

## See also

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
#> Restricted PMX schedule-uniqueness screen: 0 of 30 patients
#> have a UNIQUE OBSERVATION SCHEDULE (0%): no other patient shares their
#> list of observation times, so the schedule works as an identifier.
#> 
#>   unique observation schedule:    0  <- of which:
#>     unique observation time:        0  <- sampled when nobody else was; the grid's job
#>     unique set of visits:      0  <- every time shared; dropout, the screen's job
#>   unique observation count:    0  <- the residual that leaves, for the screen
#>   unique dose signature:      30  <- dose structure/amount; coarsening cannot change it
#> 
#> Twelve most exposed:
#>  subject_id n_obs n_doses signature_class obs_time_class min_time_share
#>           1    14       2               1             30             30
#>          10    14       2               1             30             30
#>          11    14       2               1             30             30
#>          12    14       2               1             30             30
#>          13    14       2               1             30             30
#>          14    14       2               1             30             30
#>          15    14       2               1             30             30
#>          16    14       2               1             30             30
#>          17    14       2               1             30             30
#>          18    14       2               1             30             30
#>          19    14       2               1             30             30
#>           2    14       2               1             30             30
#>  n_obs_class unique_schedule
#>           30           FALSE
#>           30           FALSE
#>           30           FALSE
#>           30           FALSE
#>           30           FALSE
#>           30           FALSE
#>           30           FALSE
#>           30           FALSE
#>           30           FALSE
#>           30           FALSE
#>           30           FALSE
#>           30           FALSE
#> ... 18 more row(s) in the returned table.
#> 
#> Source-derived; not releasable unless separately public or privately budgeted.
```
