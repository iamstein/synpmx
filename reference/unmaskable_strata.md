# Which strata can mask their own avatars

One row per declared stratum level, counting the patients in it that no
avatar can safely be built on. An avatar is anchored on one real patient
and copies that patient's dose times verbatim, so an arm has to contain
somebody who can be masked; re-anchoring picks a different anchor, but
never one from another arm, because the anchor carries its `strata`
values into the output and moving it would silently rewrite the arm
sizes.

## Usage

``` r
unmaskable_strata(data, roles, min_pattern_share = 2L, coarsen_time = TRUE)
```

## Arguments

- data:

  A PMX dataset – the source, before generating.

- roles:

  Explicit roles from
  [`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md).

- min_pattern_share:

  The floor
  [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
  will be run with: how many patients must share a pattern before it may
  be reused.

- coarsen_time:

  Score the coarsened visit grid (`TRUE`, the default, and what
  [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
  does) or the recorded times as given.

## Value

A `pmx_unmaskable_strata` data frame, worst arm first, with `stratum`,
`patients`, `unmaskable_dosing`, `unmaskable_visits` and `safe_anchors`.
One row named `all` when the roles declare no strata.

## Details

`safe_anchors` is therefore the number that matters. An arm with none of
them will emit avatars carrying a real patient's schedule whatever the
generator does, which is what
[`synpmx_scorecard()`](https://iamstein.github.io/synpmx/reference/synpmx_scorecard.md)
reports as B1a and B1b, and the fix is to the study description rather
than to the run – see the alert
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
raises.

Two ways a patient cannot be masked, counted separately because the
remedies differ:

- `unmaskable_dosing` – their set of dose times is shared by fewer than
  `min_pattern_share` patients, and no prefix of it can be given to an
  avatar either: every shorter opening is one that either a single
  patient stopped at, or that fewer than `min_pattern_share` patients
  passed through. Dose events are copied from the anchor as they stand,
  since resampling them would emit regimens the protocol never
  permitted, so there is nothing to substitute.

- `unmaskable_visits` – the set of visits they have observations at is
  shared by fewer than `min_pattern_share` patients, and their schedule
  group holds no set that is shared widely enough to put in its place.
  Usually fixable: a wider pool, or `strata = NULL` so the arms share
  one.

This reads real patient data and is marked
`"restricted_not_releasable"`.

## See also

[`skeleton_uniqueness()`](https://iamstein.github.io/synpmx/reference/skeleton_uniqueness.md)
for the same question per patient,
[`synpmx_scorecard()`](https://iamstein.github.io/synpmx/reference/synpmx_scorecard.md),
which reports this as rows B1a and B1b.

## Examples

``` r
data <- pmx_simulated_fixture(20)
data$ARM <- ifelse(as.integer(data$ID) %% 2L == 0L, "A", "B")
roles <- pmx_roles(
  id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
  cmt = "CMT", dvid = "DVID", covariates = "WT", strata = "ARM"
)
unmaskable_strata(data, roles)
#> Patients no avatar can safely be built on, by arm
#> 
#>  stratum patients unmaskable_dosing unmaskable_visits safe_anchors
#>        B       10                 0                 0           10
#>        A       10                 0                 0           10
#> 
#> Every arm can mask its own avatars.
```
