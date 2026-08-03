# Report what each masking mechanism did, and what it cost

[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
records everything it removed on the `"pmx_settings"` attribute of its
result. This turns that flat list into the table to read after a run:
who was left to build on, how many real patients reach one avatar, what
the visit grid managed to collapse, which visit sets were too rare to
reuse, and whether dose amounts were recomputed.

## Usage

``` r
pmx_masking_report(synthetic, source = NULL, roles = NULL)
```

## Arguments

- synthetic:

  A dataset from
  [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md),
  carrying its `"pmx_settings"` attribute.

- source:

  Optionally the source dataset. Supplying it (with `roles`) adds the
  before-coarsening schedule count, so the table shows what coarsening
  removed rather than only what was left.

- roles:

  Explicit roles from
  [`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md).
  Required with `source`.

## Value

A `pmx_masking_report` data frame with columns `Quantity`, `Value`, and
`What it means`. Section headers appear as rows whose `Quantity` is bold
and whose other cells are empty.

## Details

Every row carries a sentence saying what the number means, because none
of them mean anything on their own. The rows worth looking at hardest:

- **Unique observation schedules, after coarsening** – patients whose
  list of observation times nobody else shares. An avatar anchored on
  one wears a schedule belonging to one real person. Its two sub-rows
  have opposite remedies: a one-off observation time is what declaring
  `nominal_time` fixes, and a unique set of *attended* visits is
  dropout, which no grid touches.

- **Shared by too few patients, so not reused** – real dropout and
  dose-interruption patterns that will not appear in the synthetic data.
  Discarding them is what stops an avatar carrying a schedule traceable
  to one person. If this study's interruptions matter, lower
  `min_pattern_share` (2 is the lowest value that still guarantees no
  synthetic patient has a schedule unique to a real one).

- **Avatars keeping their anchor's own visit set** – the fallback when a
  schedule group had nothing shareable to draw from. Those avatars carry
  one real patient's pattern of absences, which is the thing
  `min_pattern_share` exists to prevent, so a high percentage here
  undoes the row above it.

- **Amounts recomputed from a covariate** – says outright whether
  weight-based or body-surface-area dosing was detected, and when it was
  not, why not. Detection is deliberately conservative: it fails closed
  and leaves amounts alone rather than rewriting a study that is not
  dose-proportional.

At the default `min_pattern_share = 2`, "shared by too few patients" and
"real patients holding those" are necessarily equal – a set is discarded
exactly when fewer than two patients share it, so every discarded set
has one holder. They diverge only at a floor of 3 or more.

Marked `"restricted_not_releasable"` when `source` is supplied, since
the before-coarsening row then reads the source.

## See also

[`skeleton_uniqueness()`](https://iamstein.github.io/synpmx/reference/skeleton_uniqueness.md),
[`compare_pmx_proximity()`](https://iamstein.github.io/synpmx/reference/compare_pmx_proximity.md),
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md).

## Examples

``` r
data <- pmx_simulated_fixture(30)
roles <- pmx_roles(
  id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
  cmt = "CMT", dvid = "DVID", covariates = "WT"
)
synthetic <- suppressWarnings(synpmx_avatar(data, roles, seed = 1))
#> synpmx_avatar(): dropped 9 undeclared column(s): NTIME, TAD, OCC, RATE, MDV, CENS, LIMIT, AGE, SEX.
#>   Declare a column in `keep` to carry it through verbatim.
pmx_masking_report(synthetic, data, roles)
#> What the masking mechanisms did
#> 
#> Who was available to build on
#>   Patients in the source                           30
#>     excluded as structurally extreme               0
#>       `screen`: follow-up or dose count over twice the cohort's 90th
#>       percentile
#>     excluded, route arm too small                  0
#>       `on_donor_shortfall`: a route arm holding fewer than k + 1 patients
#>     left to anchor avatars on                      30
#>       an excluded patient still contributes as a donor
#>   Avatars built                                    30
#>       cohort size is unaffected by the exclusions above
#> 
#> Donor pools: who may be blended with whom
#>   Administration routes                            1
#>       oral, infusion, and so on. Donors are NEVER blended across a route,
#>       so each is a separate pool
#>   Dose/schedule groups                             30
#>       patients with an identical dose pattern and endpoint set. Donors are
#>       looked for here first; many small groups means the search falls back
#>       to the wider route pool
#> 
#> How much of one real patient reaches one avatar
#>   Donor floor, k                                   5
#>       real patients blended into each avatar
#>   Largest share one donor may hold                 0.5
#>       `max_donor_weight`
#>     that cap actually bound on                     70%
#>       of avatars. Near 100% means the cap, not distance, is setting the
#>       weights
#>   Effective donors per avatar, mean                2.86
#>       1 / sum(w^2). This, not k, is how many patients an avatar is really
#>       made of
#> 
#> Visit schedule: WHEN patients were observed
#>   Visit grid used                                  derived
#>       no usable `nominal_time`, so a grid was inferred from the recorded
#>       times themselves. Declaring `nominal_time` is better
#>   Unique observation schedules, before coarsening  0
#>       patients whose list of observation times nobody else shares
#>   Unique observation schedules, after coarsening   0
#>       the count that matters: an avatar copies its anchor's times verbatim
#>     because of a one-off observation time          0
#>       sampled when nobody else was. Declaring `nominal_time` is the fix
#>     because of which visits they attended          0
#>       every time is shared; this is dropout, and no grid can fix it
#> 
#> Visit sets: WHICH of those visits each patient attended
#>   Distinct visit sets in the source                1
#>       a visit set is which of the shared grid visits one patient actually
#>       had
#>     held by fewer than 2 patients, so not reused   0
#>       `min_pattern_share` is that threshold. These visit sets are lost, not
#>       approximated
#>     real patients holding those                    0
#>       those patients are NOT removed -- they still anchor avatars and still
#>       act as donors. Only their particular pattern of absences stops being
#>       copied
#>   Avatars given a visit set from the pool          100%
#>       drawn from the sets that cleared the threshold, or built from their
#>       shape -- never from their own anchor alone
#>     of those, misses placed fresh                  0%
#>       the kind of missingness was reused; exactly which visits were missed
#>       was invented
#>     of those, miss count moved                     0%
#>       no placement at the wanted number of misses was free, so the count
#>       moved by a visit or two. Dropout is the usual reason: a
#>       discontinuation at a given depth has only one possible placement
#>   Avatars keeping their anchor's own visit set     0%
#>       the last-resort fallback, and the one row you want at 0%: these
#>       avatars carry one real patient's absences exactly. The run alerts
#>       past 10%
#> 
#> Dose
#>   Amounts recomputed from a covariate              yes, from `WT` (inferred)
#>       the 30 distinct dose amounts are a fixed multiple of `WT`, at 9
#>       protocol level(s)
#>     protocol levels found                          1.17, 1.205, 1.281, 1.376, 1.43, 1.492, 1.622, 1.741, 1.788
#>       dose per unit of `WT`; every amount was snapped to the nearest of
#>       these
```
