# Flag structurally unusual – and so easily identifiable – subjects

A post-generation screen for subjects that stand out from the cohort and
are therefore easy to single out and re-identify: the per-subject
counterpart to
[`compare_pmx_distributions()`](https://iamstein.github.io/synpmx/reference/compare_pmx_distributions.md),
which compares whole distributions. Each subject is scored, one axis at
a time, on a robust median/MAD statistic across four structural
features:

## Usage

``` r
flag_identifiable_subjects(data, roles, threshold = 3.5)
```

## Arguments

- data:

  A PMX dataset – typically the synthetic output, or the source.

- roles:

  Explicit roles from
  [`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md).

- threshold:

  Absolute modified-z cutoff above which a subject is an outlier on an
  axis. Default 3.5, the Iglewicz–Hoaglin value.

## Value

A `pmx_identifiability` data frame, most-unusual first, one row per
subject: `subject_id`, the four axis values (`follow_up_time`,
`n_doses`, `max_dose`, `max_dv`), `outlier_axes` (a comma-separated list
of the axes on which it is unusual, empty if none), and `flagged`.

## Details

- **follow-up time** – the last observation time (catches the lone
  long-followed subject);

- **number of doses** – an unusual dosing-history length;

- **dose magnitude** – a rare dose level (needs an `amt` role); and

- **DV value** – an extreme peak measurement.

A subject is flagged when it is an outlier on any axis. This matters
because
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
copies each avatar's event skeleton from a single anchor, so a
structurally unique source subject yields a structurally unique – and
identifiable – avatar even though its measurements are blended. Run it
on the synthetic data before the data leaves the source's access
controls and drop or regenerate the flagged subjects; it can also be run
on the source itself to see which real subjects are hardest to hide. It
is a heuristic screen, not a privacy guarantee, and is marked
`"restricted_not_releasable"`.

Scores are computed **within each declared stratum**
([`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md)
`strata`), because "does this patient stand out?" needs a comparison
group and the whole cohort is the wrong one as soon as a study assigns
anything. On a six-arm dose-ranging study the top arm sits far from the
cohort median dose purely by protocol: scored cohort-wide,
[`xgxr::case1_pkpd`](https://rdrr.io/pkg/xgxr/man/case1_pkpd.html) flags
59 of 180 avatars, 31 of them for receiving the dose their arm was
assigned. Scored within arm it flags 1, and a patient given twice their
arm's dose is still flagged. Strata holding fewer than five subjects are
scored against the whole cohort instead, since a scale estimated from
four patients describes the four rather than the one being screened.
With no `strata` declared, every subject is scored against the cohort,
as before.

## See also

[`compare_pmx_distributions()`](https://iamstein.github.io/synpmx/reference/compare_pmx_distributions.md),
[`compare_pmx()`](https://iamstein.github.io/synpmx/reference/compare_pmx.md).

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
flag_identifiable_subjects(synthetic, roles)
#> PMX outlier / identifiability check: 1 of 30 subjects flagged
#> Flag = a robust outlier in follow-up time, dose count, dose magnitude, or DV value.
#> 
#> Twelve most unusual:
#>  subject_id follow_up_time n_doses max_dose max_dv outlier_axes flagged
#>          32             20       2       92    133     DV value    TRUE
#>          31             20       2       87     77                FALSE
#>          33             20       2      103   91.6                FALSE
#>          34             20       2      107    101                FALSE
#>          35             20       2     93.6     73                FALSE
#>          36             20       2     88.6   79.8                FALSE
#>          37             20       2     91.9   72.8                FALSE
#>          38             20       2      103   76.8                FALSE
#>          39             20       2      124   67.7                FALSE
#>          40             20       2      114   51.4                FALSE
#>          41             20       2     83.3   89.8                FALSE
#>          42             20       2      111   84.6                FALSE
#> ... 18 more row(s) in the returned table.
#> 
#> Source-derived; not releasable unless separately public or privately budgeted.
```
