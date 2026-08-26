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
flag_identifiable_subjects(data, roles, threshold = 3.5, min_relative_gap = 1)
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

- min_relative_gap:

  How far a subject must sit from the nearest other subject in its
  comparison group before that outlier counts, as a fraction of the
  group's median on that axis. Default 1, chosen against ordinary
  between-subject variability: a subject 40 minutes from a 216-hour
  follow-up is not a finding, and one followed to 1730 hours where
  everybody else has finished by 672 is. Lower it to widen the net – 0.5
  adds the merely striking, 0 leaves the z alone.

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

A subject is flagged on an axis when **both** hold: it is a robust
outlier there, and its gap to the nearest other subject in the
comparison group is at least `min_relative_gap` of that group's median.
The second condition is what makes the answer readable. A modified z
says a subject is unusual for its cohort; it cannot say the difference
is big enough to single anybody out, and on protocol-driven data that is
the question that decides. Every subject in
[`xgxr::mad`](https://rdrr.io/pkg/xgxr/man/mad.html) completed the
study, so follow-up times sit within two hours of each other and the
median absolute deviation of that is six minutes: a subject forty
minutes from its arm's median scores 4.7 and is identifiable to nobody.
Requiring the gap as well takes that study from 6 flagged of 60 to 0,
`pheno_sd` from 25 of 59 to 0 and `mavoglurant` from 41 of 120 to 0,
while leaving a subject given twice its arm's dose, or followed three
times as long as anyone else, flagged. The gap is measured to the
nearest other subject rather than to the median, because two subjects
sharing an extreme value single out neither – the same reasoning as
`min_pattern_share = 2` on visit sets.

The default of 1 is set against ordinary between-subject variability
rather than against any cohort: a 30-50% coefficient of variation is
unremarkable on a pharmacokinetic parameter, so a subject separated by
less than that is inside the noise of the study and could not be picked
out of it. At 1 the separation must exceed the group's whole median
value, which is the glaring case and nothing smaller –
[`nlmixr2data::wbcSim`](https://nlmixr2.github.io/nlmixr2data/reference/wbcSim.html)
reports the two patients followed to 1730 and 4580 hours in a cohort
that otherwise ends by 672, and not the one at 1130.

On that setting no synthetic dataset in
[`vignette("avatar-public-data-examples")`](https://iamstein.github.io/synpmx/articles/avatar-public-data-examples.md)
reports anybody, which is the screen working rather than idling: an
avatar is a blend of several donors, so producing a subject that extreme
takes a defect, and this is the row that would say so.

This matters because
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
anything. The clearest case is a patient given a *higher* arm's dose:
cohort-wide that dose is thirty other patients' dose, so nothing is
reported, and it is only unusual next to the arm the patient was
actually allocated to. Strata holding fewer than five subjects are
scored against the whole cohort instead, since a scale estimated from
four patients describes the four rather than the one being screened.
With no `strata` declared, every subject is scored against the cohort.

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
#> PMX outlier / identifiability check: 0 of 30 subjects flagged
#> Flag = a robust outlier in follow-up time, dose count, dose magnitude, or DV value,
#> that is also at least 100% of the group median away from the nearest
#> other subject.
#> 
#> Twelve most unusual:
#>  subject_id follow_up_time n_doses max_dose max_dv outlier_axes flagged
#>          31             20       2       87     77                FALSE
#>          32             20       2       92    133                FALSE
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
