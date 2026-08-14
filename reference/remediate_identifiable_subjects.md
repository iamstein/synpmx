# Remove or shorten the subjects `flag_identifiable_subjects()` flags

Applies a remediation policy to the outliers found by
[`flag_identifiable_subjects()`](https://iamstein.github.io/synpmx/reference/flag_identifiable_subjects.md).
By default a subject flagged **only** for an unusually long follow-up is
*truncated* back to the cohort's longest ordinary follow-up (its late
rows dropped), and a subject flagged for **any other** reason – an
extreme DV, a rare dose level, or an unusual dose count – is *dropped*
entirely.

## Usage

``` r
remediate_identifiable_subjects(
  data,
  roles,
  source = NULL,
  time = c("truncate", "drop", "keep"),
  other = c("drop", "keep"),
  threshold = 3.5,
  seed = NULL,
  max_tries = 20L
)
```

## Arguments

- data:

  A PMX dataset, typically the synthetic output.

- roles:

  Explicit roles from
  [`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md).

- source:

  Optional source PMX data. When given, dropped subjects are replaced by
  fresh avatars generated from it, so the cohort size is preserved. When
  `NULL` (default), dropped subjects are simply removed.

- time:

  Action for a subject whose *only* outlier axis is follow-up time:
  `"truncate"` (default) to shorten it to the longest ordinary
  follow-up, `"drop"` to remove it, or `"keep"` to leave it.

- other:

  Action for a subject flagged for any non-time reason: `"drop"`
  (default) or `"keep"`.

- threshold:

  Passed to
  [`flag_identifiable_subjects()`](https://iamstein.github.io/synpmx/reference/flag_identifiable_subjects.md).

- seed:

  Reproducibility seed for the replacement generation. The caller's
  random-number state is restored by
  [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md).

- max_tries:

  Maximum regeneration batches when refilling dropped subjects.

## Value

`data` with the policy applied, carrying attributes `dropped`,
`truncated` (affected subject ids), `replaced` (count refilled), and
`horizon` (the follow-up truncation used, or `NA`).

## Details

The split is deliberate. Truncation is offered only for an unusually
*long* follow-up, the one structural outlier a value-level edit can
genuinely fix: shortening a long timeline leaves a shorter but ordinary
subject. Everything else is dropped – an unusually *short* follow-up has
nothing to trim, an extreme-DV subject is elevated across its whole
trajectory so removing points would only mangle it, and a rare dose
cannot be trimmed without breaking the regimen. A subject flagged for
both a long follow-up and another reason is also dropped, since
truncation would not resolve the other reason.

When `source` is supplied, each dropped subject is **replaced**: fresh
avatars are generated from `source`, screened by the same policy, and
appended (with new ids) until the cohort is back to its original size.
So the output keeps the same number of subjects, minus any it could not
refill within `max_tries`. Truncation keeps its subject, so it never
triggers a replacement.

Detection is per subject, so one long-followed patient is truncated once
and one extreme patient dropped-and-replaced once – there is no
row-level outlier spray. With replacement, this is a self-contained
alternative to preventing structural outliers at generation time
(skeleton sampling).

## See also

[`flag_identifiable_subjects()`](https://iamstein.github.io/synpmx/reference/flag_identifiable_subjects.md).

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
cleaned <- remediate_identifiable_subjects(synthetic, roles, source = data)
#> remediate_identifiable_subjects(): dropped 1, truncated 0, replaced 1.
```
