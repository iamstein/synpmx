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
  time = c("truncate", "drop", "keep"),
  other = c("drop", "keep"),
  threshold = 3.5
)
```

## Arguments

- data:

  A PMX dataset, typically the synthetic output.

- roles:

  Explicit roles from
  [`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md).

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

## Value

`data` with the policy applied, carrying attributes `dropped` and
`truncated` (the affected subject ids) and `horizon` (the follow-up time
truncation used, or `NA`).

## Details

The split is deliberate. Truncation is offered only on the time axis
because it is the one structural outlier a value-level edit can
genuinely fix: shortening a long timeline leaves a shorter but ordinary
subject. An extreme-DV subject is elevated across its whole trajectory,
so removing points would only mangle it, and a rare dose cannot be
trimmed without breaking the regimen – so those subjects are dropped
rather than edited. A subject flagged for both a long follow-up and
another reason is dropped, since truncation would not resolve the other
reason.

This is a stop-gap. The durable fix is to sample each avatar's event
skeleton from the cohort so structural outliers are not generated in the
first place (`REV-026`); see
[`vignette("synpmx-method")`](https://iamstein.github.io/synpmx/articles/synpmx-method.md).

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
cleaned <- remediate_identifiable_subjects(synthetic, roles)
#> remediate_identifiable_subjects(): dropped 1 subject(s); truncated 0 subject(s).
```
