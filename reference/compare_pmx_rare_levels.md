# Which rare source levels reached the synthetic output

Censuses every categorical axis – `strata` and each non-numeric
covariate – on both sides, and marks the levels **too few source
patients held** to be safely copied out. That floor is
`min_pattern_share`, the same rule the generator applies to visit sets
and the scorecard applies to copied vectors, so no new threshold is
introduced.

## Usage

``` r
compare_pmx_rare_levels(source, synthetic, roles, floor = NULL)
```

## Arguments

- source:

  Source PMX data.

- synthetic:

  Generated synthetic PMX data. When it carries a `"pmx_settings"`
  attribute, its `min_pattern_share` is used as the floor.

- roles:

  Explicit roles from
  [`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md).

- floor:

  Levels held by fewer than this many source patients are `exposed`.
  Left `NULL` it is taken from the run, and `2` otherwise – the lowest
  value that means "more than one real patient".

## Value

A `pmx_rare_levels` data frame, one row per categorical column and
level, with `source_patients`, `synthetic_patients`, `exposed` and
`reached`. Zero rows when the roles declare no categorical axis.

## Details

Read the `exposed` rows, and among them the ones that `reached` the
output. A level held by two real patients, appearing in a released
table, says that someone with that attribute was in this study; for a
named trial with public inclusion criteria that can be close to
identifying on its own, and no cohort size helps. The remedies are
upstream of generation: drop the covariate from `covariates`, or
collapse its rare levels before generating.

This reads real patient data on both sides and is marked
`"restricted_not_releasable"`.

**What it cannot see.** Rarity *in the world*. If every living carrier
of a mutation is in this study, the source count is the whole population
and looks unremarkable. And it censuses each column on its own: with `d`
covariates there are `2^d` combinations that could single a patient out,
and enumerating them is not something this does.

## See also

[`synpmx_scorecard()`](https://iamstein.github.io/synpmx/reference/synpmx_scorecard.md),
which reports this as row B5b,
[`compare_pmx_distributions()`](https://iamstein.github.io/synpmx/reference/compare_pmx_distributions.md),
[`vignette("scorecard-synthetic-data-checks")`](https://iamstein.github.io/synpmx/articles/scorecard-synthetic-data-checks.md).

## Examples

``` r
data <- pmx_simulated_fixture(20)
roles <- pmx_roles(
  id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
  cmt = "CMT", dvid = "DVID", covariates = c("WT", "SEX")
)
synthetic <- suppressWarnings(synpmx_avatar(data, roles, seed = 1))
#> synpmx_avatar(): dropped 8 undeclared column(s): NTIME, TAD, OCC, RATE, MDV, CENS, LIMIT, AGE.
#>   Declare a column in `keep` to carry it through verbatim.
compare_pmx_rare_levels(data, synthetic, roles)
#> PMX rare-level census (source against synthetic)
#> 
#> 0 level(s) held by fewer than 2 source patients; 0 of them reached the output.
#> 
#> Every level in the output is one that at least 2 source patients held.
#> 
#> Source-derived; not releasable.
```
