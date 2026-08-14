# Are synthetic subjects sitting too close to real ones?

The measurement for donor blending, the one masking mechanism
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
applies to the *values* rather than the structure.
[`skeleton_uniqueness()`](https://iamstein.github.io/synpmx/reference/skeleton_uniqueness.md)
answers the structural question by counting who shares which schedule;
this answers the geometric one, by asking whether each subject's nearest
neighbour lies in its own dataset or the other one.

## Usage

``` r
compare_pmx_proximity(
  source,
  synthetic,
  roles,
  replicates = 50L,
  seed = 1L,
  pca_variance = 0.9
)
```

## Arguments

- source:

  Source PMX data.

- synthetic:

  Generated synthetic PMX data.

- roles:

  Explicit roles from
  [`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md).

- replicates:

  Split-half replicates used to build the null. Default 50.

- seed:

  Seed for the subsampling and splits. The caller's RNG is left
  untouched.

- pca_variance:

  Variance retained when both datasets are projected into a common
  profile space. Default 0.90, matching
  [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md).

## Value

A one-row `pmx_proximity` data frame: `adversarial_accuracy`,
`null_lower` / `null_upper` (the central 95% of the split-half null),
`verdict`, `n_compared` (patients per side, the same on both arms and in
the null), `n_null_replicates`, and the 5th-percentile nearest-neighbour
distances `synthetic_to_source_q05` and `source_to_source_q05`.

## Details

The reported statistic is a nearest-neighbour adversarial accuracy in
\\\[0, 1\]\\:

- **near 0.5** — a synthetic subject is no more like a real subject than
  one real subject is like another. This is the target.

- **toward 0** — synthetic subjects sit closer to real subjects than to
  each other. That is memorisation, and it is the privacy failure.

- **toward 1** — the two sets have separated. Privacy is fine and
  utility is not.

Raw distance to the closest real record is deliberately not the
headline. It has no natural scale, and measured against zero it mostly
tracks cohort size — the nearest of `N` points gets closer as `N` grows,
so a larger source would score worse while blending across more donors
actually makes it safer. The quantiles are still returned for context,
alongside the real-to-real quantiles they should be read against.

The null interval comes from running the identical statistic on two
halves of the **source** cohort, so every small-sample artefact is
present in the null and the observed value alike and cancels. At the
cohort sizes pharmacometrics works with, that interval is wide: this
will catch a blatant leak, not a subtle one. Treat a value inside the
interval as "nothing detected", never as "nothing there".

Marked `"restricted_not_releasable"`: it reads the source.

## See also

[`skeleton_uniqueness()`](https://iamstein.github.io/synpmx/reference/skeleton_uniqueness.md),
[`compare_pmx_distributions()`](https://iamstein.github.io/synpmx/reference/compare_pmx_distributions.md).

## Examples

``` r
data <- pmx_simulated_fixture(40)
roles <- pmx_roles(
  id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
  cmt = "CMT", dvid = "DVID", covariates = "WT"
)
synthetic <- suppressWarnings(synpmx_avatar(data, roles, seed = 1))
#> synpmx_avatar(): dropped 9 undeclared column(s): NTIME, TAD, OCC, RATE, MDV, CENS, LIMIT, AGE, SEX.
#>   Declare a column in `keep` to carry it through verbatim.
compare_pmx_proximity(data, synthetic, roles, replicates = 10)
#> PMX nearest-neighbour proximity check
#> 
#>   Question: is a synthetic patient closer to a real patient than real
#>     patients are to each other?
#>   Measured 0.800, on a scale where 0.5 means 'no closer' and is the target;
#>     0 would mean every synthetic patient is glued to a real one.
#>   Expected 0.356 to 0.666 if nothing were wrong. That interval is not
#>     assumed -- it is the same statistic run 10 times on two halves of the
#>     real cohort, 20 patients per half, which is also how many synthetic
#>     patients were compared.
#>   Verdict: Too far apart. The two sets have separated, so a classifier
#>     could tell them apart. That is a utility problem, not a privacy one.
#>   For context, distance to the nearest neighbour (5th percentile, so the
#>     closest pairs): synthetic-to-real 0.588 versus real-to-real 0.106.
#>     These are only comparable to each other; the units are PCA profile
#>     space.
#> 
#> Source-derived; not releasable unless separately public or privately
#> budgeted.
```
