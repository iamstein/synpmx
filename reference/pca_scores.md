# The score model each arm is generated from

One row per arm and retained component, giving the mean score that arm
is centred on and the residual standard deviation it is scattered by. A
generated subject's scores are their arm's `mean` plus a draw whose
spread is `sd`, so these two columns are the whole of the
between-subject variability the synthetic data will have.

## Usage

``` r
pca_scores(x)
```

## Arguments

- x:

  A dataset from
  [`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md),
  or its trial summary.

## Value

A data frame with `arm`, `component`, `mean` and `sd`. Marked
`"restricted_not_releasable"`.

## Details

The `sd` values differ by arm on purpose. An arm sitting on an assay
limit is genuinely tighter than one well above it, and giving every arm
the pooled spread would smear the low arms upward.

Both `dose_term` settings report the same shape. Under `"factor"` the
mean is the arm's own; under `"log"` it is the regression evaluated at
that arm's planned total dose, and `sd` is the shared residual.

## See also

[`synpmx_pca_summarize()`](https://iamstein.github.io/synpmx/reference/synpmx_pca_summarize.md),
[`pca_components()`](https://iamstein.github.io/synpmx/reference/pca_components.md),
[`pca_report()`](https://iamstein.github.io/synpmx/reference/pca_report.md).

## Examples

``` r
data <- pmx_simulated_fixture(60)
roles <- pmx_roles(
  id = "ID", time = "TIME", nominal_time = "NTIME", dv = "DV", amt = "AMT",
  evid = "EVID", cmt = "CMT", dvid = "DVID", mdv = "MDV"
)
pca_scores(synpmx_pca_summarize(data, roles))
#>   arm component        mean       sd
#> 1 all       PC1 6.29589e-16 3.741645
```
