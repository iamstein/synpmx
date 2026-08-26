# Every feature the components are built on

One row per column of the matrix
[`synpmx_pca_summarize()`](https://iamstein.github.io/synpmx/reference/synpmx_pca_summarize.md)
decomposed: one per baseline covariate, and one per endpoint per
retained nominal time. This is the grid the whole method sits on, so it
is where to look first when a generated dataset is missing a visit or a
covariate.

## Usage

``` r
pca_features(x)
```

## Arguments

- x:

  A dataset from
  [`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md),
  or its trial summary.

## Value

A data frame with `feature`, `kind`, `endpoint`, `time`, `covariate`,
`level`, `patients`, `center`, `scale` and `transform`. Marked
`"restricted_not_releasable"`: the counts and moments are read from real
data.

## Details

`center` and `scale` are the column's mean and standard deviation on the
modelling scale, which is the log scale for a positive endpoint — see
`transform` in the same row. `patients` is how many subjects hold an
observation in that cell; cells held by fewer than `min_column_patients`
were dropped rather than modelled, so they do not appear here at all.

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
head(pca_features(synpmx_pca_summarize(data, roles)))
#>    feature          kind endpoint  time covariate level patients    center
#> 1 dv_cp__1 endpoint_cell       cp  0.25      <NA>  <NA>       60 0.3629025
#> 2 dv_cp__2 endpoint_cell       cp  1.00      <NA>  <NA>       60 2.1297321
#> 3 dv_cp__3 endpoint_cell       cp  2.00      <NA>  <NA>       60 1.4877196
#> 4 dv_cp__4 endpoint_cell       cp  6.00      <NA>  <NA>       60 0.6605288
#> 5 dv_cp__5 endpoint_cell       cp 12.25      <NA>  <NA>       60 0.4927632
#> 6 dv_cp__6 endpoint_cell       cp 13.00      <NA>  <NA>       60 2.1872652
#>        scale  transform
#> 1 0.05953933 log_offset
#> 2 0.08140515 log_offset
#> 3 0.07734409 log_offset
#> 4 0.06632245 log_offset
#> 5 0.06274771 log_offset
#> 6 0.08165741 log_offset
```
