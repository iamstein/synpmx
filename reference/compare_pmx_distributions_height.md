# How tall the distribution figure should be drawn

The figure gets one panel per endpoint and per baseline covariate, so a
study with five endpoints needs several times the height of a study with
one. `fig.height` is fixed before a chunk runs and cannot be read off
the plot, so pass this to the chunk option instead of letting every
figure inherit one default and arrive squashed or stranded in white
space.

## Usage

``` r
compare_pmx_distributions_height(source, roles, per_row = 2.1, minimum = 3)
```

## Arguments

- source:

  Source PMX data.

- roles:

  Explicit roles from
  [`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md).

- per_row:

  Inches per row of panels.

- minimum:

  Inches below which the figure is never drawn, so a single-panel study
  does not come out as a letterbox.

## Value

One number, in inches.

## Details

    ```{r, fig.height = compare_pmx_distributions_height(raw, roles)}
    compare_pmx_distributions(raw, synthetic, roles)
    ```

## See also

[`compare_pmx_distributions()`](https://iamstein.github.io/synpmx/reference/compare_pmx_distributions.md).

## Examples

``` r
data <- pmx_simulated_fixture(20)
roles <- pmx_roles(
  id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
  cmt = "CMT", dvid = "DVID", covariates = c("WT", "SEX")
)
compare_pmx_distributions_height(data, roles)
#> [1] 4.8
```
