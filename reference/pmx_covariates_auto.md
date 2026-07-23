# Declare bootstrap-resampled covariates by column name

A low-ceremony alternative to
[`pmx_covariates()`](https://iamstein.github.io/synpmx/reference/pmx_covariates.md)
for a long list of covariates whose fidelity does not matter. Instead of
a public range or level set per column, the columns are named and their
values are drawn directly from the source data: a uniform draw over the
(clipped) observed range for continuous columns, and a proportional
resample for categorical ones. Column type is detected from the data at
fit time.

## Usage

``` r
pmx_covariates_auto(names, clip = c(0.01, 0.99))
```

## Arguments

- names:

  Character vector of covariate column names.

- clip:

  Two probabilities giving the quantiles a continuous column is clipped
  to before its range is taken, so the exact minimum and maximum are not
  exposed. Defaults to the 1st and 99th percentiles. Pass `NULL` to use
  the raw observed minimum and maximum, matching `synadam` exactly.

## Value

A `pmx_covariates` object of bootstrap covariates.

## Details

This is the approach used by Novartis's `synadam`, and it is **not**
differentially private: it exposes the data-derived support of each
column. A model that uses it is marked as having non-private covariates,
and its privacy report says so. Use it only inside a trusted
environment, and never when the covariate columns may cross a trust
boundary.
