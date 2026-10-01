# What leaves the study: the release of a fitted model

The part of a fitted model that synthetic data is generated from, and
nothing else.
[`synpmx_model_estimate()`](https://iamstein.github.io/synpmx/reference/synpmx_model_estimate.md)
returns a fit that also carries diagnostics for the person who ran it:
the candidate comparison and its AIC values, the starting values, the
correlations between covariates and the individual random effects, the
design notes and timings. None of those is needed to generate, so none
of them is in the release.
[`synpmx_model_generate()`](https://iamstein.github.io/synpmx/reference/synpmx_model_generate.md)
attaches the release, never the fit, to the data it returns.

## Usage

``` r
model_release(fitted_model)
```

## Arguments

- fitted_model:

  A `pmx_fitted_model` from
  [`synpmx_model_estimate()`](https://iamstein.github.io/synpmx/reference/synpmx_model_estimate.md),
  or a release, which is returned unchanged.

## Value

A `pmx_model_release`, which is also a `pmx_fitted_model`.

## Details

Save the release, not the fit, when a fingerprint has to be carried out
of the environment that holds the study:
[`synpmx_model_generate()`](https://iamstein.github.io/synpmx/reference/synpmx_model_generate.md)
accepts either and produces the same data from both.

## See also

[`model_privacy_checks()`](https://iamstein.github.io/synpmx/reference/model_privacy_checks.md),
[`synpmx_model_generate()`](https://iamstein.github.io/synpmx/reference/synpmx_model_generate.md),
[`model_report()`](https://iamstein.github.io/synpmx/reference/model_report.md).
