# Check whether a private release is worth its budget, before spending it

Reports `f = d / (epsilon * N)`, the fraction of each prior's width that
survives as noise, and the resulting fold-error. This consumes no
privacy budget and reads no data: it depends only on the configuration.

## Usage

``` r
pmx_preflight(priors, epsilon, n_subjects, covariates = NULL)
```

## Arguments

- priors:

  A
  [`pmx_priors()`](https://iamstein.github.io/synpmx/reference/pmx_priors.md)
  object.

- epsilon:

  The privacy budget under consideration.

- n_subjects:

  Number of subjects in the fit.

- covariates:

  Optional
  [`pmx_covariates()`](https://iamstein.github.io/synpmx/reference/pmx_covariates.md),
  so the reported `d` matches a fit that also releases covariate
  summaries.

## Value

A `pmx_preflight` report. The arithmetic behind it is worked through at
<https://iamstein.github.io/synpmx/articles/privacy-background.html>.

## Details

The fold-error is `exp(f * span)`, capped at the prior's half-width
because clipping prevents a release from landing outside the prior. The
uncapped form is accurate for `f` below roughly 0.25 and increasingly
pessimistic above it; see the feasibility article at
<https://iamstein.github.io/synpmx/articles/feasibility.html>.
