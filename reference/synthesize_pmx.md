# Synthesize a structurally faithful PMX dataset (AVATAR-style)

Samples complete subject event templates and fills them with
AVATAR-like, endpoint-specific blends of compatible subjects' baseline
covariates and longitudinal measurements. Event-control fields such as
EVID, AMT, RATE, CMT, and DVID are never averaged or independently
generated.

## Usage

``` r
synthesize_pmx(
  data,
  roles,
  n_subjects = NULL,
  seed = 123,
  event_method = "template",
  dv_method = "avatar_blend",
  k = 5,
  pca_variance = 0.9,
  subject_noise_sd = 0.15,
  residual_noise_sd = 0.05,
  residual_phi = 0.6,
  time_jitter = 0
)
```

## Arguments

- data:

  A source PMX data frame or tibble.

- roles:

  Explicit roles from
  [`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md).
  Columns listed in `roles$exclude` are omitted from the generated
  output.

- n_subjects:

  Number of synthetic subjects. `NULL` retains the source count.

- seed:

  Reproducibility seed. The caller's random-number state is restored on
  exit.

- event_method:

  Event generation method. The prototype supports `"template"`.

- dv_method:

  Measurement method. The prototype supports `"avatar_blend"`.

- k:

  Maximum number of compatible non-anchor donors.

- pca_variance:

  Fraction of usable profile variance retained for neighborhood
  distances.

- subject_noise_sd:

  Nonnegative subject perturbation multiplier.

- residual_noise_sd:

  Nonnegative within-trajectory noise multiplier.

- residual_phi:

  AR(1) correlation in observation order, strictly between -1 and 1.

- time_jitter:

  Standard deviation for coherent tied-time jitter. Zero, the default,
  leaves the event template's times unchanged.

## Value

An ordinary data frame or tibble with retained source columns, order,
and practical classes. A lightweight `pmx_settings` attribute records
the generator choices and endpoint transformations.

## Details

This is an AVATAR-inspired adaptation, not an exact implementation of
published AVATAR software. It creates synthetic data for model-workflow
exploration. It does not provide formal anonymization or preserve
scientific parameter or covariate-response relationships.

For selected compatible donors, randomized raw weights are
`Exp(1) / max(distance, epsilon) * 2^(-randomized_rank)`. They are
normalized and, when multiple donors are available, a dominant weight is
capped at 0.80 with its excess redistributed. The same subject weights
are used for covariates and all endpoints; weights are renormalized
locally when a donor lacks a requested endpoint/time value.

Positive-like endpoints use an offset log scale and are constrained to
be nonnegative after back-transformation. Other endpoints use the
identity scale. Transform choices and interpolation alignment are
recorded in the returned `pmx_settings` attribute.

## Examples

``` r
source <- data.frame(
  ID = rep(1:3, each = 4),
  TIME = rep(c(0, 0, 1, 2), 3),
  DV = c(0, 0.2, 2, 1, 0, 0.3, 3, 1.5, 0, 0.4, 4, 2),
  AMT = rep(c(100, 0, 0, 0), 3),
  EVID = rep(c(1L, 0L, 0L, 0L), 3),
  CMT = rep(c(1L, 2L, 2L, 2L), 3),
  WT = rep(c(60, 70, 80), each = 4)
)
roles <- pmx_roles("ID", "TIME", "DV", "AMT", "EVID", "CMT", NULL,
                   NULL, NULL, "WT")
synthetic <- synthesize_pmx(source, roles, n_subjects = 2, seed = 123)
#> Warning: Synthetic generation used documented small-group/profile fallbacks:
#> - `k` was reduced to 2 in at least one compatible event-pattern group.
validate_pmx(synthetic, roles)$valid
#> [1] TRUE
```
