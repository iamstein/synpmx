# Fit a subject-level differentially private PMX population generator

This is the only stage that reads source patient data. It
deterministically bounds each complete subject contribution, invokes a
validated DP backend for every released source-dependent summary, and
returns no raw patient records. There is deliberately no fitting seed.

## Usage

``` r
fit_private_pmx(
  data,
  roles,
  endpoints,
  epsilon,
  delta,
  bounds,
  public_design,
  contribution_limits,
  budget_allocation,
  delta_justification = NULL,
  backend = "opendp",
  public_source = FALSE
)
```

## Arguments

- data:

  Confidential PMX event data (or a public fixture when
  `public_source = TRUE`).

- roles:

  Explicit column roles from
  [`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md).

- endpoints:

  Named endpoint declarations from
  [`pmx_endpoint()`](https://iamstein.github.io/synpmx/reference/pmx_endpoint.md).

- epsilon, delta:

  Explicit requested subject-level privacy parameters.

- bounds:

  Public clipping domains from
  [`pmx_bounds()`](https://iamstein.github.io/synpmx/reference/pmx_bounds.md).

- public_design:

  Public schema and protocol metadata from
  [`pmx_public_design()`](https://iamstein.github.io/synpmx/reference/pmx_public_design.md).

- contribution_limits:

  Public contribution limits from
  [`pmx_contribution_limits()`](https://iamstein.github.io/synpmx/reference/pmx_contribution_limits.md).

- budget_allocation:

  Explicit fractions from
  [`pmx_budget_allocation()`](https://iamstein.github.io/synpmx/reference/pmx_budget_allocation.md).

- delta_justification:

  Required justification when `delta > 0`.

- backend:

  Production fitting supports `"opendp"`. `"public"` is a noiseless
  structural backend allowed only with `public_source = TRUE` for public
  examples; it makes no privacy claim.

- public_source:

  Logical assertion that the complete input is already public. Never set
  this for confidential or patient data.

## Value

A `private_pmx_model`. It contains public configuration, noisy
population summaries, privacy accounting, and one release-ledger entry;
it contains no raw IDs, rows, profiles, templates, or residuals.
