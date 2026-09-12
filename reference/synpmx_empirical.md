# Generate a dataset from a dense differentially private release

Rather than asserting the curve shape, this mode measures it: it
releases noised summaries for the subject count, event and regimen
structure, observation timing, endpoint trajectories, baseline
covariates, and censoring, then rebuilds subjects from those summaries
alone.

## Usage

``` r
synpmx_empirical(
  data,
  roles,
  endpoints,
  epsilon,
  delta,
  bounds,
  public_design,
  contribution_limits,
  budget_allocation,
  n_subjects = NULL,
  seed = 123,
  n_datasets = 1L,
  delta_justification = NULL,
  backend = "opendp",
  public_source = FALSE
)
```

## Arguments

- data:

  The confidential dataset.

- roles:

  A
  [`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md)
  declaration for `data`. It names the columns of the output too, so the
  synthetic table wears the source study's schema.

- endpoints:

  Named list of
  [`pmx_endpoint()`](https://iamstein.github.io/synpmx/reference/pmx_endpoint.md)
  declarations.

- epsilon:

  The privacy budget. A governance decision, not a default.

- delta:

  Additive slack in the probability bound. The implemented Laplace
  releases spend none, so realized accounting reports `delta = 0`.

- bounds:

  Public clipping domains from
  [`pmx_bounds()`](https://iamstein.github.io/synpmx/reference/pmx_bounds.md).

- public_design:

  A
  [`pmx_public_design()`](https://iamstein.github.io/synpmx/reference/pmx_public_design.md).

- contribution_limits:

  A
  [`pmx_contribution_limits()`](https://iamstein.github.io/synpmx/reference/pmx_contribution_limits.md).

- budget_allocation:

  A
  [`pmx_budget_allocation()`](https://iamstein.github.io/synpmx/reference/pmx_budget_allocation.md)
  splitting `epsilon` across release groups.

- n_subjects:

  Number of subjects to generate. Defaults to the released noisy count.

- seed:

  Ordinary generation seed. Unrelated to privacy noise.

- n_datasets:

  Number of datasets to draw from the single release. One dataset is
  returned directly; several are returned as a list.

- delta_justification:

  Required when `delta > 0`.

- backend:

  Privacy backend. Defaults to the validated OpenDP adapter and fails
  closed if it is unavailable.

- public_source:

  Assert that `data` is genuinely public. Required by, and only
  meaningful for, `backend = "public"`, which makes no DP claim.

## Value

A data frame in the source event-table schema, carrying its release so
that
[`privacy_report()`](https://iamstein.github.io/synpmx/reference/privacy_report.md)
and
[`synpmx_generate()`](https://iamstein.github.io/synpmx/reference/synpmx_generate.md)
can read it. A list of such data frames when `n_datasets > 1`.

## Not developed further

This engine exists so both ends of the privacy and utility axis are in
the package and the tradeoff can be checked. It is not being developed,
it has no documents of its own, and the one section that explains why it
is rarely the right choice is "Why not release more quantities?" in the
[methods
article](https://iamstein.github.io/synpmx/articles/synpmx-methods.html).
Every clipping range, contribution limit and budget share is an explicit
public input, and getting one wrong breaks the guarantee without any
error being raised, so do not reach for it on the strength of this help
page alone.

It asserts less than
[`synpmx_calibrated()`](https://iamstein.github.io/synpmx/reference/synpmx_calibrated.md)
but releases far more numbers, so one epsilon is split many ways.
Utility therefore collapses below a few hundred subjects; this mode
earns its keep on large pooled corpora.

**This function spends privacy budget.** Calling it again spends the
budget again. To draw further datasets from the release you have already
paid for, use
[`synpmx_generate()`](https://iamstein.github.io/synpmx/reference/synpmx_generate.md)
or ask for several at once with `n_datasets`.

## Maintenance status

A secondary, provided-as-is path.
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
is the primary, maintained method. The differentially private modes are
complete and tested but not under active development, carry known open
findings, and have not been independently privacy-audited. Use them to
demonstrate the privacy/utility tradeoff, not as a production release
mechanism; a real regulated release needs specialist review and the
external OpenDP backend. Requires a session-level acknowledgment: call
[`synpmx_enable_dp_engines()`](https://iamstein.github.io/synpmx/reference/synpmx_enable_dp_engines.md)
once before this function will run (unless `backend = "public"`, which
makes no DP claim).

## See also

[`synpmx_generate()`](https://iamstein.github.io/synpmx/reference/synpmx_generate.md)
to draw more datasets for free,
[`privacy_report()`](https://iamstein.github.io/synpmx/reference/privacy_report.md)
for the realized accounting.
