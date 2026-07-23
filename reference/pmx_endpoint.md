# Declare endpoint scientific-clock behavior

Endpoint alignment is public scientific metadata. It controls how a
small fixed-dimensional trajectory summary is built and how new
trajectories are generated; it is not a PK or PD model.

## Usage

``` r
pmx_endpoint(
  dvid = NULL,
  alignment,
  transform = c("auto", "log", "identity"),
  shape,
  units = NULL,
  grid = NULL,
  cmt = NULL,
  subject_sd = 0.2,
  residual_sd = 0.08,
  censoring = NULL
)
```

## Arguments

- dvid:

  Public DVID value for this endpoint, or `NULL` when no DVID column is
  used.

- alignment:

  One of `"dose_relative"`, `"study_time"`, `"occasion"`, or `"hybrid"`.

- transform:

  One of `"log"`, `"identity"`, or `"auto"`. `"auto"` uses only the
  public DV bounds: a nonnegative domain uses an offset log scale.

- shape:

  One of `"occasion"` or `"global"`; this is a broad public shape
  expectation, not a fitted structural model.

- units:

  Optional public unit label.

- grid:

  Optional strictly increasing public grid on the declared clock. This
  is a discretization basis, not a sampling schedule. When omitted, a
  generic basis is constructed from public bounds and contribution
  limits; sampling-cell occupancy is then learned from the fitted data.

- cmt:

  Optional public observation compartment value.

- subject_sd, residual_sd:

  Public generation variability multipliers.

- censoring:

  Optional public list with `left`, `right`, or a two-value `interval`.
  Source-dependent censoring frequencies are separately private.

## Value

A `pmx_endpoint` declaration.
