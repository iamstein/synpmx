# Declare public event-design information

Every supplied value is treated as public protocol or schema metadata
and is not charged to the privacy budget. Omitted regimen and timing
quantities are estimated through budgeted private summaries.

## Usage

``` r
pmx_public_design(
  schema,
  dose_times = NULL,
  dose_interval = NULL,
  n_doses = NULL,
  dose_amount = NULL,
  dose_rate = NULL,
  infusion_duration = NULL,
  dose_evid = 1,
  dose_cmt = 1,
  endpoint_grids = NULL,
  endpoint_occasion_grids = NULL,
  endpoint_cmt = NULL,
  category_levels = NULL,
  defaults = NULL,
  time_jitter_sd = 0.02,
  subject_count = NULL
)
```

## Arguments

- schema:

  A schema from
  [`pmx_schema()`](https://iamstein.github.io/synpmx/reference/pmx_schema.md).
  It is required.

- dose_times, dose_interval, n_doses, dose_amount, dose_rate:

  Optional public regimen overrides. Normally omit these so the budgeted
  event summaries infer the regimen from the input data. Supply them
  only when they are independently public protocol facts, never by
  inspecting confidential records outside the private fit.

- infusion_duration:

  Optional public infusion duration.

- dose_evid, dose_cmt:

  Public event and dose-compartment values.

- endpoint_grids:

  Optional named list of fixed discretization bases on each endpoint's
  declared scientific clock. These are not sampling schedules. When
  omitted, the package constructs generic bases from public bounds and
  contribution limits and privately learns their occupancy.

- endpoint_occasion_grids:

  Optional named list of endpoint-specific, public dose-occasion
  sampling schedules. This exceptional override should be used only for
  an independently public protocol. Each endpoint entry is a list named
  by positive occasion number; omitted occasions generate no
  observations.

- endpoint_cmt:

  Named list or vector of public observation compartments.

- category_levels:

  Named lists of allowed values for categorical covariates and subject
  properties. Supplying levels forces even a numeric covariate to be
  treated categorically. Factor levels may instead come from the public
  schema. These levels are public domains, not values discovered by
  inspecting confidential records outside the private fit.

- defaults:

  Named values for otherwise unmodeled public columns.

- time_jitter_sd:

  Nonnegative generation-time jitter scale, expressed as a fraction of
  the closest public grid spacing.

- subject_count:

  Optional independently public source subject count. When omitted,
  generation uses the privacy-accounted fitted count.

## Value

A `pmx_public_design` object.
