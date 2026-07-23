# Concentration-time profile for a public structural model

Concentration-time profile for a public structural model

## Usage

``` r
.pk_profile(model, time, doses, dose_times, params = NULL, duration = 0)
```

## Arguments

- model:

  A `pmx_structural_model`.

- time:

  Numeric times, measured from the first dose.

- doses:

  Numeric dose amounts.

- dose_times:

  Times at which those doses are given.

- params:

  Named parameter vector overriding the model's typical values.

- duration:

  Infusion duration, recycled over doses.

## Value

Numeric concentrations, one per `time`.
