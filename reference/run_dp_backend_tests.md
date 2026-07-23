# Run canonical checks against the configured DP backend

The checks verify construction, invocation, finite output, and OpenDP's
own privacy-map relation. They are implementation checks, not an
empirical proof of differential privacy.

## Usage

``` r
run_dp_backend_tests()
```

## Value

A structured check result. The function fails closed if OpenDP is
unavailable.
