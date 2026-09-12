# Declare one public baseline covariate

A covariate is either continuous, with a public plausible `range`, or
categorical, with public `levels`. The range or level set must be chosen
without inspecting the confidential data.

## Usage

``` r
pmx_covariate(
  range = NULL,
  levels = NULL,
  source,
  median = NULL,
  cv = NULL,
  distribution = NULL,
  prob = NULL,
  integer = FALSE
)
```

## Arguments

- range:

  Two increasing numbers bracketing a continuous covariate.

- levels:

  Character levels of a categorical covariate.

- source:

  Required provenance string.

- median:

  Typical value of a continuous covariate. For the default lognormal
  this is the median, equivalently the geometric mean, following the
  same convention as a population parameter in
  [`pmx_structural_model()`](https://iamstein.github.io/synpmx/reference/pmx_structural_model.md).

- cv:

  Coefficient of variation of a continuous covariate, as a proportion:
  `0.18` for 18 percent. Requires `median`.

- distribution:

  `"lognormal"` (the default when `median` and `cv` are given) or
  `"normal"`. Lognormal is right-skewed and cannot go negative, which
  suits body size and clearance; normal suits age and a laboratory value
  that is roughly symmetric.

- prob:

  Probability per level of a categorical covariate, in the order of
  `levels`. Normalized if it does not sum to one.

- integer:

  Round draws to whole numbers. For a covariate recorded as a count of
  years or a score.

## Value

A `pmx_covariate`.

## The range is a bound, and the distribution is separate

`range` is the clipping bound that makes a differentially private
release of this covariate possible, so it should be generous: a value
outside it is pulled to the edge. `median` and `cv` say where the
population actually sits inside that bound. Supplying neither leaves the
generator with only the range to work from, and it falls back to a
normal draw centred on the midpoint with the range spanning six standard
deviations. That couples two unrelated things – widening the bound for
safety also widens the distribution – and for a right-skewed covariate
such as body weight the shape is wrong as well. State `median` and `cv`
for anything whose distribution matters, and see
[`pmx_covariates_reference()`](https://iamstein.github.io/synpmx/reference/pmx_covariates_reference.md)
for defensible starting values.

For a categorical covariate, `prob` is the same point: without it every
level is equally likely, which is rarely what a trial looked like.

## See also

[`pmx_covariates()`](https://iamstein.github.io/synpmx/reference/pmx_covariates.md)
to collect them,
[`pmx_covariates_reference()`](https://iamstein.github.io/synpmx/reference/pmx_covariates_reference.md)
for reference-population values,
[`pmx_covariates_auto()`](https://iamstein.github.io/synpmx/reference/pmx_covariates_auto.md)
to resample from the data instead.

## Examples

``` r
# Body weight: a generous bound, and the distribution stated separately.
pmx_covariate(range = c(35, 160), median = 75, cv = 0.18,
              source = "protocol inclusion criteria")
#> $type
#> [1] "continuous"
#> 
#> $range
#> [1]  35 160
#> 
#> $median
#> [1] 75
#> 
#> $cv
#> [1] 0.18
#> 
#> $distribution
#> [1] "lognormal"
#> 
#> $integer
#> [1] FALSE
#> 
#> $source
#> [1] "protocol inclusion criteria"
#> 
#> attr(,"class")
#> [1] "pmx_covariate"

# A trial that enrolled seven men for every three women.
pmx_covariate(levels = c("M", "F"), prob = c(0.7, 0.3),
              source = "protocol enrollment targets")
#> $type
#> [1] "categorical"
#> 
#> $levels
#> [1] "M" "F"
#> 
#> $prob
#> [1] 0.7 0.3
#> 
#> $source
#> [1] "protocol enrollment targets"
#> 
#> attr(,"class")
#> [1] "pmx_covariate"
```
