# Summarize a trial into the quantities a synthetic copy is built from

The only stage that reads patient data. Reduces each subject's
trajectories and baseline covariates to principal-component scores,
models those scores against the arm, and fits a dosing model and a visit
model per arm. The returned object holds nothing but summaries: no
patient row survives it.

## Usage

``` r
synpmx_pca_summarize(
  data,
  roles,
  seed = NULL,
  dose_term = c("factor", "log"),
  pca_variance = 0.9,
  n_components = NULL,
  min_column_patients = NULL,
  min_arm_patients = 3L
)
```

## Arguments

- data:

  Source PMX event data.

- roles:

  Explicit column roles from
  [`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md),
  including `nominal_time`.

- seed:

  Seed for the one random step in summarizing: censored values are
  replaced by a draw inside the censoring region before the basis is
  fitted, so that a column where most subjects sit at the assay limit
  describes the patients rather than the assay. The boundary is
  reapplied at generation.

- dose_term:

  How dose enters the score model. `"factor"` gives each arm its own
  mean score vector and its own residual covariance. `"log"` regresses
  the scores on [`log1p()`](https://rdrr.io/r/base/Log.html) of the
  total dose, which spends one coefficient rather than one mean per arm
  and extrapolates to doses the study did not run, at the cost of
  assuming the dose-response is log-linear and that every arm has the
  same between-subject spread. A lower limit of quantification breaks
  both assumptions, which is why `"factor"` is the default.

- pca_variance:

  Cumulative variance the retained components must reach.

- n_components:

  Number of components, overriding `pca_variance`.

- min_column_patients:

  Minimum distinct patients holding an observation in a grid cell for
  that cell to be modelled. Defaults to the larger of 3 and a tenth of
  the cohort.

- min_arm_patients:

  Minimum patients in every arm. An arm below it has no spread of its
  own to model, so its mean score vector and its covariance would
  describe the one or two patients in it. The function refuses rather
  than summarizing them; pool the arm, drop the column from `strata`, or
  exclude those patients before calling.

## Value

A `pmx_trial_summary`.

## Details

Run this, look at what it produced with
[`pca_report()`](https://iamstein.github.io/synpmx/reference/pca_report.md),
[`pca_dosing()`](https://iamstein.github.io/synpmx/reference/pca_dosing.md),
[`pca_visits()`](https://iamstein.github.io/synpmx/reference/pca_visits.md)
and
[`pca_components()`](https://iamstein.github.io/synpmx/reference/pca_components.md),
then pass it to
[`synpmx_pca_generate()`](https://iamstein.github.io/synpmx/reference/synpmx_pca_generate.md).
Generation reads the model and nothing else, so what those four
functions show is the whole of what the synthetic data was built from.

`nominal_time` is required. The grid it names is the axis every feature
sits on, and dose rows and observation rows are placed on it together,
so it is the caller's statement about the protocol rather than something
inferred from recorded times.

No formal privacy claim is made.

## See also

[`synpmx_pca_generate()`](https://iamstein.github.io/synpmx/reference/synpmx_pca_generate.md),
[`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md),
[`pca_report()`](https://iamstein.github.io/synpmx/reference/pca_report.md).

## Examples

``` r
data <- pmx_simulated_fixture(60)
roles <- pmx_roles(
  id = "ID", time = "TIME", nominal_time = "NTIME", dv = "DV", amt = "AMT",
  evid = "EVID", cmt = "CMT", dvid = "DVID", mdv = "MDV"
)
trial_summary <- synpmx_pca_summarize(data, roles)
trial_summary
#> A trial summary, from synpmx_pca_summarize()
#> 
#>   fitted on    60 patients, 1 arm(s): all (60) 
#>   endpoints    cp (8 visits modelled), pd (6 visits modelled) 
#>   covariates   none 
#>   components   1 (100% of variance) 
#>   dose term    factor 
#>   dosing       2 planned cycle(s) per arm | no reductions, interruptions or early stops 
#> 
#> synpmx_pca_generate() reads this object and nothing else. To look inside it:
#>   pca_report()      what it read out of the source data
#>   pca_dosing()      the planned dose schedule, per arm
#>   pca_dose_rates()  reduction, interruption and discontinuation
#>   pca_visits()      the probability of a visit, per arm
#>   pca_components()  the loadings, over time
pca_report(trial_summary)
#> What the PCA fit read out of the source data
#> 
#>   subjects: 60  components retained: 1 
#> 
#>             quantity                                                      what
#>           visit grid                      Nominal times modelled, per endpoint
#>      feature centers                      Mean of each grid cell and covariate
#>       feature scales                            Standard deviation of the same
#>             loadings                        Component loadings on each feature
#>          score means                                Mean score vector, per arm
#>     score covariance                    Residual covariance between components
#>  endpoint transforms                             Log or identity, per endpoint
#>         assay limits                          Censoring boundary, per endpoint
#>         dosing model Planned cycles, the dose ladder, and three rates, per arm
#>          visit model        Probability of a visit, per arm, endpoint and time
#>        arm constants                Strata and kept columns, one value per arm
#>  numbers min_patients
#>       14           60
#>       14           60
#>       14           60
#>       14           60
#>        1           60
#>        1           60
#>        2           60
#>        0           60
#>        8           60
#>       14           60
#>        0           60
```
