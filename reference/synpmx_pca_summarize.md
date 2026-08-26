# Summarize a PMX dataset into a principal-component model

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

A `pmx_pca_model`.

## Details

Run this, look at what it produced with
[`pmx_pca_report()`](https://iamstein.github.io/synpmx/reference/pmx_pca_report.md),
[`pmx_pca_dosing()`](https://iamstein.github.io/synpmx/reference/pmx_pca_dosing.md),
[`pmx_pca_visits()`](https://iamstein.github.io/synpmx/reference/pmx_pca_visits.md)
and
[`pmx_pca_components()`](https://iamstein.github.io/synpmx/reference/pmx_pca_components.md),
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
[`pmx_pca_report()`](https://iamstein.github.io/synpmx/reference/pmx_pca_report.md).

## Examples

``` r
data <- pmx_simulated_fixture(60)
roles <- pmx_roles(
  id = "ID", time = "TIME", nominal_time = "NTIME", dv = "DV", amt = "AMT",
  evid = "EVID", cmt = "CMT", dvid = "DVID", mdv = "MDV"
)
model <- synpmx_pca_summarize(data, roles)
model
#> A synpmx PCA model
#> 
#>   fitted on    60 patients, 1 arm(s): all (60) 
#>   endpoints    cp (8 visits modelled), pd (6 visits modelled) 
#>   covariates   none 
#>   components   1 (100% of variance) 
#>   dose term    factor 
#>   dosing       2 dose(s) per arm | shared by 2%-2% of each arm 
#> 
#> Generation reads this object and nothing else. To look inside it:
#>   pmx_pca_report(model)      what it read out of the source data
#>   pmx_pca_dosing(model)      the dose schedule each arm shares
#>   pmx_pca_visits(model)      the probability of a visit, per arm
#>   pmx_pca_components(model)  the loadings, over time
pmx_pca_report(model)
#> What the PCA fit read out of the source data
#> 
#>   subjects: 60  components retained: 1 
#> 
#>             quantity                                               what numbers
#>           visit grid               Nominal times modelled, per endpoint      14
#>      feature centers               Mean of each grid cell and covariate      14
#>       feature scales                     Standard deviation of the same      14
#>             loadings                 Component loadings on each feature      14
#>          score means                         Mean score vector, per arm       1
#>     score covariance             Residual covariance between components       1
#>  endpoint transforms                      Log or identity, per endpoint       2
#>         dosing model             Dose times and amounts each arm shares       4
#>          visit model Probability of a visit, per arm, endpoint and time      14
#>        arm constants         Strata and kept columns, one value per arm       0
#>  min_patients
#>            60
#>            60
#>            60
#>            60
#>            60
#>            60
#>            60
#>            60
#>            60
#>            60
```
