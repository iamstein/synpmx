# A simulated study dosed intravenously, subcutaneously, and both

Ninety patients in three arms of thirty: intravenous only, subcutaneous
only, and an intravenous loading dose followed by subcutaneous
maintenance. The third arm is the reason the dataset exists — a patient
who receives both routes is what a mixed model is for, and no public
dataset in the [package
survey](https://iamstein.github.io/synpmx/articles/avatar-public-data-examples.html)
carries an administration column at all.

## Usage

``` r
mixroute_sim
```

## Format

A data frame with 1800 rows and 11 columns:

- ID:

  Patient identifier, 1 to 90.

- TIME:

  Recorded time in days, drifting off the plan the way visits do.

- NTIME:

  Nominal (protocol) time in days.

- DV:

  Concentration in mg/L; `NA` on dose records.

- AMT:

  Dose in mg, 100 at each of days 0, 7 and 14.

- EVID:

  1 on a dose record, 0 on an observation.

- CMT:

  1 for a dose entering the depot, 2 for central and for every
  observation.

- ADM:

  Administration identifier: 1 intravenous, 2 subcutaneous. Declare it
  with
  `pmx_roles(adm = "ADM", routes = c("1" = "iv", "2" = "extravascular"))`.

- CENS:

  1 where the value is at the quantification limit.

- ARM:

  `"IV only"`, `"SC only"` or `"IV then SC"`.

- WT:

  Baseline weight in kg.

## Source

Simulated by `scripts/build-example-data.R`. Not real patient data.

## Details

Simulated from a one-compartment model with first-order absorption, so
the truth is known: clearance 2 L/day, volume 10 L, absorption 0.5 /day,
and bioavailability 0.7, with clearance and volume scaled allometrically
on weight. Bioavailability is identifiable here precisely because the
same study doses both ways.

Three weekly doses, sampled richly after the first and last and at
troughs in between. Concentrations below 0.05 mg/L are reported at the
limit with `CENS` set, which is the convention `nlmixr2` and Monolix
share.

## See also

[onc_sim](https://iamstein.github.io/synpmx/reference/onc_sim.md)
