# A simulated oncology study with crossover, dose reduction and compressed doses

Two hundred patients shaped like the phase 3 RECORD-1 trial of
everolimus in metastatic renal cell carcinoma, following the
tumour-growth model of Stein et al. (2012). Patients are randomised to
everolimus 10 mg daily or to placebo; most placebo patients cross over
to everolimus when their disease progresses, and about a quarter of
everolimus patients reduce to 5 mg, some after a short interruption.

## Usage

``` r
onc_sim
```

## Format

A data frame with 3379 rows and 15 columns:

- ID:

  Patient identifier, 1 to 200.

- TIME:

  Study day.

- NTIME:

  Nominal (protocol) day; scans are planned, so it equals `TIME`.

- DV:

  Tumour size in cm or trough concentration in ng/mL, by `NAME`.

- AMT:

  Daily dose in mg: 10, 5, or 0 while off treatment or on placebo.

- EVID:

  1 on a dose record, 0 on an observation.

- CMT:

  1 dosing, 2 concentration, 3 tumour size.

- ADDL:

  Additional doses implied by this record.

- II:

  Interdose interval in days, 1 on every dose record.

- NAME:

  `"SLD"` or `"Everolimus trough"`; `NA` on dose records.

- CENS:

  1 where a concentration is at the quantification limit.

- ARM:

  Randomised arm, `"Everolimus 10 mg"` or `"Placebo"`.

- CROSSOVER:

  `TRUE` where a placebo patient crossed over to everolimus.

- BSLD:

  Baseline tumour size in cm, the covariate the model scales on.

- AGE, SEX:

  Baseline age in years and sex.

## Source

Simulated by `scripts/build-example-data.R`, following Stein A, Wang W,
Carter AA, Chiparus O, Hollaender N, Kim H, Motzer RJ, Sarr C. Dynamic
tumor modeling of the dose-response relationship for everolimus in
metastatic renal cell carcinoma using data from the phase 3 RECORD-1
trial. *BMC Cancer* 2012;12:311.
[doi:10.1186/1471-2407-12-311](https://doi.org/10.1186/1471-2407-12-311)
. Not real patient data.

## Details

It carries three shapes the public-data survey has no other example of.
The sum of longest tumour diameters is a slow endpoint anchored on a
per-patient baseline and measured over a year, where every other
endpoint in the survey is a concentration or a fast pharmacodynamic
signal. The dose changes *within* a patient for reasons that patient's
own data explains. And the daily regimen is written with `ADDL` and `II`
rather than one row per dose — 3379 rows instead of the 71096 the
expanded form would need — which is the compressed dose encoding the
survey states outright that it has no example of.

The tumour model is the paper's model 2 with its published parameters:
\$\$dy/dt = r_i - E\_{dose,i}\\y_i\$\$ where \\E\_{dose}\\ is
\\E\_{10}\\ at 10 mg, \\E_5\\ at 5 mg and zero off treatment, and each
patient's growth rate and drug effect scale with baseline tumour size.
Simulated one-year change from baseline is +139.8% on placebo, +21.6% at
5 mg and −13.8% at 10 mg, against the paper's reported +142.1%, +22.4%
and −15.7%.

Pharmacokinetics do not drive the tumour, exactly as in the paper: the
effect is indexed by the dose in force. Trough concentrations ride
alongside as a second endpoint, simulated from a one-compartment oral
model at everolimus's published disposition, so a pharmacokinetic fit of
them means something.

## See also

[mixroute_sim](https://iamstein.github.io/synpmx/reference/mixroute_sim.md)
