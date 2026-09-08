#' A simulated study dosed intravenously, subcutaneously, and both
#'
#' Ninety patients in three arms of thirty: intravenous only, subcutaneous only,
#' and an intravenous loading dose followed by subcutaneous maintenance. The
#' third arm is the reason the dataset exists — a patient who receives both
#' routes is what a mixed model is for, and no public dataset in the
#' [package survey](https://iamstein.github.io/synpmx/articles/avatar-public-data-examples.html)
#' carries an administration column at all.
#'
#' Simulated from a one-compartment model with first-order absorption, so the
#' truth is known: clearance 2 L/day, volume 10 L, absorption 0.5 /day, and
#' bioavailability 0.7, with clearance and volume scaled allometrically on
#' weight. Bioavailability is identifiable here precisely because the same study
#' doses both ways.
#'
#' Three weekly doses, sampled richly after the first and last and at troughs in
#' between. Concentrations below 0.05 mg/L are reported at the limit with `CENS`
#' set, which is the convention `nlmixr2` and Monolix share.
#'
#' @format A data frame with 1800 rows and 11 columns:
#' \describe{
#'   \item{ID}{Patient identifier, 1 to 90.}
#'   \item{TIME}{Recorded time in days, drifting off the plan the way visits do.}
#'   \item{NTIME}{Nominal (protocol) time in days.}
#'   \item{DV}{Concentration in mg/L; `NA` on dose records.}
#'   \item{AMT}{Dose in mg, 100 at each of days 0, 7 and 14.}
#'   \item{EVID}{1 on a dose record, 0 on an observation.}
#'   \item{CMT}{1 for a dose entering the depot, 2 for central and for every
#'     observation.}
#'   \item{ADM}{Administration identifier: 1 intravenous, 2 subcutaneous.
#'     Declare it with `pmx_roles(adm = "ADM", routes = c("1" = "iv",
#'     "2" = "extravascular"))`.}
#'   \item{CENS}{1 where the value is at the quantification limit.}
#'   \item{ARM}{`"IV only"`, `"SC only"` or `"IV then SC"`.}
#'   \item{WT}{Baseline weight in kg.}
#' }
#' @source Simulated by `scripts/build-example-data.R`. Not real patient data.
#' @seealso [onc_sim]
"mixroute_sim"

#' A simulated oncology study with crossover, dose reduction and compressed doses
#'
#' Two hundred patients shaped like the phase 3 RECORD-1 trial of everolimus in
#' metastatic renal cell carcinoma, following the tumour-growth model of Stein
#' et al. (2012). Patients are randomised to everolimus 10 mg daily or to
#' placebo; most placebo patients cross over to everolimus when their disease
#' progresses, and about a quarter of everolimus patients reduce to 5 mg, some
#' after a short interruption.
#'
#' It carries three shapes the public-data survey has no other example of. The
#' sum of longest tumour diameters is a slow endpoint anchored on a per-patient
#' baseline and measured over a year, where every other endpoint in the survey
#' is a concentration or a fast pharmacodynamic signal. The dose changes *within*
#' a patient for reasons that patient's own data explains. And the daily regimen
#' is written with `ADDL` and `II` rather than one row per dose — 3379 rows
#' instead of the 71096 the expanded form would need — which is the compressed
#' dose encoding the survey states outright that it has no example of.
#'
#' The tumour model is the paper's model 2 with its published parameters:
#' \deqn{dy/dt = r_i - E_{dose,i}\,y_i}
#' where \eqn{E_{dose}} is \eqn{E_{10}} at 10 mg, \eqn{E_5} at 5 mg and zero off
#' treatment, and each patient's growth rate and drug effect scale with baseline
#' tumour size. Simulated one-year change from baseline is +139.8% on placebo,
#' +21.6% at 5 mg and −13.8% at 10 mg, against the paper's reported +142.1%,
#' +22.4% and −15.7%.
#'
#' Pharmacokinetics do not drive the tumour, exactly as in the paper: the effect
#' is indexed by the dose in force. Trough concentrations ride alongside as a
#' second endpoint, simulated from a one-compartment oral model at everolimus's
#' published disposition, so a pharmacokinetic fit of them means something.
#'
#' @format A data frame with 3379 rows and 15 columns:
#' \describe{
#'   \item{ID}{Patient identifier, 1 to 200.}
#'   \item{TIME}{Study day.}
#'   \item{NTIME}{Nominal (protocol) day; scans are planned, so it equals `TIME`.}
#'   \item{DV}{Tumour size in cm or trough concentration in ng/mL, by `NAME`.}
#'   \item{AMT}{Daily dose in mg: 10, 5, or 0 while off treatment or on placebo.}
#'   \item{EVID}{1 on a dose record, 0 on an observation.}
#'   \item{CMT}{1 dosing, 2 concentration, 3 tumour size.}
#'   \item{ADDL}{Additional doses implied by this record.}
#'   \item{II}{Interdose interval in days, 1 on every dose record.}
#'   \item{NAME}{`"SLD"` or `"Everolimus trough"`; `NA` on dose records.}
#'   \item{CENS}{1 where a concentration is at the quantification limit.}
#'   \item{ARM}{Randomised arm, `"Everolimus 10 mg"` or `"Placebo"`.}
#'   \item{CROSSOVER}{`TRUE` where a placebo patient crossed over to everolimus.}
#'   \item{BSLD}{Baseline tumour size in cm, the covariate the model scales on.}
#'   \item{AGE, SEX}{Baseline age in years and sex.}
#' }
#' @source Simulated by `scripts/build-example-data.R`, following Stein A,
#'   Wang W, Carter AA, Chiparus O, Hollaender N, Kim H, Motzer RJ, Sarr C.
#'   Dynamic tumor modeling of the dose-response relationship for everolimus in
#'   metastatic renal cell carcinoma using data from the phase 3 RECORD-1 trial.
#'   *BMC Cancer* 2012;12:311. \doi{10.1186/1471-2407-12-311}. Not real patient
#'   data.
#' @seealso [mixroute_sim]
"onc_sim"
