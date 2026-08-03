# Synthesize a structurally faithful PMX dataset (AVATAR-style)

Samples complete subject event templates and fills them with
AVATAR-like, endpoint-specific blends of compatible subjects' baseline
covariates and longitudinal measurements. Event-control fields such as
EVID, AMT, RATE, CMT, and DVID are never averaged or independently
generated.

## Usage

``` r
synpmx_avatar(
  data,
  roles,
  n_subjects = NULL,
  seed = 123,
  event_method = "template",
  dv_method = "avatar_blend",
  k = 5,
  pca_variance = 0.9,
  subject_noise_sd = 0.15,
  residual_noise_sd = 0.05,
  residual_phi = 0.6,
  time_jitter = 0,
  screen = TRUE,
  coarsen_time = TRUE,
  min_pattern_share = 2L,
  max_donor_weight = 0.5,
  on_donor_shortfall = c("drop", "noise", "error")
)
```

## Arguments

- data:

  A source PMX data frame or tibble.

- roles:

  Explicit roles from
  [`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md).
  Columns listed in `roles$exclude` are omitted from the generated
  output.

- n_subjects:

  Number of synthetic subjects. `NULL` retains the source count.

- seed:

  Reproducibility seed. The caller's random-number state is restored on
  exit.

- event_method:

  Event generation method. The prototype supports `"template"`.

- dv_method:

  Measurement method. The prototype supports `"avatar_blend"`.

- k:

  Number of real patients blended into each synthetic subject (default
  5). Same-schedule donors are used first; when a subject's
  dose/schedule group holds fewer than `k`, the nearest subjects from
  other groups *on the same administration route* are borrowed to reach
  `k`, blending measurements across doses. Route is never crossed, so a
  route arm holding fewer than `k + 1` subjects cannot reach the floor
  at all; those subjects are dropped from the anchor pool with a loud
  alert, and the synthetic cohort does not represent that arm.

- pca_variance:

  Fraction of usable profile variance retained for neighborhood
  distances.

- subject_noise_sd:

  Nonnegative subject perturbation multiplier.

- residual_noise_sd:

  Nonnegative within-trajectory noise multiplier.

- residual_phi:

  AR(1) correlation in observation order, strictly between -1 and 1.

- time_jitter:

  Standard deviation for coherent tied-time jitter. Zero, the default,
  leaves the event template's times unchanged. This is a realism
  control, **not** a privacy control: every jittered time is clamped
  inside its own Voronoi cell, so no value of `time_jitter` moves a
  visit more than half a gap from where the source subject's visit was,
  and the source schedule stays recoverable. Use `coarsen_time` for
  that.

- screen:

  When `TRUE` (default), a source subject whose follow-up length or dose
  count is more than twice the cohort's 90th percentile is not used as
  an anchor, so no avatar inherits an extreme skeleton (the long tail a
  reader notices). Anchoring the cut on the 90th percentile, not the
  median, means ordinary spread is kept: only a subject well beyond the
  high end of normal is dropped. Only these structural axes are
  screened; dose magnitude (which weight-based dosing makes noisy) and
  DV (which is blended, not copied) are not. A source with no extreme
  subject is unaffected. Set `FALSE` to anchor on every subject. For a
  fuller, tunable screen of the generated output, see
  [`flag_identifiable_subjects()`](https://iamstein.github.io/synpmx/reference/flag_identifiable_subjects.md)
  and
  [`remediate_identifiable_subjects()`](https://iamstein.github.io/synpmx/reference/remediate_identifiable_subjects.md).

- coarsen_time:

  When `TRUE` (default), source times are collapsed onto a shared visit
  grid before generation, and per-visit deviations are pooled across the
  cohort and resampled independently onto each avatar. The grid is the
  `nominal_time` role where one is declared, and K-means centres of the
  pooled times otherwise. This is the mechanism that stops an avatar
  from carrying one real subject's exact visit schedule: the event
  skeleton is copied verbatim from a single anchor, and under actual
  recorded times almost every subject is alone in its event-signature
  class, so the copy is identifying. Snapping is many-to-one and
  *destroys* the deviation rather than perturbing it, which is what
  distinguishes this from `time_jitter`. A source already on nominal
  time has no deviation to remove or restore, so its output is
  unchanged. Run
  [`skeleton_uniqueness()`](https://iamstein.github.io/synpmx/reference/skeleton_uniqueness.md)
  on the source to see how much this has to do, and what it leaves
  behind. The cost is timing fidelity: an avatar's deviation from
  nominal is drawn from the cohort, not inherited, so `TIME` no longer
  pairs with its `DV` as precisely as the source did. Set `FALSE` to
  keep exact source timing.

- min_pattern_share:

  How many source subjects must share an attendance pattern before an
  avatar may be given it. Default 2; `1` restores copying the anchor's
  own pattern.

  Two is the smallest value that means something, and what it means is
  precise: **no synthetic patient carries a schedule unique to a real
  patient.** An attacker who links a reproduced pattern to a participant
  gets at least two candidates, never one. Higher values hide more and
  are harder to state — "shared by at least three" is more conservative
  but no more defensible — and they cost sharply more, because
  attendance patterns are distributed with one common pattern and a long
  tail of singletons rather than a populated middle.

  Once `coarsen_time` has put every subject on a shared visit grid, what
  remains of a schedule is which of those visits each subject attended.
  Every time is then shared, so no single visit is identifying — but the
  *combination of absences* is, and a patient who missed weeks 2 and 3
  can be singled out by a fingerprint made of gaps. No grid can fix that
  at any resolution, because the grid decides where the visits are and
  not which ones a subject has. So the pattern is sampled from ones at
  least `min_pattern_share` subjects hold, frequency-weighted, and whole
  patterns are drawn rather than individual visits, which keeps an
  ending monotone instead of producing implausible attend/miss/attend
  sequences.

  Unlike dropping the exposed subjects, nobody leaves the cohort: a
  subject with a rare pattern still contributes measurements as a donor,
  only their distinctive absences stop being reproduced. Dose events are
  never sampled, since that could emit a regimen no protocol permits.
  Raising this hides more and flattens the cohort's missingness further;
  where no pattern is shared widely enough, anchors keep their own and
  the run alerts loudly. Pools are formed within each
  `subject_properties` stratum and endpoint set.

  **Patterns below the floor are lost, not approximated.** An ending or
  dose-interruption pattern held by too few patients simply will not
  appear in the synthetic data, and that loss is the mechanism working —
  it is what stops an avatar carrying a schedule traceable to one
  person. Because it is a real cost to the data's realism, every run
  reports it: the number of source patterns excluded and how many
  subjects held them, both as a loud alert and as `patterns_total`,
  `patterns_dropped` and `subjects_with_dropped_pattern` in the
  settings, alongside `pattern_generated_fraction` for how often an
  arrangement had to be invented. What survives is how much missingness
  there was and of what kind; what is lost is which specific visits.
  Check those figures before deciding the default suits your study.

- max_donor_weight:

  Largest share of one synthetic subject that any one real donor may
  contribute. The default 0.50 states simply that no single real patient
  is more than half of any synthetic patient.

  The floor `k` sets how many patients are blended; this cap is what
  bounds any single patient's contribution, so it, not `k`, is the
  parameter that limits how closely an avatar can resemble one real
  person. Without a cap the randomized weights are strongly concentrated
  — a median 58% of an avatar in one donor at `k = 5`, and about 2.4
  effective donors — so the cap is what makes the floor mean anything.

  Two diagnostics in the returned `pmx_settings` say where a given value
  landed: `mean_effective_donors` is `1 / sum(w^2)`, the number of
  donors an avatar is effectively blended from, and
  `cap_binding_fraction` is how often the cap actually fired. A cap
  binding on nearly every subject is not a guardrail but the weighting
  scheme itself, with the inverse-distance term underneath it doing
  little; one that never fires is not protecting anything. At `k = 5`
  the default binds on roughly two thirds of subjects.

- on_donor_shortfall:

  What to do with a subject whose administration route holds fewer than
  `k + 1` subjects, so that no legal donor set exists for it. `"drop"`
  (default) omits those subjects from the anchor pool: no avatar is
  built on them and the synthetic cohort does not represent that arm.
  `"noise"` keeps them, blending however many same-route donors exist
  (possibly none) and relying on `subject_noise_sd` and
  `residual_noise_sd` for the rest — **not recommended**, because such a
  synthetic subject can remain close to one real patient; screen the
  result with
  [`flag_identifiable_subjects()`](https://iamstein.github.io/synpmx/reference/flag_identifiable_subjects.md)
  if you use it. `"error"` refuses to generate and names the choice.
  Every branch alerts loudly. When *every* route arm is below the floor,
  `"drop"` would leave nothing to generate, so generation proceeds as if
  `"noise"`.

## Value

An ordinary data frame or tibble with retained source columns, order,
and practical classes. A lightweight `pmx_settings` attribute records
the generator choices and endpoint transformations.

## Details

This is an AVATAR-inspired adaptation, not an exact implementation of
published AVATAR software. It creates synthetic data for model-workflow
exploration. It does not provide formal anonymization or preserve
scientific parameter or covariate-response relationships.

Donors are selected in two stages, both confined to the anchor's own
administration route, which is never crossed: same-signature donors
first, taken nearest-first by Euclidean distance between retained PCA
profile coordinates, then — if that yields fewer than `k` — the nearest
remaining route-compatible subjects regardless of dose or schedule.

For the selected donors, randomized raw weights are
`Exp(1) / max(distance, epsilon) * 2^(-randomized_rank)`. They are
normalized and then capped so that *no* donor exceeds
`max_donor_weight`, the excess being redistributed proportionally among
the donors still below the cap until none is over. A cap below `1/K` for
`K` donors cannot be satisfied and relaxes to `1/K`, i.e. uniform
weights. The same subject weights are used for covariates and all
endpoints; weights are renormalized locally when a donor lacks a
requested endpoint/time value.

Positive-like endpoints use an offset log scale and are constrained to
be nonnegative after back-transformation. Other endpoints use the
identity scale. Transform choices and interpolation alignment are
recorded in the returned `pmx_settings` attribute.

## Dose recomputed from a blended covariate

When the dose is a fixed multiple of a baseline covariate within each
assigned stratum — mg/kg, mg/m^2 — that multiplier is a protocol
property the whole stratum shares, while the covariate is individual and
is already blended across donors. `synpmx_avatar()` detects this and
recomputes each avatar's `AMT` (and any `rate`, keeping the infusion
duration) from the avatar's *own* blended covariate.

This fixes two things at once. The amount is no longer one real
patient's real dose, which under proportional dosing discloses that
patient's weight exactly. And the avatar stops violating the protocol it
claims to follow: previously `AMT` was copied from the anchor while
covariates were blended, so a cohort dosed at exactly 5 mg/kg produced
avatars anywhere from 4.4 to 5.3.

Several dose levels are handled by clustering the observed ratios rather
than averaging within a group, so a 1/2/3 mg/kg escalation is recognised
without the arm being declared — and so is **intra-patient** escalation,
where the level changes within a subject and no subject-constant
grouping could see it.

Detection is conservative and fails closed: the ratios must collapse
onto a handful of levels while the amounts they came from stay varied,
so a study whose dose is unrelated to the covariate produces about as
many ratios as amounts and is refused. The covariate and the levels
found are recorded as `dose_basis` and `dose_levels` in the settings,
`NA` when none was found.

## Examples

``` r
source <- data.frame(
  ID = rep(1:3, each = 4),
  TIME = rep(c(0, 0, 1, 2), 3),
  DV = c(0, 0.2, 2, 1, 0, 0.3, 3, 1.5, 0, 0.4, 4, 2),
  AMT = rep(c(100, 0, 0, 0), 3),
  EVID = rep(c(1L, 0L, 0L, 0L), 3),
  CMT = rep(c(1L, 2L, 2L, 2L), 3),
  WT = rep(c(60, 70, 80), each = 4)
)
roles <- pmx_roles(
  id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
  cmt = "CMT", covariates = "WT"
)
synthetic <- synpmx_avatar(source, roles, n_subjects = 2, seed = 123)
#> synpmx_avatar(): no `dvid` declared, so every observation is treated as one endpoint.
#>   Correct for a single-endpoint study; declare `dvid` if this one has more.
#> SYNPMX ALERT: source too small for the donor floor
#>   the source has 3 patients, so every avatar is blended from at most 2 real
#>   patients -- fewer than the floor of k = 5.
#>   Why it matters: blending across few patients leaves each avatar close to
#>     an individual, which markedly raises re-identifiability.
#>   Fix: use a larger source, or treat the output as individually identifying
#>     and keep it under the source's own access controls.
#> SYNPMX ALERT: every route arm is below the donor floor
#>   3 patients sit in 1 route arm holding fewer than the donor floor of k =
#>   5: 1:1:bolus (n=3).
#>   Why it matters: donors are never blended across routes, so these patients
#>     have no legal donor set. Dropping every arm would leave nothing to
#>     generate, so generation proceeded as if `on_donor_shortfall = "noise"`.
#>   Fix: treat the whole output as individually identifying, or use a larger
#>     source.
#> Warning: Synthetic generation used documented small-group/profile fallbacks:
#> - Fewer than 5 same-schedule donors were available for at least one
#>   subject; the nearest donors from other dose/schedule groups on the same
#>   route were borrowed to reach the floor, so some measurements are blended
#>   across doses.
validate_pmx(synthetic, roles)$valid
#> [1] TRUE
```
