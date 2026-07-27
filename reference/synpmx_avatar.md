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
  leaves the event template's times unchanged.

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
roles <- pmx_roles("ID", "TIME", "DV", "AMT", "EVID", "CMT", NULL,
                   NULL, NULL, "WT")
synthetic <- synpmx_avatar(source, roles, n_subjects = 2, seed = 123)
#> SYNPMX ALERT: the source has 3 subjects, so every avatar is blended from at most 2 real patients -- fewer than the floor of 5. This markedly raises re-identifiability; use a larger source or treat the output as individually identifying.
#> Warning: the source has 3 subjects, so every avatar is blended from at most 2 real patients -- fewer than the floor of 5. This markedly raises re-identifiability; use a larger source or treat the output as individually identifying.
#> SYNPMX ALERT: 3 subjects in 1 route arm below the donor floor of 5: 1:1:bolus (n=3). Donors are never blended across routes, so these subjects have no legal donor set. Dropping every arm would leave nothing to generate, so generation proceeded as if `on_donor_shortfall = "noise"`. Treat the output as individually identifying.
#> Warning: 3 subjects in 1 route arm below the donor floor of 5: 1:1:bolus (n=3). Donors are never blended across routes, so these subjects have no legal donor set. Dropping every arm would leave nothing to generate, so generation proceeded as if `on_donor_shortfall = "noise"`. Treat the output as individually identifying.
#> Warning: Synthetic generation used documented small-group/profile fallbacks:
#> - Fewer than 5 same-schedule donors were available for at least one subject; the nearest donors from other dose/schedule groups on the same route were borrowed to reach the floor, so some measurements are blended across doses.
validate_pmx(synthetic, roles)$valid
#> [1] TRUE
```
