# Version 2 mechanism-level privacy argument

Status: implementation argument for review. It has **not** received independent
differential-privacy, legal, information-security, or data-governance approval.

## Claim and privacy unit

For production fits made with the validated OpenDP adapter, the released fitted
model is subject-level, add-or-remove `(epsilon, delta)` differentially private.
The privacy unit is one person's complete bounded longitudinal contribution:
all rows, visits, actual and nominal times, occasions, doses, infusions,
covariates, endpoint values, missingness, DVIDs, and censoring states.

For neighboring datasets `D` and `D'` differing by one complete subject and any
set of possible model releases `S`, the target guarantee is

```text
Pr[M(D) in S] <= exp(epsilon) Pr[M(D') in S] + delta.
```

The implementation currently uses pure-DP Laplace releases, so realized
`delta` is zero. A positive requested delta is accepted only with an explicit
justification; it is unspent and does not alter the mechanisms.

The noiseless `backend = "public"` path is excluded from this claim. It requires
`public_source = TRUE`, records `formal_dp = FALSE`, and exists only so public
fixtures and demonstrations can test structure without pretending to provide
privacy.

## Confidential-data access inventory

`fit_private_pmx()` is the only exported synthesis function that accepts source
data. Inside the fit, source data are accessed for the following bounded
representations:

| Internal component | Source access | Released form |
|---|---|---|
| subject count | number of unique subject IDs | OpenDP scalar Laplace release, sensitivity 1 |
| event/regimen | bounded dose/event counts, interval, amount, rate, infusion and duration features; optionally a declared public subject-property stratum | OpenDP vector Laplace release |
| endpoint timing | per-subject fixed-basis presence indicators | OpenDP vector Laplace release |
| trajectories | per-subject fixed-basis presence and bounded transformed cell means | OpenDP vector Laplace release |
| covariates | bounded first/second moments or indicators over public categories | OpenDP vector Laplace release |
| censoring | bounded endpoint-level state frequencies and boundaries | OpenDP vector Laplace release |

Input schema checks, role checks, identifier/date rejection, and structural
validation also read the source inside the restricted fit. They must not be
exposed as public diagnostics or logs. They either validate public assertions
or stop the restricted process; they are not serialized into the fitted model.

`generate_pmx()`, `validate_private_model()`, `privacy_report()`, and generated-
data validation do not accept or consult source data. `compare_pmx()` does read
source data, but every source-derived component is explicitly marked
`restricted_not_releasable`; it is outside the model's DP release unless a
separate approved mechanism is added.

## Public inputs and deterministic bounding

The proof treats these caller declarations as established independently of
confidential values: column roles and exclusions; schema, classes, and category
levels (including every subject-property domain); endpoint names, DVIDs,
clocks, transforms, units, generic fixed-grid
bases (or independently public user-supplied bases), and
compartments; numeric bounds; contribution limits; approved protocol values;
generator variability; requested output size; and budget fractions.

For each complete subject, fitting deterministically:

1. retains at most `max_rows` rows;
2. retains at most `max_doses` dose starts;
3. retains at most `max_occasions` occasions;
4. retains at most the declared observations per endpoint and occasion;
5. uses at most `max_timing_cells` public cells; and
6. clips time, amount, rate, numeric covariates, DV, and censoring limits to
   declared public domains before forming features.

The implementation does not release who was truncated or clipped.

## Sensitivity and OpenDP mechanism

Each private summary group represents subject `i` as a fixed vector
`q_i` in `[0, 1]^p`, where `p` is the group's recorded number of dimensions.
The vector sum is

```text
s(D) = sum_i q_i.
```

Adding or removing one bounded subject changes each coordinate by at most one,
so the add-or-remove L1 sensitivity is at most `p`. The adapter constructs an
OpenDP vector domain over finite `f64` atoms with the L1 metric and invokes
`make_laplace(scale = p / epsilon_group)`. Scalar count uses an `f64` atom
domain, absolute distance, sensitivity 1, and scale `1 / epsilon_count`.

The property-conditioned event release uses a tighter sparse bound. For `K`
declared public strata it releases `10K` coordinates: one count indicator and
nine bounded event features per stratum. One subject belongs to exactly one
stratum, so only one 10-coordinate block can change and its L1 sensitivity is
at most 10, independent of `K`. The accounting entry records dimension `10K`
and sensitivity 10. Decoding a stratum probability/regimen, sampling it, and
reconstructing an occasion-assigned dose from generated positive AMT are
post-processing.

Before invocation, the adapter asks OpenDP's privacy map to verify the requested
`d_in`/`d_out` relation. It rejects nonfinite inputs or output and fails closed
when OpenDP is unavailable. It never substitutes an R implementation of
privacy noise. OpenDP represents floating inputs exactly as rationals and uses
discrete internal sampling; the exact installed version/platform remains an
independent-review item.

The privacy-map check uses the exact group epsilon with no acceptance slack.
If floating multiplication would make valid budget fractions sum microscopically
above the requested epsilon, the implementation reduces the largest group by
that overshoot before constructing any mechanism.

Noisy counts, means, proportions, variances, curve smoothing, interpolation,
constraint repair, event-table construction, schema casting, validation, and
all later generated datasets are post-processing.

## Budget composition

The caller explicitly allocates fractions of requested epsilon to subject
count, event/regimen, timing, covariates, endpoints, and censoring. Every active
group must receive positive budget. Each group makes at most one vector release
(or one scalar count release) at its recorded group epsilon. Basic sequential
composition gives

```text
epsilon_realized = sum_g epsilon_g <= epsilon_requested
delta_realized   = sum_g 0 = 0 <= delta_requested.
```

The accounting table records query name, mechanism, epsilon, delta,
sensitivity, and dimension. Finalization stops if composition exceeds the
request. Every fit creates one release-ledger entry. Refitting the same source
is a new release whose budget must be composed with all previous releases by
organizational governance. Repeated `generate_pmx()` calls from one fitted
model consume no additional budget.

## Serialization audit boundary

The returned model contains only:

- declared public configuration and a manifest of its asserted status;
- decoded noisy population summaries;
- backend name/version and proof assumptions;
- composed accounting and a release-ledger entry; and
- warnings derived only from released values and public dimensions.

It must not contain source IDs, rows, subject vectors, exact schedules, event
skeletons, raw residuals, unnoised sums, confidential levels or bounds, or
source-dependent caches. `validate_private_model()` checks prohibited payload
names and accounting, and tests serialize and inspect the model structure.
Factor-valued source ID levels and any public-design default for the ID column
are removed before release; generation creates a fresh synthetic-only level set.
Name-based checks are defense in depth, not a substitute for independent code
and binary serialization review.

## Assumptions and unresolved review work

The formal claim depends on all of the following:

1. every declared public input really was established independently of the
   confidential data;
2. all source-dependent releases are present in the mechanism inventory and
   ledger;
3. deterministic contribution bounds and sensitivity calculations are
   correct for complete-subject adjacency;
4. the installed OpenDP implementation and privacy map are correct for the
   exact domains, metrics, scales, and platform;
5. source-dependent errors, warnings, timing, memory, logs, and other side
   channels remain inside the restricted environment;
6. serialization and downstream systems do not expose transient unnoised
   aggregates; and
7. governance composes this fit with every other release from the source.

Before production claims, obtain independent DP-specialist review of this
argument and code, test the canonical mechanisms with the deployed OpenDP
version, audit runtime and serialization, and complete legal, privacy,
information-security, and data-governance review. Empirical attacks may find
bugs but cannot prove the guarantee.
