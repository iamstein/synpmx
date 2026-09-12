# How AVATAR, PCA and PMXmodel Are Related

Two of the three generators that read a study perform the same
operation. Both
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
and
[`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md)
return a weighted average of the observed patients, with weights that
sum to one. They differ in which patients carry weight, whether a weight
may be negative, and how many numbers the weights are drawn from.
[`synpmx_model()`](https://iamstein.github.io/synpmx/reference/synpmx_model.md)
does something else: it averages in parameter space rather than profile
space, and that change of coordinates is most of its advantage.

Nothing below is measured on a study and nothing below ranks the
generators for a dataset. The public-data surveys do that, on evidence.
What is here is the mechanism behind what those surveys find.

## Both AVATAR and PCA Average Patients

[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
is explicit about it. A synthetic subject’s trajectory is a
distance-weighted blend of its anchor’s nearest compatible neighbours,
so the output is

    y = sum_i w_i * y_i

over real subjects, with weights that are non-negative and sum to one.

[`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md)
computes the same kind of average without appearing to. Principal
components live in the span of the centred observations — that is what
an eigenvector of a sample covariance is — so each component is itself a
combination of the patients. A reconstruction `mu + sum_c s_c * phi_c`
therefore expands back into a combination of the original profiles, and
the coefficients sum to one.

The identity is checkable, so this chunk checks it. Nine profiles, six
visits, two components, one drawn score vector:

``` r

set.seed(1)
n <- 9
profiles <- matrix(stats::rlnorm(n * 6, log(5), 0.4), nrow = n)
centre <- colMeans(profiles)
centred <- sweep(profiles, 2, centre)
decomposition <- svd(centred)
components <- 2

# The weight each real subject carries, for a given score vector. Each loading
# is a combination of the centred rows with coefficients u[, c] / d[c], so the
# reconstruction unwinds into a combination of the rows themselves.
subject_weights <- function(scores) {
  contribution <- as.vector(
    decomposition$u[, seq_len(components), drop = FALSE] %*%
      (scores / decomposition$d[seq_len(components)])
  )
  (1 - sum(contribution)) / n + contribution
}

scores <- c(1.3, -0.7)
weights <- subject_weights(scores)
reconstruction <- as.vector(
  centre + decomposition$v[, seq_len(components), drop = FALSE] %*% scores
)

c(weights_sum_to = sum(weights),
  max_difference = max(abs(as.vector(t(profiles) %*% weights) -
                             reconstruction)))
#> weights_sum_to max_difference 
#>   1.000000e+00   8.881784e-16
```

Generating from a principal-component model re-weights the study rather
than departing from it.

## How the Weights Differ

Three properties of the weights separate the two generators.

|  | [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md) | [`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md) |
|----|----|----|
| Which patients carry weight | The anchor’s nearest neighbours | Every patient |
| Sign of a weight | Non-negative | Signed |
| How many numbers set the weights | One per donor | One per component |

The sign decides most of the behaviour. Non-negative weights summing to
one make a blend a **convex** combination, and a convex combination
cannot leave the convex hull of its donors.
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
can interpolate between real patients and can do nothing else. Nothing
extreme is invented, the spread of the output is narrower than the
study’s, and a relationship measured on the output is flatter than the
one measured on the source.

A signed combination carries no such restriction. Weights go negative on
ordinary draws, not just extreme ones:

``` r

score_sd <- apply(
  decomposition$u[, seq_len(components), drop = FALSE] %*%
    diag(decomposition$d[seq_len(components)]),
  2, stats::sd
)
do.call(rbind, lapply(c(0.5, 1, 2, 3), function(k) {
  weights <- subject_weights(k * score_sd)
  data.frame(draw_at_sd = k, negative_weights = sum(weights < 0),
             smallest_weight = round(min(weights), 3),
             weights_sum = round(sum(weights), 6))
}))
#>   draw_at_sd negative_weights smallest_weight weights_sum
#> 1        0.5                2          -0.020           1
#> 2        1.0                3          -0.151           1
#> 3        2.0                3          -0.414           1
#> 4        3.0                3          -0.676           1
```

A draw one standard deviation from the centre already subtracts several
patients. That is what lets a principal-component model produce a
profile outside everything observed — useful when the study
under-samples its own tails, and the reason the output can also violate
a constraint every real profile satisfies, such as staying positive.

## Curvature in Profile Space

Concentration-time profiles do not lie on a flat surface. Change
clearance and a profile moves along a curve, not a line, so the average
of a fast profile and a slow one is not the profile of intermediate
clearance. Average two peaks that occur at different times and the
result is one broad low bump.

That single fact explains a good deal:

- A global linear basis has to cover a curved surface with a flat one,
  so it needs either many components or it smooths. Visit-to-visit
  texture is the first thing to go.
- Local blending mostly evades the problem, because a curve is locally
  straight and near neighbours are already roughly aligned.
- Treating a patient as a fixed-length vector requires a common
  coordinate system, which is why
  [`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md)
  needs `nominal_time` and refuses without it.
  [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
  needs only a distance, so it can derive a grid or ride on the event
  skeleton it copies.

The grid requirement is a consequence of the representation rather than
a choice in the implementation.

## Resampling Against Fitting

The weights differ because the randomness enters differently, and the
standard name for each is worth having.

[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
is a **smoothed bootstrap**: resample real subjects, perturb them, and
the generating distribution is a mixture of one kernel per patient.
[`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md)
is a **parametric density fit**: estimate a low-dimensional law for the
scores and draw from it. This is the same dichotomy as estimating a
density by kernels or by fitting a Gaussian, and it carries the same
consequences.

**A smoothed bootstrap assumes nothing about shape.** It reproduces the
empirical distribution convolved with its noise, so skew, multimodality
and nonlinear dependence survive without being declared. A cohort with
fast and slow metabolisers keeps both modes because the method never
pools across the gap between them.

**A parametric fit assumes a shape and buys efficiency for it.**
Gaussian scores through a linear basis give elliptical, unimodal,
symmetric contours in profile space. The same two subpopulations are
smeared into one elongated blob unless the score model is told to expect
them. In exchange, every subject in the study contributes to estimating
a handful of numbers rather than only to its own neighbourhood.

### The trade reverses with cohort size

Shape-adaptivity and statistical efficiency trade against each other,
and the exchange rate depends on how many subjects there are.

At a phase 1 cohort size a kernel method is working against itself. A
neighbourhood holds few donors, so each blend is over a handful of
patients, which is both high-variance and the case where a synthetic
record sits closest to a real one. A parametric fit pools every subject
to estimate its few parameters, which is the better use of a small
sample.

As the cohort grows the ranking swaps. Neighbourhoods fill, local
structure becomes estimable, and assuming an elliptical law stops being
a bargain and starts being an error.

Theory says the parametric approach should be the better bet at the
sizes early-phase pharmacometrics runs, and the resampling approach
should win on large pooled data. Whether that is what happens depends on
how much structure a low-dimensional basis discards, which is a property
of the study rather than of the method. The public-data surveys are
where that gets measured.

## What Each Lets Out of the Study

The same dichotomy describes the disclosure question, without appealing
to any detail of either implementation.

A smoothed bootstrap places a synthetic record near a real one by
construction. How near is a bandwidth question, and as the noise goes to
zero the method reproduces its input. That is why
[`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
carries explicit masking mechanisms and why asking whether a record was
copied is a live question for it.

A parametric draw admits individuals only through the estimated
parameters, so one subject’s influence is spread across every synthetic
record rather than concentrated in a few. At a small cohort that is
weaker than it sounds: a mean profile over twelve patients carries a
twelfth of each, and the basis is itself a function of the data.

Neither bounds the influence of one patient. Bounding it is what
[`synpmx_prior()`](https://iamstein.github.io/synpmx/reference/synpmx_prior.md)
and
[`synpmx_calibrated()`](https://iamstein.github.io/synpmx/reference/synpmx_calibrated.md)
are for, and the argument for them is in [the privacy
article](https://iamstein.github.io/synpmx/articles/synpmx-privacy.html).

## PMXmodel Changes the Coordinates

[`synpmx_model()`](https://iamstein.github.io/synpmx/reference/synpmx_model.md)
also fits a law and draws from it, so it sits on the parametric side of
the previous section. What it changes is where the fitting happens.

Clearance is a number that can be averaged. A concentration-time curve
is not. Drawing a clearance from a lognormal and simulating the profile
it implies gives a curve that is a real curve; drawing a point in a
linear profile basis gives a curve that is a weighted sum of curves,
with the misalignment problem still attached. The structural model is
the change of coordinates that makes averaging safe, and it is the
reason a population model does not need many components to represent
what two parameters already say.

So
[`synpmx_model()`](https://iamstein.github.io/synpmx/reference/synpmx_model.md)
belongs on the same side of the previous section as
[`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md),
running the parametric strategy in coordinates where the parametric
assumption is close to true. That is also why its fitted quantities are
ones a pharmacometrician already reasons about. The cost is that the
structural model has to be right enough: where a principal-component
basis adapts to whatever shape the profiles have, a one-compartment
model asserts the shape and cannot represent a study that does not have
it.

## The Family These Sit In

None of this is particular to this package, and the general versions are
worth knowing when a limitation here becomes binding.

- **Smoothed bootstrap and kernel density estimation** are the general
  form of what
  [`synpmx_avatar()`](https://iamstein.github.io/synpmx/reference/synpmx_avatar.md)
  does. The literature covers bandwidth selection and the behaviour of
  neighbourhoods in high dimension, neither of which this package’s
  blending exposes as a tunable.
- **Functional principal component analysis** is the general form of
  what
  [`synpmx_pca()`](https://iamstein.github.io/synpmx/reference/synpmx_pca.md)
  does, and treats a curve as the object rather than a vector on a grid.
  It handles the alignment problem directly through registration, where
  this package requires a declared grid instead.
- **Archetypal analysis** is the convex-constrained version of a
  component model: components are themselves convex combinations of
  observed points, so the reconstruction cannot subtract a patient. It
  sits exactly between the two generators here and neither implements
  it.
- **Mixtures of probabilistic principal-component analysers** are the
  local version, fitting a separate low-dimensional law per region,
  which recovers multimodality that one global basis smears.
- **Nonlinear mixed-effects modelling** is the general form of what
  [`synpmx_model()`](https://iamstein.github.io/synpmx/reference/synpmx_model.md)
  does, with covariate models, inter-occasion variability and arbitrary
  structural models that the built-in catalogue here does not attempt.

What this package covers is the narrow case: a small catalogue of each,
usable without tuning, on an event table declared through
[`pmx_roles()`](https://iamstein.github.io/synpmx/reference/pmx_roles.md).
Everything above is where to go when that stops being enough.

## Where to go next

- [The synpmx data generation
  algorithms](https://iamstein.github.io/synpmx/articles/synpmx-methods.html)
  — the five modes side by side on one study, and which to use when.
- [`vignette("avatar-algorithm")`](https://iamstein.github.io/synpmx/articles/avatar-algorithm.md),
  [`vignette("pca-algorithm")`](https://iamstein.github.io/synpmx/articles/pca-algorithm.md)
  and
  [`vignette("pmxmodel-algorithm")`](https://iamstein.github.io/synpmx/articles/pmxmodel-algorithm.md)
  — each mechanism step by step.
- [Evaluating AVATAR on public
  data](https://iamstein.github.io/synpmx/articles/avatar-public-data-examples.html)
  and [Evaluating PCA on public
  data](https://iamstein.github.io/synpmx/articles/pca-public-data-examples.html)
  — what the two do to the same studies, measured.
- [Are relationships
  preserved?](https://iamstein.github.io/synpmx/articles/example-avatar-PKPD-covariate-treatment-effect.html)
  — the compression a convex blend produces, measured against a known
  truth.
