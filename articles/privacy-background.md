# Privacy background: where d, f, and the error law come from

Tutorial background for the quantities used throughout the privacy
documentation: the release dimension `d`, the noise-to-signal ratio `f`,
and the error law relating them. Built from first principles, with the
arithmetic worked out.

Related documents, so you land in the right one:

- **This article** — intuition and derivation. How to think about the
  tradeoff.
- [Mechanism-level privacy
  argument](https://iamstein.github.io/synpmx/articles/privacy-argument.md)
  — the formal argument for the implemented engine, written for a
  reviewer.
- [Feasibility by cohort
  size](https://iamstein.github.io/synpmx/articles/feasibility.md) — the
  measured evidence for what these quantities cost in practice.
- [`vignette("synpmx-privacy")`](https://iamstein.github.io/synpmx/articles/synpmx-privacy.md)
  — the decision rule: whether your release needs formal privacy at all.

------------------------------------------------------------------------

## 1. The mechanism, in one line

To release a number computed from patient data, add random noise. Make
the noise big enough that one person’s presence or absence is hidden
inside it.

Everything else is bookkeeping about **how big is big enough**.

------------------------------------------------------------------------

## 2. Sensitivity: how much can one person move the answer?

Call the released quantity `q(D)` for dataset `D`. **Sensitivity** is
the most one person can change it:

``` math
\Delta \;=\; \max_{D, D'} \bigl| q(D) - q(D') \bigr|
```

over all pairs `D`, `D'` differing by adding or removing one complete
subject.

Two examples:

- **A count.** “How many subjects?” Adding one person changes it by
  exactly 1, so `Δ = 1`.
- **A sum of per-subject values.** “Sum of everyone’s clearance.” One
  person could have an enormous CL, so `Δ` is unbounded — *unless we
  clip*.

Clipping is why the package insists on public ranges. If every subject’s
value is forced into `[0, 1]` before summing, then adding or removing a
person changes the sum by at most 1, so `Δ = 1`. **The clipping range
creates the sensitivity bound; without it there is no bound and no
mechanism.**

This is also why the range must be chosen without looking at the data. A
range derived from the data is itself a leak, and one the accounting
will not catch.

------------------------------------------------------------------------

## 3. The Laplace mechanism: noise scale

To release `q(D)` with `epsilon`-differential privacy, add noise drawn
from a Laplace distribution with scale

``` math
b \;=\; \frac{\Delta}{\varepsilon}
```

Larger sensitivity means more noise. Smaller epsilon — a stronger
promise — means more noise. That single fraction is the whole
privacy-utility tradeoff.

Useful facts about `Laplace(b)`:

| Quantity                  | Value               |
|---------------------------|---------------------|
| Mean absolute error       | `b`                 |
| **Median** absolute error | `b · ln 2 ≈ 0.69 b` |
| Standard deviation        | `b · √2 ≈ 1.41 b`   |

The 0.69 factor is why measured median errors come in slightly better
than the planning formula predicts. The formula uses `b`; the tables
report medians.

------------------------------------------------------------------------

## 4. `d`: releasing more than one number

**`d` is simply how many separate quantities you release.**

Privacy budget is a finite resource that is *spent*. Releasing `d`
quantities under basic sequential composition means splitting epsilon
`d` ways — each release gets `epsilon/d`, and the total loss is still
`epsilon`.

So for `d` clipped per-subject scalars, each with `Δ = 1`:

``` math
b \;=\; \frac{\Delta}{\varepsilon/d} \;=\; \frac{d}{\varepsilon}
```

**Releasing twice as many numbers doubles the noise on every one of
them.** This is the single most important consequence to internalize,
and it is why the Version 3 design works so hard to shrink `d` from
about 106 numbers to 3.

(Equivalently: release all `d` as one vector with L1 sensitivity
`Δ₁ = d`, since one subject could sit at 1 in every coordinate. Same
answer.)

------------------------------------------------------------------------

## 5. From sums to means, and the arrival of `N`

We never actually want a sum. We want a mean — typical clearance,
typical baseline. So we release the sum, release the count, and divide:

``` math
\widehat{\text{mean}} \;=\; \frac{\widetilde{\sum_i x_i}}{\widetilde{N}}
```

Dividing is **post-processing** and costs nothing extra: once a quantity
is differentially private, any function of it that does not revisit the
source is also differentially private.

The noise on the sum is `b`. Dividing by `N` divides the error by `N`
too:

``` math
\text{error on the mean} \;\approx\; \frac{b}{N} \;=\; \frac{d}{\varepsilon N}
```

**This is `f`.**

``` math
\boxed{\;f \;=\; \frac{d}{\varepsilon N}\;}
```

The count in the denominator is itself noisy, which adds a second error
term. When `N` is comfortably larger than `b` that term is second-order
and the approximation above holds. At very small `N` it is not
negligible, which is part of why the small-cohort measurements are
messier than the formula.

------------------------------------------------------------------------

## 6. What `f` actually means

Here is the part worth being precise about, because it is what makes `f`
interpretable.

Every subject’s value was **scaled into `[0, 1]`** before summing — that
is what clipping to the public range does. So the released mean is also
on a `[0, 1]` scale, where 0 is the bottom of the prior range and 1 is
the top.

Therefore:

> **`f` is the error expressed as a fraction of the prior range’s
> width.**

That gives it a natural reading:

- `f = 1.0` — the noise is as wide as the entire prior. The release has
  told you nothing you did not already assume. **You spent privacy
  budget for nothing.**
- `f = 0.5` — the noise covers half the prior. You have roughly halved
  your uncertainty.
- `f = 0.1` — the noise covers a tenth. The prior is now contributing
  little; the data is doing the work.

This is why the decision rule is “does the release beat the prior?”
rather than “is the error small?” `f` measures exactly that, and it is
dimensionless.

------------------------------------------------------------------------

## 7. Converting `f` to a fold-error

PK parameters are worked on the log scale, so the prior range is
expressed as a span in log units. If the prior spans a factor of `k`:

``` math
S \;=\; \ln k
```

and the error in log units is `f · S`, so the multiplicative (fold)
error is

``` math
\text{fold-error} \;=\; \exp(f \cdot S) \;=\; \exp\!\left(\frac{d}{\varepsilon N} \ln k\right)
```

| Prior                                | `k` | `S = ln k` |
|--------------------------------------|----:|-----------:|
| 2-fold                               |   2 |       0.69 |
| 5-fold                               |   5 |       1.61 |
| **8-fold (correction factor)**       |   8 |   **2.08** |
| 100-fold (absolute CL, new compound) | 100 |       4.61 |

------------------------------------------------------------------------

## 8. Worked example

**Setup.** 20 subjects. Release 3 quantities: cohort size, a PK
correction factor, a PD correction factor. Approved epsilon 0.5. The PK
correction is clipped to `[1/4, 4]`, an 8-fold prior.

| Step                    | Computation                 | Result           |
|-------------------------|-----------------------------|------------------|
| Budget per release      | `0.5 / 3`                   | `ε = 0.167` each |
| Sensitivity per release | one clipped subject         | `Δ = 1`          |
| Noise scale on the sum  | `b = Δ/ε = 1/0.167`         | `b = 6`          |
| Error on the mean       | `b/N = 6/20`                | `f = 0.30`       |
| Cross-check             | `f = d/(εN) = 3/(0.5 × 20)` | `f = 0.30` ✓     |
| Prior span              | `ln 8`                      | `S = 2.08`       |
| Error in log units      | `0.30 × 2.08`               | `0.62`           |
| **Fold-error**          | `exp(0.62)`                 | **1.87-fold**    |

**Reading it.** Exposures land within roughly a factor of two.
`f = 0.30` means the release cut prior uncertainty by about 70%, so it
clearly beat the prior and the budget was well spent. Against the
accuracy bar in the package specification, this passes comfortably.

**Now change one thing at a time:**

| Change | New `f` | Fold-error | Comment |
|----|---:|---:|----|
| Baseline above | 0.30 | 1.87 |  |
| Release 6 quantities instead of 3 | 0.60 | 3.5 | Doubling `d` doubles the error |
| Epsilon 1.0 instead of 0.5 | 0.15 | 1.37 | Twice the privacy cost for accuracy not needed |
| 60 subjects instead of 20 | 0.10 | 1.23 | `N` is the cheapest lever, if you have it |
| 100-fold prior instead of 8-fold | 0.30 | **4.0** | Same `f`, far worse outcome |

The last row is the one that repays study. **`f` did not change**,
because `f` knows nothing about the prior’s width — it is a fraction
*of* that width. A wide prior and a narrow prior with identical `f` give
wildly different answers. This is why the correction-factor
parameterization matters so much: it does not improve `f`, it shrinks
`S`.

------------------------------------------------------------------------

## 9. The four levers

| Lever | Effect | Cost |
|----|----|----|
| **Narrow the prior** (`S`) | Direct, exponential in the fold-error | Free if from public knowledge; otherwise a budget slice for private range-finding |
| **Release fewer things** (`d`) | Proportional | Give up quantities you wanted |
| **More subjects** (`N`) | Proportional | Not usually available |
| **Larger epsilon** | Proportional | A weaker guarantee — the thing we are trying to protect |

They are listed in the order to try them. Epsilon is last on purpose:
raising it is always available and is almost always the wrong answer.
Anything achievable by raising epsilon is usually also achievable by
tightening the prior, which costs nothing.

------------------------------------------------------------------------

## 10. Caveats and refinements

**Composition is pessimistic.** Basic sequential composition simply adds
up epsilon. Advanced composition and zCDP do better, scaling roughly
with `√d` instead of `d`. But those gains only appear at large `d` — see
[the feasibility
article](https://iamstein.github.io/synpmx/articles/feasibility.md)
section 8. At `d = 3` pure-DP Laplace is genuinely the better mechanism,
so the simple `d/(εN)` law is the right one here.

**Gaussian is not automatically better.** Gaussian noise trades an L2
sensitivity of `√d` against a constant of roughly `√(2 ln(1.25/δ))`,
which is about 5 at `δ = 1e-6`. It wins only once `√d · 5 < d`, meaning
`d > 25`.

**The formula is a planning bound, not a prediction.** It uses the
Laplace scale, while measured tables report medians, which run about
0.69 times smaller. Treat `f` as conservative. Measured confirmation is
in [the feasibility
article](https://iamstein.github.io/synpmx/articles/feasibility.md)
section 8: predicted 1.49-fold, measured 1.52.

**Clipping bias is separate from noise.** If the true value lies outside
the prior range, clipping pulls the release toward the boundary and no
amount of epsilon fixes it. This is why a correction factor pressed
against its clipping boundary must be reported as a diagnostic — it is
the signal that the prior was wrong, and it is not visible in `f`.

**Sensitivity depends on the adjacency definition.** All of the above
assumes add-or-remove one *complete subject*, and that each person
appears exactly once. If a subject can appear in more than one pooled
study, sensitivity scales with the number of appearances and the
arithmetic above must be redone.
