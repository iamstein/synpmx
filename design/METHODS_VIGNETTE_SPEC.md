# Version 2 vignette maintenance notes

The current vignette requirements are part of `design/PROTOTYPE_SPEC.md`; that
file is the source of truth. Version 2 has three documents with deliberately
separate jobs:

- `vignettes/pmxSynthData-demo.Rmd`: runnable use and public examples;
- `vignettes/pmxSynthData-privacy-intro.Rmd`: beginner-to-technical privacy
  guarantee and assumptions; and
- `vignettes/pmxSynthData-simulation-method.Rmd`: the implemented patient
  simulation algorithm.

Future edits must keep those purposes separate and describe the Version 2
fit-once/generate-many private population generator exactly as implemented. In
particular:

- confidential data may be read only by `fit_private_pmx()`;
- subject contribution bounding precedes every source-dependent aggregate;
- every released source-dependent value must pass through the validated DP
  adapter and appear in composed accounting;
- `generate_pmx()` is post-processing of the released model and public inputs;
- source rows, identifiers, event skeletons, residuals, and unnoised aggregates
  must not be serialized; and
- the document must distinguish a mathematical DP guarantee from legal
  anonymity, release authorization, and scientific fidelity.

The simulation-method vignette must state that trajectories use fixed-grid
cell summaries, a one-pass 1--2--1 smoother, and piecewise-linear interpolation
rather than splines, ODEs, or an NLME likelihood. Sampling-time generation must
be described separately from endpoint-value generation. Its recommended and
demonstrated workflow must omit public dose regimens and sampling schedules:
the private fit learns dose count, interval, amount, infusion behavior,
occasion activation, conditional sample count, and timing-cell occupancy from
the input. A source-independent automatic grid is a discretization basis, not
a disclosed visit schedule. Public regimen or occasion-schedule overrides are
exceptional and must be described as appropriate only when independently
public.

Do not reintroduce the retired Version 1 synthesis algorithm. The mechanism
inventory and proof argument are maintained in `design/PRIVACY_ARGUMENT.md`.
