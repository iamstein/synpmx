# Version 2 methods-vignette maintenance notes

The current methods-vignette requirements are part of
`design/PROTOTYPE_SPEC.md`; that file is the source of truth. The implemented
document is `vignettes/pmxSynthData-method.Rmd`, titled **“How pmxSynthData
Works”**.

Future edits must describe the Version 2 fit-once/generate-many private
population generator exactly as implemented. In particular:

- confidential data may be read only by `fit_private_pmx()`;
- subject contribution bounding precedes every source-dependent aggregate;
- every released source-dependent value must pass through the validated DP
  adapter and appear in composed accounting;
- `generate_pmx()` is post-processing of the released model and public inputs;
- source rows, identifiers, event skeletons, residuals, and unnoised aggregates
  must not be serialized; and
- the document must distinguish a mathematical DP guarantee from legal
  anonymity, release authorization, and scientific fidelity.

Do not reintroduce the retired Version 1 synthesis algorithm. The mechanism
inventory and proof argument are maintained in `design/PRIVACY_ARGUMENT.md`.

