# Instructions for the `pmxSynthData` Methods Vignette

Create a second package vignette at `vignettes/pmxSynthData-method.Rmd` titled
**“How pmxSynthData Works”**.

The vignette should be a rigorous but approachable teaching document for
pharmacometricians who want to understand exactly how the package generates
synthetic mock data. Write in the voice of an experienced PMX professor:
precise, intuitive, candid about assumptions, and able to connect mathematical
details to practical modeling workflows.

Before writing, inspect the current implementation and tests so every
explanation, equation, default, and edge case matches the code. Do not describe
an idealized algorithm that differs from the package.

Begin with a flow chart showing the complete pipeline from source PMX event
data to the final mock dataset. Create it without adding an unnecessary package
dependency.

Explain the method progressively:

1. The problem the package solves and its intended model-development use case.
2. The required input structure and explicit column roles.
3. Separation of event templates, baseline covariates, and longitudinal
   endpoints.
4. Construction and preprocessing of subject profiles.
5. Dose-relative time alignment and handling of irregular observation windows.
5b. Also describe handling of time/tad and how it works when different subjects have different times.
6. Scaling, missing-value treatment, and the rank-safe PCA representation.
6b. Handling of BLOQ data.  Here, there would be a CENS column in the dataset, as in Monolix
7. Schedule compatibility, donor selection, and distance calculation.
8. Randomized donor weighting, normalization, and the dominant-weight cap.
9. Endpoint-specific transformation, interpolation, blending, and
   back-transformation.
10. Subject shifts, AR(1) residual perturbations, and endpoint constraints.
11. Preservation of coherent event records and deterministic MDV behavior.
12. Restoration of names, column order, classes, factors, and schema.
13. Reproducibility, seeds, returned settings, validation, and diagnostics.
14. Important edge cases, failure modes, and limitations.

Use LaTeX for all important equations. Define every symbol immediately and
connect each equation to the relevant package operation. Include equations for
standardized profiles, PCA representation, distances, raw donor weights,
normalized and capped weights, interpolation/blending, transformations,
subject-level shifts, and AR(1) noise where applicable.

Use one small, transparent worked example that follows a subject through the
full algorithm. Supplement it with concise examples from the package API where
useful. Clearly distinguish:

- quantities copied from an anchor event template;
- quantities blended from donors;
- quantities perturbed stochastically; and
- metadata restored after synthesis.

Emphasize that the method is AVATAR-inspired but is not an exact implementation
of published AVATAR software. Cite the AVATAR paper stored in references (Destere26.pdf) 
State clearly that the output is mock data for
model-workflow exploration—not anonymous data—and offers no formal privacy
guarantee or basis for scientific inference.

Favor explanatory diagrams, compact tables, equations, and annotated examples
over long blocks of prose. Introduce technical concepts intuitively before
giving their mathematical form. End with a concise algorithm summary and a
practical checklist for users evaluating whether a generated dataset is
structurally suitable.

Add valid vignette metadata, render and visually inspect the result, regenerate
package documentation if needed, run the complete test suite, and run
`R CMD check`. Do not change package behavior unless an implementation
discrepancy must be corrected and explicitly documented.
