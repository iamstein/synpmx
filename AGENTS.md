# Repository instructions

This repository contains the `pmxSynthData` R package prototype for simulating
structurally faithful mock pharmacometric datasets.

- Put package functions in `R/`, tests in `tests/testthat/`, and runnable
  demonstrations in `scripts/`.
- Document public functions with roxygen2 and regenerate documentation after
  API changes.
- Keep simulation assumptions, units, schemas, and seeds explicit.
- Do not commit sensitive, proprietary, or patient-level data.
- Treat `data/` and `output/` as local/generated unless told otherwise.
- Preserve unrelated changes and avoid adding dependencies unnecessarily.
- Run the full tests and `R CMD check` after behavioral changes.
