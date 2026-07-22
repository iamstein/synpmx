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

## Testing

- A description of a simulation testing scirpt should be made available in TEST_SIM.md. 
- This will use a number of publicly available datasets, starting with the three in the demo.
- Checks should be written to see if there are issues with the data.
- This scirpt should be maintained.  Every time you find an issue with teh simulator, add it to the list of issues in the design/TEST_SIM.md file.