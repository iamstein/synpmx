# Repository instructions

This is a lightweight R research/development project for simulating synthetic
pharmacometric datasets. It may become a package later; do not assume package
structure or add package tooling unless requested.

- Put reusable R functions in `R/` and runnable work in `scripts/`.
- Keep simulation assumptions, units, schemas, and seeds explicit.
- Do not commit sensitive, proprietary, or patient-level data.
- Treat `data/` and `output/` as local/generated unless told otherwise.
- Preserve unrelated changes and avoid adding dependencies unnecessarily.
