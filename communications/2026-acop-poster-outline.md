# Deliverable

We've developed a synthetic data generation R package synpmx (https://iamstein.github.io/synpmx/) for sharing realistic-looking synthetic data outside the GxP computing environment.  This synthetic data can be used in generating code for data exploration, model building, and diagnostics.  When the latest GenAI tools are not available on the systems hosting the clinical data, synthetic data can be generated and transferred to other machines that run the latest GenAI tools.

# Background: Four synthetic data generation approaches were considered:

refer to table in 2026-acop-poster-table-synthetic_pkpd_methods.md

# What is available in `synpmx`

note from author: here there should be a table of the methods.  the approach (blending, mechanistic model, etc.), whether it offers formal privacy guarantees.  prespecification required.

Methods based fully off the data

1. `synpmx_model()` — Fits simple PK and PD models to the observation data.
2. `synpmx_avatar()` — Blended values from real patients
3. `synpmx_pca()` — Principal component analysis from real patient data.
4. `synpmx_prior()` — Simulates a public model and a public protocol.  Reads no data at all.
5. `synpmx_calibrated()` — Keeps the public model's shape, and spends a small privacy budget
   correcting its magnitude.
6. `synpmx_empirical()`

# Primary synthetic data method: mechanistic model of PK, PD, dose changes, and missed observations

Author note: Describe the method.  It will involve
- dose model
- missed observation model
- Pooled-PD model
- Pop-PK model


Check these things still happen
- BLOQ handling
- Weight-based dosing
- Proper event handling
- Removal of extreme outliers (in time, dose, observations)

# Example

probably mad from xgxr

- show simulated data
- show show distributions
- show scorecard (full scorecard if room)

# Evaluated datasets

Author note: fill in M and N

- N public datasets
- M simulated datasets
- 3 complex internal datasets (mulitple PK, mixed IV + SC, intrapatient escalation)

# Current Status

- R package available in github: https://iamstein.github.io/synpmx/
- Internal discussions ongoing.
