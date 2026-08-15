## Features added beyond AVATAR paper
- BLOQ handling
- Weight-based dosing
- Proper event handling
- Removal of extreme outliers
- Removal of dose and observation timing fingerprint

## Nice to haves
- Keeps correlation structure to some degree
- Seems to handle multiple dosing regimens (of all kinds)

## Gaps
- no time varying covariates (but would be another DV with long dataset)
- baseline levels won't match if you have a baseline column, that'd need to be rederived.

## Methods

- **1. Prior-Based Simulation**
  - **Information source:** Published/literature models, prior knowledge, protocol assumptions.
  - **Study-specific information released about true dataset:** None.
  - **Privacy:** Independent of the protected study; no patient information is communicated.
  - **Strengths:**
    - Strongest privacy.
    - Well suited for software development, AI-assisted coding, teaching, and prototyping.
    - Familiar workflow for pharmacometricians.
  - **Limitations:**
    - Does not reflect characteristics of a specific study population.
    - Requires a reasonable prior model and protocol assumptions.

- **2. Differentially Private Synthetic Data**
  - **Information source:** Learns from the protected study under differential privacy constraints.
  - **Study-specific information released:** Yes, but with mathematically bounded disclosure.
  - **Privacy:** Formal differential privacy guarantee (ε, δ).
  - **Strengths:**
    - Preserves useful study-specific information.
    - Provable protection against disclosure of individual participants.
  - **Limitations:**
    - More technically complex.
    - Privacy-utility tradeoff may reduce fidelity.
    - Requires specialized methods.

- **3. Study-Derived Synthetic Data (e.g., AVATAR)**
  - **Information source:** Directly derived from the protected study.
  - **Study-specific information released:** Yes; designed to preserve much of the original study structure.
  - **Privacy:** Reduced re-identification risk, but often no formal differential privacy guarantee.
  - **Strengths:**
    - High fidelity to the original study.
    - Preserves complex relationships that may be difficult to model explicitly.
  - **Limitations:**
    - Greater dependence on the original dataset.
    - Privacy protection is empirical rather than mathematically guaranteed.

## synpmx avatar

## synthetic data requirements

- no patient identifiable from various features like
unique doses or dose times or observation times etc

## synthetic data checks

- distance distribution
- identifiability check

## Possible Survey

- How are you permitted to use agentic coding at your company.  Direct access to GxP, not at all, something in between.
- What is your role.
- What size company are you.