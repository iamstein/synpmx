# Objective

We've developed a synthetic data generation R package synpmx (https://iamstein.github.io/synpmx/) for sharing realistic-looking data outside the GxP computing environment.  This synthetic data can be used in generating code for data exploration, model building, and diagnostics.  When the latest GenAI tools are not available on the systems hosting the clinical data, synthetic data can be generated and transferred to other machines that do run the latest GenAI tools.  

# Background: Four approaches to Synthetic Data Generation

Provide strengths, limitations, and applications.
Mention differential privacy.

1. Trial Simulation (rxode2)
2. Deep Generative Model (GenAI)
3. Sequential Conditional Models (synthpop)
4. Blending (AVATAR)

# synpmx-avatar method

AVATAR-inspired blending method that extends past work (2023, 2026).
Include ✅/⚠️/❌ checks

## The algorithm

13 steps
6 sources of masking

## Applied to example (result)

probably mad from xgxr

- show simulated data
- show show distributions
- show scorecard

## What it adds

- BLOQ handling
- Weight-based dosing
- Proper event handling
- Removal of extreme outliers (in time, dose, observations)
- Removal of dose and observation timing fingerprint

## Use at our company (not sure yet what this section will read)

- We've tried with X datasets and it gives something reasonable
- Status of where we are with it in terms of approval?

## Gaps

- no time varying covariates (but could be treated as another another DV with long dataset)
- baseline levels won't match if you have a baseline column, but that's an easy fix, just rederive baseline levels
- no formal privacy protection
- not for scientific discovery (relationships between variables not fully preserved)

# Survey

- How are you permitted to use agentic coding at your company.  Direct access to GxP, not at all, something in between.
- What is your role.
- What size company are you.