# Synthetic Data Generation Approaches for PK/PD Studies

| Approach | Example methods | How it works | Source patients | Longitudinal data | Reportable fingerprint | Main tradeoff |
|---|---|---|---:|---|---|---|
| **Mechanistic simulation** | `rxode2`, `mrgsolve`, Simulx | Simulate from explicit trial + PK/PD model | **10s** | **Yes** | **Yes, Simple** | Highly interpretable, but requires model assumptions |
| **Patient blending** | AVATAR, PCA/KNN, SMOTE-like methods | Recombine or perturb observed patient profiles | **10s–100s** | **Yes** | **No** | Preserves realistic profiles, but constrained by observed patients |
| **Sequential conditional** | `synthpop` | Generate variables sequentially from fitted conditional models | **100s** | **Challenging** | **No** | Simple and inspectable, but order-dependent and awkward for trajectories |
| **Deep generative model** | CTGAN, TVAE, medGAN, CopulaGAN | Learn a high-dimensional joint/latent distribution | **100s–1000s** | **Challenging*** | **No** | Flexible, but data-hungry and difficult to interpret |

\* Longitudinal structure is possible with specialized architectures, but is not automatic.

**Mechanistic simulation uniquely provides an explicit study fingerprint:** design + population + dosing + PK + PD + variability + observation model.
