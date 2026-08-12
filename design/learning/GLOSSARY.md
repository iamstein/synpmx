# Glossary

Terms defined **in your own words, from memory, before checking the source.**
Writing the definition is the retrieval exercise; the file is a by-product.

Block 1 of the routine in `EFFECTIVE_LEARNING.md` ends with glossary entries.
Re-test an entry by covering the definition and writing it again, then compare.

Format, and one worked example:

---

---

## Queue

Terms met in the reviews and not yet defined. Delete a line when it becomes an
entry above.

**Bold** marks the terms that bear directly on `synpmx_avatar()`: either it names
a mechanism the generator actually applies, or it names a measure that only
makes sense for a record-based generator where every avatar has an anchor. The
rest are the surrounding literature, worth recognizing rather than mastering.

- **Adversarial accuracy** — implemented as `compare_pmx_proximity()`
- Privacy loss (the delta in AA_test − AA_train)
- **Local cloaking** — per patient, and needs the anchor map
- **Hidden rate** — the aggregate of local cloaking
- **Distance to closest record (DCR)**
- **Nearest-neighbour distance ratio (NNDR)**
- Key / quasi-identifier
- Target (in the disclosure sense)
- RepU
- DiSCO
- DiO
- **Singling out** — the criterion the checks cover thoroughly
- **Linkability** — the criterion a record-based generator fails by construction
- Inference (WP29 sense)
- Membership inference
- **k-anonymity** — `min_pattern_share` = 2 is this, on visit patterns
- l-diversity
- t-closeness
- Composition attack
- pMSE
- Confidence interval overlap
- alpha-precision
- beta-recall
- **Authenticity** — the machine-learning name for local cloaking
- **Mode drop** — losing the tails, which masking causes on purpose
- Epsilon, and (epsilon, delta)
- **Sensitivity** — `max_donor_weight` = 0.50 is a bound of this kind
- Composition (DP sense)
