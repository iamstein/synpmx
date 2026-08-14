# Practice

Terms defined **in your own words, from memory, before checking the source.**
Writing the definition is the retrieval exercise; the file is a by-product.

- **Adversarial accuracy** - 
---

---

## Glossary

Terms met in the reviews and not yet defined. Delete a line when it becomes an
entry above.

**Bold** marks the terms that bear directly on `synpmx_avatar()`: either it names
a mechanism the generator actually applies, or it names a measure that only
makes sense for a record-based generator where every avatar has an anchor. The
rest are the surrounding literature, worth recognizing rather than mastering.

- **Adversarial accuracy** — How often the closest patient in your dataset (trial/synthetic) comes from the other dataset vs your own.  Target around 50%  `compare_pmx_proximity()`
- _Privacy loss_ - How closely the synthetic patients match individuals from the training set vs matching the overall population.  How much they learned from the data (the difference in $AA_test$ − $AA_train$).  
- **Local cloaking** — For patient $i$ and synthetic avatar $a(i)$, how many other avatars lie closer to $a(i)$ than $i$.  
- **Hidden rate** — Summary of local cloaking - the percentage of patients whose nearest avatar is not their own avatar.
- **Distance to closest record (DCR)** - Distance to closest real record, for each synthetic record
- **Nearest-neighbour distance ratio (NNDR)** - ratio of nearest distance (between real and synthetic) to second nearest distance record .  0 means one real record is much closer than everything else, 1 means there's a crowd.   
- Key / quasi-identifier - combination of variables that an attacker can use to identify a patient
- Target (in the disclosure sense) - variable an attacker may want to identify about a patient
- RepU - "replicated uniques" = percentage of records that are unique in the original data and also the synthetic data.  
- DiSCO - attributed disclosure - 
- DiO
- WP29 - three part test of whether dataset is anonymized.  It covercs the three items below
- 1. **Singling out** — isolate records that identify an individual
- 2. **Linkability** — ability to link multiple records as belonging to teh same individual.
- 3. Inference (WP29 sense) - deduce the value of an attribute for an individual, based on other values
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
