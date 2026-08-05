# Private company-data work

Scratch space for trying the calibrated generator on a real modeling dataset,
inside the safe computing environment.

## Be very careful in this folder what goes into git

The root `.gitignore` ignores everything here **except** `README.md` and some files that are
explicitly declared as acceptable.  Both real data and any outputs from the synthetic data generator
are ignored.  The only thing acceptable here is code (.R, .Rmd files) which can be run in the GxP 
environment.  

## Inventory of Studies thate are explored:

### PIT565 A1
- **Template** — `try_avatar_pit565a1.qmd`
- **Design shape** — Phase 1 dose escalation, with subsequent weekly dosing.
  Dose interruptions exist.
- **Why it is here** — first real schema
- **Checked** — roles, validation

### PIT565 B1
- **Template** — `try_avatar_pit565b1.qmd`
- **Design shape** — Phase 1, dose escalation
- **Why it is here** — fixed 3 doses per subject with intrapatient dose
  escalation
- **Checked** — roles, validation

### ECI - Oncology study *(to add)*
- **Template** — not yet written
- **Design shape** — repeated dosing with intra-patient escalation
- **Why it is here** — Many different regimens
- **Checked** — not yet run