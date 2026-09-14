| Dataset | Patients | Design feature | Main finding |
|---|---:|---|---|
| `case1_pkpd` | 180 | Multiple arms; censored concentrations | Censoring represented; pooled PD loses arm differences. |
| `mad` | 60 | Multiple doses; mixed endpoint types | Endpoint types retained; early PK peak underrepresented. |
| `warfarin` | 32 | Single oral dose; delayed response | PD decline represented; recovery absent. |
| `wbcSim` | 45 | Infusions; white-cell nadir and recovery | Response fitted as PK; nadir and recovery missed. |
| `mavoglurant` | 120 | Repeated occasions; resetting clock | Second occasions lost; fewer observations generated. |
| `theo_md` | 12 | Repeated oral doses; small cohort | Repeated dosing retained; fit rests on a small cohort. |
| `nimoData` | 12 | Weekly infusions; small cohort | Infusion schedule retained; concentration spread inflated. |
| `pheno_sd` | 59 | Neonatal care; sparse, irregular sampling | Nominal grid required; fewer doses generated. |
| `mixroute_sim` | 90 | Intravenous and subcutaneous dosing | Both routes retained; bioavailability close to simulation truth. |
| `onc_sim` | 200 | Trough sampling; tumour response; dose changes | Pooled tumour curve loses arm differences; fewer doses generated. |
