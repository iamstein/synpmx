# Scripts

- `demo_nlmixr2data.R` runs the formal-privacy (Version 2) fit-once/generate-many demonstrations for
  the public `theo_md`, `warfarin`, `wbcSim`, `nimoData`, and `mavoglurant`
  datasets after `synpmx` is installed. It deliberately uses the guarded
  public-fixture backend because
  those example sources are already public, and labels the resulting models as
  non-DP. Confidential fitting uses the default OpenDP backend and fails closed
  when it is unavailable. Comparison figures put source facets above synthetic
  facets, connect each subject's observations on the appropriate scientific
  clock, and use
  endpoint-specific linear DV scales. The script fails if patient counts,
  observations per patient, or endpoint time coverage differ materially.
- `test-avatar.R` exercises the primary AVATAR-style `synthesize_pmx()` path
  over the same public datasets as a fast structural smoke test.
