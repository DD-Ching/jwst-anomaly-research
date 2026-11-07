# MIRI re-rank with band-restricted NIRCam confirmation (follow-up to `PLAN_bands.md`)

Status: plan, written before measurement (2026-11-07, about 02:20Z).

- **Change.** A position-only join (no colour replacement) marks a MIRI source NIRCam-confirmed when a DJA v7.4
  SMACS source with S/N ≥ 5 in at least one NIRCam band (0.5″ aperture) lies within 0.3″ (nearest neighbour, as
  `scripts/miri_confirm.py`). `quality.assess_sources` lifts `single_band` only for sources detected in a
  `confirm_bands` band: [F770W, F1000W] (D-014 addendum 2). Other gate tests unchanged.
- **Hypothesis.** The re-run gate passes 113 + 55 = 168 sources (113 as before, plus the 45 F770W and 10 F1000W
  single-band-only matches); the top 20 changes, with some newly admitted F770W/F1000W singles entering it,
  because single-band sources have extreme (missing-colour) feature vectors.
- **Ordinary explanations for a new top-20 entry.** A chance NIRCam match (expected ≈ 0.046 × 66 + 0.020 × 20 ≈
  3.4 of the 55, from the per-band chance rates); a real faint galaxy whose MIRI colours are missing (feature
  imputation, not physics); spike or halo residuals that happen to coincide with a DJA source.
- **Checks (frozen).** (1) the gate count equals 168 on the same catalogs (else the join differs from the
  script: stop and find out why); (2) top-20 overlap with the previous gating (Jaccard of source_uid sets) is
  reported; (3) every new top-20 entry gets a cutout look on the contact sheet before any statement about it.
- **What changes the plan.** If more than half of the new top 20 are newly admitted singles, the ranking is
  dominated by missing-colour imputation, and the next step is to rank single-band MIRI sources as their own
  stratum (not done here) rather than to vet them as anomalies.
