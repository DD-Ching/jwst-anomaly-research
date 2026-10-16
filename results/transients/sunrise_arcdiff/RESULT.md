# Sunrise arc difference imaging — result

Box: 24″ × 24″ centred on Earendel (Scofield+2025 1.1e, 24.346822, −8.464511), which holds the whole Sunrise arc.
`arc_diff.py --kernel-px 4 --kernel-stamps --inject 100 --min-bands 2`, run with `OMP_NUM_THREADS=1` (without
it, two runs on 4 cores made no progress in 20 min). Every number is `derived`; thresholds are ASSUMPTIONs.
Sheets: `scripts/arc_diff_sheet.py` (epoch-1 image and S at ±10, 2″ stamps), looked at before writing this.

## o010 (2022-07-30, jwst 2.0.1) → VENUS o052 (2025-07-08, jwst 3.0.0), F150W + F444W

- Noise scales 1.20 / 1.29; flux scales 1.024 / 1.026 (the 2.0.1 → 3.0.0 calibration change); 18 min wall.
- Injected point brightenings on bright pixels (> 90th percentile of epoch 1: arc and galaxies), per band:
  F150W 0.86 / 0.14 / 0.01 at 27.5 / 28.0 / 28.5 AB; F444W 0.83 / 0.28 / 0.04. The 27.5-mag figure is close to the
  PLAN's 0.8 bar (an injection-only figure; the morphology step is not injection-calibrated here).
- **9 hits** (6 near a bright core). `arc_diff_morph.py` (defaults): 2 dipole, 2 sign_flip, 1 edge, 4 pass.
- **The 4 passes, by eye** (`o010_o052/sheet.jpg`):
  - h001, h002: mixed-sign residuals on the cores of a compact galaxy and a bright face-on galaxy (PSF/kernel
    residuals of bright cores, `near_bright_core`);
  - h005: F444W-only dipole next to a galaxy near the F444W map border; nothing in F150W;
  - h007: elongated residual along an edge-on disk (PSF mismatch along the disk).
  - None lies on the Sunrise arc. **0 candidates.**
- **Earendel (control, not a hit):** S = −5.2 (F150W), −1.5 (F444W) at the reference position; the 2-band rule
  does not fire and the stamp shows no compact residual at ±10. The forced photometry (recentred, docs/fields/sunrise.md)
  gave Δm +0.05 ± 0.16 in F150W; the unrecentred 0.15″ offset and the 2.0.1 → 3.0.0 change are the likely cause of the
  −5.2. Not followed up here.

## o010 → o120 (2023-01-10, jwst 2.0.1), F090W F115W F277W F356W
Pending (see the PR handoff).
