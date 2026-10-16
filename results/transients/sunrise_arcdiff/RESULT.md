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
  - h005: negative residual on a compact galaxy's core in both bands (S −10 F150W, −19 F444W), 0.5–0.8″ from the
    box border; in F444W an equal positive lobe sits beside it (dipole ratio 1.05; 0.36 in F150W, so the two-band
    dipole rule does not fire): a registration/PSF residual of a compact core at the box edge;
  - h007: elongated residual along an edge-on disk (PSF mismatch along the disk).
  - All 9 hits are ≥ 4.0″ from every Scofield+2025 multiple image (`scripts/sunrise_scofield.py` IMAGES; nearest:
    h001 and h003 at 4.0″), so none lies on the Sunrise arc's catalogued images. **0 candidates.**
- **Earendel (control, not a hit):** S = −5.2 (F150W), −1.5 (F444W) at the reference position; the 2-band rule
  does not fire and the stamp shows no compact residual at ±10. The forced photometry (recentred, docs/fields/sunrise.md)
  gave Δm +0.05 ± 0.16 in F150W; the unrecentred 0.15″ offset and the 2.0.1 → 3.0.0 change are the likely cause of the
  −5.2. Not followed up here.

## o010 → o120 (2023-01-10, jwst 2.0.1), F090W F115W F277W F356W
Not finished. With `--procs 2` next to the o052 run, the first two SW bands had not finished after 45 CPU-min
each, and the run was stopped at the end of the session (no output). Next: run it alone with `--procs 4`
(command in `PLAN.md`'s box and bands; `OMP_NUM_THREADS=1`). This pair has the same pipeline version in both
epochs and a ~180° PSF rotation (docs/fields/sunrise.md), so rosettes around compact sources are expected.
