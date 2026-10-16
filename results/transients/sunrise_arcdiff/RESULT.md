# Sunrise arc difference imaging — result

Box: 24″ × 24″ centred on Earendel (Scofield+2025 1.1e, 24.346822, −8.464511), which holds the Sunrise arc's
catalogued images. `arc_diff.py --kernel-px 4 --kernel-stamps --inject 100 --min-bands 2` (full command in
`o010_o052/summary.json`), run with `OMP_NUM_THREADS=1`. Every number is `derived`; thresholds are ASSUMPTIONs.
Sheets: `scripts/arc_diff_sheet.py` (epoch-1 image and S at ±10, 2″ stamps), looked at before writing this.

## Deviations from PLAN.md
- **Deviation:** validation used `--inject 100` (per-band point injections on pixels above the box's 90th
  percentile), not `--inject-grid`; the PLAN's "≥ 0.8 at 5× the box σ" bar was not evaluated as written. The
  injection figures below are per band and before the morphology step.
- **Deviation (postponed):** the o010 → o120 pair did not finish in this session (see below).
- **Post hoc:** "on the arc" is made computable as the distance to the nearest Scofield+2025 multiple image
  (`o010_o052/image_distances.json`); the PLAN did not define it.

## o010 (2022-07-30, jwst 2.0.1) → VENUS o052 (2025-07-08, jwst 3.0.0), F150W + F444W
- Noise scales 1.20 / 1.29; flux scales 1.024 / 1.026 (`summary.json`).
- Injection recovery per band (`summary.json` `injection_recovered`): F150W 0.86 / 0.14 / 0.01 at
  27.5 / 28.0 / 28.5 AB; F444W 0.83 / 0.28 / 0.04.
- **9 hits** (6 `near_bright_core`; `hits.ecsv`). `arc_diff_morph.py` (defaults; `morph.json`): 2 dipole,
  2 sign_flip, 1 edge, 4 pass.
- **The 4 passes, by eye** (`o010_o052/sheet.jpg`):
  - h001, h002: mixed-sign residuals on the cores of a compact galaxy and a bright face-on galaxy (`near_bright_core`);
  - h005: negative residual on a compact galaxy's core in both bands (S −10.0 F150W, −19.0 F444W), 0.50–0.78″ from
    the map edge (`edge_arcsec`); in F444W a positive lobe of +19.9 sits beside it (dipole ratio 1.05; 0.36 in
    F150W, so the two-band dipole rule does not fire);
  - h007: elongated negative residual along an edge-on disk (F150W extent 1.59″).
- All 9 hits lie 3.98″ or more from every catalogued multiple image (`image_distances.json`; nearest h001 to 1.7a).
- **Outcome:** 0 hits within 3.98″ of a catalogued Sunrise-arc image, so 0 candidates on the arc's images in this
  pair. This is not a limit on the whole box below the injection depth.
- **Earendel (control row, not a hit):** S = −5.2 (F150W), −1.5 (F444W) at the reference position (`summary.json`
  `known`); the 2-band rule does not fire and the stamp shows no compact residual at ±10. The catalog forced
  photometry (recentred) gave Δm = +0.05 mag in F150W, 0.3σ calibrated (docs/fields/sunrise.md). The −5.2 is not
  explained here.

## o010 → o120 (2023-01-10, jwst 2.0.1), F090W F115W F277W F356W
Not finished (not stored: the run was stopped with no output). With `--procs 2` beside the o052 run, it ran for
≈ 45 min without finishing its first two SW bands. Next: run it alone with `--procs 4` and `OMP_NUM_THREADS=1`. This
pair has the same pipeline version in both epochs and a ~180° PSF rotation (docs/fields/sunrise.md), so rosettes
around compact sources are expected.
