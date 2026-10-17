# `arcs` on WHL0137−08 (Sunrise): result (2026-10-17)

Plan: `PLAN.md` (df3bcd7, committed before any count). Script `scripts/arcs_whl0137.py` (the RXJ2129 driver
`arcs_rxj2129.decide`, outcome rule `arcs_macs0647.outcome`, both unchanged in logic). Labels: shapes `observed`;
predicted axes `model_prediction`; classes, counts and P values `derived`; thresholds ASSUMPTIONs (D-024 defaults).

## Inputs (pinned)
- F200W catalogue `jw02282-o010_t001_nircam_clear-f200w_cat.ecsv`, sha256 `b0693065…` (manifest `sunrise`).
- DJA v7.5 `sunrise-grizli-v7.5-fix.eazypy.zout.fits`, sha256 `672d4bdf31926209633c1a7feeb436bab3e5b08e847a079521820f3409f548e0`,
  streamed out of the 51.2 MB photo-z tarball (tarball not stored).
- RELICS v1 maps: `whl0137-relics-lenstool` (primary) and `whl0137-relics-glafic` (second, registered in this run;
  sha256 in `lensmodel.WHL0137_RELICS_GLAFIC`). Both pass `validate` (Lenstool κ median relative difference
  3.5 × 10⁻⁵; GLAFIC in `whl0137-relics-glafic/validate.json`).

## Screen (`summary.json`, `<model>/arcs.json`)
| Model | strong-shear background | aligned | anti | expected anti (shape noise) | P(≥ n_anti) | control (≥ 70 %) |
|---|---|---|---|---|---|---|
| Lenstool (primary) | 18 | 13 (72 %) | 1 | 0.150 | 0.139 | passes |
| GLAFIC (second) | 19 | 14 (74 %) | 1 | 0.232 | 0.207 | passes |

**Frozen outcome: "does not replicate"** (deciding model Lenstool, P 0.139 ≥ 0.05). The control passes narrowly
(13 / 18; one fewer aligned row would fail it).

## The anti row
- **2194:** m 24.36, e 0.79, g_min 0.24 (weak end of the strong-shear cut), r 45″, z_phot 1.17 (no grid-edge
  value), offset 63.1° (just past the 60° anti cut); anti in GLAFIC too.
- Looked at (`background.png`, all 18 rows, anti first): a thin, edge-on-looking disc with the same axis in F200W
  and F444W, isolated. Leading ordinary explanation: intrinsic shape in weak shear, which the shape-noise null
  models. Not sent to `/vet-candidate`; not a candidate.
- Pre-registered secondary checks: no z_phot > 15 row among the anti; F444W axis read from the cutout only
  (`arcs_bands.py` not run).

## Post hoc (not in the plan)
- Unlike S1063's four anti (m 27.4–28.2, g 1.5–1.9 for the open two), 2194 is bright and in weak shear: the S1063
  faint, high-shear pattern does not appear here. Background rows at g_min ≥ 1: 4 / 4 aligned.
- On the sheet, 1587 (m 27.5, `mixed`) is barely visible at this stretch.

## Limits
- 18 rows: the test can only detect an excess of a few anti; RELICS models are HST-era with few constraints.
- One epoch (no second F200W epoch here); single photo-z catalogue.
