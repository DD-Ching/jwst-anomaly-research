# `arcs` on MACS1423: screen result (2026-10-16)

Plan: `PLAN.md` (committed 18:13Z, before the run). Labels: shapes `observed`; predicted axes `model_prediction`;
classes, counts and P values `derived`. All thresholds are the D-024 defaults (ASSUMPTIONs).

## Inputs
- F200W pipeline catalogue `jw01208-o010_t005_nircam_clear-f200w_cat.ecsv` (CANUCS, 6,399 s; MAST level 3;
  4,987,129 B, sha256 `1a968574…b5bfa3`).
- CANUCS DR1 `macs1423-clu` photometry (25,278,098 B, sha256 `c257d8df…ac622b0`), converted with
  `inject_radial.canucs_zout`.
- Model `macs1423-canucs-maps` (pinned in `lensmodel.py`): CANUCS DR1 Lenstool maps (Desprez, preliminary v1,
  5 spectroscopic systems), 100 mas, z_lens 0.545. `validate`: κ from the deflections vs the published κ, median
  relative difference 3.0 × 10⁻⁵ (`macs1423-canucs-maps/validate.json`).

## Result (`summary.json`, `macs1423-canucs-maps/arcs.json`)
| Strong-shear background | Aligned (control ≥ 70 %) | Anti | Expected (shape noise) | P(≥ anti) | Frozen outcome |
|---|---|---|---|---|---|
| 55 | 44 (80 %) | 4 | 0.949 | 0.016 | **replicates** (screen level) |

- By g_min: g < 1 48 rows, 38 aligned, 3 anti; g ≥ 1 7 rows, 6 aligned, 1 anti.
- Anti rows (`summary.json` `anti_background`): 1818 (m 28.46, z_phot 3.11, g_min 1.03), 2197 (m 28.34, z 3.42,
  g 0.27), 2322 (m 28.08, z 3.42, g 0.27), 4234 (m 24.19, z 1.01, g 0.20).

## Vetting
Step 1 (ICL-subtracted shapes, measurement noise): `vetting/RESULT.md`. 1818 becomes `mixed`; the other three
stay anti and are not noise-explicable; post hoc P(≥ 3) 0.071 vs shape noise.

## Looked at (F200W 5″ cutouts with observed and predicted axes, 4 anti + 4 random aligned controls; not stored)
- 1818: not visible at the sheet's stretch, on the BCG's intracluster light next to a long thin arc;
- 2197: one knot of a clumpy chain whose overall axis roughly follows the prediction (a fragment candidate);
- 2322: a tiny, faint blob;
- 4234: a bright, clearly elongated galaxy at the g_min threshold (0.20), where intrinsic shape dominates;
- the four aligned controls follow their predicted axes.

## Reading (screen level only; no anti row is vetted)
- Under the pre-registered rule this is a replication of the S1063 pattern in a field never used to tune `arcs`.
  It is **not** a candidate and not evidence of anything beyond the screen: per PLAN.md the anti rows go through
  the S1063 vetting chain (#273) before any reading.
- **Post hoc observation (not tested):** 3 of the 4 anti rows are at m ≥ 28.0, as all four S1063 anti rows are at
  m ≥ 27.4 (#273). The shape-noise null (σ_ε 0.26) contains no measurement noise; this is the first ordinary
  explanation to test.
- **Limits:** one model, preliminary and constrained by 5 spectroscopic systems; photo-z background cut; no frame
  offset measured; shape-noise null only.

## Data handling
The three maps (87–89 MB each) and the photometry (25 MB) were cached in the gitignored data root and deleted after
the run.
