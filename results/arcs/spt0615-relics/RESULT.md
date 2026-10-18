# `arcs` on SPT-CL J0615−5746: result (2026-10-18)

Plan: `PLAN.md` (commit 8210bef, 18:21:39Z, before any `arcs` count; the "18:25Z" in its header was an estimate, the
commit time is the record). Script `scripts/arcs_spt0615.py` (BPZ conversion and astrometric shift, unit-tested in
`tests/test_arcs_spt0615.py`; driver `arcs_rxj2129.decide` and outcome rule `arcs_macs0647.outcome` unchanged).
Sheets `scripts/arcs_spt0615_sheet.py`. Labels: shapes `observed`; photo-z `derived` (RELICS BPZ); predicted axes
`model_prediction`; classes, counts and P values `derived`; thresholds ASSUMPTIONs (D-024 defaults). Exploratory: the
photo-z source (RELICS BPZ, 95 % interval as z160–z840) differs from the DJA eazy tables of the earlier fields.

## Inputs (pinned)
- F200W catalogue `jw04212-o001_t001_nircam_clear-f200w_cat.ecsv` (MAST, PUBLIC; the only public F200W level-3
  catalogue within 2′), sha256 in `summary.json`.
- RELICS BPZ `hlsp_relics_hst_acs-wfc3ir_spt0615m57_multi_v1_cat.txt` → 9,782 galaxy rows (stel < 0.8); converted
  table sha256 in `summary.json`. Astrometry (`photoz_offset.json`): median JWST − BPZ offset +0.306″ (RA) and
  −0.273″ (Dec) over 2,392 pairs, above the 0.1″ rule, so it was removed before the 0.3″ match (as planned).
- RELICS maps `spt0615-relics-glafic` (v2, primary) and `spt0615-relics-lenstool` (v1, second); sha256 in
  `lensmodel.SPT0615_RELICS_*`. Both pass `validate`: κ median relative difference 1.9 × 10⁻⁵ (GLAFIC) and 1.5 ×
  10⁻⁵ (Lenstool); z = 2 μ 9.7 × 10⁻⁶ and 3.1 × 10⁻⁶ (`<model>/validate.json`).

## Screen (`summary.json`, `<model>/arcs.json`)
| Model | strong-shear background | aligned | anti | expected anti (shape noise) | P(≥ n_anti) | control (≥ 70 %) |
|---|---|---|---|---|---|---|
| GLAFIC v2 (primary) | 7 | 3 (43 %) | 2 | 0.119 | 0.0065 | fails |
| Lenstool v1 (second) | 4 | 1 (25 %) | 3 | 0.132 | 0.00034 | fails |

**Frozen outcome: "no result"** (no run passes its control). Under the frozen rule neither P value decides anything:
when the predicted axes fail to match most strong-shear background shapes, the anti count cannot be separated from
model error at this field.

## Strong-shear background rows, looked at (`background_glafic.png`, `background_lenstool.png`; 4″ F200W / F444W)
- **2681** (GLAFIC `mixed`, g_min 0.67; Lenstool `anti`, g_min 4.84; z_phot 4.47, m 26.5): in Lenstool it lies at
  g_min ≥ 1 for its own redshift range, where a critical-curve flip makes the drawn z = 2 axis and the classifying
  axis differ (in the sheet the observed axis lies close to the z = 2 prediction). Ordinary explanation: model
  critical-curve position (the two models differ by 7× in g_min here).
- **2669** (anti in both; m 24.4, z_phot 3.0, g_min 0.42 / 0.27, r 46.5″): compact, with an elongated faint
  extension; the F200W stamp shows horizontal striping across it (by eye; not measured), the F444W axis is similar.
- **2773** (anti in both; m 26.4, z_phot 3.4, g_min 0.27 / 0.35): faint and noisy at this stretch.
- No row is sent to `/vet-candidate`; nothing here is a candidate.

## Limits
- Very few strong-shear background rows (7 and 4): the lens redshift 0.972 shrinks the region with g_min ≥ 0.2 for
  z ≲ 3 sources, and the 95 % BPZ interval used as the background test removes more rows than a 68 % one would.
- Both models are HST-era; their predicted axes do not match the JWST shapes here (control fails), so this field
  tests the models more than the S1063 hypothesis.
- One F200W epoch (program 4212, 2,920 s).
- The RELICS model rests on three multiple-image systems (two with spectroscopic redshifts, z 1.358 and 4.013;
  Paterno-Mahler et al. 2018, arXiv:1805.09834, abstract read 2026-10-18), so a weak shear-direction prediction
  away from those images is not surprising. The plan also cites Bleem et al. 2015 for z_lens; that citation was not
  checked (deviation note: the code cites only Paterno-Mahler et al. 2018).
