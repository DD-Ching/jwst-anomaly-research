# `arcs` on RXC J0600−2007 (MACS J0600.1−2008): result (2026-10-18)

Plan: `PLAN.md` (commit c764d6b, 17:20:09Z, before any `arcs` count; the "17:30Z" in its header is a typo, the commit
time is the record). Script `scripts/arcs_rxc0600.py` (driver `arcs_rxj2129.decide`, outcome rule
`arcs_macs0647.outcome`, both unchanged). Sheets `scripts/arcs_rxc0600_sheet.py`. Labels: shapes `observed`;
predicted axes `model_prediction`; classes, counts and P values `derived`; thresholds ASSUMPTIONs (D-024 defaults).

## Inputs (pinned)
- F200W catalogue `jw06796-o010_t001_nircam_clear-f200w_cat.ecsv` (MAST, PUBLIC), sha256 in `summary.json`.
  Choice-rule correction (pre-run, no count seen): the footprint centre is ≈ 20″ from the map reference point
  (90.0340, −20.1358), not ≈ 13″ as the plan says; still inside the 30″ rule, and program 4573 (≈ 35″) stays out,
  so the choice is unchanged (offsets computed in the session, not stored).
- DJA v7.2 `rxcj0600-grizli-v7.2-fix.eazypy.zout.fits` (sha256 in `summary.json`), streamed out of the 33.8 MB
  photo-z tarball (tarball not stored).
- RELICS maps registered as `rxc0600-relics-glafic` (v3, primary) and `rxc0600-relics-lenstool` (v1, second); sha256
  in `lensmodel.RXC0600_RELICS_*`. Both pass `validate` (κ median relative difference 1.7 × 10⁻⁵ and 2.1 × 10⁻⁵;
  Lenstool z = 2 μ map 3.7 × 10⁻² (median ratio 1.028), GLAFIC 1.9 × 10⁻⁵; `<model>/validate.json`).

## Screen (`summary.json`, `<model>/arcs.json`)
| Model | strong-shear background | aligned | anti | expected anti (shape noise) | P(≥ n_anti) | control (≥ 70 %) |
|---|---|---|---|---|---|---|
| GLAFIC v3 (primary) | 12 | 11 (92 %) | 1 | 0.203 | 0.184 | passes |
| Lenstool v1 (second) | 18 | 12 (67 %) | 3 | 0.499 | 0.014 | fails |

**Frozen outcome: "does not replicate"** (deciding model GLAFIC, P 0.184 ≥ 0.05). The Lenstool run fails its
control, so under the frozen rule its P 0.014 decides nothing.

## Anti rows, looked at (`background_glafic.png`, `background_lenstool.png`; 4″ F200W / F444W)
- **581** (anti in both): m 21.2, e 0.69, z_phot 0.70 (Δz 0.27 behind the lens); g_min 0.39 and offset 63.6° in
  GLAFIC (just past the 60° anti cut), g_min 1.19 and offset 74.1° in Lenstool (g_min ≥ 1, where the two models
  disagree on the shear by 3× and a critical-curve flip is possible; `<model>/arcs.ecsv`). A bright, isolated,
  edge-on-looking disc with the same axis in F200W and F444W. Leading ordinary explanation: intrinsic shape (the
  shape-noise null models it) plus model-dependent shear at this position. Not sent to `/vet-candidate`; not a
  candidate. Like WHL0137 2194 (bright, near the 60° cut), unlike the faint, high-shear S1063 anti.
- Lenstool-only (control-failing run, listed for completeness): **1143** (m 24.4, z 1.0) a compact source with a
  crossing feature; **952** (m 27.6, z 3.4) barely above the noise at this stretch.
- No z_phot > 15 row among the anti (pre-registered check).
- Post hoc, from the sheet: Lenstool row 1224 lies on the mosaic edge (half the stamp blank); it is `aligned`.

## Limits
- 12 (GLAFIC) strong-shear background rows: only an excess of a few anti is detectable. Both models are HST-era;
  Fujimoto et al. 2024 (arXiv:2404.03286) give a newer model, not used (maps not found public on 2026-10-18).
- One F200W epoch; one photo-z catalogue (DJA v7.2, which predates program 6796).
