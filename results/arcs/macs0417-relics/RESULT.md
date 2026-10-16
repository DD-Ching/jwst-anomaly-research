# `arcs` on MACS0417: result (2026-10-16)

Plan: `PLAN.md` (written before run 1; addendum A written after run 1's control failure and before run 2).
Labels: shapes `observed`; predicted directions `model_prediction`; classes, counts and P values `derived`. All
thresholds are the D-024 defaults (ASSUMPTIONs).

## Inputs
- F200W pipeline catalogue `jw01208-o006_t003_nircam_clear-f200w_cat.ecsv` (CANUCS, 6,399 s; MAST level 3;
  5,030,481 B, sha256 `e55dffb9…6ee57d90dd`). The MAST position query (64.3945, −11.9089, r 0.02°) found no other
  public F200W level-3 pointing on the cluster.
- CANUCS DR1 `macs0417-clu` photometry (27,820,129 B, sha256 `ca6b48aa…6203be73d7e`, 14,272 rows), converted to
  zout columns with `inject_radial.canucs_zout` (Z_SPEC where > 0, else Z_ML / Z160 / Z840).
- Models (pinned in `lensmodel.py`, `lens_consistency.py` `MODELS`):
  - run 1 `macs0417-relics-lenstool`: RELICS Lenstool v1 maps (Mahler et al. 2019), z_lens 0.443. `validate`:
    κ median relative difference 2.6 × 10⁻⁵, z = 2 |μ| 1.4 × 10⁻⁵ (both reproduce the published maps);
  - run 2 `macs0417-canucs-maps`: CANUCS DR1 Lenstool maps (Desprez, preliminary v1, 30 spectroscopic systems),
    100 mas. `validate`: κ from the deflections vs the published κ, median relative difference 3.1 × 10⁻⁵.

## Results (`run1-relics-v1/arcs.json`, `run2-canucs-v1/arcs.json`)
| Run | Strong-shear background | Aligned (control ≥ 70 %) | Anti | Expected (shape noise) | P(≥ anti) | Frozen outcome |
|---|---|---|---|---|---|---|
| 1 RELICS v1 (HST-era) | 61 | 39 (64 %) | 8 | 1.42 | 1.2 × 10⁻⁴ | **control fails: not interpreted** |
| 2 CANUCS v1 (JWST-era) | 49 | 37 (76 %) | 3 | 1.52 | 0.195 | **does not replicate** the S1063 excess |

- Run 1 by g_min: g < 1 26 / 36 aligned, 2 anti; g ≥ 1 13 / 25 aligned, 6 anti. Run 2: g < 1 30 / 33 aligned,
  1 anti; g ≥ 1 7 / 16 aligned, 2 anti. The RELICS v1 misalignment sits mostly near its critical curves, where an
  HST-era model's errors turn the predicted axis most (an interpretation, not tested).
- Run 2 anti rows (cutouts looked at, F200W 4″ with observed and predicted axes; not stored):
  - 1993 (g_min 3.7, m 27.2, z_phot 3.0): a faint, low-surface-brightness smudge;
  - 2341 (g_min 0.26, m 27.1, z_phot 3.2): a compact knot with a faint extension;
  - 4071 (g_min 1.2, m 27.4, z_phot 8.1): on a vertical detector stripe.
  The three aligned controls on the same sheet follow their predicted axes. The anti rows are not vetted further:
  the frozen rule sends anti rows to vetting only on replication.

## Reading
- In the one independent field where the model passes the control, the background anti count is consistent with
  shape noise (3 vs 1.52). The S1063 excess (4 vs 0.52, #273; P 0.031 clump-aware, #279) stays a single-field,
  exploratory result.
- **Deviation:** run 2's model was chosen after run 1 failed its control (addendum A). Its outcome is
  "does not replicate", so the choice cannot have produced a false replication.
- **Limits:** one field; the CANUCS model is a preliminary release without a published image list, so there is no
  known-case check beyond the background alignment; the anti counts use photo-z for the background cut
  (CANUCS blends with members are not tested here); shape-noise null only (no clump-aware row, since the outcome
  did not call for one); no frame offset between either model and the JWST catalogue was measured (offsets of
  0.1–0.3″ matter only near critical curves, D-034).

## Data handling
The three CANUCS maps (85–89 MB each, 263.5 MB in total; under 200 MB per file) were read into the session's
gitignored cache and deleted after the run. Nothing is stored in git except these tables.
