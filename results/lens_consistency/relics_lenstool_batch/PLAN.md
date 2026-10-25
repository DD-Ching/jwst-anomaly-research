# RELICS Lenstool v1 2017-12 map batch: does the amplitude offset follow it? (H_K; frozen 2026-10-25)

Frozen before any deflection map of this run was read. Follows D-034 addendum bi
(`../relics_lenstool_zref/RESULT.md`, "Prediction for the next test").

## Hypothesis and ordinary explanations
- **H_K (hypothesis):** the RELICS Lenstool v1 deflection amplitude offset (rule-7 factor 0.665 MACS J0553, 0.692
  Abell 3192; SMACS 0723 |α| 1.14–1.56 × GLAFIC) belongs to the map batch with header `MAPMKR` of 2017-12-03…06.
  Then the other maps of that batch (RXC J0600, MACS J0417; DATE 2017-12-05/06) are offset too, and the
  maps of other makers (WHL0137, SPT0615) are not.
- **Ordinary alternatives:** (a) the offset is a property of each cluster's mass model (H_M, per cluster, no batch
  pattern); (b) the second model, not Lenstool, is off; (c) the statistic below measures mass-profile differences,
  not a normalisation.
- **What changes the plan:** H_K supported → a RELICS-wide rule "check the 2017-12 batch maps by default" goes into
  the addendum-bd check, and the four untested batch maps are flagged before any use. H_K rejected → the per-cluster
  check stays the only rule and the batch line closes.

## Data (observed: published maps; derived: everything computed)
Deflection maps at D_LS/D_S = 1 (arcsec), registered in `scripts/lens_consistency.py` / `jwst_anomaly.lensmodel`:

| Cluster | Role | Lenstool v1 | Second model |
|---|---|---|---|
| MACS J0553 | validation, offset (f 0.665) | RELICS | RELICS GLAFIC v2 |
| Abell 3192 | validation, offset (f 0.692) | RELICS | RELICS GLAFIC v2 |
| RXC J2211 | validation, near 1 (f 1.017) | RELICS | RELICS GLAFIC v2 |
| RX J0152 | validation, near 1 (f 1.045) | RELICS | RELICS GLAFIC v2 |
| RXC J0600 | test, batch | RELICS | RELICS GLAFIC v3 |
| MACS J0417 | test, batch | RELICS | CANUCS DR1 Lenstool (JWST-era, 30 spec-z systems) |
| WHL0137 | test, other maker | RELICS | RELICS GLAFIC v2 |
| SPT0615 | test, other maker | RELICS | RELICS GLAFIC v2 |

RXC J0018 (batch-maker header, 2020) has no second model on MAST: untestable here, stated as such. Maps are streamed
into memory (no disk copy; transfer ≈ 0.8 GB, every single file < 200 MB) and discarded.

## Statistic (frozen)
1. Sample both deflection fields on a 1″ sky grid over the overlap of the two map footprints (bilinear in each map's
   own WCS); magnitudes |α| only, so sign and axis conventions cannot enter.
2. Strong-lensing region S: grid points where |α_ref| ≥ 0.5 × max(|α_ref|) over the overlap (ASSUMPTION: selects
   the annulus around the Einstein radius where multiple images form).
3. R = median over S of |α_LT| / |α_ref|. Also reported: ratio of means over S, the 16–84 % range of the pointwise
   ratio, and |S| in arcsec².
4. Classes (ASSUMPTION): **offset** R ≥ 1.25; **near 1** 0.85 ≤ R ≤ 1.15; else **intermediate**.

## Gates and outcomes
- **Gate V (validation on known cases):** all four validation clusters fall in their known class, and |1/R − f_obs|
  ≤ 0.12 for each. If V fails, the statistic is not trusted and H_K is **undetermined** (no reading of the test
  clusters beyond reporting the numbers).
- **If V passes:**
  - H_K **supported**: RXC J0600 and MACS J0417 both offset, and WHL0137 and SPT0615 both near 1.
  - H_K **rejected**: any batch map near 1, or any other-maker map offset.
  - otherwise (an intermediate class): **undetermined** for that cluster; the rest read as above.
- Caveat fixed in advance: MACS J0417's second model is JWST-era with more constraints, so a mass-profile difference
  there is likelier than for the RELICS-era pairs; its R is read with that label.
