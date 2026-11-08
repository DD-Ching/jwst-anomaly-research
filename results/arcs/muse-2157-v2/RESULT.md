# `arcs` anti row MACS0308 2157: MUSE redshift, second plan — gate M2 passed; 2157 `undetermined`; no [O II] F90 within the grid

Plan with amendment 1 (written before any data read): `PLAN.md`. Script: `scripts/arcs_muse_2157_v2.py` (reuses the
reader and line functions of `scripts/arcs_muse_2157.py`; new functions unit-tested in
`tests/test_arcs_muse_2157_v2.py` before the run). Numbers: `summary.json`; spectra, both noise estimates and scores:
`spectra.json.gz`; contact sheet `contact_sheet.png` (looked at). Cube: ESO MUSE `ADP.2020-01-04T23:56:39.902`, DATA
and STAT extensions streamed once, two 80 × 80-spaxel boxes kept in memory, nothing stored (Cloud disk decision).

## Gate M2 (validation, central galaxy): **passed**
- Aperture recentred by 0.47″; z_best 0.3568 (within 0.01 of 0.356), peak contrast 1.875 ≥ 1.5 (ASSUMPTION).
- Features at z_best (derived): Mg b 5.5, Na D 4.9, G band 4.5, Ca K 3.9 σ; continuum S/N 10.6 per pixel.
- By eye (not stored): the recentred aperture sits on the edge of the brighter of two close knots in the white
  light, not on its peak; the score curve has one clear spike at 0.357.

## Target 2157 (observed / derived)
- Recentred by 0.22″; z_best 0.2046, peak contrast 1.041 < 1.5 → frozen class **`undetermined`**.
- No feature above |S/N| 2.6 at the best z, and none above 1.8 at z 0.356 (`features_at_z_cluster`); continuum S/N
  2.2 per pixel; 94 sky-masked pixels. Noise used is max(empirical, STAT); empirical / STAT median ratio 1.34
  (control 2.06), so STAT alone under-states the aperture noise.
- By eye (not stored): no distinct source peak at the aperture in the white light, which shows a faint vertical
  stripe through that column; the 0.22″ recentre may follow that stripe rather than the galaxy.

## [O II] injection (simulated; amendment 1 criterion: [O II] window S/N ≥ 5 at z_inj)
- 40 random z in 0.40–1.45 per level; detected fraction 0 / 0 / 0 / 0 / 0.15 / 0.825 at 50 / 100 / 200 / 400 / 800 /
  1600 × 10⁻²⁰ erg s⁻¹ cm⁻² (aperture flux, no aperture correction).
- **F90 not reached within the grid**: the narrowest statement is that an [O II] doublet of 1.6 × 10⁻¹⁷ erg s⁻¹ cm⁻²
  in the 0.6″ aperture is detected in 82.5 % of the trials. A fainter [O II] line at 2157 is not excluded.
- Diagnostic: the scan identifies the injected z in at most 5 % of trials (single lines are z-degenerate here).

## Reading
- The tool now passes its validation on the cluster galaxy, so the cube and scan can return an absorption redshift
  at continuum S/N ≈ 10; at 2157's S/N ≈ 2 they return none. Membership stays undecided; 2157 stays an open anti row.
- The MUSE route for 2157 is closed at this depth: neither the absorption scan nor an [O II] line (to the level
  above) gives a redshift. No flag, no candidate.

## Limitations
- Aperture flux only (0.6″); the [O II] limit is not total flux. Gaussian σ 1.1 Å, 1:1 doublet, window and S/N 5
  are ASSUMPTIONs. 40 trials per level (binomial error about ±6 % at 82.5 %, derived, not stored).
- The contrast threshold 1.5 was chosen with the first run's readings known (control 1.79, target 1.07): not blind.
- The flux grid stopped at 1600; a higher level would locate F90 (deviation from nothing; the grid was frozen).

Not a lensing anomaly claim and not a candidate.
