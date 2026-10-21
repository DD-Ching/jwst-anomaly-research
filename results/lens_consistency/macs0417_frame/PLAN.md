# MACS0417: Jauzac et al. 2019 image list → JWST frame offset, then the extra-image forced photometry again (plan, frozen)

Frozen 2026-10-21 23:20Z before any offset was measured. Follows the "Next (optional)" of
`../macs0417_extra/RESULT.md`: there gate C failed (55 / 55 `no_reference`) and the F200W peaks sat a median
0.30″ from the published positions (post hoc, dy −0.27″ coherent). Script: `scripts/m0417_frame.py`.

## Hypothesis and ordinary explanations
- H (ASSUMPTION to test): the Jauzac et al. 2019 list is on the RELICS HST frame, which is offset from the JWST
  1208 pipeline mosaic frame; correcting it moves the catalogued positions onto their F200W peaks.
- Ordinary alternatives: the list is on another frame (then the RELICS offset does not shrink the residuals); the
  residuals come from extended arcs whose peaks need not sit at the catalogued knot (then the offset helps little and
  gate C still fails on the reference rule).
- What would change the plan: F2 fails → no shifted run; the MACS0417 test stays waiting on CANUCS v2. Gate C still
  fails after the shift → the D-034 compact-reference route is not valid for MACS0417 arcs; next would be an
  arc-aware reference (`halo_phot` totals), not more frame work. A flag with gate C passing → `/vet-candidate`.

## Inputs
- observed: RELICS macs0417m11 HST catalogue v1 (`hlsp_relics_hst_acs-wfc3ir_macs0417m11_multi_v1_cat.txt`, MAST
  HLSP); JWST 1208 o006_t003 F200W pipeline `_cat.ecsv` and `_i2d.fits` (as `../macs0417_extra`).
- unchanged: `macs0417-canucs-maps`, `../macs0417_images/jauzac19_images.ecsv`, `m0417_extra.py` rules.

## Method and frozen rules (ASSUMPTIONs)
1. Offset O = (ΔRA cos δ, ΔDec), JWST minus RELICS, from RELICS sources with `fwhm` ≤ 0.5″ and `nf5sig` ≥ 3
   matched to JWST `is_extended` False or True sources: pass 1 nearest match ≤ 1.0″, median; pass 2 re-match
   ≤ 0.3″ around the pass-1-shifted position, median = O. Error per axis 1.2533 × 1.4826 MAD / √N.
2. Gate F1 (measurement): N ≥ 30, error ≤ 0.02″ per axis, and the two halves (even / odd RELICS id) agree within
   0.03″ per axis.
3. Gate F2 (does O apply to the Jauzac list?): with the list shifted by O, the median |F200W peak − position| of the
   56 images (same diagnostic as `m0417_extra_offsets.py`) is ≤ 0.18″ (0.6 × the unshifted 0.300″).
4. If F1 and F2 pass: `m0417_extra.run_field("macs0417")` with the shifted list, every other rule unchanged (gate C ≥ 0.8
   of testable controls recovered or confused; flag = targeted predicted image `absent`). Outputs in this directory.
