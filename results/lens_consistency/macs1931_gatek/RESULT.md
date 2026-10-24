# MACS J1931: unit check U and gate K pass on the VENUS AstroLensPy v2.1 maps (2026-10-24)

Plan `PLAN.md` (frozen in b7ac41b, before any map value was read). Script `scripts/macs1931_gatek.py`; tests
`tests/test_macs1931_gatek.py`. Numbers from `summary.json` (plan) and `posthoc.json` (post hoc). Every image-plane
and source-plane number is a model_prediction; the image list is observed.

## Result (frozen rules)
- **α normalisation check passes** (before the maps were read): (D_LS/D_S)(z) / (D_LS/D_S)(5.079) is 0.890432 /
  0.998854 / 1.032566 at z_s 2 / 5 / 9, against the authors' `k_dlds` 0.89043 / 0.99885 / 1.03257. So the maps are
  normalised at z_s = 5.079 and are divided by (D_LS/D_S)(5.079) to reach D_LS/D_S = 1.
- **Unit check U passes; the header unit is the right one.** Median distance of our mean back-traced sources to the
  authors' `sourcepos.txt` Bary sources over the 9 spec-z systems: **0.110″** with the map values read as pixels
  (× 0.04″; header `BUNIT = 'Pixels'`), against 363″ read as arcsec (README). Per-system distances with the pixel
  reading: 0.023–0.265″. The README's "in arcsec" is wrong for these files.
- **Gate K passes:** 26 spec-z images at their system z, 0.25″ grid, mean family source: rms **0.919″** (rule
  ≤ 1.0″), **100 %** within 3″ (rule ≥ 80 %). Largest per-system rms: system 2 (1.76″, 3 images) and system 4
  (1.17″).
- No early deviation: the only code change after the freeze and before the first gate output is that the maps are
  read from their image extension (the primary HDU is empty), committed before the run.

## Reported, not gates (fixed in the plan)
- With the D-041 WLS family source: rms **0.630″**, 26 / 26 within 3″.
- Over the 44 used images with a non-zero catalogue z: rms 4.83″, 84.1 % within 3″ (it would fail gate K). Systems
  8 (rms 15.8″) and 10 (5.08″) get one predicted image for three and two catalogued; 15 (9.16″) is one of the five
  systems with different z per image (11, 13, 15, 16, 17), which our solver runs at the first image's z.
- Spec-z systems with more predicted than catalogued images: 1, 5, 6, 7, 100, 101 (none with fewer). They are the
  targets of a later counter-image test; none is read here.
- Authors' `chires.txt` Bary rms: 0.444″ over all 51 used images; 0.606″ over the 26 spec-z images (post hoc, below).

## Post hoc (`posthoc.json`; not in the plan; does not change the gate outcome)
- U distances for the photo-z / model-optimised systems with the pixel reading: system 8 **3.18″**, 11 2.41″,
  others ≤ 0.88″ (system 10: 0.10″). So for system 8 our source lies 3″ from the authors' source at the catalogue z
  1.08688. A redshift that differs from the one the authors used for this map is one possible cause; it is not
  tested here. System 10 (z 0.66) matches the authors' source yet gets one predicted image, so its miss is on the
  image-plane solve side (cause not isolated).
- Per-image comparison: the authors' Bary residuals on the 26 spec-z images have rms 0.606″, ours 0.919″ with the
  mean source and 0.630″ with the WLS source (summary.json). Their residual definition was not reproduced.

## Reading
- The VENUS maps, read in pixels and normalised at z_s 5.079, reproduce the spec-z images the model was fitted to
  (own-constraint validation only). The maps can be used with the spec-z systems; the photo-z systems 8, 10 and the
  mixed-z systems are not reproduced as catalogued and are excluded from any later test until their z is settled.
- No deblended reference catalogue is public for this field (D-034 addendum aa), so gate C2 and the extra-image test
  still cannot run. Next: a counter-image test for the six spec-z systems with extra predicted images needs
  imaging; the VENUS mosaics on MAST could not be checked from this session before (connection failure).
- Data: the two maps (67.7 / 67.0 MB zipped) were held in the session scratchpad only and deleted after the run.
