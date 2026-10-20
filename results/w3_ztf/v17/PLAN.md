# W3 ZTF v17 — evidence check of the two v16 final survivors (post hoc)

Plan fixed 2026-10-20 before any recomputation. Owner authorization: after #385, 2 h, ≤ 200 MB. **Post hoc:**
the v16 image results have been seen; this is a follow-up of a known measurement conflict, not a blind test. Scope:
708205400009336 (Npix 1624, r) and 709212400017577 (Npix 1625, g) only. No search, no threshold change, no new veto
rule, no new PSF solver. Nothing here changes the frozen v13 rule or the v16 counts.

## 1. Evidence file (`evidence.json`, from committed v16 files only)
- Member objectids, DR24 position, run window and every vetting test of the final run (`search_holdout_v16_time.json`).
- Per checked exposure: the v16 automated class, 2″ aperture flux and error, field-star ratio, infobits, maglimit,
  seeing, mask bits at the source, dark-column distance (`image_check.json`, `mask_check.json`).
- The visual reading of `image_check.png` is kept as a separate field and never replaces an automated class.
- 1624's "≈ 2,450 µJy" is cited with its 4 bracket detections (2,446–2,560 µJy); its v16 `reference_uJy` is NaN.

## 2. Checks, in this order (`check_script.py.txt` → `check.json`; criteria are code, unit-asserted at start)
New data: IBE science and official mask cutouts (180″) of the same exposures plus 1624's unchecked ND, and the
official per-exposure PSF catalogue (`_psfcat.fits`) of each exposure. Caps: ≤ 40 cutouts (science + mask), ≤ 20
psfcat files, ≤ 200 MB; no new download once 180 MB is passed. No light-curve partition is read.

1. **Consistency:** the files are `_sciimg` (not difference images); flux is µJy from each header's `MAGZP` (AB 23.9);
   the DR24 position falls inside each image; the flux-weighted centroid within 3″ of the DR24 position in the bracket
   detections lies ≤ 1″ from it.
2. **Aperture, background, mask:** the existing `hit_lightcurve.aperture_photometry`, radius r = max(2″, seeing FWHM),
   annulus 2r to 2r + 5″. Per exposure: the target centroid offset, the annulus median, and the official-mask bits
   in the target aperture. **Bad bits:** 2 LOW RESPONSIVITY, 8 SATURATED, 9 DEAD, 10 NAN. For an exposure with bad
   bits in the aperture, a **column-interpolated** flux replaces the bad pixels by the mean of the nearest good pixels
   in the same row. The same pixel offsets are masked and interpolated in the bracket detections, so the ratio is
   self-calibrated.
3. **Comparison sources (same exposure):**
   - stars found in the first bracket image (the existing peak finder), with bracket flux 0.25–4× the target's;
   - ≥ 15″ from the target and inside the cutout margin, with no bad bits in their aperture in either image;
   - scale s_e = median of control flux ratios (exposure / median bracket). Stable only with ≥ 3 controls and a ratio
     MAD / median ≤ 0.15.
   - Target ratio R_e = (target flux / s_e) / target median bracket flux. σ_R combines the random-aperture noise
     (existing 200-aperture MAD), the bracket scatter, and the control MAD / √n.
4. **Depth:** an exposure can distinguish the frozen dimming (DIM_MAG 1.0, i.e. ratio 0.398) only if
   (1 − 0.398) / σ_R ≥ 3.
5. **psfcat (dependency stated):** the nearest psfcat source within 1.5″ (`psfcat_match`), with mag + MAGZP, sigmag and
   flags. DR24 light curves are built from these same PSF catalogues, so this explains why DR24 has no record; it is
   **not** an independent observation. A missing psfcat row is never zero flux. PSF magnitudes are compared only with
   PSF magnitudes, never with aperture fluxes.

## 3. Per-exposure classes (computable)
- **Judgeable:** infobits 0, stable controls, no bad bits in the target aperture, and the depth condition. Otherwise
  the class is **unknown**, with every failing reason listed.
- **normal:** judgeable, R ≥ 0.75 and (1 − R) / σ_R < 3.
- **dim:** judgeable, R ≤ 0.6 and (1 − R) / σ_R ≥ 3.
- **undetermined:** judgeable, but neither normal nor dim.
- **For an unknown exposure with bad bits only:**
  - **data problem supported** when the column-interpolated R_i ≥ 0.75 and the psfcat row is missing or flagged;
  - **undetermined** when R_i < 0.75.

## 4. Per-source answer
- **Real dimming evidence:** ≥ 2 dim exposures on ≥ 2 nights inside the run, and ≥ 1 normal exposure before and
  ≥ 1 after.
- **Data problem supported:** every checked ND is normal or data problem supported, and none is dim or
  undetermined.
- **Otherwise undetermined:** name the missing data.
- **Each source also answers:**
  - Is a brightness drop confirmed in the images?
  - Does a specific data problem explain the flag?
  - Which exposures are judgeable, and which are unknown?
- **If a real dimming is confirmed:** the remaining time goes to ordinary variables, occultation or known events.
  No exotic language.
