# W3 ZTF v18 — PSF photometry check of the two undetermined v17 sources (post hoc)

Plan fixed 2026-10-20 after the data fetch and before any target PSF fit. Owner authorization: after #387, 3 h,
≤ 200 MB. Only 708205400009336 (Npix 1624, r) and 709212400017577 (Npix 1625, g).

**Post hoc:** these are known items whose v16 and v17 results have been seen. This is not a blind test. No search,
no change to the v13 rule, no new veto.

## 0. Tools and products (checked; `manifest.json`)
- **Products:** the official per-exposure DAOPHOT PSF `*_sciimgdao.psf` (GAUSSIAN plus 47×47 lookup tables with
  linear spatial terms) exists for all 20 exposures.
  - `daopsf.py.txt` is a reader adapter only. It evaluates the stored model as DAOPHOT II USEPSF does, at the
    target's own quadrant position, never by assuming the centre PSF.
  - **Reader check:** rendered at the quadrant centre, the model matches the official `*_sciimgdaopsfcent.fits` to
    ≤ 5.24e-4 (peak ≈ 0.19) in all 20 exposures. A 0.1 px shift gives ≈ 0.013, so the centring convention is
    confirmed.
- **Fitter:** photutils 3.1.1 `PSFPhotometry` with an `ImagePSF` built from that model (oversampling 4), the
  official mask, a pixel error array and fixed positions. No new solver.
- **Data:** 600″ science and mask cutouts, psfcat and PSF files for each exposure, 48.4 MB.

## 1. Common measurement (both sources, controls and targets)
- **Background:** photutils `Background2D` (box 50 px, 3σ clip) on the cutout, with bad pixels and psfcat sources
  masked. Pixel error = sqrt(bkg_rms² + max(data − bkg, 0) / GAIN).
- **Bad pixels (official mask):** bits 0, 2–10 and 12 are excluded from every fit. Bits 1 and 11 only mark
  detections.
- **Positions:**
  - fixed from RA/Dec: the mean psfcat position over the bracket detections;
  - each exposure is registered by the median WCS offset of its isolated psfcat stars.
  - Fluxes are free.
- **Fit box:** 9×9 px (1624) and 11×11 px per source (1625). 1625's two stars form one group and are fitted jointly.
- **PSF:**
  - **Targets:** the exposure's own model at the target's quadrant position.
  - **Controls:** the model of a *different* exposure of the same quadrant (the nearest in time), so no control is
    fitted with a PSF built from its own pixels.
  - **214647008** (the only exposure in field 1748) uses its own model; flagged.
- **Information condition:** the PSF fraction on valid pixels in the fit box must be ≥ 0.30. Otherwise the exposure
  is **unknown**, and no extrapolated flux is given.
- **Scale:** each flux is divided by the v17 D1 psfcat comparison scale of its exposure (`../v17/followup.json`).

## 2. Validation (must pass before a target exposure is classified)

### V1, bad columns (1624), per ND geometry g
Geometry g is the bad pixels of that ND relative to the target position.

- **V1a, real stars.**
  - **Stars:** isolated psfcat stars (no other source within 15″, flags 0, off their own bad pixels).
    - Class A: |m − 15.14| ≤ 0.5.
    - Class B: |m − 16.14| ≤ 0.5, i.e. 1 mag fainter.
  - **Measurement:** in all 13 exposures, the masked flux with geometry g removed, over the full flux. Ratio
    r = F_mask / F_full.
  - **Passes when**, for A and B separately: |median r − 1| ≤ 0.05, robust σ(r) ≤ 0.20, and ≥ 95 % of
    |r − 1| ≤ 3 σ_tot,rel.
  - **σ_tot,rel:** the empirical robust σ(r) combined with each fit's formal error.
- **V1b, real column injection.**
  - **Placement:** stars of target flux (the median bracket flux) and 0.398× that flux, injected into each ND image
    with a *different* exposure's PSF model. Each goes at 6 positions along the same bad columns, at the target's
    offset from the column, ≥ 15″ from any psfcat source.
  - **Fit:** with geometry g masked.
  - **Passes when:** median recovered / true is within 1 ± 0.05 for both levels. The 1.0 level is classed normal in
    ≥ 90 % of trials and the 0.398 level dim in ≥ 90 % (class rules in section 3, σ from V1a).

### V2, blend (1625), per exposure (3 NDs and 4 detections)
- **Stamps:** real isolated psfcat stars with g 17.7–18.7 in the same exposure. Each is cut out after background
  subtraction, scaled to the target flux and to 0.398× it.
- **Injection:** at the target's offset from the real 6.6″ neighbour, rotated about it by 90°, 180° and 270°
  (positions with no psfcat source within 4″).
- **Fit:** jointly with the neighbour. Truth = scale × the stamp's own isolated fit in that exposure.
- **Passes when**, for both levels:
  - |median(recovered / true) − 1| ≤ 0.10 and robust σ ≤ 0.20;
  - normal trials classed normal and dimmed trials classed dim in ≥ 90 %;
  - the neighbour's fitted flux changes by ≤ 5 % between the two levels (degeneracy check).
- **Also stored:** the joint-fit residual in a ring around the neighbour (systematic wing residual).

### Not allowed
Choosing a model, box, mask or reference by whether the target comes back normal. Every setting above is fixed now.

## 3. Target classes per exposure (as v17, now with PSF fluxes)
- **Ratio:** R = F_e / F_ref, with F_ref the median of the 4 bracket fluxes.
- **σ_R combines:**
  - the fit error (background and Poisson);
  - the validation scatter of that exposure (V1a or V2, relative);
  - the bracket scatter;
  - a PSF-mismatch term |F_own − F_other| / F_ref (the target fitted with a different exposure's model).
- **Classes:**
  - **unknown:** the exposure's validation failed, or the information condition failed;
  - **normal:** R ≥ 0.75 and (1 − R) / σ_R < 3;
  - **dim:** R ≤ 0.6 and (1 − R) / σ_R ≥ 3;
  - **undetermined:** otherwise.
- **The 1-mag drop:**
  - **excluded** in an exposure when (R − 0.398) / σ_R ≥ 3;
  - **supported** when the exposure is dim.
  - Excluding a 1-mag drop does not mean no variation.
- **Per source:** as v17 PLAN section 4. Real dimming evidence needs ≥ 2 dim exposures on ≥ 2 nights inside the
  run, plus normal exposures before and after. Judgeable and unknown exposures are reported separately.
- **If real dimming is confirmed:** the rest of the budget goes to ordinary explanations. No exotic language.

## 4. End of this follow-up
- **If a validation fails:** the failing geometry or exposure stays unknown and the case is saved. No other method
  is tried on these two sources in this round.
- **Reopening needs** new data or an already-validated method.
- **Kept for reuse** (no re-download): the manifest with URLs, small stamps (sci, mask) around each source, PSF
  files and control lists (`stamps.npz`, ≤ 1 MB).
