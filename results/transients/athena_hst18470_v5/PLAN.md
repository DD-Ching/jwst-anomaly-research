# SN Athena reappearance in HST 18470 visit 5: plan (written before measuring, 2026-11-08)

## Question
Has the next image of SN Athena (MACS J0417 h002, D-034 / #323) appeared in the first HST 18470 visit not yet
checked (visit 5, `ifu905010`, WFC3/IR F105W, MJD 61341.56 = 2026-10-28)?

## Predictions (model_prediction)
- Discovery team (JWST DD 12781, three JWST-era lens models): reappearance mid-Aug 2026 to early Mar 2027, at the
  12781 target "Athena-2" (64.394375, −11.90758).
- Our single pre-JWST RELICS glafic model (#309): the 7.6″ image (64.394334, −11.907649; 0.29″ from Athena-2) lags
  708–778 d → 2027-09-14 to 2027-11-22. A reappearance now would put the RELICS glafic delay ≳ 1 yr too long
  (an ordinary model error, not an anomaly).

## Ordinary explanations for any change
Lensed-galaxy light ~1″ away leaking into the annulus, persistence / cosmic rays / detector artefacts, registration
differences between drizzled visits, and the visit-to-visit scatter already seen in visits 1–4 (+1.4σ maxima).

## Method (frozen; same as `results/transients/confirm_r1/macs0417_hst18470/hst_check.py.txt`)
- F105W `_drz.fits` of visits 1 and 5; 0.2″-radius aperture, 0.6–1.0″ annulus; errors from 300 random apertures
  (robust σ); both positions above.
- Difference error σ_d = √(σ₁² + σ₅²).
- **Rule (ASSUMPTION):** `appeared` if Δ = F₅ − F₁ ≥ 5 σ_d at either position **and** a point source is visible at
  the position in the visit-5 stamp (looked at); `not seen` if Δ < 3 σ_d; otherwise `undetermined`.
  Expected brightness if it has appeared: ~200–330 nJy (F150W-based estimate from #323; F105W may differ), i.e.
  ~5–8 σ_d.

## What would change the plan
`appeared` → record that the discovery team's window holds and our RELICS glafic window is too late (D-034
addendum); no further action from us (followed by its discoverers). `not seen` → the window stays open to early Mar
2027; repeat at each new public visit.
