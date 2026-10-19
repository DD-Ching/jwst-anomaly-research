# W3 ZTF v9 — official single-exposure PSF catalogues fill the DR24 gap for 3 calibration stars

- Plan: `PLAN.md`, committed in 601e927 before any psfcat download.
- Scope: no new photometry code, no relaxed quality rule, no search, no injection re-run. Labels are `observed` for
  psfcat, the IRSA light-curve API and TAP, and `derived` for counts and comparisons.
- Conclusions are limited to field 672, CCD 8, quadrant 2 in Npix 7134 and to these 3 stars.

## 1. The product exists and is readable (`psfcat_overlap.json`, `psfcat_fill.json`)
- All 72 psfcat files requested were readable: 12 for the test and 60 for the fill, on the official IBE path
  built from the TAP metadata.
- Each file holds 13 columns: sourceid, x/y, ra/dec, flux and sigflux in DN, mag, sigmag, snr, chi, sharp and flags.
  The header gives MAGZP, CLRCOEFF and PCOLOR (g−R) and INFOBITS.
- **Magnitude definition, checked in the data:** the `mag` column is instrumental. In every file, mag − (−2.5 log₁₀
  flux + MAGZP) equals −MAGZP, so the calibrated value is mag + MAGZP.
  - **Deviation (labelled):** the first overlap run omitted MAGZP and gave a median Δ of −26.10 mag; its output is
    kept as `psfcat_overlap_run1_no_zeropoint.json`. The plan already required the zero point, so the fix applies the
    plan's own definition. No coefficient was fitted.
  - No colour term was applied. DR24 light curves also assume zero colour (release notes).

## 2. Overlap test (fixed criteria, PLAN §3): **passed**
- 13 pairs (star × overlap exposure, both clean).
- Δ = psfcat − DR24: median **−0.015 mag** (limit ±0.05); robust scatter **0.0013 mag**, against a limit of 0.072 mag
  (2 × the median combined error).
- The −0.015 mag offset is within the light-curve bias corrections the release notes describe (up to about 0.025
  mag). It was not removed.

## 3. Gap fill (PLAN §4)
- **Candidates:** gap exposures of the quadrant (2022-10 to 2023-10, infobits 0, at least one star inside its
  footprint and without a clean DR24 detection). There were more than 60; the fill took 60 evenly spaced ones, so
  some remain unfetched. Download: 12.8 MB.

| star | clean DR24 detections in the gap | usable psfcat added (unique match, flags 0) | total in the gap after fill | missing | flagged |
|---|---|---|---|---|---|
| 672208200000073 | 2 | 28 | 30 | 30 | 2 |
| 672208200000239 | 2 | 60 | 62 | 0 | 0 |
| 672208200002085 | 1 | 31 | 32 | 25 | 4 |

- **All 55 missing targets** fall outside the exposure footprint or inside the 10″ margin (`psfcat_missing.json`: 0
  inside the footprint). They are coverage gaps, not non-detections, and none is read as a disappearance or limit.
- **0 ambiguous matches.** Matched separations are ≤ 0.68″ in the fill (≤ 0.36″ in the 6 test gap exposures).
- **Baseline:** each star now has ≥ 20 usable observations inside the gap (30, 62 and 32). This is indicative only. The
  frozen pre-screen counts its 20-epoch bracket per event, and psfcat lists detections only.

## 4. ZFPS (PLAN §5)
- Not needed for these 3 stars, because psfcat filled the gap. No request package was prepared.
- Difference images, which the service would need, exist for 5 of the 6 test gap exposures (HTTP 200, 7.8–8.5 MB).
  The sixth check hit a connection reset and was not retried (`diffimg_check.json`).
- No account was used or requested.

## 5. Cost (`cost_estimate.json`)
- **Measured:** 214 kB and 1.5 s per psfcat file, including the cross-match. One file serves every star in its
  quadrant.
- **Estimate for the 49 #350 bracket cases (one cell):** at most 533 files, about 114 MB and 13 min. That is an upper
  bound; carriers that share quadrants need fewer files.
- **Downloads this round:**
  - psfcat: 2.44 MB (test) + 12.83 MB (fill);
  - TAP: 0.15 + 0.15 + 0.11 MB;
  - LC API: 0.33 + 0.33 MB.

## 6. What this does and does not settle
- It shows that an official product holds calibrated detections for exposures missing from DR24 light curves, at
  least in this quadrant. Those detections agree with DR24 on shared exposures to 0.015 mag (offset) and 0.0013 mag
  (scatter).
- It does **not** re-run the calibration. Using psfcat in the W3 chain would need a new data path: series and
  comparison counts built from psfcat for the gap exposures. That needs a plan and the owner's authorization. Whether
  the 10 cells would then reach 0.5 is unknown. 2024 coverage stays sparse, and pre-event shortfalls remain.
