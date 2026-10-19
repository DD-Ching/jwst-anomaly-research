# W3 ZTF v9 — official single-exposure PSF catalogues fill the DR24 gap for 3 calibration stars

- Plan: `PLAN.md`, committed in 601e927 before any psfcat download.
- Scope: no new photometry code, no relaxed quality rule, no search, no injection re-run. Labels are `observed` for
  psfcat, the IRSA light-curve API and TAP, and `derived` for counts and comparisons.
- Conclusions are limited to field 672, CCD 8, quadrant 2 in Npix 7134 and to these 3 stars.

## 1. The product exists and is readable (`psfcat_overlap.json`, `psfcat_fill.json`)
- All 66 distinct psfcat files requested were readable, on the official IBE path built from the TAP metadata. Twelve
  were for the test and 60 for the fill; 6 of the fill files were reused from the test cache.
- Each file holds 13 columns: sourceid, x/y, ra/dec, flux and sigflux in DN, mag, sigmag, snr, chi, sharp and flags.
  The header gives MAGZP, CLRCOEFF and PCOLOR (g−R) and INFOBITS.
- **Magnitude definition, checked in the data:** the `mag` column is instrumental. In every file, mag − (−2.5 log₁₀
  flux + MAGZP) equals −MAGZP, so the calibrated value is mag + MAGZP.
  - **Deviation (labelled):** the first overlap run omitted MAGZP and gave a median Δ of −26.10 mag; its output is
    kept as `psfcat_overlap_run1_no_zeropoint.json`. The plan already required the zero point, so the fix applies the
    plan's own definition. No coefficient was fitted.
  - **Deviation (labelled):** no colour term was applied, although the plan says "any colour term stated" and the
    header gives CLRCOEFF 0.12–0.15 with PCOLOR g−R. DR24 light curves also assume zero colour (release notes), and
    no star colours were used.

## 2. Overlap test (fixed criteria, PLAN §3): **passed**
- 13 pairs (star × overlap exposure, both clean).
- **Pooled:** Δ = psfcat − DR24 has a median of −0.015 mag (limit ±0.05) and a robust scatter of 0.0013 mag (limit
  0.072 mag, 2 × the median combined error).
- **The pooled scatter hides a per-star offset** (`audit.json`):
  - 672208200000073: −0.014 to −0.015 (3 pairs);
  - 672208200000239: −0.014 to −0.016 (6 pairs);
  - 672208200002085: +0.013 to +0.030 (4 pairs; median +0.020).
  - Full range −0.016 to +0.030 mag, so the stars differ by up to about 0.045 mag.
- A star-dependent offset fits a colour term or a magnitude-dependent bias correction. It is not tested here and was
  not removed. It is small against the 1-mag dimming threshold.
- **Not an independent validation.** DR24 light curves are themselves built from single-exposure PSF-fit catalogues
  (release notes, advisory 7), so this agreement shows consistency between two entries to the same measurements.

## 3. Gap fill (PLAN §4)
- **Candidates:** gap exposures of the quadrant (MJD 59853–60249, i.e. 2022-09-29 to 2023-10-30; infobits 0; at least
  one star inside its footprint without a clean DR24 detection).
  - There were 63 candidates (`audit.json`). The fill took 60 evenly spaced ones, so 3 remain unfetched.
  - 12.8 MB of psfcat was read, including 6 cached files.

| star | clean DR24 detections in the gap | usable psfcat added (unique match, flags 0) | total in the gap after fill | missing | flagged |
|---|---|---|---|---|---|
| 672208200000073 | 2 | 28 | 30 | 30 | 2 |
| 672208200000239 | 2 | 60 | 62 | 0 | 0 |
| 672208200002085 | 1 | 31 | 32 | 25 | 4 |

- **All 55 missing targets** are coverage gaps, not non-detections (`psfcat_missing.json`). None is read as a
  disappearance or a limit.
  - 52 lie outside the exposure footprint.
  - 3 lie inside the footprint but within its 10″ margin.
  - 0 lie at least 10″ inside the edge.
- **0 ambiguous matches.** Matched separations are ≤ 0.68″ in the fill (≤ 0.36″ in the 6 test gap exposures).
- **Strict recount** (`audit.json`), applying infobits 0 to both columns and dropping one exposure (223128439) where DR24
  has a flagged detection:
  - clean DR24 observations in the gap: 1 per star;
  - psfcat added: 27, 59 and 30;
  - totals: 28, 60 and 31.
- **Baseline:** each star has ≥ 20 usable observations inside the gap under both counts. This is indicative only: the
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
- **What it shows:** for these 3 stars in field 672 / CCD 8 / q2, an official product holds calibrated detections for
  exposures missing from DR24 light curves. On shared exposures these detections agree with DR24 to within −0.016
  to +0.030 mag, with per-star offsets.
- It does **not** re-run the calibration. Using psfcat in the W3 chain would need a new data path: series and
  comparison counts built from psfcat for the gap exposures. That needs a plan and the owner's authorization. Whether
  the 10 cells would then reach 0.5 is unknown. 2024 coverage stays sparse, and pre-event shortfalls remain.
