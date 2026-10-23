# W2 J2244−0550 HST test: this lenstronomy two-point-source + Sérsic implementation fails validation in both bands (v3, after review fixes: F814W detects 2 of 6 close known lenses; F160W has 2 usable known lenses, below the 3 needed), so the target is not measured; no new evidence for or against an ordinary lens

- **Plan:** `PLAN.md`, committed in three versions, each before the runs it governs:
  - v1 (926367d), before any control run;
  - v2 (00716b1), with development settings and thresholds, before any validation run;
  - v3 (2b4b058), with the implementation fixes from the review, before any rerun. The v3 thresholds were committed
    (7f0065d) before the v3 validation run.
- **Primary result:** v3 (§3–§6). The v2 run is kept as the "before" record (§8).
- **Labels:**
  - HST pixels and Gaia stars: observed;
  - fits, Δχ², status and gate: derived;
  - injections: simulated;
  - thresholds and rules: ASSUMPTIONs (`scripts/w2_j2244_hst_rules.py` `Params`).
- **Scope:** a test of whether the public HST imaging, analysed with this implementation, can reliably detect or limit
  extended light next to two bright images ~1″ apart. Nothing here bears on exotic physics.

## 1. Prior work and tool
- **D-060 addendum:** the earlier Moffat chain failed its negative control (S/N 13 from PSF mismatch alone) and missed
  the ~1″ lens SDSSJ0246−0825. Its revisit condition named a lenstronomy joint model.
- **Earlier analyses of these J2244 HST data:** none found this time (MAST, arXiv / web searches 2026-10-23).
- **Tool:** lenstronomy 1.14.2 `FittingSequence` (PSO → `psf_iteration` → PSO; `UNLENSED` point sources,
  `SERSIC_ELLIPSE`, `UNIFORM`) with a photutils 3.1.1 `EPSFBuilder` initial PSF.
  - Installed in an isolated uv environment; versions in `cost.json`. The project environment is unchanged.
  - Only data interface, set-up and result handling were written (`scripts/w2_j2244_hst.py`,
    `w2_j2244_hst_rules.py`, `w2_j2244_hst_jobs.py`, `w2_j2244_hst_reflag.py`).

## 2. Data (`products.json` holds the S3 ETags of the 19 products used; `PLAN.md` §1)
- **Target products:**
  - F160W: `if8a5k030_drz` (GO 17455, 4 × 105.5 s, IRSUB512, 2024-07-05);
  - F814W: `jf56c1010_drc` (GO 17308, 2 × 337 s, 2024-08-31).
  - One product per band; the HAP SVM/MVM resamplings of the same exposures were not used.
- **F814W saturation:** the brighter image's peak pixel is DQ-saturated in both exposures.
- **Star libraries** (Gaia DR3 stars in programme fields):
  - 91 stars from 144 GO 17455 fields (`library_17455_F160W.json`) and 14 from 8 GO 17199 fields
    (`library_17199_F160W.json`);
  - 375 stars from 80 GO 17308 fields (`library_17308_F814W.json`). v3 re-derives the saturation flags at 337 s per
    exposure; 364 stars remain (`library_17308_F814W_v3flags.json`).

## 3. Development (development controls only; `runs/dev_v3.json`, `png/dev_*`)
- **Development controls:** SDSSJ0818+0601 and ULASJ0743+2457 in both bands, plus 3 negative pairs per band.
- **Thresholds** by the unchanged rule T = max(25, 2 × the largest development-negative Δχ²):

| band | development negatives Δχ² (reduced χ²) | T | development positives Δχ² → status |
|---|---|---|---|
| F160W | 13601.5 (7.8), 649.3 (9.6), 1750.7 (2.3) | **27202.9** | J0743 2742.9 → none; J0818 957.6 → none |
| F814W | 64.9 (7.9), 5465.4 (2.6), 805.3 (1.0) | **10930.8** | J0743 2249.0 → none; J0818 3119.8 → none |

- **Stated before validation (PLAN v3):** all four development positives are below T, so both bands were expected to
  fail. The rule was not changed.

## 4. Validation (held-out controls; `runs/val.json`, `gate.json`, `png/val_*`)
- **Positive-pool geometry:** the frozen rule keeps systems with a fitted image separation of 0.7–1.8″.
  - Excluded: F160W SDSSJ1620+1203 (2.27″); F814W ULASJ0820+0812 (2.23″), SDSSJ0946+1835 (3.06″), SDSSJ1640+1932
    (2.46″) and WGD2021−4115 (2.74″).
  - Deviation: the plan's "third peak" exclusion was not implemented. Only the separation was applied.
    VDESJ2325−5229 shows three compact sources (`png/val_F160W_VDESJ2325-5229.png`; one point source sits on the
    central galaxy), so that rule would have excluded it and left F160W with 1 eligible positive. The F160W gate
    fails either way.

| band | T | eligible positives: Δχ² → status | negatives (pairs 3–8): Δχ² → status | gate |
|---|---|---|---|---|
| F160W | 27202.9 | SDSSJ0806+2006 35254.2 → extended light; VDESJ2325−5229 18190.8 → none | 1947.7, 2033.2, 2549.7, 436.4, 7613.2, 1338.3 → none | **fail**: 2 eligible positives < 3 (1 detected) |
| F814W | 10930.8 | J0140−1152 184500.3 and DESI-011.9124+03.3756 19529.6 → extended light; SDSSJ0114+0722 24738.4 and 022418+024703 40135.9 → offset (0.65″, 0.85″ > 0.6″); SDSSJ2146−0047 1130.0 and 144320−012537 828.8 → none | 962.3, 4104.7, 2469.1, 6468.4, 2174.0, 1478.5 → none | **fail**: 2 / 6 detected (Wilson 0.10–0.70) < 80 % |

- **Fit quality:** reduced χ²(M1) 0.7–4.5 for the eligible positives and 1.0–9.6 for all negative pairs in both bands
  (development and validation).
- **No false detections:** 0 / 6 negatives per band (Wilson upper 0.39). Six negatives without a false detection do
  not show a zero false-positive rate.
- **Contact sheet** (looked at, `png/val_*`):
  - Most positives show lens light in the data − M0 residual. Several are wide, ring- or arc-like systems with bright
    deflectors, easier than the J2244 case.
  - In SDSSJ2146−0047 the data − M0 residual is dominated by core dipoles. Its M1 model shows the square edge of the
    41 px PSF kernel.
  - The negatives' residuals are compact dipoles at the star cores.

## 5. Injections (simulated; `runs/inj.json`, `gate.json` → `injections`)
- **Set-up:** Sérsic n = 4, q = 0.8, r_eff 0.2″ and 0.5″, placed 0.4 of the way from the fainter image, into
  validation negative pairs 3 and 4 of each band.
- **F160W:** no injection is recovered at any magnitude from 19 to 23 (all "none"; Δχ² 129–4533, all below T).
- **F814W:** recovered at 20 mag only, for both r_eff and both pairs (Δχ² 11951–37787). At 21 mag, Δχ² is 2829–7438,
  below T. The recovery limit is 20.0 mag in all four series.
- **Absorption check** (fitted / injected Sérsic flux at the brightest injection ≥ 0.5): passes in both bands, but not
  because the galaxy flux is measured.
  - Raw ratios are 2.01–3.92 (F160W, 19 mag) and 2.00–3.56 (F814W, 20 mag): the Sérsic takes up more light than was
    injected.
  - After subtracting the pair's own no-injection Sérsic flux, the ratios are −5.93 to 0.33 (F160W) and 1.32–2.24
    (F814W). The fits are not nested (§7), so this subtraction is not stable.
  - The fitted flux is not a measurement of galaxy light in either band.

## 6. Target
- **Not measured.** PLAN §4 and §6 allow target residuals only in a band whose gate passes; neither did.
  - No target model or residual was computed or looked at.
  - No second method was tried and no threshold was changed.

## 7. What the failure shows (descriptive; numbers in `runs/*.json`)
- **Real-star mismatch and lens light give overlapping Δχ².**
  - Real-star pairs matched to the target's Gaia brightness (G 17.58–18.42) reach Δχ² 13601.5 in F160W (development
    pair 0) and 6468.4 in F814W (validation pair 6).
  - Of the 8 fitted F814W lenses at 0.93–1.60″ (6 validation, 2 development), 4 fall below 6468.4: SDSSJ2146−0047
    1130.0, 144320−012537 828.8, ULASJ0743+2457 2249.0 and SDSSJ0818+0601 3119.8.
  - In F160W, 2 of the 4 lenses within the geometry rule fall below 13601.5 (the two development positives).
  - So in this implementation Δχ² does not separate the lens light of these lenses from star-to-star PSF mismatch. This is
    consistent with the D-060 failure mode, now measured with real stars instead of a PSF model.
- **Scope:** this is a failure of this implementation and settings, not of lenstronomy or of the data in general.
- **Possible reasons, none tested:**
  - **PSF halo beyond the kernel** (41 px = 2.05″ in F814W, 21 px = 2.69″ in F160W). F814W pairs 7 and 8 put the
    Sérsic at its R = 2″ bound; SDSSJ2146−0047's M1 model shows the kernel edge.
  - **Two PSFs in one pair.** All 9 F160W pairs and 7 of 9 F814W pairs combine stars from two fields, while the
    target's images share one PSF. But F814W pair 6, both stars from one field (`library_17308_F814W.json`), has the
    largest validation Δχ² (6468.4).
  - **Non-nested fits.** M0 and M1 run separate PSF iterations, so Δχ² is not a likelihood ratio and can be negative
    (v2 F160W pair 8: −2666.6, `runs/val_v2run.json`).
  - **Brightness matched in Gaia G only.** The pairs' F160W brightness relative to the target's images was not
    checked.
- **The point-source position rule (0.15″) matches the fit's bound** (a ±0.15″ box per image), so it can fire only
  near the box edge. It fired once (SDSSJ1640+1932, excluded by geometry anyway).

## 8. Before: the v2 run (kept; `runs/*_v2run*.json`, `gate_v2run.json`, `png/v2run/`)
- **Two implementation errors found by the review** (PLAN v3), fixed with this record kept:
  1. The PSF error term was zero (`keep_psf_variance_map: True` restored the absent initial variance map).
  2. The F814W saturation cut used 674 s instead of 337 s per exposure, so every star at the target's brightness was
     flagged saturated and the F814W PSFs came from G 19–21 stars only.
  - **Caveat on fix 1 (from the re-check):** lenstronomy keeps the estimated variance map only when a PSF-iteration
    step is accepted. SDSSJ1620+1203's M1 fit is identical in v2 and v3 (χ²(M1) 49944.0), so it kept the zero map.
    The outputs do not record which fits accepted a step. SDSSJ1620+1203 is excluded by geometry, so no gate
    changes.
- **v2 thresholds:** F160W 7774.3, F814W 9184.5 (`runs/dev_v2run.json`; iteration 3 in `runs/dev_v2_iter3.json`).
- **v2 gate:**
  - F160W failed: 2 eligible positives (SDSSJ0806+2006 extended light, VDESJ2325−5229 offset), and the absorption
    check failed.
  - F814W failed: 2 / 6 detected (J0140−1152, DESI-011.9124+03.3756), absorption check passed.
  - 0 false detections in 5 decided F160W and 6 F814W negatives. The largest validation-negative Δχ² was 7375 in
    F814W (pair 4) and 11780 in F160W (pair 7, "unresolved").
- **v2-only caveat:** F160W pair 3, an injection base, was itself "fit unreliable" (reduced χ² 12.8).
- **Sérsic flux units** (found after the v2 runs): lenstronomy amplitudes are per arcsec², so the stored flux lacked
  the pixel area.
  - Statuses, Δχ² and T are unaffected.
  - The v2 F160W absorption check changed from pass to fail.
  - Earlier outputs: `runs/*_v2run_units_v1.json`.
- **Same outcome after the fixes:** both bands fail in v2 and v3, and in F814W the same two lenses are detected.

## 9. Deviations and corrections
- **Traffic budget exceeded (estimate ≈ 28.5 GB against 3 GB; `traffic.json`).**
  - Cause: fsspec's 5 MB read-ahead block made each ACS star candidate cost ≈ 12.5 MB (21 MB for one full cutout),
    and the F814W star library was built four times (the first build crashed at field 3).
  - No image data was kept on disk. The code now uses 64 kB blocks (one cutout 21 → 4.1 MB).
  - A worker restart reset the network counter, so the figure is estimated from Gaia candidate counts
    (`traffic_candidates.json`) times the measured traffic per candidate.
- **F814W star G range 16.5–21.5 in the code; the plan says 17.5–21.5.** The 72 stars with G < 17.5 are all flagged
  saturated (`library_17308_F814W_v3flags.json`), so none enters a PSF or a negative pair (pair G 17.90–18.11).
- **Negative-pair separations:** the plan fixes 1.02″, but integer-pixel shifts of re-centred stars give fitted
  separations 0.875–1.232″ (F160W) and 1.000–1.117″ (F814W).
- **Worker restart (~07:19 UTC):** the first v2 validation run was killed. It was rerun identically with the same
  seeds, and the runner now appends results as they finish.
- **The tool test on SDSSJ0818+0601 F160W** before v1 is disclosed. J0818 is a development control.
- **Validation positives:** all 10 F814W pool members were run. The geometry rule then left 6, the plan's "at most 6".

## 10. Cost (`cost.json`)
- **Wall clock:** 06:23 UTC to the end time in `cost.json`.
- **CPU:** 1.5 CPU-h of stored fits on 4 cores (v2 0.8, v3 0.73), plus unstored development iterations and library
  builds (estimated < 0.6 CPU-h, not stored).
- **Traffic:** ≈ 28.5 GB estimated (§9).
