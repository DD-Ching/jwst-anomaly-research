# W2 J2244−0550 HST test: the frozen lenstronomy two-point-source + Sérsic chain fails validation in both bands (F814W finds 2 of 6 close known lenses; F160W has 2 usable known lenses and fails the absorption check), so the target is not measured; no new evidence for or against an ordinary lens

- **Plan:** `PLAN.md`. v1 was committed at 926367d before any control run; v2 (development settings and
  thresholds) at 00716b1, before any validation run.
- **Labels:**
  - HST pixels and Gaia stars: observed;
  - fits, Δχ², status and gate: derived;
  - injections: simulated;
  - thresholds and rules: ASSUMPTIONs (`scripts/w2_j2244_hst_rules.py` `Params`).
- **Scope:** a test of whether the public HST imaging can reliably detect or limit extended light next to two bright
  images ~1″ apart. Nothing here bears on exotic physics.

## 1. Prior work and tool
- **D-060 addendum:** the earlier Moffat chain failed its negative control (S/N 13 from PSF mismatch alone) and missed
  the ~1″ lens SDSSJ0246−0825. Its revisit condition named a lenstronomy joint model.
- **Earlier analyses of these J2244 HST data:** none found this time (MAST, arXiv / web searches 2026-10-23).
- **Tool:** lenstronomy 1.14.2 `FittingSequence` (PSO → `psf_iteration` → PSO; `UNLENSED` point sources,
  `SERSIC_ELLIPSE`, `UNIFORM`) with a photutils 3.1.1 `EPSFBuilder` initial PSF.
  - Installed in an isolated uv environment; versions in `cost.json`. The project environment is unchanged.
  - Only data interface, set-up and result handling were written (`scripts/w2_j2244_hst.py`,
    `w2_j2244_hst_rules.py`, `w2_j2244_hst_jobs.py`).

## 2. Data (`products.json` holds the S3 ETags of the 19 products used; `PLAN.md` §1)
- **Target products:**
  - F160W: `if8a5k030_drz` (GO 17455, 4 × 105.5 s, IRSUB512, 2024-07-05);
  - F814W: `jf56c1010_drc` (GO 17308, 2 × 337 s, 2024-08-31).
  - One product per band; the HAP SVM/MVM resamplings of the same exposures were not used.
- **F814W saturation:** the brighter image's peak pixel is DQ-saturated in both exposures.
- **Star libraries** (Gaia DR3 stars in programme fields; `library_*.json`): 91 stars from 144 GO 17455 fields, 14
  from 8 GO 17199 fields, 375 from 80 GO 17308 fields.

## 3. Development (development controls only)
- **Development controls:** SDSSJ0818+0601 and ULASJ0743+2457 in both bands, plus 3 negative pairs per band.
- **Iterations:**
  - Iteration 1 used the v1 settings (numbers in `PLAN.md` v2; not stored).
  - Iteration 3 is stored in `runs/dev_iter3.json`.
  - Iteration 4 is the frozen setting (`runs/dev.json`); the changes are listed in PLAN v2.
- **Iteration-4 thresholds:**
  - F160W: the negative pairs give Δχ² −838, −899 (reduced χ² 18.7–26.5) and 3887, so T = 7774.3.
  - F814W: the negative pairs give 350, 4592 and 930, so T = 9184.5.
- **Stated before validation:** both F160W development positives (7707, 5750) fall below T, so F160W was expected to
  fail. The rule was not changed.

## 4. Validation (held-out controls; `runs/val.json`, `gate.json`, `png/val_*`)
- **Positive-pool geometry:** the frozen rule keeps systems with a fitted image separation of 0.7–1.8″.
  - Excluded: F160W SDSSJ1620+1203 (2.27″); F814W ULASJ0820+0812 (2.24″), SDSSJ0946+1835 (3.05″), SDSSJ1640+1932
    (2.36″) and WGD2021−4115 (2.75″).
  - Deviation: the plan's "third peak" exclusion was not implemented. Only the separation was applied.

| band | T | eligible positives: Δχ² → status | negatives: Δχ² → status | gate |
|---|---|---|---|---|
| F160W | 7774.3 | SDSSJ0806+2006 152134 → extended light; VDESJ2325−5229 25040 → offset (0.82″ > 0.6″) | 915 (unreliable, red. χ² 12.8), −40, 4517, −226, 11780 (**unresolved**, R = 0.04″), −2667 | **fail**: 2 eligible positives < 3; absorption check fails (§5) |
| F814W | 9184.5 | J0140−1152 366033 and DESI-011.9124+03.3756 33924 → extended light; SDSSJ0114+0722 45255 and 022418+024703 48884 → offset (0.64″, 0.87″); SDSSJ2146−0047 3825 and 144320−012537 718 → none | −353, 7375, 1556, 4141, 2277, 2060 → none | **fail**: 2 / 6 detected (Wilson 0.10–0.70) < 80 % |

- **No false detections:** 0 / 5 decided F160W negatives (Wilson upper 0.43) and 0 / 6 F814W (0.39). Few negatives
  without a false detection do not show a zero false-positive rate.
  - F160W pair 7 exceeded T (Δχ² 11780) and escaped only because its Sérsic was unresolved.
- **Contact sheet** (looked at): most positives show lens light in the data − M0 residual. Several are wide, ring- or
  arc-like systems with bright deflectors, easier than the J2244 case.
  - SDSSJ2146−0047 shows lens light by eye, yet its Δχ² (3825) is below T.
  - The negatives' residuals are dipoles at the star cores.

## 5. Injections (simulated; `runs/inj.json`, `gate.json` → `injections`)
- **Set-up:** Sérsic n = 4, q = 0.8, r_eff 0.2″ and 0.5″, placed 0.4 of the way from the fainter image, into
  validation negative pairs 3 and 4.
- **F160W:** no injection is recovered as "extended light" at any magnitude from 19 to 23. At 19 mag the fits are
  "unresolved" (R < 0.15″) or "fit unreliable".
  - Fitted / injected flux at 19 mag: 0.62–0.69 (r_eff 0.2″) and 0.39–0.66 (0.5″). Pair 4 at r_eff 0.5″ is below
    0.5, so the PSF-absorption check fails.
- **F814W:** recovered to 21 mag (r_eff 0.2″, both pairs) and to 20 / 21 mag (r_eff 0.5″), i.e. only bright
  extended light.
  - Fitted / injected flux at 20 mag is 1.08–2.71, so the Sérsic also takes up PSF-mismatch light.
  - After subtracting the pair's own no-injection Sérsic flux, the ratio is 0.38–1.46.
  - The flux is not a reliable measurement of galaxy light.

## 6. Target
- **Not measured.** PLAN §4 and §6 allow target residuals only in a band whose gate passes; neither did.
  - No target model or residual was computed or looked at.
  - No second method was tried and no threshold was changed.

## 7. Why it fails (descriptive; numbers in `runs/*.json`)
- **PSF mismatch rivals faint lens light.** For real-star pairs at the target's brightness (Gaia G 17.6–18.4), the
  mismatch gives Δχ² up to 4592 (F814W) and 11780 (F160W). That is as large as a faint lens galaxy's (SDSSJ2146−0047:
  3825).
- **The threshold built from negatives therefore misses faint lenses.** This is the D-060 failure mode again, now
  measured with real stars instead of a PSF model.
- **F160W core residuals.** The cores leave dipole residuals that a single empirical PSF plus PSF iteration does not
  remove: 0.128″ pixels, a 4-exposure snapshot, and stars from many field positions.
- **The centre rule is too strict.** The frozen 0.6″ rule rejects real lens galaxies that sit toward the fainter image
  (SDSSJ0114+0722 at 0.64″). This is a rule-design flaw, recorded and not changed after the fact.

## 8. Deviations and corrections
- **Traffic budget exceeded (estimate ≈ 28 GB against 3 GB; `traffic.json`).**
  - Cause: fsspec's 5 MB read-ahead block made each ACS star cutout cost about 12 MB (measured), and the F814W star
    library was built three times.
  - No image data was kept on disk. The code now uses 64 kB blocks (one cutout 21 → 4.1 MB).
  - A worker restart reset the network counter, so the figure is estimated from Gaia candidate counts
    (`traffic_candidates.json`) times the measured traffic per candidate.
- **Implementation fixes, before / after kept:**
  - Star re-centring and weight rule, before any control run (`PLAN.md` v2).
  - Sérsic flux units, after the validation and injection runs. lenstronomy amplitudes are per arcsec², so the flux
    lacked the pixel area.
    - Statuses, Δχ² and T are unaffected (only the flux sign is used).
    - The F160W absorption check changed from pass to fail.
    - Earlier outputs: `runs/*_units_v1.json`.
- **Worker restart (~07:19 UTC):** the first validation run was killed. It was rerun identically with the same seeds,
  and the runner now appends results as they finish.
- **The tool test on SDSSJ0818+0601 F160W** before v1 is disclosed. J0818 is a development control.
- **Validation positives:** all 10 F814W pool members were run. The geometry rule then left 6, exactly the plan's
  "at most 6".

## 9. Cost
- **Wall clock:** 06:23 UTC to the end in `cost.json`.
- **CPU:** 0.8 CPU-h of stored fits on 4 cores, plus unstored development iterations and library builds (estimated
  < 0.6 CPU-h, not stored).
- **Traffic:** ≈ 28 GB estimated (§8).
