# AT 2024jfz: NIRISS 4368 F150W check of the z_s 1.35–1.40 counter-image (result, 2026-10-19)

Plan: `PLAN.md` (frozen in `a779a53` before any pixel read). Every number is in `result.json`. Positions are
`model_prediction` (Mahler+22 ICLv2, `../posthoc.json`), fluxes `derived` (aperture r = 0.15″, annulus 0.35–0.6″,
S3 byte-range boxes, nothing stored), thresholds ASSUMPTIONs. Stamps were looked at before this was written.

## Frozen outcome: **absent**
- Epoch E: `jw04368-o001_t001_niriss_clear-f150w` (MJD 60557.96, 0.0656″ px); template T: ERO
  `jw02736-o001_t001_nircam_clear-f150w` (2022, 0.0312″ px).
- Coverage: 119 grid points within 0.5″ of the z_s 1.35–1.40 track; none NaN in E or T.
- Injection control: 542.5 nJy Gaussians at the 6 track points give Δ = 536.6–548.3 nJy (6 / 6 pass; σ_Δ = 9.36 nJy).
- Grid: max Δ = 9.5 nJy (1.01 σ_Δ), min −12.6 nJy; 0 points ≥ T_ap = 271 nJy, 0 ≥ 5 σ_Δ. Track fluxes are
  −3.7 to 6.4 nJy in E and −1.3 to 2.2 nJy in T: blank sky in both (`stamps_track.png`).
- Cross-check: T at the TNS position gives 31.3 ± 3.2 nJy, as `../../smacs_arcdiff/o4043_o12774/vet/lightcurves.ecsv`
  has 35.0 ± 4.9 nJy for the same image and aperture.

Reading: a repeat of AT 2024jfz at ≥ 1/8 of the predicted 2170 nJy (aperture units) is excluded at the predicted
z_s 1.35–1.40 counter-image within 0.5″. The z_s 1.35–1.40 branch of the lensed-host hypothesis is disfavoured under
this model and the #358 plateau assumption.

## Also measured (in the plan, not in the rule): AT 2024jfz itself at +116.1 d
E at the TNS position: 254.7 ± 5.6 nJy (T: 31.3 ± 3.2). The transient is still visible in NIRISS F150W 116.1 d after the
2024 F115W detection (2688 nJy, different instrument; aperture fractions differ, so no decline rate is derived here).
`stamps_tns.png` shows a point source at the TNS position in E with no counterpart in T. This NIRISS point is not in
the NIRCam-only light curve of `hit_lightcurve.py`.

## Post hoc (not in PLAN.md; exploratory)
Because the transient is still bright in E, every counter-image whose delay lies in [0, 116.1 d] shows a phase at
which (ASSUMPTION: single-peaked light curve) the event was at least as bright as at the TNS position in E.
- Rule `posthoc_outcome` (tested), per counter-image over a 0.5″ disc: "present" if any Δ ≥ 0.5 × |μ_k/μ_0| ×
  223.4 nJy (E − T at the TNS position; same instrument, so no aperture correction) and ≥ 5 σ_Δ.
- Injection gate: the expected counter-image (aperture flux |μ_k/μ_0| × 223.4 nJy) is injected at the predicted
  position and must pass the same cuts, else "untestable". All 25 rows pass (Δ_inj 181.2–228.5 nJy).
- Rule `classify_flag` (tested, written after the flags were seen) handles a "present" within 0.3″ of a source whose
  E and T centroids agree to 0.05″. That flag is the source's cross-instrument residual. It becomes "absent" if the
  source's r = 0.8″ flux change |E − T| is below the threshold, and "undetermined" otherwise.

| z_s rows | Outcome | Numbers |
|---|---|---|
| 1.05 (far image, 22.7″), 1.35–1.50, 1.75–2.50 (21 rows) | absent | max Δ ≤ 40.9 nJy vs thresholds 93.5–106.5 nJy |
| 1.55, 1.60, 1.65, 1.70 | present → **undetermined** | max Δ 188.6–1285.1 nJy, each 0.197–0.210″ from one compact source |

The four flags sit on the wing of a source present in both epochs: centroids agree to 0.014″; r = 0.8″ flux 7950.6 nJy
in E and 7744.9 nJy in T (`static_check`; `stamps_posthoc_flags.png`). So the flag is the PSF and pixel-scale
difference at 0.2″ from a bright core, not a new source. The predicted positions lie 0.574–0.698″ from that source.
Its +205.7 nJy change (`flux_change_frac` 0.0266) is above each threshold (96.9–99.7 nJy), and the cross-instrument
calibration is not known to that level here. A counter-image blended with it is therefore neither seen nor excluded.

Reading (post hoc): with a single-peaked light curve, nearby counter-images are not seen for z_s 1.35–1.50 and
1.75–2.50, nor the far image at z_s 1.05. The 1.75–2.50 range includes 1.75–1.80, next to system 6's model z 1.71.
z_s 1.55–1.70 is undetermined, and z_s < 1.35 or > 2.5 is not covered. Rows step by 0.05 in z_s; values between rows
are not tested. Positions come from the best model only (no chain).

## Limitations
- Model positions: best ICLv2 model, no position uncertainty beyond the assumed 0.5″.
- E and T are different instruments: the static-source flags show the limit of a cross-instrument difference
  within ≈ 0.2″ of bright cores.
- The injected Gaussian (FWHM 0.08″) is narrower in its wings than the real NIRISS PSF, so injections slightly
  overstate recovery; T_ap assumes an aperture fraction ≥ 0.5.
- No candidate can come from this step: the frozen outcome is "absent", and the post hoc flags are a static source.
