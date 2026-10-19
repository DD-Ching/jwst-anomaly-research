# AT 2024jfz z_s 1.55–1.70: PSF-matched NIRISS 4368 vs ERO NIRCam F150W (plan, frozen before any pixel read)

**Why:** handoff of #361. In `../niriss4368/RESULT.md` the post hoc rows z_s 1.55, 1.60, 1.65 and 1.70 are
"undetermined": the plain cross-instrument difference E − T flags the wing of a static compact source
(T centroid 110.8355656, −73.4530888; 0.574–0.698″ from the predictions). A pixel-kernel PSF match removes the
instrument PSF difference; the same rows are then re-tested. Not a candidate search; an ordinary-explanation step.

**Hypothesis under test (hypothesis):** a counter-image of AT 2024jfz at z_s 1.55–1.70 (Mahler+22 ICLv2, best
model, `../posthoc.json`) was bright in E (phase 61.8–83.8 d after t_0, single-peaked light curve ASSUMPTION).
**Ordinary alternatives:** the earlier flags are PSF / pixel-scale residuals of the static source; the model z_s
or position is wrong; the event is not a lensed-host event.

## Data (unchanged from #360)
E `jw04368-o001_t001_niriss_clear-f150w`, T `jw02736-o001_t001_nircam_clear-f150w` (`../niriss4368/result.json`
`uris`); one S3 byte-range box of half-size 6.5″ centred on the mean of the four predictions (it contains the TNS
position, the positive control). Nothing stored.

## Method (`scripts/smacs_at2024jfz_psfmatch.py`)
1. T (SCI and ERR²) resampled onto E's pixel grid (`arc_diff.resample`).
2. Pixel kernel K (half-width 4 E pixels) fitted with `arc_diff.kernel_fit` so that K * T ≈ E, weights 1 / (σ_E² +
   σ_T²), 4σ clipping; pixels within 0.5″ of each of the four predictions and of the TNS position are excluded from
   the fit (so neither a counter-image nor the transient trains K).
3. D = E − K * T. Aperture photometry on D exactly as #360 (r 0.15″, annulus `hit_lightcurve.ANNULUS`, grid
   0.1″ within 0.5″ of each prediction); σ_D = `hit_lightcurve.empirical_error` on D.

## Frozen rule (function `outcome_psfmatch`, unit-tested before the run), per row z_s, then combined
- Thresholds as #360's post hoc rows: thr = `threshold_njy` of that row (0.5 × |μ_k/μ_0| × 223.4 nJy).
- **Positive control (global):** D at the TNS position ≥ 5 σ_D and ≥ 0.5 × 223.4 nJy. If it fails, every row is
  "untestable".
- **Subtraction gate (per row):** |D| aperture flux at the static source centroid < thr; else "untestable".
- **Injection gate (per row):** a Gaussian (FWHM 0.08″) of aperture flux thr / 0.5 injected into E at the
  prediction, the whole chain (fit, convolution, difference) re-run: D at the prediction ≥ thr and ≥ 5 σ_D, else
  "untestable".
- NaN fraction of the disc > 0.2 → "untestable".
- Otherwise "present" if any disc point has D ≥ thr and D ≥ 5 σ_D, else "absent".
- Combined: the four rows are reported separately; the window 1.55–1.70 is "absent" only if all four are "absent".

## What changes the plan
- All "absent": the AT 2024jfz counter-image line is closed for z_s 1.35–2.50 under this model (with #360).
- Any "present": stamps looked at, `/vet-candidate`-style ordinary tests (static source, artefact) before any reading.
- "untestable": record which gate failed; no further tuning in this run.

ASSUMPTIONs: kernel half-width, the 0.5″ masks, the 5σ and threshold factors (inherited from #360).
