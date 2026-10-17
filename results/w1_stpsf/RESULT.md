# W1 shear screen: model (STPSF) PSF anisotropy vs the star χ_PSF and Abell 2744's B-mode floor (2026-10-17)

Plan and rules: `PLAN.md` (S1, S2 frozen before any model PSF; the addendum lists the deviations). Script:
`scripts/w1_stpsf_floor.py`. Labels: model PSFs are `model_prediction` (STPSF 2.2.0, data 2.2.0, F200W,
`DET_DIST`, 10″ stamps, 3 positions on each of the 8 SW detectors); everything else `derived` from the committed W1
inputs (F200W catalogue sha256 `9953c869…`, DJA v7.2 zout `43662686…`, Bergamini+23 model). Thresholds and
isophote levels are ASSUMPTIONs.

## Results
- **Direction (`model_chi.json`):** at isophote levels 3e-2 to 3e-4 × peak, the median sky χ of the 24 grid PSFs lies
  at PA 125.4–128.3° (star median 133°, `results/w1_psf`). The rotation uses the mosaic's PC matrix (pixel x axis at
  PA 130.8°; ASSUMPTION: detector and mosaic pixel axes agree). In the pixel frame, the stars' catalogue
  `orientation` is 175.8° / 179.8° / 15.5°; the area-matched model gives 176.9° / 0.6° / 178.8°
  (`stars_check.json`). Two stars agree within 1.2°; the third is 16.7° off. In `model_psf_isophotes.png` (looked
  at) the faint isophotes follow the horizontal spike pair. The figure was drawn in the session from the first run's 6″ PSFs (NRCA1, NRCB3, NRCA3 centres); the
  plotting code is not in the script.
- **S1 (per component ≤ 0.05; deviation: median over the 24 grid PSFs, not over the screened sources):**
  consistent at 2 of 6 levels (3e-2: |χ| 0.183; 3e-4: |χ| 0.209). It fails at 1e-1 (0.008), 1e-2 (0.077),
  3e-3 (0.107) and 1e-3 (0.091). Model |χ| is not monotonic in the level. No segment touches the stamp edge at any
  level (`n_edge` 0).
- **Area-matched stars (post hoc, `stars_check.json`):** at the level where the model segment area equals each
  star's `isophotal_area` (464 / 1,591 / 1,034 px²; levels 1.3e-4 to 5e-4), model |χ| is 0.106 / 0.229 / 0.251 and
  star |χ| is 0.053 / 0.158 / 0.193. |χ| grows with area in both, but the model is larger by 0.05–0.07 for each
  star (component differences up to 0.109). The implied pixel threshold (`level_times_flux`) agrees across the
  three stars within a factor 1.8.
- **S2 (deviation: constant χ per level instead of per source, `screens_abell2744.json`):** b_max 4.21
  (isotropic) → 3.97 / 3.74 / 3.90 / 3.89 / 3.89 / 3.57 for levels 1e-1 … 3e-4. Only the two S1-consistent levels
  reach ≤ 3.79. They cost 11–12 % of the shear sources (4,951 → 4,423 / 4,336). At 3e-2 the peak moves to
  (−72, 58). At 3e-4 the E-mode maximum rises to 3.87, above that run's rotation-null 95th percentile (3.82).
  **Not carried to injections:** the frozen S2 is per source. The level that passes for stars (≈ 3e-4) would
  over-correct fainter sources, whose segments sit at higher levels relative to their peak, where the model |χ|
  is 0.008–0.107. The correction also creates E-mode power. W1 limits unchanged.

## Conclusion (narrow)
The star-measured PSF anisotropy in Abell 2744's F200W catalogue is close to the axis of the modelled NIRCam PSF's
faint isophotes (2 of 3 stars within 1.2°, one 16.7° off), so it is plausibly optical. Its size depends on each
source's isophote level, and at the stars' own isophotal areas STPSF predicts |χ| larger by 0.05–0.07. A
constant correction at the stars' level lowers b_max to 3.57 but over-corrects fainter sources and adds E-mode
power. No correction is adopted; the 4.21 floor is not explained.

## Limitations
- STPSF at detector sampling, no drizzle, jitter or segmentation smoothing kernel (each rounds the PSF; a likely
  part of the over-prediction, not tested). Median over 24 grid PSFs; each mosaic point mixes dithers.
- Three stars; isophote level ↔ area matching on a 36-level grid.

## Next
A per-source forward model: convolve each source's deconvolved Gaussian with the STPSF PSF (drizzled/jittered,
validated on the three stars' χ first), apply the catalogue isophote and remeasure moments; or move the screen to
PSF-weighted moments that do not reach the spike wings.

## Regenerate
```
curl -sSL <stpsf-data-LATEST.tar.gz URL, SOURCES.md> | tar xz -C <dir>   # 88 MB, streamed, deleted after
STPSF_PATH=<dir>/stpsf-data python scripts/w1_stpsf_floor.py psf --out results/w1_stpsf/model_chi.json
python scripts/w1_stpsf_floor.py stars --chi results/w1_stpsf/model_chi.json \
    --catalog <data>/cache/mast/jw02561-o001_t003_nircam_clear-f200w/jw02561-o001_t003_nircam_clear-f200w_cat.ecsv \
    --out results/w1_stpsf/stars_check.json
python scripts/w1_stpsf_floor.py screen --chi results/w1_stpsf/model_chi.json --work <dir> \
    --out results/w1_stpsf/screens_abell2744.json
```
Time: PSF grid ≈ 1.9 min on 4 processes (10″); screens ≈ 1 min (not stored).
