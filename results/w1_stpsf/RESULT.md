# W1 shear screen: model (STPSF) PSF anisotropy vs the star χ_PSF and Abell 2744's B-mode floor (2026-10-17)

Plan and rules: `PLAN.md` (S1, S2 frozen before any model PSF; the addendum lists the deviations). Script:
`scripts/w1_stpsf_floor.py`. Labels: model PSFs are `model_prediction` (STPSF 2.2.0, data 2.2.0, F200W,
`DET_DIST`, 6″, 3 positions on each of the 8 SW detectors); everything else `derived` from the committed W1
inputs (F200W catalogue sha256 `9953c869…`, DJA v7.2 zout `43662686…`, Bergamini+23 model). Thresholds and
isophote levels are ASSUMPTIONs.

## Results
- **Direction (`model_chi.json`):** at isophote levels 3e-2 to 3e-4 × peak the median model χ lies at sky PA
  124.6–129.4° (stars: median 133°, `results/w1_psf`). Rotation uses the mosaic's PC matrix (pixel x axis at PA
  130.8°); the assumption that detector and mosaic pixel axes agree is checked by the catalogue itself: the three
  stars' `orientation` is 175.8° / 179.8° / 15.5° from +x, the model's 178.0° / 1.4° / 178.9° (`stars_check.json`).
  The figure `model_psf_isophotes.png` (looked at) shows the faint isophotes following the horizontal spike pair.
- **S1 (frozen, per component ≤ 0.05): fails at all six levels.** The model |χ| depends on the level and is not
  monotonic: 0.016 (1e-1), 0.230 (3e-2), 0.075 (1e-2), 0.096 (3e-3), 0.090 (1e-3), 0.234 (3e-4); the star median
  0.157 lies inside that range.
- **Area-matched stars (post hoc, `stars_check.json`):** at the level where the model segment area equals each
  star's `isophotal_area` (464 / 1,591 / 1,034 px²), model |χ| 0.105 / 0.229 / 0.252 vs star 0.053 / 0.158 /
  0.193: same axis and the same ordering with area, but the model is larger by 0.05–0.09 per star (component
  differences up to 0.109). The implied pixel threshold (level × isophotal flux) agrees across the three stars
  within a factor 1.8.
- **S2 (deviation: constant χ per level, `screens_abell2744.json`):** b_max 4.21 (isotropic) → 3.98 / 3.71 / 3.90
  / 3.89 / 3.89 / 3.49 for levels 1e-1 … 3e-4. Only the two levels with |χ| ≈ 0.23 reach ≤ 3.79; they cost 14 % of
  the shear sources (4,951 → 4,260 / 4,249), move the 3.71 peak to (−64, 55), and raise the E-mode maximum to 3.78
  / 3.77 (rotation-null 95th percentile 3.78 / 3.74). Picking a level after seeing b_max would be tuning, and the
  stars show the model |χ| is too large, so **not adopted**; W1 limits unchanged.

## Conclusion (narrow)
The star-measured PSF anisotropy in Abell 2744's F200W catalogue has the axis of the modelled NIRCam PSF's faint
isophotes (the horizontal spike pair), so it is plausibly optical; its size depends on each source's isophote
level and STPSF over-predicts it for the three stars by 0.05–0.09. A constant correction cannot be fixed from the
model; none of the tested corrections is adopted, and the 4.21 floor stays unexplained.

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
Time: PSF grid ≈ 1.7 min on 4 processes; screens ≈ 1 min (not stored).
