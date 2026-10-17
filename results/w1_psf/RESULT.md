# W1 shear screen: PSF anisotropy from stars vs Abell 2744's B-mode floor (2026-10-17)

Plan and rules: `PLAN.md` (R1, R2 frozen before the run; R3 added before injection results; a post-review model
correction: all deviations, labelled there). Script: `scripts/w1_psf_floor.py`.
Labels: all numbers `derived` from the committed W1 inputs (F200W GLIMPSE catalogue sha256 `9953c869…`, DJA v7.2
zout `43662686…`, Bergamini+23 model) and Gaia DR3 (VizieR I/355, G < 21, 5′ cone, retrieved 2026-10-17);
thresholds are ASSUMPTIONs.

## Results (`psf_abell2744.json`)
- **PSF stars:** 3 catalogue rows pass the frozen star rule (Gaia match within 0.3″, S/N > 50, mag ≥ 19.5):
  |χ★| 0.053 / 0.158 / 0.193 at PA 127° / 131° / 146°. Median χ_PSF (−0.015, −0.156): |χ| 0.157 at PA ≈ 133°.
  Of all 14 Gaia-matched rows (any magnitude or S/N), 11 have PA 114–149° (session check, not stored).
- **R1 map feasibility: fails.** 3 stars (rule ≥ 10); median distance from a used source to its nearest star 92″
  (rule ≤ 10″; 90th percentile 157″). A star-based PSF map cannot resolve the 10″ aperture scale in this catalogue.
- **R2 (constant correction, corrected PSF model): fails.** b_max = max |S_×| falls from 4.21 to 3.85 at the same
  centre (88, −42); the rule needs ≤ 3.79. Not carried to injections. Side numbers: used sources 4,951 → 4,496
  (the larger PSF leaves more rows unresolved); responsivity 0.451 → 0.476; field-mean ε₂ of used sources −0.026 →
  −0.006; rotation-null 95th percentile of max S 3.76 → 3.78. E-mode still null: S_max 3.39 → 3.62 (< 3.78).
- **Figure** `sxmap_abell2744.png` (looked at): S_× before and after, |S_×| = 3.79 contour, the floor centre (+) and
  the stars as sticks (East left). The pattern is visibly the same; the floor peak stays at (88, −42); two of the
  three stars lie outside the screened ±120″ grid.
- **Withdrawn (post-review):** a first run with a PSF model whose minor axis was smaller than measured passed R2
  (b_max 3.66, not stored) and went to injections that painted images with a round PSF; both defects were found by the code
  review, so those numbers are not used (not stored).

## Conclusion (narrow)
A constant PSF distortion of |χ| ≈ 0.16 at PA ≈ 133°, measured from 3 stars, lowers Abell 2744's B-mode floor from
4.21 to 3.85 and removes most of the field-mean ε₂ offset, so an uncorrected PSF anisotropy is a plausible part of
the floor, not all of it. It misses the frozen 3.79 rule and costs 9 % of the shear sources; not adopted, W1 limits
unchanged (5.8 × 10³ / 1.1 × 10³ / 5.7 × 10² deg⁻²).

## Limitations
- 3 stars; Gaia positions (epoch 2016.0) are matched without proper-motion propagation and without a stellar
  classification check (G < 21 at high latitude; a fast-moving star could be missed, a compact galaxy kept).
- One Gaussian PSF for isophotal catalogue moments of all sizes (ASSUMPTION).

## Next
A model PSF per detector and dither (STPSF, resampled like the mosaic; reuse-check first) gives χ_PSF at each
source without stars, and a measured size instead of the 1st-percentile σ.

## Regenerate
```
python scripts/gaia_stars.py 3.586 -30.40 --radius-arcmin 5 --gmax 21 --out <dir>/gaia21.ecsv
python scripts/w1_psf_floor.py --stars <dir>/gaia21.ecsv --work <dir> --out results/w1_psf/psf_abell2744.json \
    --png results/w1_psf/sxmap_abell2744.png
```
`inject_shear.py --psf-chi CHI1 CHI2` (one field per run) runs injections with the same PSF model when a rule calls
for them.

## Time
Screens about 10 s each; two (withdrawn) injection runs of ≈ 4.5 min in parallel; ≈ 40 min in all (not stored).
