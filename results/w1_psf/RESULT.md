# W1 shear screen: PSF anisotropy from stars vs Abell 2744's B-mode floor (2026-10-17)

Plan and rules: `PLAN.md` (R1, R2 frozen before the run; R3 added after R1/R2, before any injection result:
deviation, labelled there). Script: `scripts/w1_psf_floor.py`; injections: `inject_shear.py --psf-chi`.
Labels: all numbers `derived` from the committed W1 inputs (F200W GLIMPSE catalogue sha256 `9953c869…`, DJA v7.2
zout `43662686…`, Bergamini+23 model) and Gaia DR3 (VizieR I/355, G < 21, 5′ cone, retrieved 2026-10-17);
thresholds are ASSUMPTIONs.

## Results
- **PSF stars (`psf_abell2744.json`):** 3 catalogue rows pass the frozen star rule (Gaia match, S/N > 50, mag ≥
  19.5): |χ★| 0.053 / 0.158 / 0.193 at PA 127° / 131° / 146°. Median χ_PSF (−0.015, −0.156), i.e. |χ| 0.157 at
  PA ≈ 133°. Of all 14 Gaia-matched rows (any magnitude or S/N), 11 have PA 114–149° (session check, not stored).
- **R1 map feasibility: fails.** 3 stars (rule ≥ 10); median distance from a used source to its nearest star 92″
  (rule ≤ 10″; 90th percentile 157″). A star-based PSF map cannot resolve the 10″ aperture scale in this catalogue.
- **R2: passes for the constant correction.** b_max = max |S_×| falls from 4.21 to 3.66 at the same centre
  (88, −42); the field-mean ε₂ of used sources moves from −0.026 to −0.004; the rotation-null 95th percentile of
  max S is 3.84 with the correction (3.76 without). The E-mode stays null: S_max 3.49 (3.39 before) < 3.84.
  Map figure `sxmap_abell2744.png` (looked at): the S_× pattern is visibly the same; two of the three stars lie
  outside the screened ±120″ grid.
- **R3 (injections, 200 lenses per mass, seed 50): fails; not adopted.** The isotropic run reproduces the committed
  recoveries exactly (0/0/14/72/146 of 200; `inject_summary_abell2744_iso.json`). With the correction
  (`inject_summary_abell2744_psf.json`): 0/0/15/50/97. Recombined five-field 95 % limits
  (`combined_5fields_psf.json`): 5,765 → 5,593 / 1,082 → 1,240 / 567 → 666 deg⁻² at 2 × 10¹² / 8 × 10¹² /
  2 × 10¹³ M☉. Lower at one mass only, so the rule fails.
- **Why the efficiency falls (post hoc, from the stored summaries):** the median best injected S drops from 3.70 to
  3.29 (8 × 10¹²) and 4.96 to 3.85 (2 × 10¹³), and the mean number of painted images used from 130 to 113 and 311
  to 268. That fall (11–22 %) is larger than the floor's (13 %). Which step of the correction removes them (more
  rows unresolved after the anisotropic subtraction, the responsivity rising from 0.451 to 0.485, or both) is not
  tested.

## Conclusion (narrow)
A constant PSF distortion of |χ| ≈ 0.16 at PA ≈ 133°, measured from 3 stars, lowers Abell 2744's B-mode floor to
3.66, below the rotation null's 95th percentile, so an uncorrected PSF anisotropy is a plausible contributor to the
floor (supported by the field-mean ε₂ shift; 3 stars only). The same correction lowers injected-lens recovery more,
so it is not adopted and W1 limits stay as committed (5.8 × 10³ / 1.1 × 10³ / 5.7 × 10² deg⁻²).

## Next
- A model PSF per detector and dither (STPSF, resampled like the mosaic; reuse-check first) gives χ_PSF at each
  source without stars, and a moment-matched trace instead of the 1st-percentile σ.
- Before that, find which part of the correction removes injected images (post-hoc question above): if it is the
  unresolved-row cut, a correction that keeps small rows (e.g. KSB-style PSF-smear terms) could keep both gains.

## Regenerate
```
python scripts/gaia_stars.py 3.586 -30.40 --radius-arcmin 5 --gmax 21 --out <dir>/gaia21.ecsv
python scripts/w1_psf_floor.py --stars <dir>/gaia21.ecsv --work <dir> --out results/w1_psf/psf_abell2744.json \
    --png results/w1_psf/sxmap_abell2744.png
python scripts/inject_shear.py --fields abell2744 --out <dir>/inj_iso
python scripts/inject_shear.py --fields abell2744 --psf-chi -0.01519762273987317 -0.15611039032845844 \
    --out <dir>/inj_psf
```
The recombination (`combined_5fields_psf.json`) used `inject_radial.combine` on the Abell 2744 summary plus the
committed `results/w1_shear_s1063/rerun_4fields/*_summary.json` and `inject_summary_abells1063.json` (session
script, stored output).

## Time
About 10 min of computation (screens 10 s each; two injection runs of ≈ 4.5 min in parallel), ≈ 35 min in all
(not stored).
