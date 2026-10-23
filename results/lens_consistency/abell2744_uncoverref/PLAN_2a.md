# Abell 2744: UNCOVER-reference extra-image test, stage 2a (reference inputs; plan, frozen)

Frozen 2026-10-23 11:25Z, after the DR3 SUPER catalogue was downloaded and its column list and value ranges were read
(no flux was compared with any prediction; no predicted or control table exists for this field yet). Follows
`PLAN_1b.md` (gate K-WLS passed). Stage 2a only fixes the reference inputs; the extra-image stage (predicted /
control tables, gate C2, flags) needs its own frozen plan (stage 2b).

## Why this stage
The column mapping of the reference catalogue decides the shape and flux inputs of the stage-3 rules
(`m0416_canucsref.stage3`). The README of the release (`README_DR3.txt`) says `theta_J2000` is in degrees, clockwise
from the first image axis, but the column holds values in [−1.5708, 1.5708], i.e. radians. A wrong position-angle
convention would mis-set the arc-referenced controls. So the convention is measured on the release's own
segmentation map before any stage-2b rule uses it (validate the tool on known cases first).

## Hypotheses (method)
- H_a: `theta_J2000` is SEP's `theta` (radians, counter-clockwise from +x of the 40 mas mosaic pixel grid).
- H_b: `theta_J2000` is radians, measured as the README says (clockwise from +y).
Both readings are tested against the second moments of each source's segment in the release segmentation map
(pixel grid of the catalogue `x`, `y`). Ordinary explanations of a failure of both: segment shapes differ from the
windowed SEP ellipse for blended or irregular sources; the segmap grid differs from the catalogue grid.

## Frozen rules (computable; function `a2744_uncoverref.theta_check`, unit-tested before the run)
1. Sample: rows with `use_phot` = 1, `flag_kron` = 0, `flag_nearbcg` = 0, `a_image` ≥ 4 px and
   1 − `b_image`/`a_image` ≥ 0.4 (ASSUMPTIONs), and a segment of ≥ 50 pixels. At most the 500 with the largest
   `a_image` (ties by `id`).
2. Per source: segment orientation φ_seg = ½ atan2(2 μ11, μ20 − μ02) (radians, counter-clockwise from +x, unweighted
   second moments of the segment's pixel coordinates). Under H_a the prediction is θ; under H_b it is π/2 − θ. The
   residual is the angle difference folded into [0°, 90°].
3. A reading is **adopted** if its median residual over the sample is ≤ 10° **and** ≥ 70 % of the sample has a
   residual ≤ 20° (ASSUMPTIONs). If both readings pass, the one with the smaller median is adopted; if neither
   passes, stop: no stage 2b with a shape rule until the convention is resolved.
4. Also frozen (from the README, not tested here): fluxes `f_f200w`, `e_f200w` are total fluxes in units of 10 nJy
   (`PHOT_ZP` 28.9); `ext` = not `flag_star`; rows with non-finite `f_f200w` are dropped; `a_image`, `b_image` are
   pixels of 0.04″ (README: all mosaics measured at 0.04″/pixel; checked against the catalogue's `x`, `y`, `ra`,
   `dec` with `m0417_canucsref.pixel_scale`).

## Inputs (observed; not committed; deleted after the run)
- `UNCOVER_DR3_LW_SUPER_catalog.fits`, 83,684,160 bytes, sha256
  `b8f265311d041b815e1511d74eadec2d88bf8c80ef2fee2394d7f78dbe4adb4c`, Google Drive id
  `19Jd3Xg5BTrHX5icb81txSn_P0LqxmLsz` linked from `https://jwst-uncover.github.io/DR3.html` (accessed 2026-10-23).
- `README_DR3.txt`, sha256 `93cfa19494b077a60541a34fc17db6125c5667f3a17a9dcd2c10d2915828c397`, id
  `1Xbfkl2ZxIm-XnGj1ResBSKWCDjkgULqc`.
- `UNCOVER_DR3_SEGMAP.fits.gz` (the "LW detection" segmap of DR3.html), id `1s7DnoX8G73sOLiMiDttlsyZLjwSB2WlX`; sha256 pinned in the script at first
  download, before the run.

## Noise source for stage 2b (decision recorded here, used only by the stage-2b plan)
The bCG-subtracted F200W mosaic (1.84 GB gzip) has no byte ranges. Stage 2b will use the JWST pipeline level-3 F200W
mosaic of GO 2561 by S3 byte-range cutouts for the local-noise rule (`m1423_canucsref.local_noise`). That mosaic is
not bCG-subtracted, so near cluster members the random-aperture σ is larger; with the frozen classes
(`m1423_canucsref.local_class`) a larger σ moves a prediction to `undetectable_local` or `ambiguous`, never to
`absent`. This is conservative against false flags and costs completeness near bright members (stated in 2b).
