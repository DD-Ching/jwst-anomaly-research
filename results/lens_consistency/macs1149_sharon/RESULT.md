# MACS1149: Sharon v4cor as a second map model — G1 passes, G2 fails as written (system 9 off the map)

Plan: `PLAN.md` (frozen before the run). Script: `scripts/m1149_sharon.py`. Numbers: `summary.json`,
`imageplane_*.ecsv`. All model outputs are `model_prediction`; the image positions and redshifts are `observed`
(readme Table 5, transcribed into `sharon_gold_images.ecsv`).

**G1 (maps): passes.** κ from the deflection maps vs the published κ: median relative difference 0.00025
(11,623 pixels). z = 2 |μ| vs the published map: 0.0164 (95th percentile 0.075; 11,135 pixels). Rule ≤ 0.05.

**G2 (gold images): fails as written.** Image-plane rms 0.226″ over the 19 solved images (rule ≤ 1.0″), and 7 of the
8 systems reproduced with every image within 3″. System 9 (z 0.981) is not solved: 9.2, 9.3 and 9.4 lie outside the
published 200″ × 200″ maps (north edge; `summary.json` `post_hoc.images_outside_maps`), so its source cannot be
traced. This is a map-footprint limit, not a test of the model there.

**Post hoc (deviation).** Restricted to systems inside the maps, G2 passes (7 / 7 systems, rms 0.226″;
`post_hoc.G2_inside_footprint`). Under the frozen rule the model is **not** pinned for use; the `MODELS` entry
`macs1149-sharon` exists so it can be loaded, and its use is limited to positions inside the map footprint
(a post hoc condition).

**Comparison: the same 23 images through CATS v4.1 (no gate).**
- At the Sharon spectroscopic redshifts: all 23 images solved, rms 4.35″, carried by system 13 at z 1.24: CATS
  predicts one image for it (offsets 16.8″, 11.8″, 1.95″).
- System 13 is free in CATS v4.1 (`params.txt` z_m_limit 1.301). At z 1.301 (post hoc), CATS reproduces it
  (offsets 0.74″, 1.73″, 0.63″) and the 23-image rms becomes 0.853″. So the system 13 difference is a redshift
  input, not a lensing residual.
- Largest remaining CATS offsets at the gold redshifts: 3.2 (1.71″), 2.3 (1.66″, a shared match with 2.2, as in
  D-037's system 2 topology note), 4.1 / 4.2 (1.30″ / 1.23″). Sharon puts every solved image within 0.53″ of the
  gold positions; the largest are 3.2 and 4.1 (0.52″ and 0.40″, read from `imageplane_macs1149-sharon.ecsv`).

**Limits.** The gold list is Sharon's own constraint set, so G2 checks that the model loads correctly, not that it
predicts independent data. Gold systems only (no photometric systems). CATS was fitted to a different but
overlapping constraint set.

**Next.** Use `macs1149-sharon` as the second model for D-037 MACS1149 rows inside the footprint (the `forced_check`
D-037 rules, cross-model classification as for Abell 2744 UNCOVER v2.0).
