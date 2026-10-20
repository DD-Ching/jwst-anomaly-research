# MACS1149: Sharon v4cor as a second map model (frozen plan, 2026-10-20 23:15Z)

**Question.** Does the HFF Sharon & Johnson v4cor Lenstool model of MACS J1149.5+2223 load correctly as a
`MapLensModel`, so it can serve as the second, independent model for the D-037 MACS1149 rows (TASKS
"second models as pinned `MapLensModel` entries")?

**Hypothesis / ordinary explanations.** A failure would most likely be a loading error (WCS orientation, units,
D_LS/D_S scaling, cosmology), not a lensing anomaly.

**Known cases (validation, before any unknown).** The model's own gold spectroscopic constraints, readme Table 5
(`hlsp_frontier_model_macs1149_sharon_v4cor_readme.pdf`): systems 1, 2, 3, 4, 5, 13, 110, 9 (23 images, z_spec from
the table). Positions transcribed into `sharon_gold_images.ecsv`.

**Frozen gates (computable; `scripts/m1149_sharon.py`, function `gates`):**
- G1 (maps): the published κ map and the z = 2 magnification map agree with the values our model derives from
  the deflection maps: `map_check`'s median relative difference ≤ 0.05 for κ (pixels 0.05 < κ < 2) and for |μ|
  (pixels |μ| < 10) (ASSUMPTION).
- G2 (images): image-plane rms of the 23 gold images, back-traced to the per-system mean source and re-solved, is
  ≤ 1.0″ (ASSUMPTION: Sharon's rms is not quoted in the readme; HFF v4 models quote 0.3–0.8″; CATS v4.1 0.63″), and
  every one of the 8 systems is reproduced with each catalogued image matched within 3″ (ASSUMPTION).
- Both pass → the model is pinned as `macs1149-sharon` and usable for cross-model checks. Either fails → not
  pinned for use; record the failure.

**Comparison (no gate).** The same 23 images and redshifts through `macs1149-cats` (CATS v4.1), per-image
image-plane offsets in both models.

**What changes the plan.** A G2 failure with G1 passing would point at the image list or redshifts (transcription)
before the model.

**Amended before the first run (23:13Z), against the pushed 23:11Z version:** G1 uses `map_check`'s median
relative difference (the pushed text had an absolute |Δκ| ≤ 0.02 and |μ| < 20); the gold list has 23 images, not
24 (counting error); the CATS comparison was added. No result was seen before these edits.
