# MACS1149: Sharon v4cor as a second map model (frozen plan, 2026-10-20 23:15Z)

**Question.** Does the HFF Sharon & Johnson v4cor Lenstool model of MACS J1149.5+2223 load correctly as a
`MapLensModel`, so it can serve as the second, independent model for the D-037 MACS1149 rows (TASKS
"second models as pinned `MapLensModel` entries")?

**Hypothesis / ordinary explanations.** A failure would most likely be a loading error (WCS orientation, units,
D_LS/D_S scaling, cosmology), not a lensing anomaly.

**Known cases (validation, before any unknown).** The model's own gold spectroscopic constraints, readme Table 5
(`hlsp_frontier_model_macs1149_sharon_v4cor_readme.pdf`): systems 1, 2, 3, 4, 5, 13, 110, 9 (24 images, z_spec from
the table). Positions transcribed into `sharon_gold_images.ecsv`.

**Frozen gates (computable; `scripts/m1149_sharon.py`, function `gates`):**
- G1 (maps): the published κ map and the z = 2 magnification map agree with the values our model derives from
  the deflection maps: median |Δκ| ≤ 0.02 and median |Δμ/μ| ≤ 0.05 over pixels with |μ| < 20 (ASSUMPTION).
- G2 (images): image-plane rms of the 24 gold images, back-traced to the per-system mean source and re-solved, is
  ≤ 1.0″ (ASSUMPTION: Sharon's rms is not quoted in the readme; HFF v4 models quote 0.3–0.8″; CATS v4.1 0.63″), and
  every one of the 8 systems is reproduced with each catalogued image matched within 3″ (ASSUMPTION).
- Both pass → the model is pinned as `macs1149-sharon` and usable for cross-model checks. Either fails → not
  pinned for use; record the failure.

**What changes the plan.** A G2 failure with G1 passing would point at the image list or redshifts (transcription)
before the model.
