# MACS0717: Sharon v4cor and Keeton v4 as pinned map models (frozen plan, 2026-10-21 01:11Z)

**Question.** Do the HFF Sharon & Johnson v4cor (Lenstool) and Keeton v4 (parametric) models of MACS J0717.5+3745
load correctly as `MapLensModel`s, so they can replace the scratch-only second models of D-036 / D-041 for this field
(TASKS "second models as pinned `MapLensModel` entries: Sharon v4cor (MACS0717), Keeton v4 (MACS0717)")?

**Hypothesis / ordinary explanations.** A failure would most likely be a loading error (WCS orientation, units,
D_LS/D_S scaling, cosmology, map footprint), not a lensing anomaly. The maps are the files whose sha256 SOURCES.md
recorded on 2026-10-08 (deflection, κ) plus the z = 2 magnification maps (new pins).

**Known cases (validation, before any unknown).** Sharon readme Table 7 (the shared HFF v4 readme,
`hlsp_frontier_model_macs0717_sharon_v4cor_readme.pdf`): the systems with a spectroscopic z used by the model,
1, 3, 4, 6, 12, 13, 14 (Gold) and 19 (z 6.387): 26 images, transcribed into `sharon_gold_images.ecsv`. Keeton's readme
(`hlsp_frontier_model_macs0717_keeton_v4_readme.txt`) lists 29 images of gold systems 1, 3, 4, 6, 12, 13, 14, 15, 19
without positions (quoted image-plane rms 0.79″); the same 26 Table 7 images serve as its known cases (ASSUMPTION:
the HFF v4 shared numbering and positions; system 15 is not used because Sharon did not use its z).

**Frozen gates, per model (computable; `scripts/m0717_models.py`, function `gates`, same thresholds as the
MACS1149 Sharon pin):**
- G1 (maps): `map_check`'s median relative difference ≤ 0.05 for κ (pixels 0.05 < κ < 2) and for |μ(z = 2)|
  (pixels |μ| < 10) (ASSUMPTION).
- G2 (images): image-plane rms of the 26 images, back-traced to the per-system mean source and re-solved, ≤ 1.0″
  (ASSUMPTION; Keeton quotes 0.79″, Sharon none), and every one of the 8 systems reproduced with each catalogued image
  matched within 3″ (ASSUMPTION).
- Both pass → the model is pinned for use (`macs0717-sharon`, `macs0717-keeton`). Either fails → the entry may stay
  loadable, but its use is limited or withdrawn as the failure dictates; record it.

**Comparison (no gate).** The same 26 images and redshifts through `macs0717-cats` (CATS v4.1, image-plane rms 3.2″
on its own 132 images, D-035), per-image offsets in all three models.

**What changes the plan.** G1 failing → loading/units problem, fix before anything else. G2 failing with G1 passing
→ check the transcription and the footprint first (as for MACS1149 system 9). Both models passing → next unit: the
D-041 cross-model rule on the MACS0717 forced rows with pinned models (replacing the scratch readings in
`docs/fields/macs0717.md`).

**Note after the run (01:20Z, wording only).** The gate function is `m1149_sharon.gates` (imported by
`scripts/m0717_models.py`, which defines no `gates` of its own); its thresholds (0.05, 1.0″, 3″) are those written
above, and `m1149_sharon` is unchanged in this PR.
