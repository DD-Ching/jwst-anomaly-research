# SMACS 0723 `f200w_1096`: ICLv2 shear and magnification at the arc — consistent (weak shear)

Plan: `PLAN.md` (frozen before evaluation). Script: `scripts/smacs_1096_shear.py` (`decide` unit-tested,
`tests/test_smacs_1096_shear.py`). Numbers: `summary.json`. Contact sheet: `contact_sheet.png` (looked at).

**Outcome (frozen rule, F200W):** `consistent`. H-1096 (orientation inconsistent with the model) is not supported.

- **Model prediction** (Mahler+2022 ICLv2, best model, r = 46.96″) at the NED z = 1.3618: κ 0.1804, γ 0.1082,
  reduced shear g 0.132, μ 1.515, predicted stretch PA 122.77°. At z_s = 1 / 2 / 4: g 0.110 / 0.152 / 0.175,
  μ 1.417 / 1.61 / 1.721 (same PA). The source is weakly sheared and mildly magnified; it is outside the
  D-024 strong-shear selection (g ≥ 0.2), which is why the 2026-10-08 `arcs` run did not test it.
- **Observed** (jwst 3.0.0 pipeline catalog moments, 2736 o001): F200W label 3298 (0.298″ from the record
  position), PA 97.38°, e 0.661; F444W label 624 (0.009″), PA 97.38°, e 0.676. Offset from the predicted axis 25.39°
  in both bands.
- **Shape-noise null** (D-024 addendum model; σ_ε 0.26, e ≥ 0.5; ASSUMPTIONs): P(offset ≥ 25.39° | g 0.132) =
  0.4661, above the frozen 0.01. With g = 0 the same offset has P 0.7176, so at this shear the orientation carries
  little information either way.
- **F150W** (post hoc reading): the closest segment (label 3043, 0.054″) has PA 154.03° and e 0.388 (below the
  e ≥ 0.5 condition of the null); the contact sheet shows it is one deblended part of the clumpy galaxy (F150W splits
  the galaxy into three segments within 0.63″, not stored), so its moments do not describe the whole galaxy.
  F200W and F444W segments cover the whole galaxy on the sheet.

**Deviations.** No posterior: ICLv2 has no pinned `bayes.dat` (`chain.note` in `summary.json`), so model
uncertainty on g and the PA is not evaluated; best model only.

**Limits.** Single model; catalog moments (no PSF correction, clumpy morphology); the g < 0.2 regime gives low
power to the orientation test (P 0.47 vs 0.72 at g = 0).

**Record closure.** The candidate record's open question (shear and magnification at 1096) is answered: weakly
sheared (g 0.132), μ 1.515, orientation consistent with ICLv2 under shape noise. Verdict unchanged (`explained`).
