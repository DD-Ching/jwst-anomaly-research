# W1 shear screen: does the model (STPSF) PSF anisotropy set Abell 2744's B-mode floor? (plan, 2026-10-17)

Claimed by a cloud routine session. TASKS "JWST focus" 1 (W1); follows `results/w1_psf/RESULT.md` "Next".
Written before any model PSF was computed (frozen at this commit; deviations get a labelled addendum).

## Hypothesis and ordinary alternatives
- **H (instrumental):** NIRCam's optical PSF anisotropy at each source's detector position, as STPSF models it,
  is what the 3 stars measured (|χ| 0.157, PA ≈ 133°) and sets part of the 4.21 B-mode floor at (88, −42).
- **Alternatives:** the star χ is not optical (isophotal-moment bias of the catalogue, drizzle kernel, a
  misidentified star), shape noise in the rotation null's tail, bright-galaxy light gradients (`results/w1_floor`).
- **What would change the plan:** a model |χ| far below 0.157 (rule S1) says the star χ is not the optical PSF,
  so the next floor test is the catalogue moments (light gradients), not a PSF model. A model χ field that
  lowers b_max (rule S2) goes to injections.

## Rules (ASSUMPTIONs, computable)
- **S1 consistency:** model χ (median over the screened sources) agrees with the star median if
  |χ_model − χ★| ≤ 0.05 per component.
- **S2 adoption:** the per-source model χ correction is carried to injections only if b_max ≤ 3.79 (as R2 of
  `results/w1_psf/PLAN.md`).

## Addendum (deviations, written after the first S1/S2 numbers; labelled in RESULT.md)
- **Deviation (replaced):** S2 ran with each isophote level's median model χ as a constant correction, not a
  per-source χ: the catalogue does not record each source's detector/dither mix or its isophote level relative to
  its peak, so a per-source model χ was not built this cycle.
- **Added (post hoc):** each star's model χ at the level whose segment area matches its catalogue
  `isophotal_area` (`w1_stpsf_floor.py stars`), to test the model where the catalogue fixes the level.
- **Deviation (statistic):** S1 compares the median over the 24 grid PSFs (each rotated to the sky first), not
  over the screened sources.
- **Changed after code review:** stamps 6″ → 10″ (faint isophotes could reach the 6″ edge; at 10″ none does), and
  the median is taken after rotating to the sky. Before this, S1 failed at all six levels; at 10″ with the
  sky-frame median it passes at two. Both runs used the same rules.
