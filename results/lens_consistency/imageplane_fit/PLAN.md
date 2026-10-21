# Image-plane source fit: route check on MACS1149, then Abell S1063 z_spec systems (frozen 2026-10-21 21:25Z)

**Question.** Abell S1063's CATS v4.1 map gives rms 18.5″ on the Chan+2020 z_spec systems with the mean-source route
(`../abells1063_d035/RESULT.md`), although CATS fit the same images and redshifts to a quoted 0.48″. Is the failure the
route (the mean back-traced source can land where the map forms fewer images) rather than the map?

**Hypothesis (H-route).** Fitting each system's source position in the image plane brings S1063's z_spec rms near
the quoted value. **Ordinary alternatives:** the map itself misfits (ROI / resolution / smoothing near critical
lines), identification or frame differences, Lenstool's per-image σ weighting.

**Route** (`scripts/imageplane_fit.py`): per system with ≥ 2 traced images, minimise Σ dθ² over the source position
(β_x, β_y), where dθ pairs catalogued with predicted images one-to-one (`linear_sum_assignment`); a source giving fewer
predicted than catalogued images costs a fixed penalty. Starts: the mean traced source and each image's own traced
source; Nelder–Mead from each; best kept. Same solver grid step (0.25″) as `validate`.

**Frozen gates (ASSUMPTIONs):**
- V1 (synthetic, on the S1063 map): 10 sources at known β with ≥ 3 predicted images, started from β + 1″ offsets:
  recovered images rms ≤ 0.05″ in ≥ 9 / 10.
- V2 (known case): MACS1149 CATS `arcs.txt` (145 images) with the fit: every image solved and rms ≤ 1.5 × 0.63″.
- R (only if V1 and V2 pass): the D-035 rules of `../abells1063_d035/PLAN.md` (R1 whole-list ≤ 0.72″ with every
  image solved on subset S; else R2 per-system ≤ 0.5″), with fitted instead of mean sources.
**What changes the plan:** an opened S subset → `images` on S1063 next; V1/V2 fail → the route is not usable; R fails
with V passing → the map (not the route) is the limit and S1063 stays map only.
