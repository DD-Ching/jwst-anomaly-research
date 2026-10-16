# Plan: PA convention of the deblend `pa_sep` versus the catalogue `pa_obs` (frozen before the run)

**Why.** The MACS1423 and RXJ2129 handoffs (CHANGELOG 2026-10-16) list the convention of the two-knot separation
axis `pa_sep` (`arcs_vet_deblend.sky_pa`, astropy `position_angle` between the two fitted knots) against the
catalogue axis `pa_obs` (pipeline `sky_orientation`) as unchecked. Statements that use both depend on it: S1063
1316 "the separation is the anti axis (78.5°)" and 1557 "separation 75° from the prediction" (#273).
The clump-aware census null draws random pair orientations, so its P values do not use `pa_sep`.

**Hypothesis (H0, convention agrees).** Both are degrees east of north, mod 180. For a source the deblend calls
`two-knot`, the moment axis of the whole segment follows the knot separation, so `pa_sep ≈ pa_obs`.
**Alternatives:** mirrored axis (`180° − pa_obs`, an x-flip as tested for DJA in #284) or a 90° rotation
(`pa_obs + 90°`, the θ_image convention).

**Data (stored on `main`, no download).** `results/arcs/abells1063-cats/vetting/census.json` (24 rows) with
`results/arcs/abells1063-cats/arcs.ecsv`; `results/arcs/macs1423-canucs/vetting/census.json` (55 rows) with
`results/arcs/macs1423-canucs/macs1423-canucs-maps/arcs.ecsv`. Join on `label`.

**Primary set.** Gate-passed rows classed `two-knot`. **Secondary set** (reported, not in the rule): every
gate-passed row with catalogue `ellipticity ≥ 0.2` (ASSUMPTION), whatever its class.

**Statistic.** For each candidate convention X, the median of `axis_offset(pa_sep, X)` (0–90°) over the set.
Random axes give 45°.

**Outcome rule (ASSUMPTION thresholds, `outcome()` in `scripts/arcs_pa_convention.py`, unit-tested):**
- fewer than 8 primary rows → `undetermined`;
- the identity has median ≤ 20° and is at least 15° below both alternatives → `convention agrees`;
- an alternative meets those conditions → `convention differs: <name>` (then every stored statement combining the
  two axes is corrected in a follow-up);
- otherwise `undetermined`.

**What would change the plan.** `convention differs` reopens the 1316 / 1557 axis readings and the MACS1423 census
sheet labels. `convention agrees` closes the open item in both handoffs.
