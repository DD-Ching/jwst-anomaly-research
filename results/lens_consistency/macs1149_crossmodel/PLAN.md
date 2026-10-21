# MACS1149: D-041 cross-model rule with the pinned `macs1149-sharon` on every CATS extra-image row (frozen plan)

Frozen 2026-10-21 before any run. Script: `scripts/m1149_crossmodel.py`.

**Why.** #417 pinned Sharon v4cor as `macs1149-sharon` (maps gate passes; use limited to the map footprint, post hoc).
The D-037 verdicts (systems 2, 16, knot 1192) and the 2026-10-14 `no_reference` addendum used Sharon loaded in scratch
code (`m1149_refs.py`, z_lens 0.542). This run applies the D-041 cross-model rule with the pinned entry to **all 16**
non-`observed` rows of the tracked CATS `images_predicted.ecsv`, so every MACS1149 forced-photometry verdict rests on
pinned inputs.

**Hypothesis.** No CATS extra image is both predicted by Sharon and `absent` in the data (D-037: 0 surviving flags).
**Ordinary explanations** for any disagreement: CATS-only image (model-dependent), a catalogued image misplaced by
CATS, a fitted source redshift (D-037 "fitted redshift"). **Would change the plan:** a row whose frozen outcome is
`lensing_residual_input` (below), which goes to `/vet-candidate`.

**Inputs.** `results/lens_consistency/macs1149-cats/images_predicted.ecsv` (CATS v4.1, `model_prediction`, forced
classes `derived`, 2026-10-14 run); CATS `arcs.txt` (`observed`); `macs1149-sharon` x/y maps (sha256-pinned).
No photometry is re-measured.

**Procedure (per system with a non-`observed` row).**
1. Footprint: every catalogued image of the system and the CATS row inside the Sharon maps (finite deflection);
   else the row is `outside_footprint`.
2. Sharon images of the system's magnification-weighted least-squares source at the CATS `z_sys`
   (`m1149_refs.second_model_images`, grid step 0.25″); each assigned to the nearest catalogued image within 1.5″
   (`split_second`), the rest are Sharon extras.
3. D-041 verdict per CATS row (`m1149_refs.d041_verdict`, loose radius 3.0″, ASSUMPTION): `both models`,
   `catalogued image misplaced by CATS`, or `model-dependent (CATS only)`.

**Known-case gate (G, computable).** Over the tested systems, the fraction of catalogued images that Sharon
reproduces (a Sharon image within 1.5″) must be ≥ 0.8 (ASSUMPTION), counting only images inside the footprint.
If G fails, the outcome is "tool not validated" and no verdict is used.

**Frozen outcome per row (function `outcome`).**
- `outside_footprint` → `untestable`;
- forced class `absent` and D-041 `both models` → `lensing_residual_input`;
- forced class `absent` otherwise → `explained_model_dependent`;
- any other forced class → `no_flag` (not absent).

**Secondary (pre-declared).** System 16 is re-solved at z 2.7 (D-037 addendum 2026-10-14: z ≈ 2.5–2.9) with the same
rule, because its CATS z 4.419 is model-fitted.

**Comparison (diagnostic).** Each `no_reference` row's verdict is compared with
`results/lens_consistency/macs1149-cats/no_reference_refs.json` (scratch Sharon, z_lens 0.542); differences are listed,
not used for the outcome.
