# `arcs`: second-band shapes of the 3 single-component pooled anti rows (plan, frozen before the run)

Continues `../photoz-contam/RESULT.md` "Next" (D-024 addendum 10). Exploratory: the rows were chosen after the pooled
counts were seen. No second-band shape of these rows has been read before this plan was committed.

## Question and hypotheses
- Rows: WHL0137 2194 (primary run `whl0137-relics-lenstool`), RXC J0600 581 (`rxc0600-relics-glafic`), MACS0308 2157
  (`macs0308-relics-glafic`); positions, F200W axes and classes from each run's stored `arcs.ecsv`.
- H_blend (ordinary): the F200W axis is set by a blend of differently coloured objects, or by colour-dependent
  structure (a red bulge plus a blue disc or clumps), so the axis changes between bands.
- H_achromatic: the axis is the same in every band. This is what a lensing stretch gives, but also an intrinsic
  shape (addendum 4: band agreement is a blend veto, not a lensing test).

## Inputs (public; read in the session scratchpad, deleted after the run, none stored in git)
- Pipeline `_cat.ecsv` of the same observation as each F200W run, from `stpubdata` S3 (sha256 recorded in
  `summary.json`):
  - SW band: F150W of `jw02282-o010_t001`, `jw06796-o010_t001`, `jw05293-o005_t001`.
  - LW band: F444W of `jw02282-o010_t001` and `jw06796-o010_t001`; F410M of `jw05293-o005_t001` (no F444W
    catalogue exists for that observation on S3, checked 2026-11-08).
- Matching, usability and re-classing reuse `scripts/arcs_bands.py` unchanged (`match_second_band`, `reclassify`,
  `usable_mask`): nearest source within 0.3″, not shared by two F200W rows, e ≥ 0.3, S/N ≥ 10 (ASSUMPTIONS, as
  in addendum 4).

## Frozen rules (functions in `scripts/arcs_anti_bands.py`, unit-tested before the run)
1. **Band gate (per band, SW and LW separately, pooled over the 3 fields):** control = F200W `aligned`,
   background, strong-shear rows of the 3 primary runs that are usable in that band. The band passes if the
   control has ≥ 10 usable rows and ≥ 0.8 of them stay `aligned` in that band (ASSUMPTION). A failed band is not
   used for any row.
2. **Per row and band:** `agree` if |Δpa(band, F200W)| ≤ 15° (ASSUMPTION; addendum 4 control medians 1.7–3.9°),
   `differs` otherwise; `unusable` if the row is not usable in that band.
3. **Per row label:** `blend_flag` if any passing band is `differs`; `same_axis` if ≥ 1 passing band is `agree`
   and none `differs`; `undetermined` otherwise.
4. **Frozen outcome:** ≥ 2 `blend_flag` → "open part of the pooled anti excess explained by colour blends or
   colour structure (exploratory)"; ≥ 2 `same_axis` and 0 `blend_flag` → "blend veto passed; the rows stay open
   (not a lensing test)"; otherwise "mixed / undetermined".

Secondary (reported, never a label): each row's class in each band; control Δpa median and 95th percentile per band.

What would change the plan: `blend_flag` rows drop out of the open anti count (next: a deblended re-measure or a
re-run of the pooled tally without them). If the rows keep their axis, the next test is the deciding models' local
shear error (model-to-model spread of the predicted axis at each row), not more shape work.
