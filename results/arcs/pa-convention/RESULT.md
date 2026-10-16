# Result: deblend `pa_sep` and catalogue `pa_obs` share one convention (derived)

Plan frozen first: `PLAN.md` (`2599b94`). Code: `scripts/arcs_pa_convention.py`; rule `outcome()` unit-tested in
`tests/test_arcs_pa_convention.py`. Outputs: `summary.json`, `rows.ecsv`. Stored inputs only, no download.

**Primary set** (gate-passed `two-knot`, S1063 5 + MACS1423 13 = 18 rows): median offset of `pa_sep` from
`pa_obs` **3.0°**; from the mirrored axis 35.2°; from the 90°-rotated axis 87.0° (`summary.json` `primary`).
**Frozen outcome: `convention agrees`.** Per field: S1063 1.3° / 15.5° / 88.7°, MACS1423 3.9° / 36.7° / 86.1°.
The S1063 per-field mirrored margin (14.2°) is under 15°, but the rule applies to the pooled set only.

**Secondary set** (gate-passed, ellipticity ≥ 0.2, any class; 66 rows): 2.5° / 46.1° / 87.5°.

**Reading.** The S1063 statements that combine the two axes stand: 1316's knot axis is 78.5° from the z_s = 2
prediction and 5.5° from its catalogue axis; 1557's is 75.2° and 3.2° (`rows.ecsv`). The open item in the MACS1423
and RXJ2129 handoffs is closed for these two fields.

**Post hoc (not in the plan; `rows.ecsv` columns `offset_sep_pred_z2`, `offset_obs_pred_z2`).** Two primary rows
have knot axes far from their catalogue axes: MACS1423 2197 (81.3°) and 691 (51.6°). 2197 is a screen-anti row whose
fitted knot pair lies 6.6° from the prediction while its catalogue axis is 74.7° from it, so the anti class of
2197 does not come from the fitted pair's axis. The MACS1423 census conclusion (#282, P 0.259) draws random pair
orientations and does not use `pa_sep`, so it is unchanged. Labels: `pa_obs` observed; `pa_sep`, offsets derived.
