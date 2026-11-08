# `arcs`: clump-pair census of the pooled anti rows (exploratory)

Plan: `PLAN.md` (committed before the run). Script: `scripts/arcs_pooled_census.py` (main's `arcs_vet_census.py`,
#273's rules, unchanged); every number below is in `summary.json`, `posthoc.json` or `census.ecsv`. Pixels: 8″ F200W
S3 cutouts of each field's screen mosaic (not stored). Contact sheet: `census_anti.png` (looked at).

## Frozen outcome: **anti excess explained by clump orientation (exploratory)** (derived)
- 104 strong-shear background rows in the 5 variant-B fields; 73 pass the own-flux gate (73 / 104 ≥ 0.5, ≥ 12); 19 of
  them are `two-knot`, f̂ 0.26 (95 % upper 0.376).
- The shape-noise sum reproduces the pooled tally (2.273); under the clump-aware null the expectation is 5.538 and
  P(≥ 6) is 0.484 (shape-noise check 0.0227) → ≥ 0.05, so the frozen label is "explained".
- Sensitivity (mixture for every row): P 0.477 at f̂, 0.707 at the 95 % upper f. Secondary: P(≥ 4) for the 4 rows
  that persist under second models 0.789 (deviation: `summary.json`'s 0.824 left the 2 non-persisting anti rows in the
  null; the review fix drops them from count and null, as #273's fix 4, `posthoc.json`); variant A with MACS1423's stored census rows, P(≥ 10) 0.293.

## What the anti rows are (observed pixels, derived classes; inputs in `posthoc.json` `anti_rows_inputs`)
| Field | Row | mag | z_phot | g_min | Class | Gate |
|---|---|---|---|---|---|---|
| MACS0417 | 1993 | 27.2 | 3.03 | 3.69 | undecided | failed |
| MACS0417 | 2341 | 27.1 | 3.15 | 0.26 | two-knot | passed |
| MACS0417 | 4071 | 27.4 | 8.12 | 1.18 | two-knot | failed |
| WHL0137 | 2194 | 24.4 | 1.17 | 0.24 | single | passed |
| RXC J0600 | 581 | 21.2 | 0.70 | 0.39 | single | passed |
| MACS0308 | 2157 | 20.6 | 0.76 | 0.64 | single | passed |

Contact sheet (by eye; deviation: it shows the 1″-filtered stamp without the fitted knot positions the plan named,
and its 33-px median filter on a 65-px stamp is a display filter, not the fitted residual): the MACS0417 rows are faint and noisy; 4071's cutout shows vertical striping (the stripe
pattern of NIRCam 1/f noise is a likely reading, not tested). 2194 is a thin faint streak, 581 an elongated disc with
a bright core, 2157 a smooth elongated galaxy; none of the three shows two knots.

## Post hoc (not in the plan; `posthoc.json`)
Anti observed vs expected within each class:
| Class | n | Anti | Expected | P(≥ anti) |
|---|---|---|---|---|
| post hoc: gate-passed two-knot | 19 | 1 | 2.866 | 0.959 |
| post hoc: gate-passed other | 54 | 3 | 0.631 | 0.024 |
| post hoc: gate failed | 31 | 2 | 2.041 | 0.620 |

The frozen "explained" comes from the null's inflation (19 two-knot rows at P_cp, 31 gate-failed rows at the f̂
mixture), not from the anti rows being clump pairs: only 1 of the 6 is a gate-passed two-knot, and the two-knot rows
hold fewer anti than P_cp predicts (post hoc 1 vs 2.866, P(≤ 1) 0.189). The 3 single-component anti rows (2194, 581, 2157, the
3 non-MACS0417 rows that persist under second models) stay above shape noise within their class (P 0.024, post hoc).
581 and 2157 are bright (mag 21.2 / 20.6) with low background z_phot (0.70 / 0.76; cluster redshifts not stored here), so photo-z
contamination (a foreground or member galaxy counted as background) is the obvious next test; not tested here.

## Reading (narrowest statement)
- Under #273's frozen clump-aware null the pooled variant-B anti count is not in excess (P 0.484; exploratory).
- That null does not explain the 3 bright single-component anti rows (post hoc P 0.024 within their class); they are
  the open part of the pooled excess. No candidate.

## Limitations
- Exploratory: the census idea is post hoc on S1063 and the pooled counts were seen before the plan.
- P_cp assumes isotropic knot separations; the two-knot rows hold fewer anti than it predicts (post hoc), so the
  clump-aware null may overstate anti for clump pairs, which would make the frozen "explained" too lenient.
- 31 / 104 rows fail the gate (MACS0417 22 / 49), mostly faint; post hoc, their f̂ mixture carries 2.041 of the 5.538 expected.
- One band (F200W); no second band or spectra.
