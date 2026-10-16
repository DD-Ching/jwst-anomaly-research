# `arcs` on RXJ2129 with the DJA v7.0 catalogue: result (2026-10-16)

Plan: `PLAN.md` (committed `ea2feb6` before the run). A re-test of #283's field (same maps and rule, different
shape source), not an independent field. Labels: shapes `observed`; conversion and classes `derived`; predicted axes
`model_prediction`; thresholds are ASSUMPTIONs (D-024 defaults).

## PA convention gate (`convention.json`)
68 rows matched (≤ 0.2″, e ≥ 0.5 and S/N ≥ 20 in both catalogues): median |ΔPA| 2.8° between PA = θ + 90° and the
pipeline `sky_orientation` (gate ≤ 10°: **passes**); 97.1 % within 20°. The alternative PA = 90° − θ gives 49.4°.
Pipeline pixel 0.03074″ (from its centroids).

## Result (`summary.json`, `*/arcs.json`)
| Model | Strong-shear background | Aligned (control ≥ 70 %) | Anti | Expected (shape noise) |
|---|---|---|---|---|
| Lenstool (primary) | 89 | 57 (64 %) | 10 | 1.513 |
| GLAFIC (second) | 89 | 60 (67 %) | 10 | 1.304 |

Both maps pass `validate`. **Frozen outcome: no result** (both controls fail). The anti counts carry no weight.

## Deviation
- `isophotal_abmag` is the catalogue `mag_auto`, not 28.9 − 2.5 log₁₀ `flux_iso` as PLAN.md said: `flux_iso` is in
  detection-image units (the planned formula gave m ≈ 33). `mag` enters no selection or class; the counts above are
  identical before and after the change (both runs, not stored for the first).

## Post hoc (not in the plan; does not change the outcome)
- `posthoc.json`: 17 of the 89 rows have z_phot > 15 (eazy fits at the grid edge; 187 of 3,295 catalogue rows); they
  hold 7 aligned and 5 of the 10 anti in both models. The other 72 rows: 50 aligned (69 %) and 5 anti (Lenstool), 53
  (74 %) and 5 (GLAFIC). No anti row of the earlier `arcs` runs on main has z_phot > 15 (checked over every stored
  `arcs.ecsv`; not stored).
- The 10 anti rows are faint (m_auto 26.6–27.7, one NaN; `*/arcs.ecsv`).
- `sheet_anti_aligned.png` (F200W program 2767 stamps; red observed axis, cyan predicted z = 2 axis, drawn through
  the WCS) was looked at. Aligned rows 111, 281, 427 and 440 show elongated sources along the cyan axis. Of the
  anti rows, 862 and 1451 show no source at the centre in F200W (a bright neighbour at about 1″ and 1.5″), and 902
  and 903 sit on a linear stripe. 1017, 1094, 1137, 1219 and 1411 show only faint or no elongated F200W emission.
  1645 sits at the end of a long streak. So the anti rows are mostly detections of the HST+JWST `ir` image without a
  clear F200W counterpart, or rows on linear structures; this is a by-eye reading.
- #283's background anti rows: 464 matches DJA 259 (0.04″), here `mixed`; 447 has no DJA row within 3.7″.

## Reading
No result on the S1063 replication question from RXJ2129 under the frozen rule, now from both catalogues. The DJA
`ir` shapes raise the row count (89 vs 12) but add faint rows without clear F200W counterparts; an `arcs` shape source
from a mixed-PSF detection image is not a good substitute for F200W shapes (by eye, above). Not a candidate.
