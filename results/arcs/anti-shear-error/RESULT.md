# `arcs`: model shear error at the single-component anti rows (exploratory, 2026-11-08)

Plan: `PLAN.md` (frozen before any realisation map or further model was read). Script `scripts/arcs_shear_error.py`;
rule functions unit-tested before the run (`tests/test_arcs_shear_error.py`). Every number below is in
`summary.json`. Maps streamed from the RELICS HLSP (DOI 10.17909/T9SP45), sha256 in `summary.json`, deleted after
use. Observed axes and photo-z from the deciding runs' stored `arcs.ecsv`; predicted axes `model_prediction`;
classes and labels `derived`; thresholds ASSUMPTIONs (`summary.json` `params`).

## Frozen outcome: **model error explains WHL0137 2194**; RXC J0600 581 `robust`; MACS0308 2157 untested
| Row | Test | Map gate (κ median rel. diff ≤ 0.05) | Control retention (≥ 0.70) | Row result | Rule 3 label |
|---|---|---|---|---|---|
| WHL0137 2194 | RELICS Zitrin LTM-Gauss v1 | 1.4 × 10⁻⁵, passes | 13 / 13 (1.0), passes | `mixed`, offset 48.3°, g_min 0.17 | |
| WHL0137 2194 | RELICS WSLAP v1 | 1.0 × 10⁻⁵, passes | 0.923 of 13, passes | `mixed`, offset 53.5°, g_min 0.24 | **`model_sensitive`** (0 / 2 anti) |
| RXC J0600 581 | RELICS Lenstool v1, 100 `range/` realisations | map000 / 049 / 099: 0.00021–0.00028, pass | median 0.727 of 11 (range 0.545–0.818), passes | `anti` in 100 / 100; offset_min 16/50/84 % 74.6 / 75.3 / 76.3°; g_min 0.902 / 1.114 / 1.379 | **`robust`** (f 1.0) |

For reference (stored runs, not re-run): 2194 is `anti` in the deciding Lenstool run (offset 63.1°, just past the
60° cut, g_min 0.236) and in GLAFIC (addendum 8); so 2 / 4 published models put it past the cut. 581 is `anti` in the
deciding GLAFIC run (g_min 0.39) and the Lenstool best fit (addendum 8).

## Reading (narrowest statement)
- WHL0137 2194: its anti class depends on the model. Two further independent models (one free-form) predict an
  offset of 48.3° and 53.5°, below the 60° anti cut, against 63.1° in the deciding Lenstool run; the row is at the
  weak-shear end (g_min 0.17–0.24). It leaves the open set as `model-dependent` (H_err supported for this row).
- RXC J0600 581: the Lenstool model's own statistical error does not move it out of `anti` (offset ≥ 74.6° in the
  16 % quantile). Lenstool predicts much stronger shear here (g_min median 1.114) than the deciding GLAFIC run
  (0.39), so |g| is model-dependent while the anti class is not, in the two models tested.
- Contact sheet of the 3 rows looked at again (`../anti-bands/anti_bands.png`): isolated single objects; this test
  changes no shape.
- Not a lensing anomaly claim and not a candidate. The open set is now 2 rows (581, 2157), with 2157 untested here.

## Deviations and limits
- MACS0308 2157 untested (stated in the plan): its realisations are 36 MB pixel-unit maps; next route byte ranges.
- Realisations cover one model's statistical error only (Lenstool for 581; the deciding GLAFIC has none); the
  Lenstool range maps are coarse (κ check on 20–21 comparable points). Control retention under the realisations
  passes narrowly (median 0.727; minimum 0.545).
- WHL0137: 2 further models, best fits only; WSLAP and Zitrin maps are taken as scaled to D_LS/D_S = 1 (RELICS readme
  v3, ASSUMPTION). Models share HST-era constraints, so they are not fully independent.

## Next
MACS0308 2157 under the Zitrin `range/` realisations via byte-range reads; membership of 581 / 2157 only with a
spectrum.
