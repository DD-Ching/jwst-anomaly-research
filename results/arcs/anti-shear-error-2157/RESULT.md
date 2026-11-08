# `arcs`: MACS0308 2157 under the Zitrin LTM-Gauss realisations (exploratory, 2026-11-08)

Plan: `PLAN.md` (frozen before any realisation value was read). Script `scripts/arcs_shear_error_2157.py`; strip
reader and gates unit-tested before the run (`tests/test_arcs_shear_error_2157.py`); rule functions are addendum 13's,
imported unchanged. Every number below is in `summary.json`. Maps read as HTTP byte ranges from the RELICS HLSP
(DOI 10.17909/T9SP45; strip sha256 in `summary.json`), nothing stored. Observed axis and photo-z from the deciding
GLAFIC run's stored `arcs.ecsv`; predicted axes `model_prediction`; classes and labels `derived`; thresholds
ASSUMPTIONs (`summary.json` `params`).

## Frozen outcome: **model error does not remove the tested rows** (2157 `robust`)
| Gate / quantity | Value | Threshold | Result |
|---|---|---|---|
| V1 pixel × 0.06″ vs arcsec best-fit maps (strip j 803–2246, i 681–2066) | max 7.6 × 10⁻⁶″ | ≤ 1 × 10⁻³″ | passes |
| V2 stored Zitrin run reproduced from strips (18 strong-shear background rows) | classes equal, max Δoffset 0.0° | equal, ≤ 0.5° | passes |
| Map gate, κ strip, map000 / 049 / 099 | median rel. diff 7.9 / 7.5 / 7.2 × 10⁻⁶ | ≤ 0.05 | passes |
| Control retention (17 GLAFIC `aligned` rows) | median 1.0 (range 0.941–1.0) | ≥ 0.70 | passes |
| 2157 `anti` share f | 1.0 (100 / 100; f_aligned 0.0) | robust ≥ 0.68 | **`robust`** |

Secondary: offset_min 16 / 50 / 84 % 78.9 / 79.1 / 79.2°; g_min 1.189 / 1.247 / 1.323. Best fit (strip) 79.2°, g_min
1.189; deciding GLAFIC run 82.0°, g_min 0.64.

## Reading (narrowest statement)
- MACS0308 2157 stays `anti` in every Zitrin LTM-Gauss realisation, with a predicted offset_min 18.9–19.2° past the 60° cut (16–84 %).
  The realisations hardly move the prediction at this position (16–84 % span 0.3°), so this samples one parametric
  model's statistical error only, not the systematic spread between models; two independent models (GLAFIC, Zitrin)
  agree on `anti` (addendum 8).
- Contact sheet looked at again (`../anti-bands/anti_bands.png`): 2157 is a single bright (mag 20.61, `../macs0308-relics/summary.json`) compact
  elongated galaxy; this test changes no shape.
- Not a lensing anomaly claim and not a candidate. The open set stays 2 rows (RXC J0600 581, MACS0308 2157), both
  `robust` to the published models' statistical error.

## Deviations and limits
- None from the plan. The strip route (byte ranges, sliced WCS) is new here and was validated by V1 / V2 before use.
- Zitrin realisation spread is narrow at the row; the published GLAFIC v2 has no realisations here. Both models
  share HST-era constraints. Membership is still open (no public spec-z, addendum 10); mag 20.61 at z_phot 0.76 is
  bright for a background source, which a spectrum would settle.

## Next
Membership of 581 / 2157 only with a spectrum; a further independent model of MACS0308 if one is published.
