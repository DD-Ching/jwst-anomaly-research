# VENUS batch 3, Abell 370 box 2: two rule-level real changes (vetting record)

Run: `results/transients/batch_venus3` (#311). Box `A370_2` (39.968731, −1.582667; 1208 o002 2022-12-29 → VENUS o053
2026-01-09) **failed calibration** (vanish 3/7), so it is not reliable coverage. Evidence: `boxes/review/review.json`,
`lightcurves.ecsv`, `post_review_checks.json`, `final_classes.json`, `crossmatch.json`; sheets and stamps looked at.
Vetted 2026-10-17. Every number is `derived` unless labelled. Neither entry is unexplained or an exotic candidate.

## B3-A370-box2-h003 (39.966687, −1.583019)
| Test | Outcome | Evidence |
|---|---|---|
| Edge, spike, bright neighbour | does not explain | inside the box; no spike on the sheet; F150W/F444W S peaks 0.015″ apart |
| Star | does not explain | no Gaia DR3 source within 15″; HFF-DeepSpace S/G 0; a galaxy spectrum (below) — not a star |
| Catalogued object | **a compact galaxy** | [LRB2019] 6401 at 0.08″, z_spec 0.3604 (MUSE, absorption lines; confidence 2) |
| Change in the search epochs | F444W 255.4 ± 3.3 → 219.4 ± 5.9 nJy (z −5.4); F150W −1.7 | `review.json` |
| Independent epoch | **does not confirm** | 5324, 18 d after VENUS: F444W 247.4 ± 7.5 nJy (consistent with 2022); F150W stable 596–646 nJy over 6 epochs |
| Programme-to-programme scatter | present | ≥ 3σ differences of both signs between mosaics (F200W −30.3, F210M +8.1, F360M +6.5) |
- **Status:** rule-level real change, not confirmed by the independent epoch; consistent with programme-to-programme
  photometric scatter of a compact bright galaxy. Variability (e.g. a weak AGN) is not excluded but has no support.
- **Cheapest discriminating next check:** matched compact-source controls in the same mosaics (`hit_lightcurve_controls.py
  pick` + `compare`) on the 1208 → VENUS F444W pair.

## B3-A370-box2-h010 (39.970466, −1.585253)
| Test | Outcome | Evidence |
|---|---|---|
| Point source in any epoch | **none visible** | stamps of every covering mosaic 2022–2026 (50 forced measurements, F090W–F480M) |
| Position | on the Dragon arc | lensed images of system 2 (z 0.7251) at 1.5–1.9″; bright cluster galaxy 05674 (z 0.3677) at 2.8″ |
| Known transients nearby | many | Palencia+2026: 33 microlensing-event detections within 1″ (nearest 0.54″), 2022-12 to 2024-09; Fudamoto+2025 |
| Flux history | ≈ 30–50 nJy (2022-12, 2024-09), ≈ 8–16 nJy (2024-12 → 2026-02), four bands | `lightcurves.ecsv`; the 5324 drop coincides with a different orientation of its first visit |
| Steep-gradient photometry | not excluded | the 0.15″ aperture sits on the arc's steep edge (cf. batch-1 h003) |
- **Status:** rule-level real change with two open ordinary explanations: a caustic-microlensing event in the Dragon arc
  (a known population there) or aperture photometry on a steep gradient. Identity not established; no catalogued event
  within 0.5″. TNS could not be queried.
- **Cheapest discriminating next check:** matched arc-pixel controls (`hit_lightcurve_controls.py pick-arc` + `compare`)
  on the 1208 → VENUS F444W pair; if the change exceeds the control scatter, compare its epochs and position with the
  Palencia+2026 event list.
