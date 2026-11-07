# Spike-geometry check for single-band elongated detections: no rule adopted; the known streak lies on a long linear feature, not a neighbour's spike

Plan: `PLAN.md` (rule and addendum B frozen before each run). Scripts: `scripts/spike_geometry.py` (rules as
`spike_flags` and `streak_verdict`, unit-tested before the runs), `scripts/spike_geometry_sheet.py`. Pipeline runs of
`configs/reference_sample.yaml` at main `1921bc2`: CEERS `20261107T081028Z-c5afef2f`, SMACS `20261107T081028Z-369d1097`
(jwst 3.0.0 catalogs). Numbers: `summary.json` (rule A), `summary_streak.json` (addendum B). Catalog columns and
pixels `observed`; flags, rates and masks `derived`; every threshold an ASSUMPTION.

## Rule A: radial alignment to a ≥ 3 mag brighter neighbour within 10″ (frozen) — not adopted
- **Positive missed.** CEERS `f200w_2842` (single band, ellipticity 0.952) is not flagged: its best-aligned bright
  neighbour (`f200w_2834`, 6.18″, 22.81 mag) is 13.8° off its major axis (limit 10°). The 19.3 mag source 4.6″ away
  (`f200w_2879`), the one the 2026-10-08 entry named, is 32.7° off (`positive.bright_neighbours`).
- **At chance.** Flags vs the 90°-rotated null: CEERS 31 vs 28, SMACS 200 vs 93 (plan needs > 3×).
- **Arcs.** 0 of the 36 ICLv2 images matched within 0.3″ (of 60) are flagged, and 0 under the geometry-only conditions
  (post hoc). But the one gate-passing SMACS flag, `f444w_1970`, lies on a thin arc segment (contact sheet).

## Addendum B: W1 `findsat_mrt` streak mask (frozen) — not adopted, both arms
| arm | positive within 1″ | matched ICLv2 sources in mask (of 36) | CEERS single / ≥3-band in mask | SMACS single / ≥3-band in mask |
|---|---|---|---|---|
| trail (W1 settings) | no | 0 | 505 / 4392 vs 111 / 2457 (2.5×) | 441 / 3314 vs 315 / 2687 (1.1×) |
| trail + row segments 64 | yes | 24 | 2050 / 4392 vs 556 / 2457 (2.1×) | 2494 / 3314 vs 1315 / 2687 (1.5×) |

The trail arm misses the line through `f200w_2842` (mask sheet `ceers_t021_trail_mask.png`: 3 diagonal or vertical
trails) and enriches single-band sources below 3×. The segment arm covers 47 % / 75 % of single-band sources and
24 of the 36 matched ICLv2 sources (40 of all 60 catalogued image positions), so it is a stripe/edge map, not a
streak flag. Code review caught the first run counting all 60 positions instead of the plan's matched images; the
committed numbers are the corrected run (verdicts unchanged).

## Looked at (`contact_sheet.png`, 16″ cutouts)
- `f200w_2842` sits on a straight horizontal line that crosses the whole 16″ cutout in F200W and, fainter, in F356W, and
  does not pass through the bright source 4.6″ away. A radial spike of that neighbour is therefore not supported; an
  off-field spike or a trail on the simultaneous SW/LW exposure is the reading (not tested further).
- CEERS `f356w_1028` (the one gate-passing CEERS flag): a faint compact source, no streak visible.
- SMACS `f444w_1970`: on a thin arc segment, a false positive of rule A.

## Conclusion
No geometric spike or streak flag is adopted. D-014's confirmation gate stays the only automatic spike control; a
single-band elongated row in a top k still needs a visual check. The 2026-10-08 reading "spike from the bright
source 4.6″ away" is replaced by "on a long linear feature" (visual).

**Revisit if** a per-observation spike-angle model (from the pointing's V3 position angle) is available to test
alignment with the actual spike directions of off-field stars, or single-band streaks enter a galaxy top 20 again.
