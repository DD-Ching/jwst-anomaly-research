# Abell 2744 unused NIRCam epochs: footprint check (frozen plan, 2026-10-26)

Written and pushed before any footprint was queried. Tool `scripts/a2744_epoch_footprint.py`; rule functions
tested in `tests/test_a2744_epoch_footprint.py`.

## Question
`results/transients/epoch_inventory` lists 78 Abell 2744 "new" pairs from eight unused epochs (1324 o001 / o006 /
o007, 2561 o002 / o006, 3073 o005 / o008 / o009). It used a 3′ cone, no footprint test, so 78 is an upper bound.
Which pairs actually image the same sky, and which of them cover the strong-lensing core, in which shared filters?
This is metadata only (`observed` footprints, `derived` areas). It decides whether a D-027 two-epoch grid on an
unused Abell 2744 pair is worth a run, and which pair.

## Ordinary explanations for "no useful pair"
GLASS 1324 NIRCam images are parallels and may lie off the cluster; 2561 / 3073 outer pointings may not overlap
the earlier mosaics; a shared filter may not exist in the overlapping part.

## Data
One batched MAST `Observations.query_criteria` for the Abell 2744 level-3 NIRCam `_i2d` observations of every epoch
in the inventory pairs (proposals 1324, 2561, 2756, 2883, 3073, 3516, 3538, 4111, 1176 and any other named in
`inventory.json`). Footprint = the `s_region` polygon(s) of each (epoch, filter) observation.

## Rules
1. Projection: gnomonic onto the tangent plane at the inventory's Abell 2744 centre (MAST name resolution of
   "Abell 2744", as stored by the inventory run, else (3.58750, −30.39750)).
2. Areas by point sampling on a 1″ grid over ±5′ (ASSUMPTION: grid step; area error ≲ 1 %). A point is in a
   footprint when it is inside any of its polygons (`matplotlib.path.Path.contains_points`).
3. Per pair and shared filter: overlap area A_ov (arcmin²) and the core overlap A_core = overlap within 1.0′ of the
   centre (ASSUMPTION: 1′ ≈ the strong-lensing region; the disc is π arcmin²).
4. Class per pair (best shared filter): **core** when A_core ≥ 0.5 π arcmin² (ASSUMPTION: half the core disc);
   **overlap** when A_ov ≥ 1 arcmin² (ASSUMPTION) but not core; **none** otherwise.
5. A filter counts only when both epochs have that filter's footprint.

## Outcomes
- **G (grid input):** ≥ 1 pair is core in ≥ 2 shared filters → name the pair with the largest A_core summed over
  filters (ties: longer Δt) as the next D-027 grid input.
- **O (overlap only):** no such pair, but ≥ 1 overlap pair → list them; a grid there is outside the core.
- **N (none):** no pair overlaps → the 78 are 3′-cone coincidences; the Abell 2744 count drops to 0.

Post hoc work (other radii, plots) is labelled post hoc. No pixel is read.
