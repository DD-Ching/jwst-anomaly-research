# W1 row/column-segment stripe veto for the CANUCS flanking fields: result (2026-10-18)

Plan: `PLAN.md` (frozen and committed before any mask). Code: `scripts/w1_streak_mask.py --seg-len`
(`segment_stripes`, `to_mask_grid`, `cost_ok`), `scripts/w1_row_veto_eval.py` (K1r, cost check, adoption).
Labels: pixels and shapes `observed`; masks, S and p `derived`; every threshold an ASSUMPTION.

## Outcome: the route fails on the cost check (both variants); the headline stays at fourteen fields
`adoption_A.json`, `adoption_B.json` (decided by `w1_row_veto_eval.adopt`):

| Variant | K1r (macs0417-ncf peak covered) | NCF passing the cost check (kept ≥ 75 %) | adopted |
|---|---|---|---|
| A (L = 80 block pixels, 10″) | yes | 1 of 5 (macs1423-ncf) | no |
| B (L = 40, 5″) | yes | 0 of 5 | no |

Rows kept by the screen (`kept_union` / `kept_d097`, `fields.*` in the adoption files):

| Field | D-097 | A | A share | B | B share |
|---|---|---|---|---|---|
| macs0416-ncf | 3,445 | 2,370 | 0.688 | 1,936 | 0.562 |
| macs1149-ncf | 3,325 | 2,386 | 0.718 | 1,930 | 0.580 |
| abell370-ncf | 3,155 | 2,196 | 0.696 | 1,873 | 0.594 |
| macs0417-ncf | 3,829 | 2,603 | 0.680 | 2,143 | 0.560 |
| macs1423-ncf | 3,348 | 2,517 | 0.752 | 2,252 | 0.673 |

Mask pixels (0.25″ grid; `masks_<V>/*.json`, `mask_pixels`): A 33,064–58,713, B 41,033–66,965, against
0–23,297 for the D-097 trail masks alone (`segment_stripes.trail_mask_pixels`). No injection run (the plan runs
them only for an adopted variant).

## What the masks and screens show (looked at; post hoc, not tests)
- Sheets (`masks_A/macs0417-ncf_streak_mask.png`, `masks_B/...`; the title still says "status-2 trail mask", but
  the red overlay is the union mask; `masks_A/macs0416-ncf` looked at in the session, not stored): in A most
  flagged segments sit on the halos and spikes of bright stars and galaxies and along detector edges; B also flags
  many short row segments across the striped upper detectors of macs0417-ncf. How the rows lost split between
  halos, edges and stripes is not measured.
- What covers the macs0417-ncf peak for K1r: in the A mask a horizontal row-segment band runs through the peak's
  mask row (plus a short column segment 4 pixels, 1″, to its right); in B a horizontal band starts 4 pixels (1″) below (session printout of the mask
  around the peak, not stored).
- With either union mask the macs0417-ncf screen no longer passes the E/B rule (S 3.15, p 0.95 under A; 3.21,
  p 0.835 under B; `screen_union`). That supports the stripe explanation, but only at a 32–44 % row cost, so it
  is not a frozen explanation.
- macs0416-ncf still passes under both masks (S 4.60 / 4.05, p 0.0); it is already explained by the frozen V2
  (spike) rule (`results/w1_shear_ncf/RESULT.md`).
- abell370-ncf passes under A (S 4.62, p 0.0, max |S_×| 4.04) though it does not under the D-097 mask (S 3.94,
  p 0.065) or B (S 3.28, p 0.61). The veto is not adopted, so this is not a flag; it shows how much S_max moves
  when a third of the rows go (exploratory; not vetted).

## Failed approach and next step
This field-wide median-segment stripe detector at 5σ (both segment lengths) costs 25–44 % of the shape rows in
these mosaics; on the sheets it also flags bright-source halos and detector edges. Before another stripe veto: mask the
segmentation footprints (or a bright-source halo model) before taking the segment medians, and require a
stripe to extend over ≥ 2 adjacent segments of the same row. That needs its own frozen plan. Alternatively
leave the two NCF out and add area from other blank fields (TASKS).
