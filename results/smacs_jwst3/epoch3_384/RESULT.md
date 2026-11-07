# Third epoch (12774) for SMACS `f200w_384`: frozen rule undetermined; post hoc, the 2026 position is reproduced

Plan: `PLAN.md` (frozen in the claim commit `4b98de3`, before any 12774 centroid was measured). Script:
`scripts/smacs_epoch3_384.py` (reuses `scripts/smacs_epochs_img.py` unchanged; `decide` unit-tested in
`tests/test_smacs_epoch3_384.py`). Output: `summary.json`. Contact sheet: `contact_sheet.png` (looked at).

## Data (observed)
- Program 12774 o002 t032 ("Taking aim at Precision Cosmology with SN Atalanta", PI Agrawal), F150W and F444W
  `_i2d.fits` by S3 byte range, observed 2026-07-06, 31 d after 6882 o057; epochs 2736 o001 (2022) and 6882
  o057 (2026-06-05) re-cut the same way.
- Contact sheet (0.6″, stored orientation of each mosaic): one compact isolated source in all 6 cutouts, the
  centroid on the peak; the 12774 F444W cutout is noisier (816 s vs 1245 s in 6882).

## Validation (simulation)
Synthetic +20 mas shifts on the 12774 peer cutouts are recovered to ≤ 0.043 mas (n = 10 per band). Passed.

## Result (derived; frozen rule)
| band | peers | 2736 → 12774 residual (mas) | peer median / p95 | angle to 2736 → 6882 residual | |
|---|---|---|---|---|---|
| F150W | 10 | 25.87 (−23.27, +11.31) | 5.32 / 21.74 | 4.3° | above p95 |
| F444W | 10 | 21.17 (−19.65, +7.89) | 14.52 / 22.07 | 12.5° | below p95 by 0.9 mas |

**Decision: undetermined** (F444W misses the p95 condition), the same split as the 2736 → 6882 run
(25.03 / 22.61 mas).

## Post hoc (not part of the rule)
- 6882 → 12774 residual of `384`, same peers: 1.69 mas (F150W) and 3.48 mas (F444W), below the peer median
  (5.43 / 9.19 mas) in both bands. The two 2026 programs (different dither, drizzle grid, exposure) place the
  source at the same position, ~23 mas from its 2022 position, in both bands. A 6882-specific systematic
  (H-6882) is therefore disfavoured.
- The four 2022 → 2026 residual vectors (two programs × two bands) point within 13° of each other.

## Reading
- Remaining explanations: real motion of ≈ 6 mas/yr (H-move; would make `384` a likely Galactic star), or a
  2736-epoch systematic at this position (shared by every comparison here). The frozen rule stays undetermined
  because the F444W peer p95 (22.07 mas) is set by 1 of 10 peers.
- Record label: "compact, unresolved; 2022 → 2026 shift ≈ 23 mas reproduced by two 2026 programs (post hoc);
  frozen rule undetermined; needs spectroscopy or a 2736-independent epoch". MAST (2″ cone, 2026-11-07)
  lists only 2736, 6882 and 12774 NIRCam imaging and one 9276 NIRSpec image-type exposure (`OPAQUE;MIRROR`,
  2026-06-05; not inspected) here, so no 2736-independent early epoch exists yet.

## Limitations
- 10 peers per band (12774 covers fewer of the 16 used before); p95 rests on one peer. Peers are not confirmed
  stars or galaxies.
- All comparisons share the 2736 epoch; a 2736 astrometric or drizzle artefact local to this position is not
  tested.
- No flux comparison here.
