# Image-level two-epoch centroids of SMACS `f200w_1816` / `f200w_384`: 1816 no motion detected, 384 undetermined

Plan: `PLAN.md` (frozen in the claim commit `c60d508`, before any centroid was measured). Script:
`scripts/smacs_epochs_img.py` (`decide` encodes the frozen rules; unit tests in `tests/test_smacs_epochs_img.py`).
Output: `summary.json`. Contact sheet: `contact_sheet.png` (both epochs' centroids marked; looked at).

## Data (observed)
- Level-3 `_i2d.fits` of 2736 o001 (2022) and VENUS 6882 o057 (2026), F150W and F444W, read by S3 byte range
  (1.5″ cutouts); CATs from `s3://stpubdata` (same products as `../epochs_1816_384/`).
- Contact sheet (0.6″ panels): one compact, isolated source at each target position in all 8 cutouts; the
  centroids sit on the source peak. The markers are raw positions before the frame tie (tie 14.5–25.5 mas, mostly
  in Dec, `tie_mas`), so their visible separation is not the residual.

## Validation (simulation)
- Shifting each peer's epoch-2 cutout by +20 mas in RA or Dec (Fourier shift) is recovered with a median error
  ≤ 0.031 mas in both bands (n = 30 per band; `summary.json` `validation`, `max_abs_error_of_median_mas`). Passed. This tests the centroid's
  response to a shift only; the noise floor comes from the peer null below, not from this test.

## Result (derived; frozen rule)
| target | band | peers | residual (mas) | peer median / p95 (mas) | peers with larger residual | decision |
|---|---|---|---|---|---|---|
| `1816` | F150W | 14 | 1.53 | 6.26 / 18.76 | 14 / 14 | **no motion detected** |
| `1816` | F444W | 14 | 2.45 | 9.22 / 26.04 | 14 / 14 | |
| `384` | F150W | 16 | 25.03 | 8.91 / 21.31 | 0 / 16 | **undetermined** |
| `384` | F444W | 16 | 22.61 | 10.15 / 27.60 | 3 / 16 | |

- `1816`: residuals below the peer median in both bands. The 95th-percentile peer residual over 4.0 yr is
  4.69 mas/yr (F150W) and 6.51 mas/yr (F444W) (`peer_p95_motion_mas_per_yr`): motion above about these values
  would likely have exceeded most peers; smaller motion is not excluded. H-star is disfavoured at ≳ 5–6.5 mas/yr
  only.
- `384`: F150W above every peer, F444W below the peer p95, so the frozen rule returns "undetermined".
- **Post hoc** (not part of the rule): the `384` residuals in the two bands point the same way, 12.3° apart
  ((−23.26, +9.23) and (−22.31, +3.68) mas, dRA cos Dec, dDec); for `1816` they are 61.3° apart, at a residual
  below the peer median. The catalog route's F150W residual for `384` was (−19.98, +6.29) mas
  (`../epochs_1816_384/summary.json`), the same image in both routes, so it is not an independent check.
  A ~23 mas shift over 4.0 yr would be ~6 mas/yr, of the order H-star predicts; it is suggestive, not established.

## Reading
- `f200w_1816`: no motion detected at the stated bound; the star route is disfavoured above ~5–6.5 mas/yr. Record
  label: "compact, unresolved; no motion detected; needs spectroscopy".
- `f200w_384`: undetermined. Per the plan the astrometric route ends here; the post hoc agreement of the two
  bands makes a later third epoch (if one exists) the cheapest decisive test, not a new analysis of these two.

## Limitations
- 14 / 16 peers per target; the p95 rests on 1–2 peers. Peers are compact sources of similar magnitude and
  size, not confirmed stars or galaxies; a moving peer would inflate the null.
- Different dither patterns and drizzle grids between programs, and an asymmetric source, could shift a
  centroid without motion; the same-band comparison limits but does not remove this.
- No variability measurement here (the catalog route found |Δm| 0.21 and 1.36 × the reference σ in F150W).
