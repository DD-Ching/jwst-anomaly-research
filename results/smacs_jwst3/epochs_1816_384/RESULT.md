# SMACS `f200w_1816` and `f200w_384` over two epochs (2736 o001, 2022 → VENUS 6882 o057, 2026): undetermined by the frozen rule

Plan: `PLAN.md` (frozen in the claim commit `ec9c112`, before any catalog was read). Script:
`scripts/smacs_epochs_1816_384.py` (`decide` encodes the plan's rules; unit tests in
`tests/test_smacs_epochs_1816_384.py`). Output: `summary.json`. Cutouts: `scripts/smacs_epochs_1816_384_sheet.py`,
`contact_sheet.png` (north-up, brightness shared per band; looked at), `cutouts.ecsv`.

## Data (observed)
- Level-3 CATs read from `s3://stpubdata` (the MAST portal API returned "Request failed" this session). sha256 in
  `summary.json` via the D-017 tool. The 2736 F150W/F200W hashes equal `data/manifests/smacs_jwst3_gaia_only.ecsv`
  (jwst 3.0.0); the VENUS F150W/F444W hashes equal SOURCES.md "Second epochs".
- Shared broad bands: F150W, F444W.

## Result (derived)
- **Frozen rule outcome:** motion `undetermined`, variability `undetermined` for both targets, because neither
  has an F444W measurement (`decide`, `summary.json`): `1816` has an epoch-1 F444W row (0.017″, 26.99 mag) but no
  epoch-2 row within 2.28″; `384` has no F444W row within 5.2″ (2022) or 6.2″ (2026).
- **Cutouts (post hoc, visual, not stored):** a compact source is visible at the target position in all four F444W
  cutouts (`contact_sheet.png`). No F444W flux or S/N was measured, so why the CATs lack these rows is not established.
- **F150W alone** (one band; the plan needs two):
  | | residual (mas) | D-017 significance | tool error (mas) | Δm − ref. median | ref. Δm σ |
  |---|---|---|---|---|---|
  | `1816` | 6.89 | 5.18 | 1.33 | +0.015 | 0.072 |
  | `384` | 20.95 | 17.41 | 1.20 | +0.092 | 0.068 |
- **Post hoc empirical null (F150W):** mutual matches within 60″ with |Δmag| ≤ 0.5 and `semimajor_sigma` within
  ±30 % of the target's, both target rows excluded, same frame tie as the D-017 tool (8 peers each; thresholds
  in `summary.json`): peer residual median 12.97 / 11.41 mas, 95th percentile 19.10 / 22.79 mas. `1816` exceeds
  2 of 8 peers; `384` exceeds 7 of 8 peers, which is not significant with n = 8. The tool error is about one
  tenth of the peers' median residual (D-017 addendum), so the D-017 significances are not evidence of motion here.
- **Variability (F150W):** |Δm| is 0.21 and 1.36 × the reference σ, below the plan's 5×.
- **Deviation:** `decide` also returns "not variable" when no measured band reaches 5×; the plan names only
  "variable". It does not affect this run, where both outcomes are "undetermined".

## Reading
- No motion or variability is established. Motions up to the peer 95th percentile, 4.8 / 5.7 mas/yr over 4.0 yr
  (post hoc, 8 peers, `p95_motion_mas_per_yr`), are not excluded, so H-star is neither supported nor excluded.
- The records stay open. Narrowest statement: "two-epoch F150W catalog positions show no residual that
  distinguishes these targets from 8 faint compact peers".

## Limitations
- 8 peers per target is a small null; the peer selection was chosen after seeing the tool's significances.
- No image-level F444W centroids were measured (the plan used catalog rows only).

## Next
- Measure F444W (and F150W) centroids from the S3 cutouts in both epochs with a peer-calibrated frame tie
  (new plan; it would complete the frozen rule's second band), or retire the two records as "compact, unresolved;
  needs spectroscopy".
