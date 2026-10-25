# SN Eos (MACS J1931 system 101) predicted images and delays: result (2026-10-25)

Plan `PLAN.md` (frozen in 89c06718 before any map value was read). Tool `scripts/sneos_predict.py`, tests
`tests/test_sneos_predict.py` (SIS delay against the analytic value, potential-scale fit, matching), committed before
the run. Numbers from `summary.json`. Every image row is a model_prediction of the VENUS AstroLensPy v2.1 model;
101.1 / 101.2 positions are observed (catalogue).

## Gates (frozen rules)
- **P passes:** `Lens_pot` is the potential of `Alpha_x/y` on the same grid: fitted s = 1.0000004, median relative
  residual of s·∇ψ vs α 7.0e-5 (rule ≤ 0.05) over 3,141,601 pixels within 40″ of 101.1.
- **S passes:** both catalogued images matched one-to-one within 0.5″.
- **M passes:** our |μ| over the `Mu_zs5.13` map value at the matched images: 0.998 (101.1), 1.000 (101.2).

## Predictions (z 5.133, mean source of 101.1 / 101.2; H0 70, Ωm 0.3)
| Image | Offset from cluster frame centre (x W, y N; ″) | μ (signed) | Δt vs 101.1 (d) |
|---|---|---|---|
| far counter-image | (−36.2, 0.4) | +2.17 | −19,656 |
| leading pair a | (17.6, −7.2) | +14.48 | −1,321 |
| leading pair b | (11.8, −11.6) | −11.92 | −1,238 |
| 101.2 | (6.9, 5.4) | −28.70 | −0.32 |
| 101.1 | (6.5, 4.5) | +24.70 | 0 |

Five images, the number Allingham et al. 2026 state (count only; their positions and delays were not read). No
extra image satisfies the plan's target rule (|Δt| ≤ 400 d and |μ| ≥ 2): every extra image **leads** 101.1, by
1,238 d, 1,321 d and 19,656 d (about 54 yr, not stored).

## Reading
- If 101.1 was near peak at the VENUS epoch (2025-09-01; the explosion date is not fixed here, an assumption), the
  leading pair arrived around 2022-01 to 2022-04 (post hoc; dates not stored). `results/transients/epoch_inventory/inventory.json`
  lists no MACS J1931 NIRCam epoch before 2025-09-01, so this model's leading images cannot be tested on NIRCam
  epochs; whether any other archival epoch covers that window was not checked.
- Per the plan, the extra-image delay test drops to "watch for new epochs"; in this model no SN image arrives after
  101.2, so new epochs can only test the 101.1 / 101.2 pair and the absence of unpredicted point sources.
- Post hoc (a next-stage plan, not tested): 101.1 and 101.2 are 0.32 d apart in this model, so their SN flux ratio
  across the three epochs is fixed by |μ_101.2 / μ_101.1| ≈ 1.16 (not stored; from the table) independent of the light curve,
  up to microlensing and the merging-pair μ sensitivity. That is a flux-ratio test the static arcs cannot give.
- The leading pair sits beside a cluster member's galaxy-scale critical curve (μ-map contact sheet, viewed in the
  session, not stored), so its μ and Δt depend on that member's mass; no uncertainty was propagated (the ±1σ maps
  were not used, as planned).

## Limits
- No Δt or μ uncertainty; one model; catalogue z 5.133 (the MUSE z 5.13). Data: four map zips (≈ 270 MB) held in the
  session scratchpad, md5-checked against the record, deleted after the run.
