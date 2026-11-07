# Plan: SMACS point-source screen v2 (faint locus window, spatial held-out check)

Written 2026-11-07 ~14:15Z, before the v2 run. Follows `results/smacs_jwst3/point_source/RESULT.md` (#766, D-022
addendum 6): the v1 screen failed V1 (held-out star recovery 0.7333 / 0.9057 < 0.90); 9 of its 14 locus outliers
had ref_mag < 20.

## What is not independent (stated before the run)
- The new window start (ref_mag 20) was chosen **after** seeing v1's outliers. By v1's stored outliers, 5 of them
  lie in [20, 24] (3 at 20.4–20.7, of which 2 sit within 0.15 px of v1's cut). A parity re-test on the same stars
  (V1a) is therefore expected to pass and is weak evidence; it is kept only as a necessary condition.
- New tests not used to choose the window: V1b (spatial halves, never run in v1) and V3 (stars fainter than the
  window, never in any locus).

## Screen (ASSUMPTION thresholds)
- Locus stars: the run's Gaia-split star-stratum sources matched to DJA v7.4 phot at 0.2″ (as v1), ref_mag in
  [20, 24]. No fallback window; fewer than 10 → undetermined.
- `r_psf` = median `flux_radius`; `s` = 1.4826 × MAD; point-like ⇔ `flux_radius ≤ r_psf + 3 s` (as v1).

## Validation, frozen before the run (all gates combined with "and")
- V1a parity held-out: cut from even DJA `id`s, recovery on odd, and the reverse; pass if ≥ 0.90 both ways.
- V1b spatial held-out: split locus stars at their median DJA `ra` (east half: `ra` > median; west: the rest);
  cut from one half, recovery on the other, both ways; pass if ≥ 0.90 both ways. Tests PSF variation across the
  mosaic.
- V2 lensed images kept: Mahler+2022 ICLv2 images matched to DJA at 0.2″; pass if ≥ 10 match and ≤ 10 % are
  point-like (as v1).
- V3 faint extrapolation (gate only if ≥ 10 such stars, otherwise descriptive and not part of the verdict):
  matched star-stratum sources with ref_mag in (24, 27]; pass if ≥ 0.80 are point-like (ASSUMPTION: lower bar, as
  the radius is noisier there).

## Decision rule (function `decide_v2`, unit-tested before the run)
- Any of V1a, V1b, V2 fails → **not validated**; record numbers; both rows stay for normal vetting.
- V1a, V1b, V2 pass and V3 passes → **validated to ref_mag 27**.
- V1a, V1b, V2 pass and V3 fails or is descriptive → **validated to ref_mag 24 only**.
- Too few locus stars → undetermined.
- Targets: `f200w_1816` (rank 14) and `f200w_384` (rank 19) have ref_mag 25.897 and 26.305 (v1 `summary.json`), so
  they get a label only under "validated to ref_mag 27": point-like → labelled point-like (not a lens-residual
  candidate); otherwise kept for vetting. Under any other verdict they stay for normal vetting.
- Descriptive only: point-like fraction of galaxy-stratum rows per ref_mag bin and of the top 20.

## What changes the plan
- Validated → propose the label as a pre-vetting screen next to D-019–D-021 (adoption is a separate PR); note that
  validation rests on one field and on V1b/V3 as the only independent star tests.
- V1b fails → PSF size varies across the mosaic; a per-position (or per-detector) locus would be needed.
- Validated to 24 only → the two targets need another route (e.g. MIRI or colour), not a size label.
