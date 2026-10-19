# AT 2024jfz (SMACS 0723): counter-images and delays in Mahler+22 ICLv2: result (2026-10-19)

Plan: `PLAN.md` (frozen in `d012032` before the run). Predicted positions, magnifications and
delays are `model_prediction` (`result.json`, `posthoc.json`); separations `derived` (`sanity.json`); epochs `observed`
(MAST metadata, `mast_epochs.json`). Thresholds are ASSUMPTIONs (PLAN.md). Best model only: no delay spread.

## Validation
- Delay integral: unit tests on a singular isothermal sphere (Δφ = 2 θ_E β) and an analytic curl-free field pass
  (`tests/test_smacs_at2024jfz.py`).
- Path independence on the real model (straight vs L-shaped path): max relative difference 1.5 × 10⁻¹³ over all
  counter-images (`path_check_max_rel_diff`; gate ≤ 1 %): **passes** for every z_s. This checks the quadrature only
  (any sum of dPIE deflections is curl-free); the delay factor is unit-tested against the independent MACS0417 #308
  value (110.66 d arcsec⁻² at z_l 0.443, z_s 0.8718) and the Fermat sign against the SIS case.
- Geometry check (`sanity.json`): at system 6's model redshift (z 1.712766) the counter-images fall 1.88″ from 6.2 and
  2.01″ from 6.3, as the TNS position lies 1.75″ from 6.1 (the 1.81″ in PLAN.md is the `arc_diff` h000 detection
  position's distance, `../smacs_arcdiff/RESULT.md`; the TNS position is 0.081″ from h000).

## Result (`result.json`, frozen rule `testable`)
| z_s | Images | Nearest counter-image: sep, μ_k/μ_0, delay | Far counter-image: sep, μ_k/μ_0, delay | Outcome |
|---|---|---|---|---|
| 0.6 | 1 | — | — | not testable |
| 1.0 | 3 | 14.6″, −2.33, +680.1 d | 21.1″, +1.19, +575.8 d | not testable |
| 1.425 | 3 | 6.2″, −0.84, +85.5 d | 30.6″, +0.33, −3,497.4 d | not testable |
| 1.712766 | 3 | 4.5″, −0.89, +31.0 d | 33.8″, +0.21, −5,463.6 d | not testable |
| 2.5 | 3 | 2.0″, −0.95, +2.6 d | 38.4″, +0.08, −8,767.9 d | not testable |
| 4.0 | 3 | 0.087″, −1.00, 0.0 d | 42.2″, +0.003, −11,646.2 d | testable (2024 epoch) |

Positive delay: the counter-image trails the 2024 event. **Frozen outcome (function `testable`):** "not testable" for
z_s 0.6, 1.0, 1.425, 1.712766 and 2.5; **"testable" for z_s 4.0**, in the 2024 discovery epoch itself. There the source
lies on the model's critical curve (μ ≈ 864) and the counter-image is a merging pair 0.087″ away. The frozen rule has
no resolvability condition, so this outcome stands; testing it needs a resolved-pair measurement of the 2024 point
source (not done). For z_s 1.0–2.5 the nearby counter-image's window (`next_window_mjd`, [arrival − 2 d, arrival
+ 30 d]) falls between the 2024 and 2026 epochs, all already past; the far image led by 9.6–24.0 yr for z_s 1.425–2.5
and 31.9 yr at z_s 4.0.

**Deviation (bug fix after review):** PLAN.md wrote the next window as [t_0 + Δt_k − 30 d, t_0 + Δt_k + 2 d], the
mirror of the frozen `testable` window; `next_window_mjd` now uses the `testable` window. The outcomes do not depend on
it. The path gate now withholds the outcome when it fails (it passed everywhere).

## Post hoc (not in the plan; does not change the outcome)
- `mast_epochs.json`: every HST and JWST imaging day within 0.02° (one row per collection, program and day; footprint
  coverage of each predicted position not checked). `posthoc.json`: z_s 1.00–2.50 in steps of 0.05 against those days.
- **NIRISS imaging, program 4368, CLEAR;F150W, MJD 60557.96** would see the nearby counter-image if z_s ≈ 1.35–1.40
  (6.4–6.8″ from the TNS position, μ_k/μ_0 ≈ −0.84, delay 94–114 d).
- HST WFPC2 (program 11103, MJD 54512.3 / 54791.02) passes the rule's timing for the far image if z_s ≈ 1.75–1.80
  (34″, μ_k/μ_0 ≈ 0.19). The flux threshold is the NIRCam one: whether WFPC2 could detect ≈ 490 nJy in its bands is
  not checked, so this route may not be usable.
- No image was examined in this step (predictions only); no candidate can come from it.

## Reading
AT 2024jfz's lensed-host hypothesis is **not tested by the existing JWST NIRCam epochs** under Mahler+22 ICLv2 for any
z_s ≤ 2.5; at z_s ≈ 1.7 (system 6's model z) the counter-image near 6.2 trailed by about a month and was not observed
then. The two post hoc routes above are the only archival timing matches found; both need a footprint and depth check
first, and both depend on a host redshift that is not known here. At z_s ≥ ≈ 4 the model puts the source on the
critical curve, an unresolved-pair case the rule counts as testable.
