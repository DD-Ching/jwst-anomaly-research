# MACS J1931 stage 4 / 4b: MUSE Lyα references; gate L fails on count (3 < 4); 0 of 10 targets testable (D-034 addendum ah)

Plans: `PLAN.md` (frozen in dea0952, before any spectrum was read) and the deviation `PLAN_4b.md` (frozen in d2da6f9,
after the stage 4 run and a post hoc look; **not blind**). Code: `scripts/m1931_muse.py`, tests
`tests/test_m1931_muse.py`. Cube: ESO MUSE `ADP.2017-03-27T15:26:07.707` (095.A-0525(A), 8,127 s), read as four
full-field SODA `BAND` cutouts of 105,131,520 B each, deleted after use (URLs and sha256 in each `summary.json`).
Image list observed; |μ| model_prediction (VENUS AstroLensPy v2.1, K-WLS as stage 3d).

| Step | Outcome (stored in) |
|---|---|
| Stage 4 (as frozen) | rule 6 fails for every system: best reference S/N 1.31–1.83, system 101 none; **no reading** (`summary.json`) |
| Post hoc frame check | MUSE − Gaia DR3 on 6 stars with proper motion (epoch 2015.5): median ≈ (+1.30″, −0.88″) in (RA cos δ, Dec), larger than the 0.8″ aperture (`posthoc_gaia_offsets.json`) |
| 4b rule 0 (registration) | applied in all four bands: n = 6, offset (+1.30 to +1.35″, −0.87 to −0.89″), per-axis MAD ≤ 0.062″ (`stage4b/summary.json`) |
| 4b line windows | peaks at 6985.0, 7390.0, 7455.0, 7712.5 Å for systems 5, 6, 100 / 101, 7 (nominal 6984.0, 7390.1, 7455.7, 7706.1) |
| 4b references | 5.2 (S/N 15.4), 6.1 (11.2), 100.4 (7.0); system 7: 7.1 at 2.30 < 5; system 101: all images `blended` with system 100 → `no_line_reference` |
| **4b gate L** | **fails**: 3 testable controls (rule ≥ 4), all 3 `recovered` (5.1 ratio 0.78, 6.2 0.78, 100.3 0.62); 6.3 `undetectable_local` (3.91 σ), 100.1 / 100.2 `blended` (`stage4b/controls.ecsv`) |
| 4b targets | 3 `off_cube` (systems 5, 6, 100: the μ ≈ 2.2 images outside the 1′ field), 1 `undetectable_local` (system 6 `no_flux_ref`, 1.36 σ), 6 `no_line_reference` (systems 7, 101); **0 testable, 0 flags** (`stage4b/targets.ecsv`) |

Reading (frozen rules): no flag and **not a null**: gate L fails on its minimum count, and no predicted extra image
could be tested. The three testable controls are all recovered within the 1/3–3 range (ratios 0.62–0.78; reported,
not a gate pass).

## Looked at
Narrow-band contact sheet of all 13 observed images and 10 targets (6″ stamps, registered frame; 50 kB, not stored):
centred line blobs at 5.1, 5.2, 6.1, 6.2, 6.3; faint emission at the 100.x positions; system 7 stamps noisy with no
clear source; 101 stamps overlap system 100 light; no visible source at any in-cube target.

## Failed approach
Stage 4 as frozen had no frame registration. The 2016 MUSE pipeline cube is offset by ≈ 1.6″ from the Gaia-tied
image list, so every reference aperture missed its image. Any MUSE reference stage needs a registration rule (as
PLAN_4b rule 0) before measurement.

## Limitations
- PLAN_4b was written after a post hoc 2″ line search around the system 6 and 100 / 101 images; its registration
  uses only Gaia stars, but the stage is not blind.
- One MUSE pointing (1′): the three μ ≈ 2.2 extra images of systems 5, 6 and 100 lie outside it.
- System 7's Lyα is at 2.3 σ in this 8.1 ks cube, and system 101 cannot be separated from system 100 at MUSE
  resolution (101.1 / 101.2 lie ≤ 0.4″ from 100.1 / 100.2).

## Next
The MUSE route cannot test the MACS J1931 extra images in this cube: its gate has too few controls and the
testable-in-principle targets lack references or coverage. A deeper or wider IFU / NIRSpec line map would be needed;
otherwise move to another field with public JWST-era maps.
