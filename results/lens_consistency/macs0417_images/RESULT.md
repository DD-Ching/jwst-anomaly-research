# MACS0417 and MACS1423: published multiple images against the pinned CANUCS v1 maps (result)

Plan: `PLAN.md` (frozen 03:20Z; a loading check on the MACS0417 spectroscopic images came first and is stated there).
Numbers: `summary.json`, `readings_*.ecsv`, `imageplane_*.ecsv`. Script: `scripts/m0417_images.py`.

## Frozen outcome: tool not validated (gate K fails for every model); no reading
| Field | Model | Spectroscopic images | rms (″) | matched ≤ 3″ | K |
|---|---|---|---|---|---|
| MACS0417 | `macs0417-canucs-maps` | 25 | 1.12 | 0.96 | fails (rms) |
| MACS0417 | `macs0417-relics-lenstool` | 25 | 9.29 | 0.12 | fails |
| MACS1423 | `macs1423-canucs-maps` | 12 | 1.59 | 0.92 | fails (rms) |

So the frozen rule produces no `lensing_residual_input` row and no flag (`summary.json` `flags` empty).

## Post hoc (deviations from PLAN.md, labelled; not frozen readings)
- **The two CANUCS failures each come from one catalogued radial / central image next to the BCG.**
  - MACS0417: without 1b.4 the other 24 spectroscopic images give rms 0.38″, all matched (K would pass). CANUCS
    predicts three images for knot 1b and no central one. Jauzac et al. 2019 found 1b.4 only through its [OII] line
    in MUSE and show that it controls the dark-matter-to-BCG offset in their model (their "Variation 1" without it
    misaligns). A published, model-core-dependent constraint that CANUCS v1 does not reproduce; ordinary.
  - MACS1423: leave-one-out over system 1 (1.1–1.5): dropping the central image 1.5 gives rms 0.37″ for 1.1–1.4
    (dropping any other image: 1.92–3.10″). For the source of 1.1–1.4, CANUCS predicts a demagnified central image
    (μ 0.43) 2.56″ from 1.5 (`held_out_offset`). System 2, the other component of the same galaxy, reproduces at
    ≤ 0.47″. A central image displaced by a few arcsec is the usual sensitivity to the BCG core and truncation;
    ordinary, not a candidate.
- **Non-spectroscopic systems at a fitted z (no frozen reading because K failed):** CANUCS reproduces all 31 MACS0417
  images (max 0.82″) and all 5 of MACS1423 system 4 (max 1.05″); no catalogued image is unmatched.
- **RELICS v1 maps load correctly:** κ and |μ(z = 2)| from its deflection maps match the published maps (median
  relative difference 2.6 × 10⁻⁵ and 1.4 × 10⁻⁵, `map_check`). Its failure is the model against this (MUSE-era)
  image set, not loading. It reproduces the non-spectroscopic systems only at fitted z ≈ 0.98–1.35 (`z_fit`),
  which is not physical for most of them; it is not used.
- `extra_predicted` (CANUCS): MACS0417 families 8, 9, 10, 11, 17 (2 catalogued, 3–5 predicted); MACS1423 system 2
  (4 vs 5, its central image). Forced photometry at these positions is the next possible step; not done.

## Limitations
- CANUCS v1 image catalogues are not public, so which of these images constrained the CANUCS models is unknown
  (ASSUMPTION that the spectroscopic ones did). The test is model-versus-catalogue consistency, not an independent
  prediction.
- One model per cluster after RELICS fails; positions are the published HST-frame coordinates (no frame offset
  measured; CANUCS is JWST-frame).
- No images were looked at in this cycle: no candidate is claimed, and the two central images are published
  identifications.
