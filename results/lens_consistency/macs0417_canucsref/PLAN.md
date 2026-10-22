# MACS0417: extra-image test with CANUCS DR1 reference fluxes (plan, frozen)

Frozen 2026-10-22 (commit time of this file), before any CANUCS DR1 flux at a predicted or catalogued image was
looked at. Only the catalogue's header, column units and its X/Y-to-RA/Dec relation were read (below).
Follows `../macs0417_arcref/RESULT.md`: the JWST pipeline catalogue's segments hold several lensed knots or other
light (gate C2 10 / 27). CANUCS DR1 detects on BCG-subtracted images and deblends, so its sources may separate the
knots. Script: `scripts/m0417_canucsref.py`; it runs `m0417_arcref.run` with every frozen rule unchanged and only
the catalogue replaced (rules 1–8 of `../macs0417_arcref/PLAN.md`, constants in its `summary.json` "rules").

## Hypothesis and ordinary explanations
- Question: with CANUCS DR1 Kron totals as the reference and search catalogue, do the catalogued images of
  MACS0417 recover their model-predicted fluxes (gate C2), and if so, are the 6 CANUCS v1 extra predictions
  (families 8, 9, 10, 17) present?
- Ordinary explanations for an `absent`: wrong model μ or position near critical curves, a low-surface-brightness
  image below detection, deblending or BCG-model residuals, a reference holding unrelated light. Exotic readings are
  not considered here.
- What would change the plan: C2 fails → the CANUCS DR1 catalogue is not a valid reference here either, and
  MACS0417 waits for the CANUCS v2 image lists or an aperture method on BCG-subtracted images. A flag with C2
  passing → contact sheet, then `/vet-candidate`.
- Known before freezing (from `../macs0417_frame`): families 9 and 17 exceed the |μ| ≤ 50 reference cap and stay
  `no_reference`; only families 8 and 10 can become testable.

## Inputs
- observed: CANUCS DR1 MACS0417 cluster-field photometry catalogue
  `hlsp_canucs_jwst-hst_multi_macs0417-clu_multi_v1_photometry-cat.fits.gz` (DOI 10.17909/18nv-np70; 14,272 rows;
  sha256 in `summary.json`). Streamed to the scratch disk, not committed.
- model_prediction: `../macs0417_frame/predicted_macs0417.ecsv` and `control_macs0417.ecsv` (unchanged); the
  Jauzac+19 list shifted by O = (+0.101″, −0.292″) as in `../macs0417_arcref`. CANUCS positions are taken as on
  the same frame as the JWST pipeline (ASSUMPTION; both are JWST 1208 F200W-era astrometry; not tested here).

## Column mapping (frozen; ASSUMPTIONs where marked)
| Rule input | CANUCS DR1 column |
|---|---|
| position | `RA`, `DEC` |
| ellipse semi-axes (1σ) | `A`, `B` (pix) × pixel scale fitted from `X`, `Y` vs `RA`, `DEC` (0.0400″, north up, east left: read before freezing) |
| ellipse orientation | `PA` (deg, −90…90) read as counter-clockwise from +x (SExtractor `THETA_IMAGE` convention; ASSUMPTION), so sky PA East of North = `PA` + 90° |
| flux | `FLUX_KRON_TOTAL_F200W` (nJy → Jy) |
| S/N | `FLUX_KRON_TOTAL_F200W` / `FLUXERR_KRON_TOTAL_F200W` |
| extended | not `FLAG_POINTSRC` |
| ellipticity | 1 − `B` / `A` |
| depth | `lens_consistency.depth_mag` definition (median AB of rows with S/N 4–6) on this catalogue |

Rows with a non-finite F200W flux are dropped before the run. No other row cut (in particular `USE_PHOT`,
`FLAG_BCG`, `FLAG_HIGH_BCGCONTAMI` are not used).

## Frozen rules
Rules 1–8 of `../macs0417_arcref/PLAN.md`, verbatim, with the mapping above. Gate C2: testable controls ≥ 8 and
≥ 0.8 recovered or confused, and the same for arc-referenced controls with n ≥ 3. Flag: a predicted extra image
classed `absent` while C2 passes. The mapping functions are unit-tested in `tests/test_m0417_canucsref.py`.
