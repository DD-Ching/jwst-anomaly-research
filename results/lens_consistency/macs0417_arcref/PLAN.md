# MACS0417: extra-image test with an arc-aware (catalogue-segment) reference flux (plan, frozen)

Frozen 2026-10-22 00:12Z (commit time of this file), before any segment flux at a predicted or catalogued image was looked at. Follows the
"Next (optional)" of `../macs0417_frame/RESULT.md`: there the 6 predicted extra images (families 8, 9, 10, 17)
were all `no_reference`, because their catalogued siblings are thin arcs that fail the D-034 compact rule.
Script: `scripts/m0417_arcref.py`; frozen rules as functions, unit-tested in `tests/test_m0417_arcref.py`.

## Hypothesis and ordinary explanations
- Question: are the CANUCS v1 extra images predicted for families 8, 9, 10, 17 present at the flux the model
  predicts, when the reference is the total (isophotal segment) flux of the catalogued sibling instead of a 0.2″
  aperture?
- Ordinary explanations for an `absent`: wrong model μ or position near critical curves, a low-surface-brightness
  image below the catalogue's detection, deblending (the image merged into a neighbour's segment), the reference
  segment including unrelated light. Exotic readings are not considered at this stage.
- What would change the plan: gate C2 fails → segment totals do not scale with the model's |μ| for MACS0417 and the
  route is dropped (next: CANUCS v2). A flag with C2 passing → contact sheet, then `/vet-candidate`.
- Known before freezing (from `../macs0417_frame/control_macs0417.ecsv`): families 9 (|μ| 77, 76) and 17 (64, 64)
  exceed the D-034 reference cap |μ| ≤ 50, so under this plan they stay `no_reference`; only families 8 and 10 can
  become testable.

## Inputs (unchanged from `../macs0417_frame`)
- observed: JWST 1208 o006_t003 F200W pipeline `_cat.ecsv` (sha256 in summary.json); the Jauzac+19 list shifted by
  O = (+0.101″, −0.292″) (`m0417_frame.shifted_loader`).
- model_prediction: `../macs0417_frame/predicted_macs0417.ecsv` (targets, their |μ|, search radii) and
  `control_macs0417.ecsv` (|μ| of each catalogued image's matched prediction). The lens model is not re-run.

## Frozen rules (ASSUMPTIONs; per item unless stated)
1. **Counterpart** of a position: the catalogue source whose 3σ ellipse (semimajor/semiminor_sigma, sky_orientation;
   pixel scale from the catalogue's own pixel and sky centroids) contains it, the smallest normalised radius winning;
   none → no counterpart.
2. **Merged images:** catalogued images of one family that share a counterpart form one unit: F = its
   `isophotal_flux`, μ = Σ|μ|.
3. **Usable reference unit:** counterpart found, catalogue S/N > 5, every member |μ| ≤ 50 (D-034 cap).
   Luminosity L = F / Σ|μ|. Several usable units with max L / min L > 3 → `inconsistent_reference`; else the
   reference is the unit with the smallest Σ|μ| (as D-034).
4. **Prediction:** F_pred = L_ref × |μ_target|; m_pred = −2.5 log10(F_pred / 3631 Jy).
5. **Search radius:** as `lens_consistency.system_search_radius` (base 1.0″) on the predicted table.
6. **Class of a target** (sources = catalogue sources with centroid within the radius, excluding the reference's
   counterpart): `undetectable` if m_pred > depth (`lc.depth_mag`, the S/N 4–6 median); else `recovered` if any
   source has F / F_pred in [1/3, 3]; else `confused` if any has F / F_pred > 3; else `absent` if
   m_pred ≤ depth − 0.753 (10σ); else `ambiguous`.
7. **Gate C2 (positive control):** every catalogued image (unit) is a target, with its reference chosen among the
   other units of its family (leave-one-out; a unit sharing the target's counterpart is excluded). Passes when the
   testable controls (`recovered`, `confused`, `absent`, `ambiguous`) number ≥ 8 and ≥ 0.8 of them are
   `recovered` or `confused`, **and** the same holds for controls whose reference is an arc (a unit whose
   counterpart has `is_extended` True and ellipticity ≥ 0.5) with n ≥ 3. If arc-referenced controls number < 3,
   C2 fails (not validated for arcs).
8. **Flag:** a predicted extra image (class not in `observed`, `demagnified`, `outside`) classed `absent` while C2
   passes.
