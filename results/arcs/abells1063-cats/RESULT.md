# `arcs` on Abell S1063: result (2026-10-16; vetting open)

Plan: `PLAN.md` (committed before the run). Model `abells1063-cats` (CATS v4.1 maps; D-035), GLIMPSE 3293 o001 F200W
pipeline catalogue (sha256 `a750a902…`, manifest `abells1063`), DJA v7.5 eazy zout (sha256 `03bc4c28…`). D-024
defaults (ASSUMPTIONs). Labels: shapes `observed`; classes and tests `derived`; predicted directions `model_prediction`.

## Screen (`arcs.json`)
- 109 selected; strong-shear background: **20 of 24 aligned** (p 6.6e-7 under random axes), **4 anti vs 0.524
  expected from shape noise, P(≥ 4) = 0.0021**. All strong-shear sources: 47 / 77 aligned, 15 anti.
- **Deviation:** the plan's known-case check (catalogued images) did not run: a map model has no parsed image list in
  `arcs` (`multiple_images: null`), as for MACS1149. The background alignment above and the checks below are the
  controls.

## Ordinary-explanation tests on the 4 background anti (1316, 1557, 2004, 2184)
| Test (file) | Control | Outcome for the anti |
|---|---|---|
| F444W axis (`results/arcs_bands/abells1063-cats/summary.json`) | 18 / 18 usable aligned stay aligned | none of the 4 has an F444W catalogue source within 0.3″ (nearest 2.1–4.2″): untestable |
| Magnitude (`faint.json`, **post hoc**) | m < 27: 17 / 17 aligned, 0 anti (0.308 expected) | all 4 at m 27.43–28.18; P = 0.0014 that 4 random rows all fall among the 6 faintest |
| Same split, five earlier fields (`results/arcs/<model>/faint.json`, **post hoc**) | — | m ≥ 27 anti there: 0 / 6, 0 / 2, 0 / 12, 1 / 7, 2 / 17 (Abell 2744, El Gordo, MACS0416, MACS1149, Abell 370). The faint end elsewhere does not show S1063's 4 / 7, so faintness alone does not explain it |
| ICL gradient: axis re-measured after a 1″ median-filter subtraction (`bgsub_shapes.json`) | 18 / 20 aligned stay aligned (median axis change 3.7°) | 1316, 1557, 2184 stay anti (e 0.25–0.36); 2004 becomes mixed (moment ellipticity undefined) |
| Model error: GLAFIC v4, Sharon v4, Keeton v4 κ/γ maps (`second_models.json`) | 20 / 20 aligned in each model | all 4 anti in all three models |
| κ at photo-z, three team models (`second_models.ecsv`; the CATS values were computed in the session, not stored) | — | 2004: κ 0.95–1.00 (at or near κ = 1, where a few-percent κ or z error turns the predicted axis by 90°); 1316 κ 0.80–0.83, g 1.6–1.9; 2184 κ 0.71–0.74, g 1.5; 1557 κ 0.67–0.69, g 0.43–0.45 |

Cutouts (F115W / F200W / F444W, 8″ and 2.4″) were looked at (not stored): all 4 are faint smudges on ICL gradients;
1557 looks like two knots side by side along its axis; 2184 lies ~2″ from the giant tangential arc, whose direction
matches the predicted one. Observed axes lie 10–15° from the radial direction to the model centre (not stored; the
re-measured offsets are in `bgsub_shapes.ecsv`).

## Classification (per item; derived)
- **2004:** explained: κ 0.95–1.00 in the three team models at its photo-z, where a few-percent κ or z error turns the predicted axis by 90°;
  the re-measured axis is mixed.
- **1557:** weak shear (g 0.43–0.45), e 0.25 after ICL removal, apparent two-knot blend: intrinsic shape or blend
  is the leading ordinary explanation; not yet tested by a deblend.
- **1316, 2184: open.** Anti in four independent models (CATS and the three team models), anti after ICL removal, g ≈ 1.5–1.9, m ≈ 27.4–27.5,
  photo-z 2.30 / 1.81 (no z_spec in DJA). **Unvetted; not candidates.** Untested ordinary explanations: photo-z
  failure (a member or foreground source is unlensed), a fragment of nearby extended light (2184 next to the giant
  arc), noise in the shapes of m ≈ 27.5 sources (the shape-noise null models intrinsic shapes only).

## Next (for the next run)
`/vet-candidate` on 1316 and 2184: MUSE redshifts (Karman+2015 / Caminha+2016 RXC J2248 catalogues), a deblend
against the giant arc for 2184, and an injection of round faint sources at their positions (shape noise at m 27.5
on the local ICL) to calibrate how often noise alone gives an anti axis.
