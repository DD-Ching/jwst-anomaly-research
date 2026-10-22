# MACS0417: extra-image test with CANUCS DR1 reference fluxes (result)

Plan: `PLAN.md` (frozen and committed before the run, c6a48c9). Numbers: `summary.json`, `control_macs0417.ecsv`,
`predicted_macs0417.ecsv`; post hoc checks `posthoc.json`. Scripts: `scripts/m0417_canucsref.py` (runs
`m0417_arcref.run` with the frozen rules and the CANUCS DR1 catalogue; mapping tested in
`tests/test_m0417_canucsref.py`), `scripts/m0417_canucsref_posthoc.py`.

## Frozen outcome: gate C2 passes; two family-10 predictions flagged `absent`
| Step | Result | Rule |
|---|---|---|
| Gate C2, all controls | 37 / 40 testable recovered or confused (0.925) | n ≥ 8 and ≥ 0.8 |
| Gate C2, arc-referenced controls | 10 / 10 (1.0) | n ≥ 3 and ≥ 0.8 |
| Control classes (54 units) | 33 recovered, 4 confused, 3 absent, 6 inconsistent_reference, 8 no_reference | — |
| Predicted extra images (6) | family 10: 1 `confused`, 2 `absent`; family 8 `inconsistent_reference`; 9 and 17 `no_reference` | flag = `absent` with C2 passing |

- Depth (S/N 4–6 median of this catalogue) 27.58 AB; pixel scale 0.0400″; 14,036 rows with a finite F200W flux.
- The reference catalogue changes the gate from 10 / 27 (JWST pipeline segments, `../macs0417_arcref`) to 37 / 40:
  the CANUCS DR1 sources (BCG-subtracted detection, deblended) scale with the model's |μ| on this field.
- Flags (frozen reading `lensing_residual_input`): family 10 `missing` at (64.394629, −11.906917), |μ| 4.85,
  predicted 131 nJy; family 10 `candidate` at (64.384714, −11.915405), |μ| 2.85, predicted 77 nJy. Reference for
  both: 10.2 (CANUCS 1116189).

## Post hoc checks of the flags (not frozen rules; `posthoc.json`)
Looked at: contact sheet of 10.1, 10.2 and the three family-10 predictions on the F200W science and BCG-subtracted
images (not stored).
- **`missing` |μ| 4.85 (on the BCG core):** the position lies on the BCG's light and its subtraction residual. The
  local 0.3″-aperture scatter is 179 nJy (σ-clipped over 200 random apertures 1.5–4″ away, centres ≥ 0.5″ from any
  CANUCS source, no NaN pixel; it includes BCG residual structure as well as noise), so the predicted total (131 nJy)
  is 0.73× the local 1σ; the aperture at the position holds +1.4σ above the local median. The frozen `absent` uses the
  field-median depth; at the local noise the image is not detectable, so this flag is untestable here (processing
  explanation: BCG light and residual). Its partner `missing` image (|μ| 6.1, 2.31″ away) is `confused`: its nearest CANUCS source (0.68″)
  is 204× brighter than predicted.
- **`candidate` |μ| 2.85:** no flux at the position (0.3″ aperture −0.24σ; the predicted total is 5.7× the local
  13.4 nJy aperture scatter, same method). CANUCS 1115938 lies 1.49″ away with F200W 1.07× the predicted flux and photo-z 1.81
  (z160–z840 1.72–1.88), matching 10.1 (z_ML 1.80) and 10.2 (1.92). Family 10 has two catalogued images, no z_spec,
  and the CANUCS map fit puts it at z 2.15 (`../macs0417_images`), so the predicted position of a further image moves
  with z. Reading: 1115938 is a plausible counter-image ≈ 1.5″ from the model position, outside the frozen 1.0″
  search radius (known-model explanation: source redshift and model position error). Not tested further.
- The target apertures are measured where they lie; the random apertures avoid catalogue sources (`posthoc.json`
  "note"), so the local median and scatter describe source-free sky near the position.
- Neither flag is a candidate; `/vet-candidate` was not run (each has an ordinary reading above). Both readings are
  post hoc.

## Limitations and deviations
- One model (CANUCS v1 maps via `../macs0417_frame`); no injections; no limit. The 1.0″ search radius is smaller than
  the 1.49″ offset found post hoc.
- CANUCS positions are assumed to be on the JWST pipeline frame (plan ASSUMPTION; not tested). The PA convention is an
  ASSUMPTION (plan). Counterpart assignment (rule 1) uses the ellipse orientation, so a post hoc sensitivity run
  (`pa_sensitivity.json`, `--pa-sensitivity`) repeated the run with the PA read as East of North, clockwise from +x,
  and clockwise from +x plus 90°: gate C2 (37 / 40, arc 10 / 10), the target classes and both flags are unchanged in
  all four conventions.
- The 3 `absent` controls (1a.1, 8.1, 15.2) were not inspected.

## Next
- MACS0417 extra images: a z-scan of family 10 in the CANUCS maps (does the third image land on 1115938 at
  z ≈ 1.8?), or the CANUCS v2 image lists. Family 8 stays `inconsistent_reference`, 9 and 17 above the |μ| cap.
- The CANUCS DR1 catalogue as the reference for the same test on MACS1423 and the other CANUCS fields (gate first).
