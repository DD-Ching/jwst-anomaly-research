# MACS J1931 stage 3: extra-image test with the own deblended reference (result: stopped at frame check F)

Plan: `../macs1931_counter/PLAN.md`, "Stage 3: frozen rules" (052ba6a), committed before the run; driver and the
footprint-rule tests committed before the run as well. Script `scripts/m1931_stage3.py`. Numbers: `summary.json`
(frozen rules), `posthoc.json` (post hoc, not rules), `frame/` (rule 3 tables).

## Frozen outcomes
| Rule | Result |
|---|---|
| 1 inputs | VENUS v2.1 maps md5-verified against the Zenodo record; 51 used images, 26 spec-z |
| 2 gate K-WLS repeated | rms 0.630″, 26 / 26 within 3″: passes, same outcome as `../macs1931_gatek` |
| 3 frame tables | written (`frame/`); pipeline depth proxy undefined (`pipeline_depth_mag` null) |
| 4 own catalogue | 653 segments → 897 after deblending, all kept; section y 557–3900, x 833–3456; depth 26.76 AB (V2-a) |
| **5 frame check F** | **18 / 51 matched within 0.5″ (rule: ≥ 20); median dx 0.045″, dy 0.004″: fails** |
| 6–9 | not run: the plan stops at a failed F, so no gate C2 and no flag reading |

Reading (frozen rule): **no stage 3 reading for MACS J1931**. No counter-image test was made; this is not a null.

## Post hoc (not rules; `posthoc.json`)
- The matched images agree in frame (medians ≤ 0.045″, median separation 0.107″); the failure is the count.
- Images with any row within 0.5 / 1.0 / 1.5″: own catalogue 18 / 20 / 27 of 51, the JWST pipeline catalogue 18 / 24 / 33.
  Spec-z only (26): own 10 / 11 / 14, pipeline 13 / 16 / 19. Only 15 / 51 (10 / 26 spec-z) positions fall inside
  any own deblended segment.
- Looked at (contact sheet of the 51 positions with the segment outlines, 6″ boxes, not stored): systems 5, 6, 7,
  100 and 101 show little or no F200W emission at the catalogued positions; several arcs (8.2, 9.1, 10.1, 12.2) are
  visible but not segmented at 3σ / 25 px after the 2″ background filter. No sign of a frame shift. Reading by eye,
  not measured: the builder does not detect most catalogued images in this field.

## Limitations and deviations
- None from the frozen plan; the posthoc diagnostic was added after the failure and is labelled as such.
- The rule 3 depth proxy (median mag at S/N 4–6 of the pipeline catalogue) is undefined here, as V1's was; it only
  enters the model-prediction classes of `frame/`, which were not read.
- F's minimum (20) was an ASSUMPTION scaled from the Bullet's 30 / 135; the rule was not changed after the run.
