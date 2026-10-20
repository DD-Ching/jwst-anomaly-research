# W3 dimming: `forced` + footprint/depth test on the 20 MACS0416 / Abell 2744 `on_sky` persistence targets

Run 2026-10-20 03:10–03:16 UTC under the frozen `PLAN.md` (one amendment before any forced flux: uid check by
detection pattern, the verdict tables carry no positions). Screens rebuilt with the production chain on `main`
(weight-map coverage term on; catalogues from S3 by sha256, Gaia from the ESA TAP), `forced` unchanged with 300
controls per field. Numbers: `<field>_summary.json` (rule outcome, depth), `<field>_targets.ecsv` (forced fluxes and
flags). Provenance: derived from observed pipeline products.

## Target match (derived)
- Without the coverage term (`--no-weight-coverage`) all 20 uids are flags vetoed by `single_epoch` alone and
  F200W-catalogued only in their detection epoch (MACS0416 e2 `jw01208-o004`, Abell 2744 e1 `jw02561-o002`), as in
  the persistence run.
- With the term on, 4 are no longer measurable flags: macs0416-d00193 (its detection-epoch pair is a weight edge)
  and abell2744-d04554, d05952, d07198. They are explained by the committed term (D-052 addendum h) and were not
  measured or looked at.

## Calibration (derived)
| field | `calibrated` | controls flagged | F200W noise scale | median control F200W 1σ (µJy) per epoch |
|---|---|---|---|---|
| MACS0416 | true | 0 / 300 | 3.26 | e0 0.0165, e1 0.0137, **e2 0.0121**, e3 0.0142 |
| Abell 2744 | true | 0 / 300 | 4.12 | e0 0.0163, **e1 0.0123**, e2 0.0123 |

Both fields meet the frozen requirement (≤ 1 %), so both count as calibrated.

## Targets (derived; frozen rule `summarise.explanation`, unit-tested before the run)
| outcome | MACS0416 | Abell 2744 |
|---|---|---|
| survivor | 0 | 0 |
| `catalogue_depth_or_deblend` | 11 | 5 |
| `footprint_or_image_quality` | 0 | 0 |
| `weight_edge` (not measured) | 1 | 3 |

- No `forced_vanish`, no `forced_confirmed`, no chromatic veto among the 16 measured targets.
- Forced F200W S/N in the detection epoch 4.4–11.8 (MACS0416) and 5.1–9.5 (Abell 2744); in the other covering
  epochs 2.4–8.4 and 3.8–5.8.
- **Post hoc (not part of the frozen rule):** inverse-variance mean flux ratio, other epoch / detection epoch:
  MACS0416 0.997 ± 0.044 (28 epoch pairs, per target 0.68–1.21), Abell 2744 0.94 ± 0.11 (5 pairs, 0.86–1.09).
  Constant flux, no dimming at the ~5–10 % level per field.
- **Sheets looked at** (`macs0416_sheet.png`, `abell2744_sheet.png`): every one of the 16 is visibly present at the
  same position in every covering F200W epoch and in F444W. abell2744-d03574 e1 has a diagonal streak a few pixels
  away; macs0416-d02264 sits near a bright galaxy halo (`cut_low_weight`); abell2744-d04698 has `cut_spike`. None
  changes the reading. Blank panels are epochs without coverage (MACS0416 e4 has F444W only).

## Reading (narrowest statement)
- None of the 20 MACS0416 / Abell 2744 D-039 `on_sky` single-epoch vanish flags survives. 16 are faint
  (0.019–0.047 µJy) real sources at constant flux that the other epochs' pipeline catalogues missed; 4 sit at
  weight-map edges that the committed coverage term already removes. Consistent with H0 (catalogue depth or
  deblending, footprint), not with a vanish.
- With the 12 NEXUS targets (`nexus_o014`), all 32 `on_sky` persistence targets are now explained, so the D-039
  swap no longer leaves unexplained unvetoed flags in the three fields.

## Next
1. Fold the swap into `inject` with the per-(source, epoch) efficiency by exposure count (0.412 / 0.961, pooled
   0.374 as sensitivity) and re-run `combine` (`expcount_eff/RESULT.md` Next 1).
2. The 11 held-out NEXUS `undetected` control stamps; NEXUS e2 controls (`expcount_eff/RESULT.md` Next 2–3).
