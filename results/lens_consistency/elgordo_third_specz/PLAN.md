# El Gordo third-image candidates at DJA NIRSpec redshifts — frozen plan

Frozen 2026-10-25 before the DJA table is read for this field. Script `scripts/elgordo_third_specz.py`.

## Hypothesis and ordinary explanations
- Addendum ay / az: the Caminha+2023 model misses the JWST-proposed third images 6c (3.05″), 13c (6.26″) and 20c
  (5.10″) of Diego+2023 (`tab_arcs`); post hoc they are isolated misses. Either the identification is wrong (the
  object is another galaxy) or the model's counter-image position is off by 3–6″.
- H_id: a candidate is the family's image, so a NIRSpec spectrum at it gives the family's MUSE z.
- Ordinary explanations for a miss that survives H_id: Caminha model error (mass near cluster members, the
  model's own rms 0.754″, its third-image predictions sit at m ≈ 29). Wrong identification is also ordinary.
- Plan change: a candidate whose spectrum agrees with the family z becomes a lensing-residual input (the model
  misplaces a confirmed image) for `/vet-candidate`; one that disagrees closes as a misidentification; no spectrum
  leaves it open.

## Inputs
- Images: `results/lens_consistency/elgordo_heldout/diego23_images.ecsv` (observed positions; families 1–23 carry
  Caminha's MUSE z_spec, observed).
- Targets: primary 6c, 13c, 20c; reported 4c, 15c, 18c, 19c (the other post hoc third-image candidates of
  addendum ay, `posthoc.json` `dag_ids`).
- DJA v4.6 NIRSpec emission-line table (the same bytes as addendum ai: sha256 `ee64bb46…c6ddf3`, 99,478 rows),
  streamed, not stored; the script stops unless every row is read.

## Rules (ASSUMPTIONs)
1. Match: DJA rows within 0.5″ of a Diego position; only `grade` 3 rows give a z. No frame offset is applied
   (rule 3 tests the frame).
2. Image z: median of its grade-3 z when they agree within 0.01 (1 + z), else `ambiguous`.
3. Control gate Z: every Diego image with a z_spec family (families 1–23), the targets excluded, that has an
   image z: ≥ 3 such images and ≥ 80 % agree with the family z within 0.01 (1 + z). If Z fails, targets are
   `not decided`.
4. Target reading (only if Z passes): `confirmed` (agrees with the family z within 0.01 (1 + z)), `z_differs`
   (disagrees), `no_spectrum` (no grade-3 row within 0.5″), `ambiguous` (rule 2).
5. A `z_differs` reading is a misidentification only after the matched source is checked to be the Diego
   position's object (separation, contact sheet if imaging is at hand); that check is post hoc and labelled so.
