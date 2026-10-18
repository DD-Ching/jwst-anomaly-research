# Batch 4: minimal discriminating check of the two rule-level real changes (fixed before any control result)

Written 04:22Z after the full-rule classes (`../boxes/review/final_classes.json`) and before any control was picked.
Per the owner instruction, no further search volume is added; this is the cheapest check that separates "the source
changed" from "bright compact sources change by this much between these two mosaics", the method of confirmation
round 1 (`../../confirm_r1/PLAN.md`), unchanged.

Entries: **B4-A370-box1-h000 (39.968857, −1.571911)** (bright point-like source; search pair 1208 o002 → VENUS o053) and
**B4-SJ0723.3-7327-box1-h005 (110.809928, −73.463102)** (bright compact galaxy; ERO 2736 o001 → VENUS o057).

- Controls: `hit_lightcurve_controls.py pick` on each entry's epoch-1 F150W mosaic, `--half-arcsec 16 --flux-range
  <0.5× and 2× the entry's epoch-1 F150W 0.15″ flux> --n 80 --seed 1`; if fewer than 20 controls, the same with
  `--half-arcsec 30` (fallback fixed now). Controls and the entry measured together by `hit_lightcurve.py` (0.15″).
- Decisive pair: the search pair in **F150W** (the band of the larger |S| peak for both entries). Secondary
  (reported, not deciding): the search pair in F444W; for h005 also F444W VENUS o057 → 12774 o002 (31 days later).
- Rule (unchanged): with ≥ 20 controls, control tail fraction < 0.05 ⇒ "change exceeds matched-control scatter"
  (supports a real change); ≥ 0.05 ⇒ "consistent with the matched-control scatter" (not proof of an instrumental
  cause). < 20 controls ⇒ not decided.
- Also recorded (no decision attached): Gaia DR3 within 1″ of h000 (a star would be an ordinary explanation for a
  point source; no Gaia match cannot by itself exclude a star).
- Stop: 05:45Z.
