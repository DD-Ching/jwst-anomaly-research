# Abell 370: the MARS held-out images under BUFFALO model A (plan, frozen)

Frozen 2026-10-24 (commit time of this file), before any BUFFALO residual of these images was computed. Follows
D-034 addendum ap (`../mars_hff_heldout/RESULT.md`): MARS fails gate H on Abell 370 (40 held-out CANUCS images,
rms 14.81″, 90 % within 3″; 3b.1, 5.3, 7.6 and 24e.2 over 3″). Script `scripts/a370_mars_misses.py`; test
`tests/test_a370_mars_misses.py` (before the run).

## Inputs
- model_prediction: BUFFALO model A Lenstool maps of Abell 370 (Niemiec et al. 2023; `abell370-buffalo-a-maps`,
  pinned by sha256; 0.12″ pixels; gate K2 passed 0.910″ on the BUFFALO gold list, `../abell370_buffalo/`).
- observed: the same CANUCS list and MARS catalogue as addendum ap. Before this plan only the nearest BUFFALO gold
  row of each of the four MARS misses was looked at (not a residual): 3b.1 → 3.2 at 4.37″, 5.3 → 5.3 at 1.07″,
  7.6 → 10.2 at 1.77″, 24e.2 → 24.2 at 1.46″; none lies within 1″, so none is a BUFFALO constraint at that radius.

## Question, ordinary explanations, what changes the plan
- Question: are the MARS misses MARS-specific (coarse 1.515″ map, free-form model away from its constraints), or
  does an independent 0.12″ parametric model also miss them?
- "BUFFALO only" for an image: the MARS miss has an ordinary model-side reading (MARS-specific); not pursued.
- "neither": the image is listed for a CANUCS family / z check before anything else (catalogue explanations first);
  still not a flag in this stage.

## Frozen rules (functions with a unit test before the run)
1. Same seen / held-out split and anchors as addendum ap (MARS catalogue, offset to MARS frame for the split), but
   CANUCS positions unshifted when solving with BUFFALO (ASSUMPTION: frame offset 0 to BUFFALO, as K2 used for
   the gold list). Grid step 0.25″, mean-source anchors (rule 3 of addendum ap).
2. `classify`: per held-out image, "both" / "MARS only" / "BUFFALO only" / "neither" by whether each model's
   residual is ≤ 3.0″ (ASSUMPTION; gate H's match radius); NaN or no prediction = not reproduced.
3. Reported, not a gate: `gate_h` with BUFFALO on all 40 held-out images.

## Stop rules
- No rule change after the first BUFFALO residual; departures labelled deviation or post hoc.
- Data: MARS zip (2.7 MB) and BUFFALO maps (3 × 36 MB, under the 200 MB rule) are deleted after the run; only
  `summary.json` and `residuals.ecsv` are stored.
