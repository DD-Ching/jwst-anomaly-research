# VENUS batch 3 with fixed `arc_diff`: result

Plan: `PLAN.md` (rules, seeds, caps and stop times before selection); `boxes.json` committed before any difference image.
Owner instruction 2026-10-17 (after #309); budget 16:06Z → 20:06Z. Every number is `derived` unless labelled; injected
sources are simulations. An ordinary transient is not wormhole evidence; zero entries is not an exclusion.

## New coverage (`boxes.json`)
- **20 new, non-overlapping boxes** of 24″ (3.2 arcmin²) in 10 cluster fields (two each): MS0451, PLCK G165, MACS0416,
  MACS1149, Abell S1063, SMACS 0723, Abell 370, Bullet, MACS1423, MACS0417. Each pairs an earlier public NIRCam epoch
  with VENUS (6882) in F150W + F444W; footprints and ≥ 0.95 finite pixels checked (metadata only).
- **New ranking (owner instruction):** among boxes of similar depth (all pairs share the VENUS 1,245 s F150W depth),
  those with more extra public F150W visits came first: MS0451 (4 extra visits), PLCK G165 (4 and 3), MACS0416 (3), MACS1149
  (2), AS1063 and SMACS (1); the remaining 8 boxes have none (`extra_f150w_epochs` in `boxes.json`).
- **Not reached:** RXCJ0018, RXJ2129, SPT0615 and Sunrise (ranked after the 20th box). El Gordo had no free box: all
  nine grid positions overlap searched boxes. Extra epochs were
  counted only among the programmes in `venus_pairs.json` and other VENUS visits.
- Pipeline versions: jwst 3.0.0 except epoch 1 of MS0451, SMACS and Bullet (2.0.1).

## Calibration (per box, unchanged rule; `boxes/batch.json`, `totals.json`)
- **Completed: 20 of 20; timed out: 0; not started: 0** (caps 900 s screen / 1,200 s injection never reached; slowest
  task 661 s).
- **Passed: 10 boxes (1.6 arcmin²):** PLCK G165 1, MACS0416 1, MACS1149 1, SMACS 1, A370 1, Bullet ×2, MACS1423 1,
  MACS0417 ×2.
- **Calibration failed: 10 boxes** (not counted as searched):
  - too many unchanged-control false hits: AS1063 1 (15/16), AS1063 2 (17/19);
  - too few vanish recoveries (26 + 27 AB): MS0451 1 (5/10), MS0451 2 (7/14), MACS0416 2 (3/8), PLCK G165 2 (6/10),
    MACS1149 2 (5/12), SMACS 2 (9/14), A370 2 (3/7), MACS1423 2 (7/14).
- **Pooled over the 10 passed boxes** (review list with sign; 68 % Wilson intervals):

  | Change | 26.0 AB | 27.0 AB | 27.5 AB |
  |---|---|---|---|
  | vanish (no flux left) | 60/64 = 0.94 [0.90, 0.96] | 51/63 = 0.81 [0.76, 0.85] | 50/54 = 0.93 [0.88, 0.95] |
  | partial fade (50 % drop) | 59/64 = 0.92 [0.88, 0.95] | 38/54 = 0.70 [0.64, 0.76] | 31/64 = 0.48 [0.42, 0.55] |
  | appear (brightening) | 48/55 = 0.87 [0.82, 0.91] | 48/64 = 0.75 [0.69, 0.80] | 51/64 = 0.80 [0.74, 0.84] |
  | unchanged controls, any hit | 3/183 = 0.02 [0.01, 0.03] | | |

- Valid only for point sources injected on galaxy and arc pixels, not within 0.8″ of bright cores, not extended
  changes, not fainter than 27.5 AB, and only for the two search epochs. The control fraction is not a real-sky
  false-alarm rate. Controls per passed box: 11–20.

## Screening and review (`boxes/review_queue.json`, `boxes/review/`)
- **136 hits → 22 not near a bright core → 9 review-list entries** (morphology pass): 3 in passed boxes, 6 in
  calibration-failed boxes. **9 of 9 reviewed** (capacity 10); pending 0.
- Each review: sheet looked at; forced photometry in every covering public mosaic (367 measurements); coded class;
  the #296 post-review checks and the #304 full-rule function (`final_classes.json`, script stored beside it).

| Rank | Entry | Box | z F150W / F444W (search epochs) | Final class |
|---|---|---|---|---|
| 1 | B3-SJ0723.3-7327-box1-h007 (110.830283, −73.460338) | passed | −1.0 / −1.1 | no significant change |
| 2 | B3-MJ0416.1-2403-box1-h004 (64.048155, −24.074940) | passed | +0.6 / +2.7 | undetermined (≥ 3σ of mixed sign between other mosaics) |
| 3 | B3-MJ0416.1-2403-box1-h005 (64.043175, −24.073775) | passed | +0.1 / −0.5 | no significant change |
| 4 | B3-MJ0416.1-2403-box2-h001 (64.031122, −24.078916) | failed | +1.4 / +1.2 | undetermined (large mixed-sign differences between other mosaics, next to a bright galaxy) |
| 5 | B3-SJ0723.3-7327-box2-h011 (110.842938, −73.454572) | failed | +0.8 / −3.7 | undetermined (opposite signs; dipole along an arc) |
| 6 | B3-A370-box2-h003 (39.966687, −1.583019) | failed | −1.7 / −5.4 | real change by rule; see vetting |
| 7 | B3-PLCKG165+67.0-box2-h002 (171.811463, +42.477850) | failed | −3.0 / −2.0 | undetermined (mixed-sign differences between other mosaics) |
| 8 | B3-AS1063-box1-h002 (342.189344, −44.530953) | failed | −1.4 / −1.3 | no significant change |
| 9 | B3-A370-box2-h010 (39.970466, −1.585253) | failed | −0.1 / −3.0 | real change by rule; see vetting |

- Names carry batch, field, box and hit: batch-1 PLCK G165 h002 (171.805899, +42.481802) and batch-2 AS1063 h002
  (342.179822, −44.533984) are different objects.

## The two rule-level real changes (both in A370 box 2, a calibration-failed box)
Vetting record: `docs/candidates/batch_venus3_A370_box2.md`.
- **B3-A370-box2-h003:** a compact galaxy ([LRB2019] 6401, z_spec 0.3604, 0.08″; not a star: no Gaia source, galaxy
  spectrum) present in all six programmes (1208, 3538, 2883, 5324, 5890,
  VENUS). The rule's change is F444W 255.4 ± 3.3 → 219.4 ± 5.9 nJy (1208 2022-12 → VENUS 2026-01). **The independent
  epoch 18 days after VENUS (5324, 2026-01-27) gives 247.4 ± 7.5 nJy**, consistent with 2022; between programmes the
  differences have mixed signs (F200W −30σ, F210M +8σ). Reading: not confirmed by the independent epoch; consistent
  with programme-to-programme photometric scatter of a compact bright galaxy; variability not excluded.
- **B3-A370-box2-h010:** a 0.15″ aperture on the steep edge of a bright extended structure; **no point source is visible
  in any of the stamps**. In F200W, F277W, F356W and F444W the fluxes are 29.9–51.5 nJy in 1208 (2022-12) and 5324 (2024-09), then
  7.7–16.2 nJy in 5324 (2024-12, 2026-01) and VENUS; F150W is noisy (1.0–23.1 nJy). By eye on the stamps (not stored),
  the first 5324 visit has a different orientation from its later ones. The position is on the
  Dragon arc (lensed galaxy, z 0.725), where Palencia+2026 list 33 microlensing-event detections within 1″ (nearest
  0.54″). Reading: two open ordinary explanations, a caustic-microlensing event or steep-gradient photometry.
- Neither is an unexplained change or an exotic candidate. TNS: one request per position, response not recorded, so
  the TNS result is unknown (`crossmatch.json`).

## Deviations and corrections (labelled)
- **Tools changed before selection (tested):** `batch_pairs.py --prefer-extra-epochs` (owner-requested ranking);
  `batch_run.py` per-task timeout and stop time; `_run` made portable without `os.wait4` after Windows CI failed on the
  new timeout test (fixed by code, not by skipping).
- **Not run:** the matched arc-pixel control check for h010 (named as the next check, not run, per the owner's
  instruction not to open new work).

## Blind spots
- The 10 calibration-failed boxes (AS1063 for false hits; 8 others for vanish recovery), the 16 held boxes of batches
  1–2, cores within 0.8″ of bright galaxies, changes fainter than 27.5 AB, extended changes, and the five fields not
  reached. The time caps removed no box this time; the per-stamp fallback (stamp pixels > `max_pix`) remains a known
  slow class.

## Cost (`boxes/*_task.json`, `totals.json`; session clock where stated)
- Box selection: 290 s (session clock, `run_params.json`).
- Search + calibration: 16:14Z → 16:46Z wall on 4 workers (session clock, `run_params.json`); 6,803 task-seconds; passed boxes 191–576 s each (median
  231 s), calibration-failed boxes 4,138 s in all; peak memory 1.41 GB per task.
- Review: 16:46Z → 16:55Z computation for 9 entries (≈ 1 min each, session clock), plus looking at 9 sheets and stamp
  sets and the catalogue checks (subagent in parallel).

## Pending
- h010: matched arc-pixel controls (`hit_lightcurve_controls.py pick-arc` + `compare`) on the 1208 → VENUS F444W pair,
  then a comparison with the Palencia+2026 event list. h003: matched compact-source controls on the same pair.
- The 10 calibration-failed boxes join the held set (26 in all); no follow-on batch is started.

## Regenerate
- `python scripts/batch_pairs.py --pairs results/transients/sfft_trial/venus_pairs.json --exclude … --prefer-extra-epochs
  --out boxes.json` (exclusions in `boxes.json` `rule`).
- `python scripts/batch_run.py run --boxes boxes.json --out boxes --scratch <dir> --workers 4 --seed0 300
  --timeout-screen 900 --timeout-inj 1200 --no-start-after 2026-10-17T18:20:00`, then `totals`.
- `python scripts/batch_review.py --batch boxes --boxes boxes.json --scratch <dir> --capacity 10 --sheet-verdicts
  boxes/review/sheet_verdicts.json`; then `post_review_checks_script.py.txt <dir> boxes/review` and
  `final_classes_script.py.txt boxes/review B3`.
