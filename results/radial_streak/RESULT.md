# Radial (W1) screen with the calibrated streak veto: result (2026-10-17)

Plan: `PLAN.md` (committed before any run). Code: `scripts/radial_more.py --streak-masks DIR` (candidate arcs within
1″ of a `findsat_mrt` status-2 trail are dropped after selection, `exotic_screens.drop_streak_arcs`; the null draws
use the same arcs). Labels: shapes `observed`; masks, arc classes, centres and null draws `derived`; every threshold an
ASSUMPTION. Table: `summary.json` (both runs per field).

## Inputs reproduced
- Masks regenerated with `scripts/w1_streak_mask.py` (the `.fits.gz` masks are not tracked): the masked-pixel count
  of every field equals the one recorded in `results/w1_streak_veto/inject_summary_<field>.json` (12 / 12).
- Without the veto, the five `results/radial_more` fields give the same `max_lines` and `p_random_max` as that run.

## Gate K (known case): **failed**
MACS0647's strongest centre (6 lines, p 0.075; five arcs on one faint linear streak, `results/radial_more`): the mask
covers 2 of its 6 contributing arcs (labels 793, 806). With the veto the centre stays (grid point 0.5″ away) with 5
lines (741, 754, 773, 789, 842), and the thinner null makes it p_random 0.03. Sheet `sheet_macs0647_veto.png`, looked
at: three of the five arcs lie on the same faint diagonal streak through the centre, which the mask does not contain
(W1 noted this missed streak at the MACS0647 S_max, `results/w1_streak_veto/RESULT.md`). The p < 0.05 centre is a
streak artefact in a field whose model fails the `arcs` control: no result, not a candidate.

So, per the plan, the streak veto is **not adopted** for `radial`; the runs below are a diagnostic only.

## Diagnostic: all twelve W1 fields, no veto → veto
| Field | model | anti (not model-radial) | vetoed | centres | max lines | p_random_max | flags |
|---|---|---|---|---|---|---|---|
| Abell 2744 | abell2744-bergamini23 | 134 | 0 | 30 → 30 | 5 → 5 | 0.505 → 0.505 | 0 → 0 |
| Abell 370 | abell370-cats | 100 | 22 | 53 → 39 | 5 → 5 | 0.495 → 0.37 | 0 → 0 |
| Abell S1063 | abells1063-cats | 32 | 2 | 6 → 6 | 3 → 3 | 0.95 → 0.96 | 0 → 0 |
| El Gordo | elgordo-caminha23 | 36 | 0 | 9 → 9 | 4 → 4 | 0.65 → 0.65 | 0 → 0 |
| MACS0416 | macs0416-cats | 120 | 7 | 81 → 81 | 7 → 7 | 0.225 → 0.215 | 0 → 0 |
| MACS0417 | macs0417-canucs-maps | 69 | 1 | 17 → 17 | 5 → 5 | 0.17 → 0.195 | 0 → 0 |
| MACS1149 | macs1149-cats | 61 | 0 | 12 → 12 | 3 → 3 | 1.0 → 1.0 | 0 → 0 |
| MACS1423 | macs1423-canucs-maps | 41 | 2 | 1 → 1 | 3 → 3 | 0.915 → 0.93 | 0 → 0 |
| RXJ2129 | rxj2129-lenstool-maps | 21 | 0 | 6 → 6 | 3 → 3 | 0.825 → 0.825 | 0 → 0 |
| SMACS 0723 | smacs0723-iclv2 | 30 | 5 | 9 → 5 | 4 → 3 | 0.98 → 1.0 | 0 → 0 |
| WHL0137 | whl0137-relics-lenstool | 26 | 3 | 2 → 2 | 3 → 3 | 0.705 → 0.735 | 0 → 0 |
| MACS0647 | macs0647-clash-nfw | 36 | 10 | 2 → 5 | 6 → 5 | 0.075 → 0.03 | 0 → 1 |

- In the ten control-passing fields the veto changes no flag and no `max_lines` except SMACS 0723 (4 → 3); Abell 370
  loses 22 of 100 arcs (its two trails cross the arc-rich core) and 14 of 53 weak centres.
- Removing arcs also thins the null, so `p_random` can fall when a centre survives (MACS0647): a partial veto can
  manufacture a flag. Any future streak veto for `radial` needs its own injections and a check on known streaks.
- The seven D-049 fields here use `inject_shear.field_spec` inputs with the Gaia spike veto, not the D-049 inputs, so
  their no-veto rows are not the D-049 numbers.

## Reading
No exotic-lens candidate; nothing for `/vet-candidate`. The `radial` conclusions of `results/radial_more` and D-049
stand. The faint-streak problem stays open: `findsat_mrt` at the W1 settings misses the faint MACS0647 streak.

## Limits of this result
- One known streak case only; the mask's completeness for faint streaks is not measured.
- MACS0417's CANUCS map failed its sha256 on the first fetch (truncated download) and passed on the retry, as before.
