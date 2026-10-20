# W3 ZTF v23 — ordinary fades pass the v13 edge test mostly through unmeasured edges; the one pre-registered fix (bounds) removes every fade pass but loses 35–99 W3 recoveries per cell and fails its gate; no validation, no search

- **Pre-registration:**
  - plan, diagnosis items and the only alternative rule: da3411b, 10:22:21 UTC, before any light curve was read;
  - rule code and unit test committed with the diagnosis: fa8f19d, 10:33:57, before any `--edge bounds` run.
- **Labels:** DR24 light curves and TAP metadata are `observed`; injections and controls are `simulated`; counts
  are `derived`. **This is a development analysis** on calibration data already used in v22, not an independent
  validation.
- **Scope:**
  - a method diagnosis of the frozen v13 rule;
  - no search, no limit; nothing here bears on wormholes;
  - 8124, 8126, 1624 and 1625 are unchanged, with no photometry.

## A. Diagnosis (4991, fade 13.8 %, failed in v22; 2344, fade 8.1 %, passed)
- **Inputs:** only the 500 calibration groups and the comparison stars of each area were re-read; no search group
  (`dev<npix>/source_usage.jsonl`).
- **A0 Reproduction:** the v22 injections were re-run with the same seeds. Every cell and control count, and every
  per-injection dump row, equals v22 in both areas (`cost.json`: `a0_counts_equal`, `a0_dump_rows_equal`).
- **The mechanism (`edge_time()`):** a side with no significant intermediate detection within one run span is
  recorded as 0 d and passes. A1 looked at each event and, for each side, how close the data actually bound the
  edge (U: the time to the nearest normal detection).

### Passing fade controls (`dev<npix>/diag_summary.json`)

| | 4991 | 2344 |
|---|---|---|
| fades passing vetting / 1,000 | 138 | 81 |
| with a zero-day side (L = 0, U > 0.25 × span) | **129** | **76** |
| post side: ambiguous epochs / unusable exposures / no exposure / intermediate measured | 53 / 52 / 23 / 10 | 25 / 36 / 13 / 7 |
| pre side: measured fast / intermediate measured / unusable / no exposure / ambiguous | 19 / 59 / 40 / 19 / 1 | 12 / 40 / 11 / 18 / 0 |
| band check passed with no other-band epoch in the run | 38 | 7 |

- **What the side categories mean:**
  - **ambiguous epochs:** series epochs that are neither normal (dm < 0.3) nor significantly intermediate, i.e.
    within 2σ of the band;
  - **unusable exposures:** official exposures of the position that did not enter the series (infobits, sparse
    comparison coverage, flags, or missing);
  - **no exposure:** neither an ambiguous series epoch nor an unusable official exposure in the interval. Dim
    epochs, or significant intermediates beyond the one-span window, may still lie there (not stored).
- **Reading:** the recovery (post) side of a passing fade is never measured as fast (0 in both areas). In 128 / 138
  and 74 / 81 it passes because the data between the last dim epoch and the next normal detection hold no
  *significant* intermediate point under the fixed thresholds (the rest have a measured intermediate within 0.25 ×
  span). This has three
  observed sources (post side, 4991 / 2344: 53 / 25, 52 / 36, 23 / 13):
  - noise-level points;
  - unusable exposures;
  - no exposures.
- **Decision rule A:** 129 / 138 = 0.935 ≥ 0.5 (`diag_a3.json`), so the mechanism counts as important.
- **The same convention carries W3 recovery.** Of the recovered W3 injections, 35–99 per cell have a zero-day side
  (`diag_summary.json`, `w3`).
- **A2, same events in both areas** (1,000 fixed fade models, fixed carrier fractions and noise; `diag_a3.json`):
  - 4991 passes 142 (Wilson 95 % 12.2–16.5 %) and 2344 passes 75 (6.0–9.3 %), with 129 and 70 of them on a zero-day
    side.
  - The area difference therefore persists with identical event models (carriers and noise draws still differ
    between areas; t1 spans the overlap of the two areas' first-to-last r-band carrier epochs). It is not an artefact of the event draws, but
    comes from the two areas' carriers (sampling and data). Which property drives it was not isolated.
- **A3, v22 fade fractions with Wilson 95 % intervals** (`diag_a3.json`):
  - the 7 failed areas: from 8.5–12.2 % (1148) to 11.8–16.1 % (4991);
  - 2344: 6.6–10.0 %.
  - The intervals of 6 of the 7 failed areas include 10 %; only 4991's lies wholly above it. No area is re-judged.

## B. The pre-registered alternative ("bounds", `edge_bounds()`)
- **Rule:**
  - each side is **fast** if U ≤ 0.25 × span, **slow** if L > 0.25 × span, and **insufficient** otherwise;
  - the run is fast only if both sides are fast;
  - no new parameter. The code reads only (t, dm, nd, limit).
- **Same events, seeds and noise as A0** (`compare.json`):

| | 4991 v13 → bounds (insufficient) | 2344 v13 → bounds (insufficient) |
|---|---|---|
| fade / 1,000 | 138 → **0** (138) | 81 → **0** (81) |
| eclipse, PSPL, common mode | 0 → 0 | 0 → 0 |
| W3 cells / 200, range | 115–164 → **43–103** | 110–164 → **43–114** |
| W3 lost per cell | 48–99 | 35–91 |
| cells below 100 / 200 | 0 → **9 of 10** | 0 → **8 of 10** |

- **Discrimination:** bounds passes no fade in either area, while it still recovers 43–114 of 200 W3 per cell.
  All its losses are injections whose edges the data do not bound ("insufficient"). They stay in the denominator.
- **Gate B (fixed in PLAN B):** every cell ≥ 100 / 200 and fade ≤ 100 / 1,000 in both areas. **Failed**
  (`compare.json`: `gate_B_passed` false). The compare script also checked eclipse, PSPL and common mode ≤ 10 %
  (not in PLAN B; all 0, so no effect).
- **Post hoc:** the time stamp in PLAN's status note was corrected after the comparison (see PLAN).
- **As the plan requires:**
  - method changes stop here, and no second rule is tried;
  - **no validation area was chosen, no new light curve was read, and nothing was searched.**
  - The bounds code stays in `scripts/w3_ztf.py` as an option (`--edge bounds`); the frozen search rule remains v13
    (`--edge time`).

## C–D. Validation, search, survivors
- None: gate B failed.
- New searched sources: 0. Final survivors: 0. Reliable real dimmings: 0.
- The v13-rule cumulative count is unchanged (37,357 groups; not stored as one number).

## E. Cost (`cost.json`)
- **Downloads:** 1,221.3 MB in total:
  - two light-curve partitions streamed by `fetch` (estimated from parquet metadata): 573.4 + 465.2 MB (parquet metadata);
  - objects columns: 95.0 + 73.4 MB;
  - TAP: 14.3 MB.
- **Time:** from the start (10:18) to the gate result (10:35:54), about 18 min wall clock including planning. Stage times were
  hand-logged from the session (`cost.json`).
- The light-curve caches were deleted after use.

## F. What this shows and what it does not
- **Shown:**
  - With the v13 test set on ZTF DR24 sampling, ordinary fades that pass do so almost entirely through edges the
    data do not measure.
  - The same convention carries a quarter to two thirds of the recovered W3 injections per cell (35–99 per
    cell; `diag_summary.json`).
  - Requiring measured edges removes every fade pass (0 / 1,000) in these two areas, but loses 35–99 W3 recoveries
    per cell, leaving fewer than half recovered in 17 of 20 cells.
- **Not shown:**
  - that W3-like events cannot be found in these data;
  - that any searched star did or did not dim. Earlier v13-rule searches keep their results and wording.
- **Limits of the bounds rule as pre-registered** (found in review after the gate; not tested, no rule changed):
  - L is the v13 edge, so a noise-level significant intermediate beyond a closer normal detection still makes a
    side "slow". Some of the W3 losses may come from this definition rather than from unmeasured edges.
  - A "normal" detection (dm < INTER[0]) needs no significance, so a near-limit point can set U and make a side
    "fast". The 0 / 1,000 fade result may depend on this asymmetry.
  - Any rule that changes either definition is post hoc and needs calibration data not used here.
  - A2's fixed-event runs both flag 857 / 1,000 (different carriers); their pass counts differ (142 vs 75), so the
    runs are not duplicates. The equal flag count is read as a coincidence; nothing beyond that was checked.
- **Next step:** a method decision for the owner, not taken here. In these events the data hold no significant
  intermediate point and no close normal detection under the fixed thresholds (EDGE_FRAC 0.25, EDGE_NSIG 2); no
  threshold was varied. The options are:
  - keep v13, with its 7-of-8 area failure rate in new fields;
  - accept a lower recovery target for a measured-edge rule;
  - add information to the edges (for example, the official psfcat route of v9–v10 for unusable exposures).
