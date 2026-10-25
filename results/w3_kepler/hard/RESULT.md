# W3 / Kepler hard tests: ordinary events fool the #500 shoulder rule; the local rule rejects them at a cost of 177 / 2,000 W3 recoveries; the 200-star search has 0 flags

The #500 rule is fooled by the deliberately hard cases: a single deep eclipse inside a real or
synthetic in-quarter brightening, and single unflagged bright cadences next to a dip. The frozen
local shoulder rule rejects all of them (0 / 100 in every family). It costs 8.9 % of W3 recoveries
(177 / 2,000 on the unread val stars), mostly at t_E = 10 d.

In 10 real deep EBs, none of the 2,826 eclipse windows with bounded edges has a shoulder under
either rule.

Validation passed. The search of 200 new stars found no flagged run.

- **Pre-registration:** `PLAN.md`, `scripts/w3_kepler_hard.py`, the local rule in
  `scripts/w3_kepler.py`, the tests and `sample_hard.json` were committed at 10:53:04 UTC (db00218),
  before any new flux or comparison.
  - A timing-wording fix followed before any data.
  - The first light-curve read was at 10:53:21 (`source_usage.jsonl`, `cost.json`).
- **Labels:** injections are simulated; light curves and their flags are observed; outcomes are
  derived. Hard cases give **conditional** false passes, not sky rates. Repeated injections on one
  star are one source (source counts given).
- **Deviations:** none in method, thresholds or seeds.
  - Phase A was run twice with an identical summary; the second run only moved the per-window runs
    into `phaseA_runs.json.gz`, because the plain file exceeded 1 MB.
  - The ledger is tracked as a gzip snapshot for the same reason.
  - §5 is post hoc (`posthoc_script.py.txt` → `posthoc.json`, `cases.png`).

## 1. Phase A: hard tests on seen backgrounds and real EBs (`phaseA.json`)
Both rules were evaluated on the same injections, 100 per family on the 23 old_bg stars.

| family (category) | old rule passed / 100 (Wilson 95 %) | local rule | sources with a pass (old) |
|---|---|---|---|
| real brightening + single eclipse (real background + deliberately placed event; 4 sources × 25) | **56** (0.46–0.65) | **0** | 4 / 4 (per source 10, 23, 22, 1 of 25) |
| single-point anomalies (+20 % cadence per side) | **72** (0.63–0.80) | **0** | 22 / 23 |
| synthetic 10–20 % brightening + eclipse | **86** (0.78–0.92) | **0** | 23 / 23 |
| long dips 0.1–40 d | 0 | 0 | — |
| quarter-boundary dips | 0 | 0 | — |
| the six #500 controls | 0 each | 0 each | — |

- **Real deep EBs** (real ordinary events; `eb_deep`, the 14 stars with a DR25 TCE depth ≥ 0.6):
  - 10 stars have flagged windows: 2,892 windows, 2,826 of them with common mode passed and both
    edges bounded "fast".
  - Of those 2,826, the old shoulder test passes **0** and the local test **0**; 137 are
    local-insufficient. So no single real eclipse window would pass either rule on its own.
  - The full chain passes 0 stars. One window (KIC 5262664, its only flagged window) passes the
    one-off test, and it fails the shoulder test (`posthoc.json`).
- **Real single events** (Kirk et al. 2016 table 9, 32 stars): none reaches the 1-mag dim threshold,
  so none is flagged. These events therefore cannot test the shoulder rule; they are rejected
  earlier, at run finding.
- **Decision:** the old rule fails 3 families, so the local rule is used (`decision.json`).
- **W3 on seen backgrounds:** the old rule recovers 77–90 / 100 per cell, the local rule 66–84.

## 2. Validation: the local rule on 25 unread stars (`validation.json`, `gate.json`)
- **W3:** all 20 cells recover **55–86 / 100**; the total is 1,494 / 2,000. Each cell is ≥ 0.50, so
  that part of the gate passes.
  - One cell's Wilson lower bound is below 0.5: ρ 0.1, blend 0.7, t_E 10 d at 55 / 100
    (0.45–0.64).
- **Controls:** all 10 families pass **0 / 100** on val, from 24–25 sources each (Wilson upper
  bound 0.037).
  - The required conditional controls also hold: real brightening 0 / 100 on 4 sources; real EB
    windows 0 / 2,892 and stars 0 / 10 (≥ 5 required).
- **Gate:** passes.
- **The old rule on the same val injections** (not gated): W3 76–91 per cell, 1,671 / 2,000;
  single-point 66 / 100 and synthetic brightening 84 / 100 pass.
- **Uninjected val stars:** 0 flagged runs.

## 3. Cost of the change (paired, same injections; `posthoc.json`)

| | val (unread) | phase A (seen) |
|---|---|---|
| W3 recovered, old → local | 1,671 → **1,494** / 2,000 (−177; gained 0) | 1,656 → 1,491 (−165) |
| lost per cell | 2–21 | 1–20 |
| lost at t_E = 10 d (4 cells) | 16–21 per cell | 13–20 |
| lost where the local test was "insufficient" | 141 of 177 | 126 of 165 |
| false passes removed (hard families) | single-point 66 → 0, synthetic brightening 84 → 0 | real brightening 56 → 0, single-point 72 → 0, synthetic 86 → 0 |

- **Where the loss comes from:** 141 of the 177 lost W3 events are local-insufficient. On at least
  one side the baseline window held fewer than 20 good cadences in the same quarter, or the
  shoulder window held none. The other 36 fail the locally-bright criterion.
  - Interpretation (not separately counted): long events have a baseline window of up to
    1.25 × span beyond the edge, which can run past the quarter end.

## 4. Search: 200 new stars (`search.json`)
- **Result:** 200 stars, 0 missing, **0 flagged runs**, so 0 survivors.
- **Vetting counts:** 0 checked real dimmings, 0 artefacts and 0 undetermined. No vetting was
  needed, and the pixel-file reserve was unused.

## 5. Post hoc observations (`posthoc.json`, `cases.png`, looked at)
- **Depth:** none of the 225 new stars (val + search) has any good cadence ≥ 1 mag below its quarter
  median. The lowest relative flux is 0.546.
- **Real brightenings:** 9 of the 225 new stars have real shoulder-like cadences (rel ≥ 1.10 at
  ≥ 5σ), similar to 4 / 100 in #500. These are the backgrounds on which the old rule was fooled.
- **Figure:** the old rule's false passes (the eclipse inside a real or synthetic brightening, and
  single bright cadences) next to a W3 event that the local rule loses at t_E = 10 d and one that
  both rules keep.

## 6. Not covered
- Flares next to dips, spotted-star dips, heartbeat stars, dippers, microlensing of blends, and
  pixel-level instrumental features (rolling band, SPSD).
- Real events shallower than 1 mag. All 32 Kirk single events are shallower.
- Real deep single eclipses: only 1 of the 2,892 real deep windows is a one-off. The single-window
  test above stands in for long-period systems.
- Single band: the achromatic test is "not tested". Efficiencies apply from the SAP product onward.
- Passing validation does not exclude every ordinary explanation. It shows that the tested
  families do not pass, under the stated conditions.
