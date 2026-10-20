# W3 ZTF v24 — existing products cannot restore the missing edge information: 8 % of undecided W3 edge sides are psfcat-fillable, optimistic gain 0–5 per cell; no download, trial closed

- **Pre-registration:** the plan was committed at 21:34:19 (6eaaaab), before any data read. It sets the categories,
  the optimistic estimate and the go / no-go rule.
- **Labels:** DR24 light curves and TAP metadata are `observed`; injections are `simulated`; counts are `derived`.
  The optimistic estimates use injection truth after classification. They are upper bounds, not recovery rates.
- **Scope:**
  - a development analysis on the v22 / v23 calibration data of 4991 and 2344 (no search group read);
  - #401's cross-area analysis is not redone;
  - v13, bounds and every threshold are unchanged;
  - 8124, 8126, 1624 and 1625 are untouched.

## 1. Per-event bookkeeping (`step12_<npix>.json`, events in `step12_events_<npix>.json.gz`)
- **Check:** replaying the v23 events gives the stored per-cell counts under both rules, and the stored
  "insufficient" counts, in both areas.
- **Bounds losses** (v13 recovered, bounds not): 773 in 4991 and 713 in 2344. **All are "edge insufficient"**, in
  every cell (`gate.json`), as expected because L ≤ U.
- **Under both rules, the W3 injections that are not recovered fail for these reasons** (per cell, `step12`):
  - the edge is slow (L > 0.25 × span);
  - another test fails: one-off or band check;
  - the event is not flagged.
- **Ordinary fades:**
  - 4991: 138 pass v13, 473 are slow, 268 fail another test (153 band, 115 one-off), 121 are not flagged;
  - 2344: 81 pass v13, 564 slow, 227 other (158 band, 69 one-off), 128 not flagged.
  - Under bounds all 138 / 81 become "insufficient"; nothing else changes.
- **The two bounds limits from v23 §F, counted per event:**
  - **(i)** A "slow" side where a normal detection lies closer to the run than the farthest significant
    intermediate. It affects 0–4 slow W3 events per cell, and 13 (4991) and 17 (2344) slow fades. These are fails
    under both rules; a corrected definition could at most turn them into insufficient or fast.
  - **(ii)** A bounds-recovered W3 run whose bounding normal detection is not significantly normal. It affects
    0–4 per cell (7 and 4 in total). No fade passes bounds, so (ii) moves no fade.
  - Neither limit is large enough to move any cell across 100 / 200 or the fade across 100 / 1,000. The limits
    stay recorded; they are not fixed.

## 2. What the missing edge information is (`gate.json`)
Each undecided side was classified within one span of the run edge, up to the bounding normal detection: a W3 side
that is insufficient under bounds, or a fade side of a v13 pass.

| class | W3 sides (both areas) | fade sides |
|---|---|---|
| A no official exposure covering the star | 1,274 (64.5 %) | 206 (64.2 %) |
| B only exposures with infobits ≠ 0 (excluded by quality) | 536 (27.2 %) | 75 (23.4 %) |
| C1 valid exposure, DR24 product gap (psfcat-fillable) | 157 (8.0 %) | 39 (12.1 %) |
| C2 valid exposure, not in the series for another reason | 3 | 0 |
| D in the series but too shallow or noisy | 4 | 1 |

- **Reading:**
  - Only C1 can be restored with existing products, and even there psfcat may lack the star (then it stays
    unknown).
  - A cannot be restored: there is no observation.
  - B is excluded by the official quality bits, which this trial does not relax.
  - C2 and D are negligible within one span. Unlike v23's interval with no span cap, ambiguous epochs mostly lie
    beyond one span.
- **Optimistic estimate (truth after classification; every C1 exposure within 0.25 × span counted as a normal
  detection when the model is at baseline and the exposure is deep enough):**
  - bounds: restorable W3 injections 0–5 per cell in 4991 and 0–4 in 2344. Optimistic bounds recovery reaches
    ≥ 100 / 200 in 1 of 10 cells in 4991 (104) and 2 of 10 in 2344;
  - v13: at most 6 of 138 (4991) and 6 of 81 (2344) passing fades could meet a significant intermediate in a C1
    exposure. The optimistic 4991 fade count is 132 / 1,000, still > 100.
- **Go / no-go (fixed in PLAN):** go only if in 4991 every cell could reach 100 under bounds, or the fade could
  reach ≤ 100 under v13. **Neither (`gate.json`: `go` false).** As planned, no psfcat product was downloaded and
  steps 3–4 did not run.

## 3. Cost (`cost.json`)
- 1,221.3 MB: the two calibration re-reads, the same files as v23. psfcat products: 0.
- About 8 min from the start (21:32:46) to the gate (21:40).
- Caches were deleted.

## 4. What this shows and what it does not
- **Shown:**
  - For these two areas, the edge information that the bounds rule finds missing is mostly absent from the archive
    (no covering exposure: 64 %) or flagged by the survey's quality bits (27 %).
  - The psfcat route of v9–v10 could fill at most 8 % of the undecided W3 sides. Even if every fill succeeded, it
    would not change any verdict under either rule.
- **Not shown:**
  - that W3-like events cannot be found;
  - that unflagged quality-B exposures hold no usable measurement (not tested; a quality relaxation is outside this
    round);
  - anything about areas other than 4991 and 2344, or edges beyond one span.
- **Restoration with existing products is closed for this rule set.** It should be reopened only with new edge
  information:
  - a denser-cadence or simultaneous survey covering the same stars (e.g. other ZTF bands or partner surveys,
    after a reuse check);
  - an owner decision on using infobits ≠ 0 exposures with their own validation;
  - an owner decision on the recovery target for a measured-edge rule.
