# W3 ZTF v24 — can existing products restore the missing edge information? (limited trial on 4991 and 2344)

Plan fixed 2026-10-20 at the start of the round (21:32 UTC), before any data read. Owner authorization: after #402,
2 h, ≤ 2 GB including re-reads, review and records; ≤ 200 new psfcat products; no new photometry tool, no account,
no new search source. #401 (cross-area fade spread, `results/w3_ztf/fade_diag/`) is read, not redone. v13, bounds and
every threshold are unchanged; the bounds limits listed in `../v23/RESULT.md` §F are kept. 8124, 8126, 1624, 1625:
no photometry. Development analysis on already-used calibration data; not an independent validation.

## 1. Per-event bookkeeping (from #402's events; `step1_script.py.txt` → `step1.json`)
- Re-read 4991 and 2344 calibration groups and controls exactly as v23 (`../v23/dev_setup_script.py.txt`), replay
  the v23 events (same seeds) and, per W3 injection with a hit and per flagged fade, record under each rule:
  recovered / edge slow / edge insufficient / another test failed (the first failing of common mode, one-off, band).
- Check: per-cell counts equal `../v23/dev<npix>/injections_v23_{time,bounds}.json`, else stop.
- Report: whether every bounds loss (v13 recovered, bounds not) is "insufficient" (expected by construction, L ≤ U).
- **Bounds limits (v23 §F), counted per event, nothing changed:**
  - (i) "slow" set by a significant intermediate lying beyond a closer normal detection;
  - (ii) a recovered or passing run whose U is set by a normal detection that is not significantly normal
    (dm + 2σ ≥ INTER[0]).

## 2. Potential of existing products (≤ 30 min, no flux read; `step2_script.py.txt` → `step2.json`)
For every side that is zero-day under v13 or insufficient under bounds, the interval (run edge, bounding normal
detection; capped at one span) is split by the official r-band exposures of the star's quadrants:
- **A no exposure:** no official exposure whose footprint (10″ margin) holds the star.
- **B excluded by quality:** exposures only with infobits ≠ 0.
- **C valid, not in the series:** infobits 0, star in the footprint, but the exposure is not `good` (sparse
  comparison coverage) or the star's entry is flagged/unknown. C1 = DR24 has no detection of any control star of the
  quadrant in it (a product gap, the psfcat-fillable case of v9–v10); C2 = the rest.
- **D in the series but too shallow / noisy:** series epochs that are neither normal nor significantly intermediate.
- A side's class is the best one present, in the order C1, D, C2, B, A (C1 first: the only one this trial can act on).

**Optimistic estimate (truth used for evaluation only, never by a classifier):** a lost (bounds) W3 injection is
"possibly restorable" if every insufficient side of a hit run that passes the other tests has a C1 exposure within
EDGE_FRAC × span of the run edge where the true model dm < INTER[0] and the exposure maglimit ≥ star baseline +
INTER[0]. A v13-passing fade is "possibly rejectable" if a zero-day side has a C1 exposure within one span where the
true dm is in [INTER[0], INTER[1]). These are upper bounds, not recovery rates.

**Go / no-go for step 3 (fixed now):** go only if, in 4991, either
- bounds: optimistic recovered (bounds + restorable) ≥ 100 / 200 in every cell, or
- v13: optimistic fade passes (138 − rejectable) ≤ 100 / 1,000.
Otherwise stop downloads and report the numbers. No threshold or quality rule is relaxed.

## 3. One small psfcat fill (only on go)
- Product list fixed before any flux read: the C1 exposures of the two areas' calibration quadrants (control-star
  product gaps, `fillset` rule of `../v10/fillset_script.py.txt`), ordered by time, ≤ 200 products in total, plus 2
  overlap exposures per quadrant within that cap.
- The v10 entry and merge (`merge_psfcat`): right band, position and zero point; DR24 never overwritten; same-exposure
  overlap test as v10 (n ≥ 8 pairs, |median Δ| and robust scatter as in v10); a missing star is unknown; a failing
  quadrant is not used. Each product downloaded once. Stop at 200 products, 2 GB or 23:00 UTC.

## 4. Before / after on the same events (only if merge validation passes)
- Same events and the original epochs' noise: noise from one generator per (cell or control seed, event index, expid),
  so added epochs do not shift other draws; the "before" run uses the same scheme (its counts are reported against
  the stored ones).
- W3 and the fade control (and the other controls) through the same fill; v13 and bounds reported separately:
  measurements added, sides turning from unknown to decided, recovered per cell, fade passes, remaining losses by
  class A/B/D or the rule. Every planned injection stays in the denominator. Partial fills are reported as partial.
- No search group is read in any case.

## 5. Close
- From 23:12 UTC only review and records. CHANGELOG, DECISIONS (D-104 addendum), TASKS.
