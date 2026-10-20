# W3 ZTF v23 — why ordinary fades pass the v13 edge test: a paired diagnosis, at most one edge-rule fix, one validation

Plan fixed 2026-10-20 at the start of the round (10:18 UTC), before any light curve is read. Owner authorization: after
#397, 4 h, ≤ 6 GB including review and records. No new inventory; the 10 % fade limit and every recovery threshold are
unchanged. No account, no external contact. 8124, 8126, 1624 and 1625 keep their status; no photometry on them.

## 0. Current behaviour to test (not assumed to be the cause)
- `edge_time()` (v13, `scripts/w3_ztf.py`): on each side of a dim run it measures the time to the farthest
  *significantly intermediate* detection (INTER[0] + 2σ ≤ dm < INTER[1]) within one run span. With no such point the
  side is recorded as **0 d** and passes (≤ EDGE_FRAC × span = 0.25 × span).
- The band check (`achromatic_ok`) passes when the other band has **no** epoch inside the run.
- The fade control (`control_model("fade")`): a linear 5–20 d decline to 3 mag, then an exponential recovery with
  τ = 30–200 d. The recovery spends about 1.2 τ between 1.0 and 0.3 mag, longer than 0.25 × span (span ≈ 1.1 τ + 0.67
  t_d). A dense post-run sampling should therefore show significant intermediate points. This is the hypothesis.

## A. Paired diagnosis (development analysis on already-used calibration data; not an independent validation)
- **Areas:** v22 4991 (fade 138 / 1,000, failed) and 2344 (fade 81 / 1,000, passed). Only their calibration groups
  and comparison stars are re-read (a calibration-only copy of each v22 sample in `v23/dev<npix>/`); no search group.
  Comparison stars are rebuilt from the full v22 sample, as in v22.
- **A0 Reproduction:** `inject --cells calibrated --n 200 --n-fade 1000 --edge time` with the v22 seed offsets
  (45000, 46000). Every cell and control count must equal the committed `v22/<npix>/injections_v22_time.json`.
  Otherwise stop and report.
- **A1 Per event** (`diag_script.py.txt` → `v23/dev<npix>/diag_events.json.gz`). The events are:
  - every fade control that passes vetting;
  - every W3 injection with a hit run, recovered or not;
  - every fade control flagged but failing.

  For each hit run and side (pre / post), from observed series information only:
  - **L** = the v13 edge (to the farthest significant intermediate within one span);
  - **U** = the time from the run edge to the nearest *normal* detection beyond the farthest intermediate. Normal
    means detected, with dm < INTER[0]; it is searched without a window. If there is none, U = ∞;
  - thr = EDGE_FRAC × span;
  - **the open interval between the run edge and that normal detection:**
    - the series epochs in it (neither normal, nor significantly intermediate, nor dim: "ambiguous");
    - the official exposures of the members' quadrants and band, with the position inside the footprint, that did
      not enter the series ("unusable": infobits, sparse comparison coverage, flagged, missing);
  - **category of a side with L = 0:**
    - `measured_fast` if U ≤ thr;
    - otherwise `ambiguous_epochs` if the interval holds series epochs;
    - otherwise `unusable_exposures` if it holds official exposures;
    - otherwise `no_exposure`;
  - **zero-day dependence:** L = 0 and U > thr on a side whose run passed the edge test;
  - **band check:** the other band's epochs inside the run; "passed by absence" if 0 and the test passed;
  - the outcome of each test.

  Injection truth (t0, t_E, τ, class) is written only after classification, for later evaluation; the
  classification code never reads it.
- **A2 Fixed-event comparison of sampling:**
  - Draw 1,000 fade models once (seed 2300) with t1 uniform over the overlap of the two areas' carrier time ranges,
    and t_d and τ as in `control_model`.
  - Pick carriers by a fixed fraction u ~ U(0, 1) (seed 2301; carrier = ⌊u × N⌋ of the area's carrier list). Each
    event's noise comes from seed 2302 + event index.
  - Apply the v13 rule in both areas. This separates sampling from event draws.
  - Reported with Wilson 95 % intervals; development analysis only.
- **A3 Uncertainty:** Wilson 95 % intervals for the fade pass fractions of all 8 v22 areas (committed decision
  files). No failed area is re-judged.
- **Decision rule (fixed now):** the mechanism counts as important if, in 4991, **≥ 50 % of the passing fade controls
  depend on a zero-day side** (A1). If it is below 50 %, report it, do not implement B, and stop method changes.

## B. The only alternative edge rule ("bounds"; implemented only if A supports it)
- **Per side:**
  - L and U as in A1; thr = EDGE_FRAC × span (unchanged 0.25; INTER and EDGE_NSIG unchanged; no new parameter);
  - **fast** if U ≤ thr; **slow** if L > thr; **insufficient** otherwise.
- **Per run:**
  - edge **fast** if both sides are fast;
  - **slow** if either side is slow;
  - **insufficient** otherwise.
  - `sharp_ok` = edge fast. Every other test is unchanged.
- **Basis:**
  - On the side's last normal detection the source was at baseline, and on the run's first dim epoch it was ≥ 1 mag
    fainter, so the transition lasted at most U. A data gap therefore bounds the edge only from above.
  - Rejecting slow edges is unchanged (L > thr fails in both rules).
- **Accounting:**
  - An injection is recovered only if a hit run passes every test.
  - An injection is "insufficient" if it is not recovered and a hit run passes every test except an insufficient
    edge. It stays in the denominator.
  - In a search, a run with an insufficient edge and every other test passed is listed as `insufficient`. It is
    neither a survivor nor excluded.
- **The code reads only the series:** (t, dm, nd, limit). Unit tests are committed with the rule, before any old/new
  comparison.
- **Comparison** (same events, seeds and noise as A0 in both areas; `inject --edge bounds`): per cell and per
  control, the v13 and bounds counts and the number insufficient.
- **Proceed to C only if bounds meets the unchanged v22 pass criteria in 2344 and 4991 together.** That means every
  cell ≥ 100 / 200 and fade ≤ 100 / 1,000 in both. Otherwise record the result and stop method changes; no second
  rule.

## C. One new validation, then search (only if B proceeds)
- **Area:** among the 17 qualifying but unqueued v22 partitions, in the frozen v22 queue order (no shared field
  first, searchable per GB), take the first whose full-sample fields (all members' objectid fields) share no field
  with the full samples of 4991 and 2344. The order is 1852, 4990, 1854, 1855, 5052, 1149, 5053, 4988, 2878, 2346,
  2345, 2876, 2813, 2812, 2879, 2347, 2877.
  - The area is chosen before any of its light curves is read. If none qualifies, stop.
- **Sample:**
  - the v22 `area_sample_script` rules, deduplicated against v1–v22 within 3″ or a shared objectid;
  - salt "v23:"; calibration = the lowest 500 by crc32;
  - search = the next ≤ 3,000 by the same order; groups beyond 3,000 are not used;
  - controls = the comparison range minus every sample member;
  - seed offset 48000.
  - Sample, rule and code are committed before any light curve.
- **Validation (once):** `inject --cells calibrated --n 200 --n-fade 1000 --edge bounds` and `search --half dev`.
  The pass criteria are the v22 ones: every cell ≥ 0.5, every control ≤ 10 %, ≤ 5 calibration survivors.
  - On a pass, search the frozen search group at once.
  - On a failure, leave it unread.
  - No seed change, retry or replacement area. A crash without a verdict may be rerun once, with the same input and
    seed, logged.
- **Stop rules:**
  - no validation start after 12:50 UTC;
  - no new download once the measured total plus the area's LC MB + 30 MB would pass 5,850 MB.
- **Results** under the new rule are a separate version ("v23 bounds"), not added to the v13-rule cumulative count.

## D. Survivor check (if searched)
- At most 5 final survivors and 50 cutouts, using the v21/v22 image and mask check.
- First list for each: the event band, that band's objectid, the ND expids, depth, mask and the same-band baseline.
- A miss is not zero flux; unvalidated photometry confirms nothing. No new PSF route.
- The last 30 min (from 13:48) are kept for review and records.

## Status before the old/new comparison (added before commit fa8f19d at 10:33:57 UTC; **post hoc correction** in b053f34, after the comparison: this
heading first read "added 10:38 UTC", which was wrong; only the time stamp changed; before any `--edge bounds` run)
- A0 reproduced v22 exactly in both areas (every cell, control and per-injection dump row).
- A decision: 129 of 138 passing fades in 4991 (0.935 ≥ 0.5) have a zero-day side, so B is implemented as written
  above, with no change (`edge_bounds()`, unit test `test_v23_edge_bounds_gap_is_insufficient_not_fast`).
- **Disclosure:** A1 also counts recovered W3 injections with a zero-day side (`dev<npix>/diag_summary.json`, 35–99
  per cell). Under B those become insufficient unless another passing run exists, so A1 already shows that B costs a
  large share of W3 recovery. B's definition was fixed at 10:22 (da3411b), before A1. It is run unchanged as planned.
