# W3 ZTF — ZTF alone vs ZTF + ATLAS on the 40 calibration stars: trial **not completed** (17 / 40 stars had ATLAS data by the deadline), but by the pre-fixed bounds it **fails whatever the missing stars show**: 6 / 63 insufficient W3 events restored (best case 27 / 63 < 0.50), 0 / 79 fade passes, 4 recovered W3 events lost

- **Plan:** `PLAN.md`. Method, rules, denominators, bounds and stop conditions were committed before any ATLAS flux of this round was read (8ab4a41). The review fixes (b768fac) came before any A / B outcome was read.
- **Labels:**
  - observed: ATLAS forced photometry (reduced mode, total flux) and the ZTF DR24 epochs, with their noise;
  - simulated: the injected W3 and fade events. Every "dim", "restored", "slow", "lost" or "pass" below is a simulated injection;
  - derived: all counts.
  - Nothing here is a real event, an astronomical anomaly or a wormhole clue.
- **Scope:** the already-used calibration events on the 40 stars (#459). The 4 #461 stars are development data. This is a feasibility trial, not a validation. No search group was read (`../v23/dev*/source_usage.jsonl`: `search_member_rows` 0).

## 1. Access and acquisition (`tasks.json`, `cost.json`)
- **Access:** the shared Astro Data Lab token documented by NOIRLab. The token is in no file.
- **Limits:** at most 2 of our tasks in flight. 19 tasks were submitted, none twice. 0 HTTP 429, 0 poll errors, none over 45 min.
- **Queue speed:** the shared queue ran our tasks one at a time:
  - queue wait 1.9–21.6 min;
  - run time 2.0–18.8 min.
- **At the 11:49 UTC deadline:**
  - 17 stars done;
  - 2 unfinished (tasks 5236004, 5236017; deleted by us after the deadline, unused);
  - 21 not submitted.
- **Which stars:** priority order put the stars with the most A-insufficient W3 events first. The 17 done stars carry 42 of the 63.
- **Data volume:** 6.96 MB of results. Compact tables are in `phot/` (17 tables of 38–129 KB, sha256 of each source in `phot/index.json`). The raw files were deleted after the tables were committed.

## 2. Measurements (`calib.json`; rules of PLAN §2, unchanged)
- **Rows:** 17 stars with 1,178–4,034 rows each; 88.7–97.0 % pass the FAQ quality cuts; 46,904 rows used.
- **Baselines:** 108 (filter, unit) groups. 88 have their own baseline; 20 fall back to the filter median (< 20 rows).
- **Error scale:** k (MAD of (rel − 1) / σ) is 1.48–8.68, median 2.19. No k < 1, so σ_cal = k × σ everywhere.
- **The value 9999:** 2 rows that pass quality have |uJy| ≥ 9999. They were excluded by the pre-fixed rule; keeping them changes nothing (§4).
- **Missing data:** a missing epoch is never zero flux. Negative and low-S/N rows keep their errors.

## 3. Reproduction
- A (ZTF-only bounds), replayed, equals the stored v24 outcome for every event of both areas: 0 mismatches over 2 × 3,000 events (`trial_<npix>.json` `checks`).
- B uses A's runs, spans and thresholds. It only adds the carrier's injected ATLAS points.

## 4. Primary result (`result.json`; both areas pooled; stars without ATLAS keep B = A and stay in every denominator)
| | A (ZTF only) | B (ZTF + ATLAS) |
|---|---|---|
| W3 recovered, of 154 | 48 | 50 |
| A-insufficient W3 restored, of 63 | — | **6** (9.5 %, Wilson 95 % 4.4–19.3 %) |
| A-insufficient W3 now "slow" (wrong class for a W3 event) | — | 4 |
| A-recovered W3 lost (→ "slow") | — | **4** |
| fade passes, of 79 | 0 | **0** (Wilson 95 % 0–4.6 %) |
| fade insufficient → slow, of 12 | — | 3 |

- **Conditional on the 17 stars with data:** 6 of their 42 A-insufficient W3 events are restored (14.3 %, Wilson 6.7–27.8 %).
- **Per area:** 4991 restores 3 / 39; 2344 restores 3 / 24. Neither area has a fade pass.
- **Without the 4 development stars:** 5 / 50 restored, 0 / 72 fade passes, 4 lost. The development stars alone restore 1 / 13.
- **Secondaries (not gating):** with the quoted σ, and with the 9999 rows kept, the outcome is the same for every event (6 / 63, 0 / 79, 4 lost).
- **Completion bounds** (PLAN §5 with the review fixes):

  | | lower | upper |
  |---|---|---|
  | restored | 6 / 63 | 27 / 63 = 0.43 |
  | fade passes | 0 / 79 | 4 / 79 |

  Success is false in both the best case and the worst case.

## 5. Why the change happens (`diag_<npix>.json`, **post hoc**; uses injection truth for evaluation only, changes no outcome)
- **The 6 restorations:** each one is decided by an ATLAS epoch closer to the ZTF run edge than any ZTF normal point.
  - Every such decisive point is where the injected model is at or above baseline (model dm −1.83 to +0.29).
  - So the gain is added sampling, not noise.
- **The 8 W3 events that became "slow"** (4 lost, 4 from insufficient):
  - **7 of the 8:** a single ATLAS epoch was read as a significant intermediate point beyond 0.25 span.
    - The injected model is at or above baseline there (dm −0.37 to −0.003).
    - The observed ATLAS flux itself is 32–79 % below its baseline (rel 0.21–0.68, 5.7–43 calibrated σ).
    - These are observed ATLAS dips that pass the FAQ cuts (artefacts or real variability; not vetted). The frozen rule cannot tell them from the event.
  - **1 of the 8** (2344, ρ 0.1, tE 10, #185): the injected W3 profile itself is intermediate there (model dm 0.61, 7.2 d after the edge vs a threshold of 6.5 d).
- **The 3 fades that became "slow":** decided by ATLAS epochs on the injected fade's true intermediate tail (model dm 0.33–0.88, observed rel 1.00–1.03). This is real added information.
- **The 32 still-insufficient W3 events on stars with data:**
  - 15 have no ATLAS epoch in any undecided side's decisive window (t_edge, t_edge + 0.25 span];
  - 17 have ATLAS epochs in at least one such window, but no undecided side gains a normal point there (S/N ≥ 3,
    dm < 0.3).

## 6. Verdict (PLAN §5)
- **Status:** **trial not completed.** 23 / 40 stars have no ATLAS data: 2 unfinished at the deadline, 21 never submitted.
- **Verdict:** **fail regardless of the missing stars.**
  - Even if every A-insufficient event on the 23 missing stars were restored, 27 / 63 = 0.43 < 0.50.
  - Fade passes stay within ≤ 0.10 at both bounds (0–4 / 79). The method fails on the restoration condition alone.
- **What follows:** this is not "completed but failed" in the sense of a full measurement. But completing the run cannot change the verdict.
- **No second try:** no seed, threshold or rule change, and no second method. An independent validation is **not** warranted by this trial.

## 7. Deviations
- **Time-hull request ranges** for the 38 stars submitted after 09:17. Added before any flux was read (PLAN "Deviation").
- **Scheduler restart** after a mistaken `pkill -f` (PLAN). Nothing was resubmitted.
- **Diagnostic:** `diag_script.py.txt` was written after `result.json` was read (post hoc). Its first run matched points with `np.isclose` (0.59 d tolerance at MJD ≈ 59,000) and was fixed to exact time matching before its numbers were used.
- **Unfinished tasks:** the 2 unfinished tasks were deleted after the deadline to free the shared queue. PLAN said "recorded for continuation". Their ids stay in `tasks.json`.

## 8. Limitations
- **Assumptions:**
  - lensing is achromatic, with the ZTF blend applied to ATLAS c / o;
  - ATLAS's wider PSF and chromatic blends are not modelled;
  - the injection rule is sky-limited at full dimming (OGLE `inject_w3`).
- **Coverage:** 17 / 40 stars; the missing 23 carry 21 of the 63.
- **Dips not vetted:** the observed ATLAS dips behind 7 W3 → "slow" moves are not vetted. A rule that needs two epochs for an intermediate point would be a new method for a new pre-registered round, not a fix to this one.
