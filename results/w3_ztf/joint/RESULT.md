# W3 ZTF — g and r together: the g epochs inside the event are real dims and tighten some edges, but the sampling allows 7 / 63 insufficient W3 events restored (noise-free replay; a generous post hoc oracle 16 / 63; needs ≥ 0.50); stopped at the worth check, no rule frozen

- **Pre-registration:**
  - PLAN §1 (question, evidence, stop conditions, the only candidate) was committed at 06:22:10 (be7b5c8). That was
    before the worth check ran (06:24:39).
  - **Deviation:** the dev-data re-read started at 06:20:26, before that commit. It read the same files as v23 / v24 /
    #459, not new flux.
  - **Deviation:** the code implementing §1 (`joint_extent`, `edge_bounds_joint`, the `edges` argument and the tests)
    was written before the run but committed with the results.
  - **Review rerun (06:34):** after the independent review, `worth_script.py.txt` was fixed and rerun once:
    - num / den are now computed exactly as §1 states, from U and L. In the first run they also moved the
      intermediate-search window. The totals did not change.
    - The reproduction checks are explicit.
    - The script now records the observed joint event count and the extension S/N.
    - It adds a **post hoc** oracle bound.
    - The committed `worth_<npix>.json` are from the rerun.
- **Labels:**
  - ZTF DR24 light curves are `observed`; injected events are `simulated`; counts are `derived`.
  - "Ideal" and "true dim" use injection truth for evaluation only. They are never given to the rule.
- **Scope:**
  - The worth check runs on the already-used calibration events of the 40 stars (#459). It is not a validation.
  - No search group was read (`../v23/dev*/source_usage.jsonl`: `search_member_rows` 0).
  - No rule was frozen, so PLAN §2–§5 were not reached: no development comparison, no new sample, no validation.
  - All "restored" or "fast" items below are simulated injections, not observed dimmings or candidates.

## 1. Reproduction (`worth_<npix>.json`: `checks`)
- **A0** is bounds on r alone. Replayed, it equals the stored v24 outcome for all 3,000 events per area (0
  mismatches).
- **B** is #459's method (g detections as edge points). It equals the stored #459 outcome for every event on the 40
  stars (0 mismatches).

## 2. Observed evidence on the 88 A0-insufficient W3 sides (`decision.json`)
- **Extended sides:** on 34 sides, g dim epochs (or deep g non-detections) continue the event past the r edge. On all
  34, a normal point (either band) follows the extension, so these sides show a normal point, dim points and the gap
  between them.
- **How much the extension moves things** (medians over those 34 sides):
  - the extension is 1.97 d;
  - the transition bracket (last dim to first normal) shrinks to 0.97 of B's;
  - the span grows to 1.06 of the r-only span.
- **The extension is real:** all 64 extending epochs fall where the injected model is ≥ 1 mag below baseline
  (evaluation). This includes the 16 on newly fast sides that condition 3 tests.
  - The g band sees the event continuing; the run is not lengthened by noise or by bad epochs.
  - The extending detections have S/N ≥ 42 (W3) and ≥ 33 (fades), so no near-limit point extends a run here.

## 3. Bracket or span? (`decision.json`)
Sides decided fast, of the 88:

| rule | fast | slow | insufficient |
|---|---|---|---|
| B (#459) | 5 | 0 | 83 |
| num: tighter bracket, r span | 8 | 0 | 80 |
| den: B's bracket, longer span | 9 | 0 | 79 |
| joint: both | 14 | 0 | 74 |
| joint on the noise-free model (ideal) | 14 | 0 | 74 |

- **Where the gain comes from:** joint decides 9 more sides than B (`new_fast_by_cause`).
  - 2 pass with the tighter bracket alone and not with the longer span alone.
  - 3 pass with the longer span alone and not with the tighter bracket alone.
  - 1 passes with either change alone.
  - 3 need both.
  - Condition 2 counts 3 sides that pass under num. This is real timing information, so the gain is **not only** a
    definition change. Most of it, though, needs the longer span.
- **What stays undecided:** no W3 side became slow. In this sample, the joint rule hid no slow-edge evidence. 74 sides
  stay undecided.
- **Fades:** of 18 A0-insufficient fade sides, joint makes 1 fast (B 0). Under the ideal, 0 fade events become fast
  at event level.

## 4. Decision (PLAN §1 stop conditions; `decision_script.py.txt`)
- **Condition 1 (pre-registered):** the noise-free replay of the joint rule at the observed epochs restores **7 / 63 =
  0.11** of the A0-insufficient W3 events (`ideal_joint_events_fast`). This is below 0.50, so condition 1 holds. The
  replay decides the same 14 sides as the noisy data, so noise is not what limits it.
- **Condition 2 does not hold:** 3 B-insufficient sides become fast under num.
- **Condition 3 does not hold:** 16 / 16 extending epochs on the newly fast sides are truly dim.
- **Post hoc, the replay is not a strict ceiling, so two more numbers** (`posthoc_bounds`):
  - The observed joint rule restores 7 / 63 events (0.11), the same as the replay.
  - **The generous oracle restores 16 / 63 (0.25).** The oracle knows every sampled epoch's true state in both bands and
    uses the true dm ≥ 1 duration as the span. It ignores intermediates, so slow evidence cannot hold it back. Even
    this oracle stays below 0.50: with these sampled epochs, no bracket-type rule at 0.25 × span can reach the #459
    condition.
- **Stop:** true. No rule is frozen, no development comparison or validation is run, and no second method is tried.

## 5. Reading
- **Does joint use add decidable information?** Yes, a little. The g epochs that fall inside the event are genuine dims
  next to the r edge (median 1.97 d). They tighten the bracket on 3 sides and lengthen the measured event on others.
- **Why it cannot be enough:**
  - Both bands share the long sampling gap after the edge. The bracket shrinks only to 0.97 of B's (median), and the
    span grows only to 1.06 of the r-only span.
  - Event level: the rule restores 7 of 63 events, and a generous oracle on the same sampling 16 of 63.
- **Limits:**
  - This checks one candidate definition (§1). A rule that also used non-adjacent g epochs or model fitting would be a
    different method, and is not tried.
  - Chromatic blending (bands crossing dm = 1 at different times) is not modelled; the injections use one blend.
  - The fade numbers are side-level only, because §3 was not reached.

## 6. Cost (`cost.json`)
- 1,221 MB: the calibration re-read of 4991 and 2344, the same files as v24 / #459.
- About 9 min from the start (06:19:32) to the decision (06:28).
