# W3 ZTF — g and r together: the g epochs inside the event are real dims and tighten some edges, but even ideal use of the sampling restores 7 / 63 insufficient W3 events (needs ≥ 0.50); stopped at the worth check, no rule frozen

- **Pre-registration:**
  - PLAN §1 (question, evidence, stop conditions, the only candidate) was committed at 06:22:10 (be7b5c8). That was
    before the worth check ran (06:24:39).
  - **Deviation:** the dev-data re-read started at 06:20:26, before that commit. It read the same files as v23 / v24 /
    #459, not new flux.
  - **Deviation:** the code implementing §1 (`joint_extent`, `edge_bounds_joint`, the `edges` argument and the tests)
    was written before the run but committed with the results.
- **Labels:**
  - ZTF DR24 light curves are `observed`; injected events are `simulated`; counts are `derived`.
  - "Ideal" and "true dim" use injection truth for evaluation only. They are never given to the rule.
- **Scope:**
  - The worth check runs on the already-used calibration events of the 40 stars (#459). It is not a validation.
  - No search group was read (`../v23/dev*/source_usage.jsonl`: `search_member_rows` 0).
  - No rule was frozen, so PLAN §2–§5 were not reached: no development comparison, no new sample, no validation.
  - All "restored" or "fast" items below are simulated injections, not observed dimmings or candidates.

## 1. Reproduction (`worth_<npix>.json`: `checks`)
- **A0** is bounds on r alone. Replayed, it equals the stored v24 outcome for all 3,000 events per area.
- **B** is #459's method (g detections as edge points). It equals the stored #459 outcome for every event on the 40
  stars.

## 2. Observed evidence on the 88 A0-insufficient W3 sides (`decision.json`)
- **Extended sides:** on 34 sides, g dim epochs (or deep g non-detections) continue the event past the r edge. On all
  34, a normal point (either band) follows the extension, so these sides show a normal point, dim points and the gap
  between them.
- **How much the extension moves things** (medians over those 34 sides):
  - the extension is 1.97 d;
  - the transition bracket (last dim to first normal) shrinks to 0.97 of B's;
  - the span grows to 1.06 of the r-only span.
- **The extension is real:** all 64 extending epochs fall where the injected model is ≥ 1 mag below baseline
  (evaluation). The g band sees the event continuing; the run is not lengthened by noise or by bad epochs.

## 3. Bracket or span? (`decision.json`)
Sides decided fast, of the 88:

| rule | fast | slow | insufficient |
|---|---|---|---|
| B (#459) | 5 | 0 | 83 |
| num: tighter bracket, r span | 8 | 0 | 80 |
| den: B's bracket, longer span | 9 | 0 | 79 |
| joint: both | 14 | 0 | 74 |
| joint on the noise-free model (ideal) | 14 | 0 | 74 |

- **Where the gain comes from:** joint decides 9 more sides than B.
  - 3 of them already pass with the tighter bracket alone. This is real timing information.
  - 4 pass with the longer span alone.
  - 2 need both.
  - So the gain is **not only** a definition change.
- **What stays undecided:** no W3 side became slow. In this sample, the joint rule hid no slow-edge evidence. 74 sides
  stay undecided.
- **Fades:** of 18 A0-insufficient fade sides, joint makes 1 fast (B 0). Under the ideal, 0 fade events become fast
  at event level.

## 4. Decision (PLAN §1 stop conditions; `decision_script.py.txt`)
- **Ideal ceiling:** the noise-free joint call at the observed epochs restores **7 / 63 = 0.11** of the A0-insufficient
  W3 events (`ideal_joint_events_fast`).
  - Condition 1 holds: the existing sampling cannot reach 0.50 even with perfect handling. The ideal decides the same
    14 sides as the noisy data, so noise is not the limit; sampling is.
- **Condition 2 does not hold:** 3 newly decided sides need no span change.
- **Condition 3 does not hold:** 64 / 64 extending epochs are truly dim.
- **Stop:** true. No rule is frozen, no development comparison or validation is run, and no second method is tried.

## 5. Reading
- **Does joint use add decidable information?** Yes, a little. The g epochs that fall inside the event are genuine dims
  next to the r edge (median 1.97 d). They tighten the bracket on 3 sides and lengthen the measured event on others.
- **Why it cannot be enough:**
  - Both bands share the long sampling gap after the edge. The bracket shrinks only to 0.97 of B's (median), and the
    span grows only to 1.06 of the r-only span.
  - Event level: even ideal handling restores 7 of 63 events.
- **Limits:**
  - This checks one candidate definition (§1). A rule that also used non-adjacent g epochs or model fitting would be a
    different method, and is not tried.
  - Chromatic blending (bands crossing dm = 1 at different times) is not modelled; the injections use one blend.
  - The fade numbers are side-level only, because §3 was not reached.

## 6. Cost (`cost.json`)
- 1,221 MB: the calibration re-read of 4991 and 2344, the same files as v24 / #459.
- About 9 min from the start (06:19:32) to the decision (06:28).
