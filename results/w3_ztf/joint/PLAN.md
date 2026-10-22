# W3 ZTF — do g and r together bound the event better? (one bounded method trial)

Owner brief 2026-10-22 after #459: 3 h, ≤ 3 GB. One method choice only. No new source inventory, no ASAS-SN retry, no new
photometry pipeline; X4 and the other paused units stay as they are. No search group is read and no search is run.
Session start 06:19 UTC.

## 0. Prior state (checked 06:19)
- `main` at 0eb3589 (#459 merged). Open PRs: #458 (MACS0416 claim), #442 (dependabot), `needs-human` backlog; nothing on
  W3 / ZTF.
- #459 (`../altdata/`): bounds on r alone (A) and r plus g detections as edge points (`edge_bounds_alt`) on the replayed
  v23 events of the 40 stars. 3 / 63 insufficient W3 events restored; fade passes 0 / 79. Post hoc: 59 of the 66 g epochs
  in the W3 windows fall inside the injected event.
- The saved #459 files hold per-side counts and states, not epoch times. So step 1 needs the v23 dev re-read (the same
  files as v23 / v24 / #459, about 1.25 GB).

## 1. Is it worth a rule? (step 1, ≤ 30 min; `worth_script.py.txt` → `worth_<npix>.json`)
**Question:** can the g epochs that fall inside the event, used with r, bound the fade / recovery times more tightly?
Or do they only lengthen the event, so that edge / span shrinks and the test passes more easily?

**Evidence used (observed series only; injection truth only where labelled "evaluation" or "ideal").** For each
bounds-insufficient side (run passing the other tests) of the 40 stars' events, A's run (r) is extended outward through
contiguous dim epochs of either band, in time order:
- **Dim:**
  - a detection with dm ≥ 1 against its own band's baseline;
  - a non-detection, only where the exposure's limit is ≥ 1 mag below that band's baseline, so that it is an upper
    limit implying dm ≥ 1.
- **Stop:** the extension ends at the first non-dim epoch of either band (normal, intermediate or ambiguous).
- **Unknown:** missing and quality-excluded exposures are never dim, never normal, and do not stop the extension.

Per side, these numbers are recorded:
- `U_r`, the bracket of the transition from A's run edge to the next r normal;
- `U_B`, as in #459 (r or g normal);
- `U_J`, from the extended edge to the next normal of either band, beyond the farthest significant intermediate of
  either band;
- `L_J`, the farthest significant intermediate of either band;
- the spans `span_r` and `span_J`;
- the evidence class: whether a g dim point lies after the r edge and a normal point (either band) follows it.

**Decomposition** (the side is fast if U ≤ 0.25 × span):
- **num:** `U_J` against 0.25 × `span_r`. The tighter bracket alone.
- **den:** `U_B` against 0.25 × `span_J`. The longer span alone.
- **joint:** `U_J` against 0.25 × `span_J`.

**Ideal upper bound (truth, evaluation only):** the same joint rule applied to the noise-free model at the observed
epoch times and limits, for the 63 A-insufficient W3 events.

**Stop (no rule is written) if any holds:**
1. The ideal joint restoration is < 0.50 of the 63, so the existing sampling cannot meet the #459 condition even with
   perfect handling.
2. The joint gain over B comes only from the span: no side that is insufficient under B becomes fast under `num`.
3. Most of the gain rests on extensions that are not real. An extension is real when the extending g dim epochs have
   true dm ≥ 1 (evaluation). The condition is: fewer than 90 % of the extending epochs on sides that become fast are
   real.

Otherwise a single rule is frozen (section 2) before any validation flux is read.

**The only candidate:** section 2 can freeze only the joint definition written here. No threshold, window or extension
choice is made after the step-1 numbers. Why a longer span is not a loosening by itself:
- The extension adds only observed dim epochs. Under the injection model's assumption of one source fraction in both
  bands, `span_J` is still a lower bound on the true dm ≥ 1 duration (EDGE_FRAC's denominator, v13), and `U_J` is
  still an upper bound on the edge duration.
- So "fast" stays a conservative call. Section 1's check 3 tests whether the extension is real.
- With chromatic blending the two bands cross dm = 1 at different times. This is a stated limit.
