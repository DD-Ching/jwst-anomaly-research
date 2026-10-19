# W3 ZTF v13 — replace the point-count edge test by an edge-duration test; validate; search only if it passes

Plan and rule frozen 2026-10-19, before any development comparison was looked at. Owner authorization: the round after
#368 (4 h, ≤ 3 GB, one W3 edge rule only).

## 0. Checks of the existing record
- **#368 sample size:** the committed `../v12/2382/sample.json` holds 3,000 groups: 500 calibration and 2,500 search,
  16,499 objectids. The "1,000 groups" figure was a slip in the chat report only; no repo file says 1,000.
- **The 80 cases:** `../v12/2382/trace_posthoc.json` at ρ 0.1, blend 0.7, t_E 30 d gives 171 flagged, 84 recovered,
  and first vetting failures sharp-edge 80, achromatic 7. "Dense cadence causes it" stays a hypothesis.
- **Reuse:**
  - The repo has a full W3 / ordinary model fitter (`scripts/w3_microlensing.py`, D-057 / D-058; ΔBIC between PSPL,
    FSPL and the exotic models). It is built for brightening events, and a finite-source fit per flagged run would
    cost seconds each over thousands of runs, so it is not used.
  - The edge shape is taken from the W3 model itself (`w3_dmag`, i.e. `exotic_sim.inject_light_curve`), not from a new
    fitter.
  - No external literature search was made for this round (time); recorded as not done.

## 1. The replacement rule (one rule; `edge_time()` in `scripts/w3_ztf.py`, selected by `--edge time`)
- **The old rule** (`sharp_ok`) counted detections with 0.3 ≤ Δm < 1.0 within one run span before or after the run,
  and allowed ≤ 1. A denser cadence puts more epochs on the same physical edge, so the count depends on sampling.
- **The new rule measures the edge in time.**
  - Within the same windows, a detection is intermediate only if INTER[0] + 2σ ≤ Δm < 1.0, where σ is the expected
    error at that epoch from its maglimit (SNR 5 at the limit, the injections' noise model).
  - **Edge extent:** before the run, the time from the run's first dim epoch back to the farthest intermediate
    detection; after the run, the time from its last dim epoch forward to the farthest one.
  - **Pass if** both extents ≤ EDGE_FRAC × run span (span = max(last − first dim epoch, 1 d), as before).
- **EDGE_FRAC = 0.25, from the model only** (`edge_frac.json`): the 99th percentile over u0 of the W3 model's edge
  duration (Δm 0.3 → 1.0) divided by its Δm ≥ 1 span, for ρ 0.1 and blend 0.7 (0.248), rounded up. It is 0.005 and
  0.008 for ρ 0.01. No injection result was used to set it.
- **Nothing is given to the vetting from the injection:** the rule uses only the observed series (times, Δm,
  non-detection flags, maglimits). Non-detections are never intermediate; a missing epoch stays unknown.
- **Unchanged:** every other test (common mode, one-off, achromatic), the pre-screen, the thresholds, the model, the
  noise and the controls (eclipse, fade, PSPL, common mode), which go through the same new rule.
- **What the rule cannot do:** an ordinary dimming with short symmetric edges (e.g. a single eclipse-like event)
  passes it as before. It separates long, asymmetric recoveries such as the fade control.

## 2. Development comparison (already-used v12 calibration group; not a validation)
- Same events and the same noise: `inject --edge count` and `--edge time` with the v12 seeds on `v13/2382dev` (= the
  v12 sample). Count must reproduce v12. Per-injection outcomes are written for both (`dump_count.json`,
  `dump_time.json`).
- **Reported:**
  - recovered per cell, old vs new;
  - how many of the old sharp-edge losses the new rule keeps;
  - controls passing vetting, old vs new;
  - real calibration light curves (survivors, new rule);
  - recovered events lost by the new rule.
- The rule is not changed after this comparison.

## 3. Re-split and validation (plan change; `split_script.py.txt`, fixed before any search-group light curve is loaded)
- **Re-split:** the 2,500 never-read v12 search groups go to 500 validation groups (seed 13) and 2,000 search groups
  (`2382val/sample.json`). The 500 development groups are left out.
  - Groups bind all objectids of one source (0 shared objectids between the sets, `split_info.json`).
  - Control stars exclude every sample member.
  - The validation set shares the area and exposures with the development set, so it is **not an independent-area
    validation**.
- **Validation:** `inject --edge time --cells calibrated --n 200 --seed-offset 2000` (seeds 3100 + cell, controls
  3200 + k) on the 500 validation groups, plus the real validation light curves (`search --half dev --edge time`).
- **Pass, all required:**
  - every one of the 10 cells recovers ≥ 0.5 (all injections in the denominator);
  - every control class passes vetting ≤ 10 %;
  - ≤ 5 survivors among the real validation light curves.
- **On any failure:** record and stop. No seed change, no threshold change, no partial pass and no other area.

## 4. Search and image check (only if section 3 passes and review leaves nothing that changes the verdict)
- One search of the 2,000 frozen search groups (`search --half holdout --edge time`).
- **Image check:**
  - **Order, fixed now:** vetted survivors first (by key), then other flagged groups by the number of
    non-detections in their longest run.
  - **Limits:** ≤ 5 events and ≤ 64 cutouts (≤ 12 each), ≤ 50 MB; `../v6/image_script.py.txt` method.
  - **Sequence:** footprint, objectid binding and same-exposure comparisons first; then images before, during and
    after the run; ordinary explanations versus W3 only after a confirmed dimming.
  - Anything beyond the limits is recorded as pending. No image-confirmed dimming means it is not called real, and
    an unexplained one is not called a wormhole.

## 5. Budget and stop
- **Downloads:** the Npix 2382 partition again (1,846.5 MB, the cache had been deleted), the objects table and TAP.
  That is ≈ 2.0 GB, plus cutouts ≤ 50 MB.
- **Small inputs kept in the repo:** samples, the split, seeds, per-injection outcomes and decisions. The large caches
  are deleted after the round.
- **If the validation fails,** the W3 / ZTF remedy is recorded as unsuccessful and paused. Reopening it needs new data
  or a new method basis.
