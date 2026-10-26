# W3 / Kepler fixed-rule search of 2,000 new stars — frozen plan

Owner authorization 2026-10-26 (about 04:45 UTC). Limits:
- 6 h;
- 18 GB of downloads in total, through `scripts/net_ledger.py` (this round's ledger `net_ledger.jsonl` here, cap
  18,000,000,000 B), with at least 2 GB kept for vetting and repairs;
- 4 cores;
- large cache ≤ 3 GB.

This file, `sample_2000.json`, `scripts/w3_kepler_2000.py` and its tests are committed before any new flux is
read. This is a fixed-rule exploratory search. A pass is not wormhole evidence, and a zero is not an exclusion.
Labels: light curves and flags are observed; injections simulated; counts, efficiencies and limits derived.

## 1. Sample (`sample_2000.json`, metadata only)
- **Pool:** the #500 / #542 FGK pool (pinned sha256 `013b8a0d…`): 18,570 stars with all of Q1–Q17.
- **Excluded:** every star an earlier Kepler W3 round read, 872 distinct (830 of them in the pool):
  - #500: development, validation, search and the interface star;
  - #529: old background, deep and single EBs, validation and search;
  - #542 / #546: the 500 stars.
  - That leaves 17,740 candidates.
- **Selection:** sort by sha256(`w3-kepler-s2000-2026-10-26:<kepid>`) and take the first 2,000.
  - Batch b holds ranks 500(b − 1) … 500b − 1.
  - The 36 stars in the Kirk et al. 2016 EB catalogue are kept.
  - No star is dropped or replaced because of a result. A star whose download fails is recorded as unread.

## 2. Data
- **Per star:** the archive listing plus every Q1–Q17 `_llc.fits` file, through the ledger. Each file's name, size
  and sha256 go to `files_b<b>.json.gz`.
- **Download change (efficiency rule 4.3):** the Q0 file is not fetched. It is recognised by its timestamp
  `2009131105131` (true of every Q0 file of #542), and #542 never used Q0.
- **Columns and cache:** as `w3_kepler._fetch_star`. A star with any failed file (after 2 retries) is unread and
  gets no cache.
- **Batch guard:** a batch is fetched only if ledger use + 3.9 GB ≤ 16 GB. The rest of the 18 GB cap is for
  vetting (pixel files) and repairs.

## 3. Search (unchanged #529 / #542 chain)
- **Screen:** `w3_kepler.screen_star(local=True)`; survivor = `pass_local`. This includes:
  - the quality mask;
  - per-quarter normalisation;
  - 1-mag dim threshold, ≥ 3 cadences, 20-cadence brackets;
  - common mode < 0.2;
  - bounded fast edges (v23 bounds rule, EDGE_FRAC 0.25, 1-d floor);
  - one-off test;
  - local shoulder rule.
- **Common-mode set of batch b:** the read stars of batch b, leave-one-out (`w3_kepler_expand.common_counts`,
  `frac_loo`). The same set and definition are used in the batch's injections.
- **Stored per batch:**
  - layers (`layers_b<b>.json`, built as #542 `layers`);
  - every star (`screen_b<b>.json.gz`);
  - every flagged run (`runs_b<b>.json.gz`);
  - the survivors in vetting order.

## 4. Injection-recovery per batch (the #546 method, unchanged)
- **Method:** `w3_kepler_effsearch.inject_star`, with 2 injections per star and W3 cell (20,000 per full batch).
  - Each injection starts from the original SAP series, before the quality mask and normalisation.
  - Event centre t0 is uniform over the star's window (first to last finite Q1–Q17 TIME, gaps included); u0 is
    uniform on [0, 2).
  - Noise, exposure average, recovery definition (matched run with `pass_local`) and loss stages are as #546.
- **Only change:** the rng key, `default_rng([54700, global rank, cell, rep])`.
- **Process pool:** 3 processes while the next batch downloads, 4 otherwise. Outcomes do not depend on this.
- **Stored per batch:**
  - `events_b<b>.json.gz` (every event);
  - `summary_b<b>.json`:
    - observing time;
    - sampling time;
    - recovery rate;
    - effective search volume E = Σ T_i k_i / 2 on the draw window (never T_obs × rate);
    - MC interval;
    - per-star cluster mean and SE;
    - loss stages;
    - limits.
- **No gate:** no "≥ 50 %" requirement. No seed change, no dropped cell, no changed threshold. A low efficiency is
  a result. Only a data-read error or a step not run as planned is repaired before reporting.

## 5. Combining batches (`combined.json`)
- **Allowed:** effective volumes are summed only over batches with the same rule, event distribution, injection
  method, recovery definition and common-mode definition:
  - the batches of this round;
  - the #546 batch (500 stars, same definitions; comparison set = its own batch).
- **Not allowed:** the earlier 60 (#500 rule) and 200 (#529 search with another comparison set and no #546-style
  injection) are not added.
- **Reported separately:** stars searched, observing time (good cadences × 29.4 min), sampling time (windows, gaps
  included) and effective star-years.

## 6. Vetting (survivors only; at most 10 distinct sources in the round)
- **Order:** batch, then rank, then run start.
- **Diagnostics per survivor run** (`w3_kepler_2000.vet_diagnostics`, computable; thresholds fixed here):
  - **V1 quality:** share of cadences within ±1 d of the run (same quarter) with any SAP_QUALITY bit set; ≥ 0.5 is
    flagged.
  - **V2 centroid:** MOM_CENTR1 / 2 median shift in the run, against the same quarter outside ± one span, in
    robust σ. |shift| ≥ 5 σ calls for the pixel check.
  - **V3 common:** same-channel, same-quarter stars of the batch at the run's cadences. Flagged if their median
    relative flux is ≤ 0.9, or if ≥ 3 of them are ≥ 1 mag dim.
  - **V4 data edge:** run edge within 1 d of the quarter start or end, or of a ≥ 0.5-d gap (context).
  - **V5 repeating:** in the Kirk 2016 catalogue, or ≥ 2 other runs (≥ 3 cadences) ≥ 0.3 mag dim outside the run
    ± one span (EB-like).
  - **V6 flare-like:** any one of —
    - the shoulder excess ratio pre / post is outside [1/3, 3];
    - the pre-side shoulder peak lies ≥ 3 cadences from the edge;
    - ≥ 3 other brightenings (≥ 2 cadences at rel − 1 ≥ max(0.01, 5σ)) outside the run.
  - **V7 baseline step:** the local-shoulder baseline medians on the two sides differ by ≥ 5 %.
  - **V8 PDCSAP:** relative flux in the run. A diagnostic only, never an independent confirmation.
  - **V9 pixel check:** only when V1 and V3 are clear, and V5 and V6 are clear or V2 is flagged.
    - Use the run quarter's official target pixel file, within the reserve.
    - Mean in-run image minus mean image of the same quarter within one span + 1 d outside the run.
    - "On target" if the most negative difference pixel is within 1 pixel of the out-of-run flux-weighted centroid;
      "off target" otherwise.
- **Classification** (`vet_class`):
  - **artifact supported:** V1, V3 or V9 off target;
  - **ordinary explanation supported:** V5 or V6 (periodic eclipses or flares);
  - **real dimming confirmed** (unexplained under V5–V7): V9 on target and no V7. That is still not exotic
    evidence, and ordinary explanations are not all excluded;
  - **undetermined:** otherwise.
  - Unvetted sources are listed as not vetted.

## 7. Limits (product level, conditional, per cell)
- **Premise:** zero survivors after vetting in the combined batches.
- **Joint limit:** `joint_limit` uses a flat prior on the rate and integrates the efficiency over its Beta(K + 1,
  n − K + 1) posterior (K = ε n_eff, n = Kish n_eff). It treats the draws as independent given the stars, which is
  conservative against the Poisson-binomial.
  - For an exactly known efficiency it equals the classical −ln 0.05 / E = 2.996 / E.
  - It is a Bayesian credible limit. It handles the count and the efficiency uncertainty in one step instead of two
    stacked 95 % steps.
- **Check of #546:** `limit_check_546.json`, from #546's stored summary; no injection re-run.
  - #546's 2.996 / E_low stacks two 95 % steps. It is conservative: 1.031–1.049 × the joint limit.
  - It is therefore relabelled "approximate, conservative" in this round's record.
- **Not done:**
  - no physical rate, abundance, wormhole number or FTL statement;
  - per-cell limits are not simultaneous.
  - Any "more data improves the limit by X" assumes a zero result and similar efficiency in the new sample.

## 8. Stop rules and pipeline
- **Next batch:** batch b + 1 is fetched only after batch b's screen.
  - If batch b has survivors, b + 1 waits until they are vetted.
  - If batch b has none, b + 1 downloads while b's injections run.
- **Too many survivors:** if the round's distinct survivor sources pass 10 after a batch's screen, that batch's
  records (layers, injections, summary) are finished, the first 10 are vetted, the rest are listed as pending, and
  no later batch is fetched.
- **Real dimming:** a reliable real dimming without an ordinary explanation (class "real dimming confirmed") stops
  later batches. The evidence and the smallest next checks are delivered; no stronger conclusion is drawn.
- **Budget:**
  - the ledger guard (section 2);
  - if vetting or repairs would pass 18 GB, stop and report;
  - no batch fetch starts after 09:15 UTC, to leave time for injections, vetting, review and merge within 6 h.
- **Cache:** a batch's column cache is deleted after its screen, injections, summary and vetting records exist
  (sources, hashes, events). Survivors' files are kept until vetted. At most two batches are on disk at once.
- **No reaction to results:** no new seed, no wider sample, no rule change, no second tool because of a result.
