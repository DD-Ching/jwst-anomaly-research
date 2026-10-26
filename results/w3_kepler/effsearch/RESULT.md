# W3 / Kepler effective search volume of the #542 sample: 1,249–1,620 effective star-years per W3 cell; product-level per-cell limits 1.9–2.5 × 10⁻³ per star-year, conditional and not a physical rate

- **What this is:** an injection-recovery measurement of search efficiency on the 500 stars #542 already searched.
  These light curves had been looked at, so it is **not a blind test**. The screen, thresholds and tests are those of
  #529 / #542, unchanged. No new star was read and no shoulder repair was made.
- **Pre-registration:**
  - `PLAN.md` was committed at 23:48:53 UTC (6e864af), before any data restore or injection.
  - The check, inject and summarise code and its tests were committed at 23:52:58 (fde145c).
  - The injection run started at 00:09:56 (`source_usage.jsonl`, `cost.json`).
- **Labels:** injections are simulated; light curves observed; recoveries, efficiencies and exposures derived.
  An efficiency is not evidence for or against any physical object.
- **Deviations:** none. Notes:
  - `examples.png` / `examples.json` (a look-at figure) and `posthoc.json` (loss diagnostics) were added after the
    run. They are **post hoc** and change no outcome.
  - After the independent review, `posthoc.json` gained the EB, grazing and u0 blocks (post hoc). Three unit tests
    were added: the match condition (c), the cascade after flagging, and the spurious flag. No stored event changed.
  - The figure was looked at. One `not_flagged` example looked findable. Checked: its dim run reaches the last
    cadence of Q17, so the 20-cadence bracket after it is empty. The classification is correct; it is not a code
    error.

## 0. Restore and check (`restore_files.json.gz`, `check.json`)
- **Files:** 8,500 Q1–Q17 files were fetched by the names in `expand/files_expand.json.gz` (3,719,520,000 B).
  Every sha256 matches, with 0 fetch errors.
  - The 245 Q0 files were not fetched (efficiency rule 4.3; #542 never used them).
  - Ledger: 3,791,889,000 B of the 6,000,000,000 B cap.
- **Re-screen:** the uninjected re-screen of all 500 stars reproduces #542 exactly: `good_cadences` per star, and
  all 86 flagged runs with every flag. Mismatches: none.

## 1. Effective search volume per cell (`summary.json`; every event in `events.json.gz`)
- **Times** (stored separately):
  - observing time T_obs = 1,732.3 star-years (good cadences × 29.4 min, as in #542);
  - sampling time ΣT = 1,997.9 star-years. Every star's window is 1,459.49 d (first to last Q1–Q17 timestamp).
    Event centres were drawn uniformly over it, gaps included.
- **Effective search volume:** E = Σ T_i × (recoveries / 2), on the same window as the draw.
  - It is **not** 1,732.3 × the recovery rate. That would remove the gaps twice.
- **Uncertainty:** the MC 95 % interval is the Wilson interval of the window-weighted recovery rate. Every star has
  the same window, so Kish's n_eff = 1,000 per cell.

| cell | ρ | blend | t_E (d) | recovery rate | E (star-yr) | MC 95 % | covered (of 1,000) | per star (yr) ± SE | Γ95 (per star-yr) |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 0.01 | 1.0 | 0.1 | 0.797 | 1,592.4 | 1,540.3–1,639.8 | 899 | 3.1847 ± 0.0522 | 0.001945 |
| 1 | 0.01 | 1.0 | 0.3 | 0.789 | 1,576.4 | 1,523.7–1,624.6 | 937 | 3.1527 ± 0.0543 | 0.001966 |
| 2 | 0.01 | 1.0 | 1.0 | 0.811 | 1,620.3 | 1,569.5–1,666.4 | 971 | 3.2406 ± 0.0489 | 0.001909 |
| 3 | 0.01 | 1.0 | 3.0 | 0.782 | 1,562.4 | 1,509.1–1,611.3 | 990 | 3.1248 ± 0.0541 | 0.001985 |
| 4 | 0.01 | 1.0 | 10.0 | 0.679 | 1,356.6 | 1,297.5–1,412.9 | 991 | 2.7132 ± 0.0595 | 0.002309 |
| 5 | 0.1 | 1.0 | 0.1 | 0.762 | 1,522.4 | 1,467.7–1,573.1 | 872 | 3.0448 ± 0.0558 | 0.002041 |
| 6 | 0.1 | 1.0 | 0.3 | 0.783 | 1,564.4 | 1,511.2–1,613.2 | 908 | 3.1288 ± 0.0507 | 0.001982 |
| 7 | 0.1 | 1.0 | 1.0 | 0.723 | 1,444.5 | 1,387.5–1,498.1 | 924 | 2.889 ± 0.0576 | 0.002159 |
| 8 | 0.1 | 1.0 | 3.0 | 0.735 | 1,468.5 | 1,412.1–1,521.3 | 958 | 2.937 ± 0.0558 | 0.002121 |
| 9 | 0.1 | 1.0 | 10.0 | 0.639 | 1,276.7 | 1,216.2–1,335.0 | 973 | 2.5534 ± 0.062 | 0.002463 |
| 10 | 0.01 | 0.7 | 0.1 | 0.77 | 1,538.4 | 1,484.3–1,588.4 | 895 | 3.0768 ± 0.0537 | 0.002018 |
| 11 | 0.01 | 0.7 | 0.3 | 0.805 | 1,608.3 | 1,557.0–1,655.0 | 940 | 3.2167 ± 0.0529 | 0.001924 |
| 12 | 0.01 | 0.7 | 1.0 | 0.809 | 1,616.3 | 1,565.3–1,662.6 | 975 | 3.2327 ± 0.0512 | 0.001914 |
| 13 | 0.01 | 0.7 | 3.0 | 0.783 | 1,564.4 | 1,511.2–1,613.2 | 984 | 3.1288 ± 0.0541 | 0.001982 |
| 14 | 0.01 | 0.7 | 10.0 | 0.693 | 1,384.6 | 1,326.1–1,440.1 | 993 | 2.7691 ± 0.0569 | 0.002259 |
| 15 | 0.1 | 0.7 | 0.1 | 0.752 | 1,502.4 | 1,447.1–1,553.9 | 872 | 3.0049 ± 0.057 | 0.00207 |
| 16 | 0.1 | 0.7 | 0.3 | 0.775 | 1,548.4 | 1,494.6–1,597.9 | 916 | 3.0968 ± 0.0539 | 0.002004 |
| 17 | 0.1 | 0.7 | 1.0 | 0.734 | 1,466.5 | 1,410.1–1,519.3 | 937 | 2.933 ± 0.0566 | 0.002124 |
| 18 | 0.1 | 0.7 | 3.0 | 0.702 | 1,402.5 | 1,344.5–1,457.6 | 940 | 2.8051 ± 0.0602 | 0.002228 |
| 19 | 0.1 | 0.7 | 10.0 | 0.625 | 1,248.7 | 1,187.9–1,307.6 | 961 | 2.4974 ± 0.0598 | 0.002522 |

- **Column notes:**
  - "Covered" means ≥ 3 good cadences at least 1 mag dim in the noiseless model.
  - "Per star ± SE" is the cluster (star-level) mean and standard error, used for extrapolation in section 4.
- **Range:** E is 1,248.7–1,620.3 star-years, 0.721–0.935 of the 1,732.3 observing star-years (ratios in
  `posthoc.json`, post hoc). The weakest cells are t_E = 10 d (1,248.7–1,384.6). The ρ = 0.01 cells with
  t_E ≤ 3 d lie at 1,538.4–1,620.3, within each other's MC intervals.
- **No false recoveries:** a new passing run away from the injection occurred in 0 of 20,000 events.
- **Per-star efficiencies are not measured:** with 2 draws per star and cell, 139–245 stars per cell have one of
  two recovered (`stars_k`).
- **Descriptive subgroups only:**
  - the 5 stars outside the #542 noise proxy recover 4–9 of 10 per cell;
  - the 10 catalogued EB stars recover 7–15 of 20 per cell.
  - Neither group was dropped.
- **Context:** #529 measured 55–86 / 100 per cell on 25 other stars. That run drew centres between the first and
  last good cadence and used the overlap-only match, so the numbers are not directly comparable.

## 2. Where events are lost (frozen loss stages, `summary.json` `stages`; post hoc groups, `posthoc.json`)
**Losses over all 20 cells:** 5,052 of 20,000 events are lost. The groups below are a post hoc grouping of the
frozen stages:

| group (frozen stages) | events lost | per cell (of 1,000) |
|---|---|---|
| insufficient data in a test window (`edge:insufficient`, `shoulder:insufficient`) | 3,050 | 79–253 |
| coverage (`gap`, `quality`): the dim part falls on missing, flagged or left-out cadences | 681 | 1–86 |
| shape tests (`edge:slow`, `shoulder:no`) | 581 | 6–72 |
| short or grazing (`short`): < 3 cadences ≥ 1 mag dim even without gaps | 483 | 1–49 |
| not flagged (covered, but no matched run) | 140 | 0–49 |
| one-off | 117 | 0–58 |

**Post hoc diagnostics** (`posthoc.json`):
- **"Near a data edge":** the farthest window any test uses reaches the data start or end, a quarter boundary or a
  ≥ 0.5-d gap.
  - At t_E ≥ 3 d almost every event is near an edge, recovered or not. The diagnostic separates the groups only at
    t_E ≤ 1 d.
  - Per t_E (near an edge / events, `per_stage`):

| t_E (d) | `edge:insufficient` | `shoulder:insufficient` | recovered |
|---|---|---|---|
| 0.1 | 81 / 102 | 214 / 225 | 205 / 3,081 |
| 0.3 | 276 / 307 | 187 / 202 | 263 / 3,152 |
| 1 | 449 / 456 | 128 / 128 | 1,095 / 3,077 |
| 3 | 435 / 435 | 238 / 238 | 2,771 / 3,002 |
| 10 | 348 / 348 | 609 / 609 | 2,633 / 2,636 |

- **Insufficient losses:** at t_E ≤ 1 d, the "insufficient" losses lie near a gap, a quarter boundary or the data
  end much more often than the recovered events. At longer t_E the test windows are long enough that nearly every
  event meets one, so the diagnostic cannot tell them apart there.
  - By definition, "insufficient" means the edge or shoulder window has too few cadences in the run edge's quarter.
- **Catalogued EB stars and shape tests:** 82 of the 136 `edge:slow` and 25 of the 445 `shoulder:no` losses are
  on the 10 catalogued EB stars (`eb_events_per_stage`), i.e. losses to the stars' own variability.
  - Per star (`eb_stars`): KIC 3954798 recovers 0 of 40, KIC 2717141 1 of 40, the other eight 15–34 of 40.
- **By construction:** KIC 3954798 has 86 uninjected flagged runs (`check.json`), so the one-off test can never
  pass for it.
  - The cascade files its losses at the first failed test, which comes before the one-off test.
  - It is 2 of the 1,000 events in each cell.
- **Short:** 452 of 483 (`short_grazing_total`) are grazing, with u0 ≥ 2 − ρ (the source never lies wholly in the
  umbra).
  - By ρ: 403 of 408 at ρ = 0.1 and 49 of 75 at ρ = 0.01.
  - So ρ = 0.1 cells lose 26–49 per 1,000 here.
- **Normalisation break (one-off 115 of 117, not flagged 82 of 140; post hoc):** in these events the model is
  ≥ 1 mag dim on at least half the good cadences of some quarter. That quarter's median is then set by the event,
  which breaks the per-quarter normalisation.
  - The one-off losses are at t_E = 10 d (57–58 per 1,000 at blend 1.0).
  - In the look-at figure, the one-off example's normal cadences come out at relative flux ~10⁵ (read off the
    figure; not stored).
  - The other 58 not-flagged events include events whose dim run reaches the end of the data (the look-at example).
- **Summary (post hoc, from the counts above):**
  - The largest frozen loss stages are the two "insufficient" stages (3,050; post hoc grouping). At short t_E they track gaps and
    boundaries.
  - Direct coverage loss is 681.
  - About 1 in 7 shape-test losses are on two catalogued EBs (KIC 3954798 and 2717141, from `eb_stars` and
    `eb_events_per_stage`; the per-star stage counts are not stored).
  - The per-quarter normalisation is a by-construction limit for long, fully dark events.

## 3. What the zero result supports
- **Premise:** #542 found 0 survivors in these 500 stars (`expand/layers.json`).
- **Product-level conditional limit per cell:** Γ95 = 2.996 / E_low, the classical Poisson 95 % upper limit for
  zero events, at the low end of E's MC interval. Values: 1.909–2.522 × 10⁻³ per star-year (table; `summary.json`
  `limit_product_level`).
  - Units: W3 signatures with that cell's ρ, blend and t_E, with t0 uniform in time and u0 uniform on [0, 2), as
    they would appear in the Kepler SAP product of these 500 stars.
  - Naive comparison (post hoc, `posthoc.json` `ratios`): 2.996 / 1,732.3 = 1.729 × 10⁻³. It is **not used**: it
    ignores the losses, and the measured limits are higher.
- **Conditions:**
  - events independent and Poisson, at one constant rate over stars and time;
  - all events in that one cell;
  - the measured efficiency;
  - rare events. At the limit, λ per star is 0.0076–0.0101. The share of events that share their star with another
    event is about λ, so ≤ 1.0 %, and the share of event-hosting stars with ≥ 2 events is 0.38–0.50 %. So the
    one-off rule removes at most about 1 % of events at these rates.
  - Repeating or periodic events are outside the model.
- **Event distribution:** u0 uniform on [0, 2) is a choice. Recovery depends on u0 (`posthoc.json`
  `recovery_by_u0`, post hoc):
  - 0.743–0.790 for u0 < 1.9;
  - 0.345 for u0 in [1.9, 2) (0.02 at ρ = 0.1, where the source never lies wholly in the umbra).
  - A population weighted towards grazing events would have a smaller E.
- **The 20 limits are separate statements,** one per hypothetical single-cell population. They are not
  simultaneous, and they are not a limit on any mixture of cells.
- **What the zero does not support:**
  - a physical lens rate, wormhole number or abundance. No event-rate model mapping physical lenses to (ρ, blend,
    t_E, u0) is adopted here;
  - any statement about the sky-to-SAP response. Aperture photometry of a real event, crowding beyond the blend
    parameter, saturation and background are untested. So is whether the archive's quality flags would have masked
    real dim cadences: the impulsive-outlier bit is set from the data, and an injection cannot set it;
  - any claim for t_E outside 0.1–10 d, u0 ≥ 2 (no dimming), other light-curve shapes, repeating events, other
    stars or other surveys;
  - "W3 excluded". The result is a product-level sensitivity with stated gaps, not an exclusion.

## 4. Would more of the same data help? (`summary.json` `expansion`; `cost.json` `expansion_estimate`)
- **Assumptions:** new stars drawn from the same pool by the same hash rule (17,740 unused after #542), with the
  same distribution of noise, variability and gaps. The expected per-star gain is the cluster mean above.
- **Gain** (expected gains; a real batch scatters around them by roughly sd(e_i) √N, not stored):
  - **+500 stars:** +1,248.7 to +1,620.3 effective star-years per cell (95 % intervals of the expectation in
    `summary.json`). If the result stays zero, each per-cell limit halves (factor 0.5).
  - **+5,000 stars:** +12,487 to +16,203 effective star-years per cell. Each limit falls to 0.091 of today's value.
- **Cost per 500 (measured):**
  - download about 3.8 GB (3,791,889,000 B this round; 3,818,397,533 B in #542 with listings and Q0);
  - fetch 1,114–1,201 s;
  - screen 31 s;
  - injection 508.8 s wall on 4 cores;
  - cache 989 MB.
- **Cost for +5,000** (linear scaling of the per-500 values, not measured):
  - about 38.2 GB streamed;
  - about 3.3 h of fetching;
  - 84.8 min of injections;
  - 9.9 GB of cache if kept at once. The cache would need batch-wise deletion in a cloud session.
  - Note: the 38.2 GB uses #542's bytes, with directory listings and Q0. PLAN §8 said "minus Q0"; new stars need the
    listings, so the larger figure is kept.
- **What it buys:** under the assumptions above, the expected gain scales with the number of stars. It buys a
  smaller product-level number, not evidence.
  - Without an event-rate model, no probability of finding an event can be given.
  - Nothing here says fainter or more crowded stars would do better: the per-star efficiency was measured only on
    this FGK sample.

## 5. Plain answers
1. **What the 1,732 star-years are worth:** 1,249–1,620 effective star-years, depending on the W3 type.
   - Events of up to 3 days on a small source keep the most (about 1,540–1,620).
   - 10-day events keep the least (about 1,250–1,385).
2. **Where the losses are:**
   - **Largest: insufficient data around the event.** The edge or shoulder test finds too few cadences in its window
     (3,050 of 5,052 lost events; post hoc grouping). For events of up to a day, these losses sit near gaps and quarter boundaries
     much more often than the recovered events (post hoc).
   - **Then:** events whose dim part falls in a gap or on flagged cadences (681); shape-test failures (581, of
     which 107 are on the 10 catalogued EBs: the stars' own variability); grazing or too-short events (483).
   - **10-day dark events:** the per-quarter normalisation breaks (one-off 115 of 117, not flagged 82 of 140; post
     hoc).
   - **By construction:** one EB (KIC 3954798) can never pass the one-off test.
3. **What the zero supports:** per W3 type, and only under the stated assumptions, fewer than about 1.9–2.5 SAP-level
   signals per 1,000 star-years (95 %). That is a product-level sensitivity. It is not a physical rate, not a joint
   limit across types, and not an exclusion of wormholes or negative masses.
4. **Is another batch worth it:**
   - Another 500 stars would double the expected effective volume and halve each limit. That costs about 3.8 GB and
     about half an hour of machine time (fetch, screen and injections; the sum is not stored).
   - 5,000 stars would cut the limits about elevenfold.
   - Whether that is worth doing depends on a rate model that says which rate matters, and none is adopted here. So
     the expansion would sharpen a number, not test a prediction.
