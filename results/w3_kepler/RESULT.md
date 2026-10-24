# W3 / Kepler dense-cadence pilot: validation passes (78–91 / 100 per cell); the 60-star search has 0 flags

Validation passes:
- **W3 recovery:** all 20 pre-declared cells recover 78–91 / 100.
- **Edge bounds:** among recovered injections, the median data bound on the edges is 0.041 d and
  the 90th percentile is 0.22–0.27 d. 164 of the 2,000 val W3 injections were lost because an edge
  was not bounded.
- **Controls:** no ordinary control passes. Each control family is rejected by a property built into
  its model (§2), so the 0 / 100 rates confirm that the chain behaves as designed. They do not
  measure contamination by real ordinary events.
- **Search:** the 60 reserved stars have no flagged run, so there are no survivors to vet.
- **Review:** the review found that 4 of the 100 real SAP light curves already show shoulder-like
  brightenings within a quarter, which weakens the W3 shoulder test (§2).

- **Pre-registration:** `PLAN.md`, `scripts/w3_kepler.py` and `tests/test_w3_kepler.py` were
  committed at 07:56:30 UTC (314ebc5). The first target file was read at 07:56:41
  (`source_usage.jsonl`, `cost.json`). Before that, the ledger holds only the NExScI query and the
  non-target interface star.
- **Hypothesis:** the existing finite-source negative-mass point lens (W3) only.
- **Labels:** injections are simulated; light curves and their flags are observed; outcomes are
  derived.
- **Scope:** efficiencies apply to processing from the DR25 SAP product onward, not to the official
  pixels-to-light-curve chain. SAP and PDCSAP are one observation; PDCSAP was not needed (no real
  flag).
- **Deviations from the plan:** none in method or thresholds.
  - Ambiguities and errata found at review are listed in `PLAN.md` "Errata after review".
  - One implementation fix: the import-time `w3_ztf.EDGE_RULE` setting leaked into
    `tests/test_w3_ztf.py`. Removing it changes no outcome; see `rerun_compare.json`.
- **Post hoc:** everything in §4 and the figure is post hoc diagnostics
  (`posthoc_script.py.txt` → `posthoc.json`, `cases.png`; `posthoc2_script.py.txt` →
  `posthoc2.json`, after the review). None of it gates.

## 1. Data and sample
- **Light curves:** 1,750 `_llc.fits` files, all `DATA_REL` 25 (`files_<group>.json`). Every one of
  the 100 stars has all 17 quarters Q1–Q17.
- **Sample:** the plan's metadata rule (`sample.json`) gives 18,570 eligible FGK dwarfs at
  Kp 12–14. Ranks 0–19 are dev, 20–39 val and 40–99 search.
- **Read order:** dev, then val after the dev phase, then search after the gate
  (`source_usage.jsonl`).
- **Download:** 766 MB of the 2 GB cap (`net_ledger.jsonl`).

## 2. Validation (fixed method, unread val stars; `injections_val.json`, `gate.json`)

| | val (gate) | dev |
|---|---|---|
| W3 recovered / 100, range over 20 cells | **78–91** (lowest Wilson 95 % lower bound 0.689) | 72–94 |
| W3 with edges not bounded ("insufficient", in the denominator) | 1–16 per cell | 2–15 per cell |
| eclipse (periodic) passed / flagged | **0** / 98 | 0 / 96 |
| single eclipse passed / flagged | **0** / 95 | 0 / 99 |
| fade passed / flagged | **0** / 75 | 0 / 78 |
| PSPL (brightening) passed / flagged | **0** / 0 | 0 / 0 |
| common mode passed / flagged | **0** / 100 (by construction) | 0 / 100 |
| quarter step passed / flagged | **0** / 43 | 0 / 44 |
| uninjected real light curves: flagged runs | **0** (20 stars) | 0 (20 stars) |

- **Gate:** every cell ≥ 0.50 and every control family ≤ 0.10. It **passes**.
  - These are pilot thresholds on simulated signals. They are not a discovery significance and
    not a sky false-positive rate.
- **Why each control family is rejected** (all by its model):
  - **PSPL** only brightens, so it is never flagged.
  - **Common mode** fails by construction.
  - **Periodic eclipses** (P ≤ 50 d over ~4 yr) always repeat, so they fail the one-off test.
  - **Quarter steps** end at a quarter roll, so their second edge is never bounded (43 / 43 flagged
    are "insufficient").
  - **Fades** (val): 27 of 75 flagged are "insufficient". The other 48 fail another test; the edge
    class of control rows is not stored (`posthoc2.json`).
  - **Single eclipses** are box dips without a brightening. Only the shoulder test stops them: 91 / 100
    pass without it in val (85 / 100 in dev).
  - **Consequence:** the control half of the gate could fail only through an implementation error. It
    is a check of the chain, not six independent tests against real ordinary events.
- **The W3 shoulder test:**
  - It costs W3 recovery 0–3 per cell.
  - In the model, the exposure-averaged brightening just outside the umbra is ≥ 0.39 above baseline
    in every cell (`model_check.json`, `model_check_script.py.txt`; window and u0 range noted there).
- **Weaknesses of the shoulder test on real data, found at review:**
  - With no injection, 4 of the 100 stars already have good cadences meeting the shoulder criterion
    (rel ≥ 1.10 at ≥ 5σ). These are SAP brightenings within a quarter up to rel 1.17: 1,797 cadences,
    1 dev star and 3 search stars, in Q1, Q4 and Q8 (`posthoc2.json`). A deep single dip next to
    such a stretch would pass.
  - The test needs only one cadence per side.
  - A deep event that covers much of a quarter lowers that quarter's median and can create a
    "shoulder". The reviewer saw this in 5 of 75 flagged val fades; that count is not stored.
  - Ordinary dips longer than 2 d (single eclipses run 0.1–2 d) were not tested.
  - The injected error for F > 1 is err √F, while the scatter is err F, so shoulder S/N is overstated
    by up to √F. This hardly matters here, because the median σ_rel per star is 1.3 × 10⁻⁴ (at most
    2.0 × 10⁻⁴; `posthoc2.json`), far below the 0.10 threshold.

## 3. Search (60 reserved stars; `search.json`)
- **Result:** 0 flagged runs, so 0 survivors.
- **Vetting counts:** 0 checked real dimmings, 0 artefacts and 0 undetermined. No survivor was vetted
  because none existed.
- **Search stars with shoulder-like stretches:** three of them have such stretches (§2). They had no
  dip ≥ 1 mag, so the weakness did not act in this search.

## 4. Post hoc diagnostics (`posthoc.json`, `posthoc2.json`)
- **Edge bounds of recovered W3 injections:**
  - The larger side's U has a median of 0.0409 d (2 cadences) in both phases.
  - The 90th percentile is 0.225 d (dev) and 0.266 d (val).
  - 51 % of recovered injections are bounded within 2 cadences.
  - These describe recovered injections only, which had to be bounded to pass. No test isolates
    the cadence as the cause.
- **ZTF comparison:**
  - Under the v23 "bounds" rule, without Kepler's per-quarter restriction, ZTF DR24 recovered 43–114
    of 200 per cell and failed its gate (`../w3_ztf/v24/gate.json`, cells_bounds_recovered).
  - Its losses relative to the earlier v13 rule were all unbounded edges
    (`bounds_losses_all_insufficient`). It also lost events to slow edges, one-off or band failures
    and missed flags (`../w3_ztf/v24/RESULT.md` §1).
  - The grids overlap only at t_E = 10 d, and the stars, survey and join differ. This compares
    capability; it is not a paired test.
- **Why 324 val W3 injections were not recovered** (labels describe the stored outcome only):

  | outcome | count | numbers behind it (`posthoc2.json`) |
  |---|---|---|
  | flagged, edge not bounded | 164 | — |
  | not flagged, umbra < 50 % covered by good cadences | 80 | — |
  | not flagged although ≥ 50 % covered | 54 | 46 with u0 ≥ 1.9; 8 in Q1 / Q17; 6 within 1 d of the data start or end |
  | flagged, shoulder missing | 12 | — |
  | flagged, one-off fails | 10 | all t_E = 10 d; 8 in Q1 / Q17 |
  | flagged, slow edge | 4 | — |

  - Before the val phase, two dev rows were replayed: an event at the end of the data and a long
    event in Q1 / Q17. Neither was a code error. This was seen in the session and is not stored.
    The review later replayed 40 stored val rows exactly.
- **Real light curves:** none of the 100 stars has any good cadence ≥ 1 mag below its quarter median
  (6,191,683 good cadences; the lowest relative flux is 0.802).
  - So the real-data check met no deep ordinary dimming.
- **Figure (`cases.png`, looked at):** replayed injections on val carriers. The W3 shoulders and the
  flat umbra are visible; a single eclipse has neither shoulder; one quarter-step edge falls at a
  quarter roll.

## 5. What this does and does not show
- **Shown:** for t_E 0.1–10 d, ρ ≤ 0.1 and blend ≥ 0.7, W3 events injected into real Kepler SAP light
  curves of bright FGK dwarfs are recovered in 78–91 % of cases. Recovered events have edges bounded
  with a median of 0.041 d.
- **Not shown:**
  - that the chain rejects real ordinary events. The controls are rejected by design; the sample
    had no deep dimming; real in-quarter brightenings can supply false shoulders.
  - efficiency of the official pixel-level chain;
  - t_E ≥ 30 d (excluded in advance);
  - heavier blending;
  - a second band (the achromatic test is "not tested", not passed);
  - any rate limit (not planned; 60 stars is a pilot).
- **Physical judgement:** unchanged. No W3-like event was found. The result is a measured injection
  capability on real light curves, not evidence for or against negative-mass lenses.
