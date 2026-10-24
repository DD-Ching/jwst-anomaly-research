# W3 / Kepler dense-cadence pilot: validation passes; the 60-star search has 0 flags

The 30-min Kepler cadence bounds both W3 edges: the median data bound is 0.041 d. All 20
pre-declared cells recover 78–91 / 100 in validation. No ordinary control passes (0 / 100 in each of
6 families). The search group (60 stars) has no flagged run, so there are no survivors to vet.

- **Pre-registration:** `PLAN.md`, `scripts/w3_kepler.py` and `tests/test_w3_kepler.py` were
  committed at 07:56:30 UTC (314ebc5). The first target file was read at 07:56:41
  (`source_usage.jsonl`, `cost.json`).
- **Hypothesis:** the existing finite-source negative-mass point lens (W3) only.
- **Labels:** injections are simulated; light curves and their flags are observed; outcomes are
  derived.
- **Scope:** efficiencies apply to processing from the DR25 SAP product onward, not to the official
  pixels-to-light-curve chain. SAP and PDCSAP are one observation; PDCSAP was not needed (no real
  flag).
- **Deviations from the plan:** none. Everything in §4 and the figure is post hoc diagnostics
  (`posthoc_script.py.txt` → `posthoc.json`, `cases.png`), not gating.

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
- **What the W3 shoulder test does:**
  - Without it (secondary), single eclipses pass 91 / 100 in val and 85 / 100 in dev.
  - With it, they pass 0. The test costs W3 recovery 0–3 per cell.
  - At Kepler precision, the brightening of the W3 lens just outside the umbra (≥ 0.39 above
    baseline after the exposure average in every cell; `model_check.json`) is what separates W3
    from a single box-shaped eclipse. The edge timing alone does not.
- **Quarter steps and fades:** these end as "insufficient" or slow, never as passes. The quarter
  rule makes an edge across a quarter roll unbounded.

## 3. Search (60 reserved stars; `search.json`)
- **Result:** 0 flagged runs, so 0 survivors.
- **Vetting counts:** 0 checked real dimmings, 0 artefacts and 0 undetermined. No survivor was
  vetted because none existed.

## 4. Post hoc diagnostics (`posthoc.json`)
- **Edge bounds of recovered W3 injections:** the larger side's U has a median of 0.0409 d (2
  cadences) in both phases. The 90th percentile is 0.225 d (dev) and 0.266 d (val). 51 % of
  recovered injections are bounded within 2 cadences.
- **ZTF comparison:** under the same v23 "bounds" rule, ZTF DR24 recovered 43–114 of 200 per cell
  and failed its gate (`../w3_ztf/v24/gate.json`, cells_bounds_recovered). Every loss there was an
  unbounded edge.
  - The grids overlap only at t_E = 10 d, and the stars, survey and join differ. This is a
    comparison of capability, not a paired test.
- **Why 324 val W3 injections were not recovered:**

  | reason | count |
  |---|---|
  | edge not bounded (gap or quarter boundary) | 164 |
  | umbra < 50 % covered by good cadences | 80 |
  | not flagged although covered (grazing u0, data end, or quarter join) | 54 |
  | shoulder missing | 12 |
  | one-off fails (umbra split by the quarter join) | 10 |
  | slow edge | 4 |

  - These fall under the limits stated in the plan (gaps, short Q1 / Q17, grazing events). The dev
    rows were checked by replay before the val phase; the cases seen were an event at the end of
    the data and long events in Q1 / Q17, not code errors.
- **Real light curves:** none of the 100 stars has any good cadence ≥ 1 mag below its quarter median
  (6,191,683 good cadences; the lowest relative flux is 0.802).
  - So the real-data check met no deep ordinary dimming. The control false-pass rates rest on
    simulations only.
- **Figure (`cases.png`, looked at):** replayed injections on val carriers. The W3 shoulders and the
  flat umbra are visible; a single eclipse has neither shoulder; one quarter-step edge falls at a
  quarter roll.

## 5. What this does and does not show
- **Shown:** for t_E 0.1–10 d, ρ ≤ 0.1 and blend ≥ 0.7, W3 events injected into real Kepler SAP light
  curves of bright FGK dwarfs are recovered in 78–91 % of cases with edges bounded to about an hour,
  and the tested ordinary shapes do not pass.
- **Not shown:**
  - efficiency of the official pixel-level chain;
  - t_E ≥ 30 d (excluded in advance);
  - heavier blending;
  - a second band (the achromatic test is "not tested", not passed);
  - the real-sky rate of ordinary mimics (the sample had no deep dimmings);
  - any rate limit (not planned; 60 stars is a pilot).
- **Physical judgement:** unchanged. No W3-like event was found. The new result is a measured
  capability, not evidence for or against negative-mass lenses.
