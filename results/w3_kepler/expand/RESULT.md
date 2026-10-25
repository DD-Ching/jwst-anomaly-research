# W3 / Kepler fixed-rule expansion search: 500 new stars, 1,732 star-years; one source triggers the pre-screen (KIC 3954798, a catalogued 0.918-d eclipsing binary, not vetted); 0 survivors

- **Pre-registration:** `PLAN.md`, `scripts/w3_kepler_expand.py`, its tests and `sample_expand.json` were
  committed at 21:52:38 UTC (927a716), before any new flux. The first light-curve read was at 21:52:49
  (`source_usage.jsonl`, `cost.json`).
- **Rule:** the #529 rule, unchanged (`w3_kepler.screen_star(local=True)`, survivor = `pass_local`): 1-mag dim
  threshold, quality mask, per-quarter join, common mode, bounded fast edges, one-off, local shoulder.
  - No threshold change, no second method, no injection.
- **Labels:** light curves and flags are observed; counts are derived. This is an exploratory search. A pass
  would not be wormhole evidence, and a zero is not an exclusion.
- **Recovery:** the numbers come from #529's validation on 25 other stars (55–86 / 100 per cell; weakest 55,
  which does not show a population recovery above 50 %). They were **not re-measured on this batch**.
- **Inherited limits:**
  - some #529 controls were rejected by construction;
  - brightenings of ≥ 2 cadences in both shoulder windows (e.g. flares bracketing a dip) are untested.
- **Deviations:** none. Notes:
  - `frac_loo` equals `w3_kepler.frac_for` on each star's own cadences, which are the only ones the screen looks
    up. The test covers stars with shared cadences, and one star missing 200 cadences (where no other star is dim,
    so the own-cadence filter is not stressed).
  - `vet_script.py.txt` was committed at 22:04, after the fetch began. It is not part of the frozen plan, and
    it never ran because there were no survivors.
  - `flagged_source.json` / `.png` are a post hoc description of the one flagged source. They are not
    vetting, because the source is not a survivor.

## 0. Previous counts (correction of a chat figure)
- **The figure:** my chat summary after #529 said "425 stars". That figure is not in the repository and
  was wrong.
- **Deduplicated counts** (`sample_expand.json`, `previous_rounds`):
  - development 20;
  - validation 45;
  - search 260;
  - ordinary-event controls 46 (EB stars);
  - interface check 1.
- **Total:** 372 distinct stars, 325 of them FGK pool stars. All 372 were excluded from this batch.

## 1. Layers (`layers.json`; per star `screen_expand.json.gz`; every run `runs_expand.json.gz`)

| layer | count |
|---|---|
| selected (frozen list) | 500 |
| read (0 download failures; 8,745 files, 3.74 GB) | 500 |
| screened | 500 |
| inside the PLAN §3 applicability range (17 quarters and median σ_rel ≤ 3.96 × 10⁻⁴; a proxy, `applicability.json`) | 495 |
| effective observing time (good cadences × 29.4 min) | 632,717 d = **1,732 star-years** |
| 2a. stars with any good cadence ≥ 1 mag below the quarter median | **1** (KIC 3954798; 2,848 cadences) |
| 2b. flagged deep runs (≥ 3 consecutive dim cadences, bracketed) | 86 runs, 1 source |
| 3. after common mode | 86 runs, 1 source |
| 3. after bounded fast edges (all 86 classed "slow") | **0** |
| 3. after one-off / local shoulder | 0 / 0 (each test on its own: one-off 0 / 86, local shoulder 0 / 86) |
| 4. survivors | **0**; nothing vetted, nothing left unvetted |

- **The 5 stars outside the applicability range:** they have median σ_rel 5.7 × 10⁻⁴ to 2.0 × 10⁻³, and all 17
  quarters. Each was screened and has no flagged run.
- **What the range is:** a proxy. It is twice the largest median σ of the #500 stars. #529's 25 validation
  stars, where recovery was measured, have no stored σ. 12 of the 495 applicable stars lie above the #500
  maximum itself (1.98 × 10⁻⁴; `applicability.json`).
  - Applicable share: 495 / 500 = 99 %, above the 90 % stop rule. The #529 recovery numbers are cited for this
    batch only with these caveats.
- **Other 499 stars:** none has any cadence ≥ 1 mag dim. Their lowest relative flux is 0.644.
- **Plain statement:** apart from one eclipsing binary, this batch did not trigger the pre-screen. The zero
  survivors is not a result of the shoulder test.

## 2. The one flagged source (post hoc description; `flagged_source.json`, `flagged_source.png`, looked at)
- **Identification:** KIC 3954798 is in the Kirk et al. 2016 EB catalogue, with P = 0.91848 d and
  morphology 0.58.
- **Light curve:**
  - deep primary eclipses every 0.918 d, down to relative flux 0.217 (1.66 mag);
  - a shallower secondary at phase 0.5;
  - 2,848 good cadences ≥ 1 mag dim.
- **Why it is rejected:** its 86 flagged runs all fail the edge test first ("slow"), and they also fail the
  one-off test (86 runs) and the local shoulder test (0 / 86).
- **What "slow" means here** (`flagged_source.json`):
  - All 86 runs are shorter than 1 d (85 span 0.041 d, one 0.94 d), so the 1-d span floor sets the edge window
    for every run: intermediate cadences are searched within one span (1 d) of each edge, and "fast" needs
    U ≤ 0.25 d.
  - The farthest significant intermediate cadence (L) is at the **neighbouring** primary eclipse, and so is the
    nearest normal cadence beyond it (U): median L 0.94 d, median U 0.96 d. That is the previous eclipse on the pre side and the next one on the post side
    (79 / 86 pre sides and 83 / 86 post sides are slow).
  - The eclipse's own edges are fast: on the fast sides, median L is 0.020 d and median U 0.041 d.
  - So the edge rejection here depends on P = 0.918 d being below the 1-d floor. With L ≈ P + 0.02 d and the
    eclipse's intermediate cadences starting about 0.07 d before its minimum, the neighbour would leave the
    window at roughly P ≳ 1.1 d. This is an estimate, not tested. An EB like this with a longer period would not
    be stopped by the edge test; it would rest on the one-off and shoulder tests.
- **Run timing:**
  - Every run starts within 0.025 in phase of the primary minimum.
  - One run (Q11, n = 3) spans 0.94 d, because two of its dim cadences are separated by a 0.919-d data gap.
    Its midpoint falls about half a period after an eclipse; that is the 0.4903 stored as
    `run_phase_at_kirk_period_abs_max`. Every other run's midpoint lies within 0.023 of the primary minimum.
  - Cadences within ±0.02 in phase of the secondary (phase 0.5) have a median relative flux of 0.852.
- **Status:** deep periodic dimming at the Kirk 2016 period of a catalogued eclipsing binary. Not vetted,
  because it is not a survivor. Not a W3 candidate.

## 3. What this round adds
- **Search coverage:**
  - 500 new stars, 1,732 star-years of 29.4-min data;
  - the Kepler W3 search total is now 760 search stars (260 earlier + 500; a sum, not stored). The 60 #500
    search stars were screened under the older rule, but had 0 flagged runs.
- **Measurement capability:** nothing new (no injection this round).
- **Physical evidence:** none. No W3-like event, and the one deep dimming is an ordinary eclipsing binary.
- **Not done:** no rate limit. The weakest cell's recovery is not shown to be above 50 %, and recovery was not
  re-measured on this batch.
