# W3 / Kepler hard tests: simulated hard cases fool the #500 shoulder rule, real deep EBs fool neither rule; local rule adopted at a cost of 177 / 2,000 W3 injections; the 200-star search has 0 flags

**What happened:**
- **Old rule fooled:** the #500 rule is fooled by deliberately placed, simulated hard cases. These
  are a single deep eclipse inside a real or synthetic in-quarter brightening, and one unflagged
  bright cadence on each side of a dip.
- **Local rule:** the frozen local shoulder rule passes none of them. Most of these 0 / 100 results
  follow from how the families were built (§2), so they are not a broad test of the new rule.
- **Real deep EBs:** none of the 2,892 eclipse windows in 10 real deep EBs has a shoulder under
  either rule. 2,826 of them have bounded edges (post hoc).
- **Cost:** on the unread val stars, W3 recovery falls from 1,671 to 1,494 of 2,000 injections. That
  is 177 injections, 8.9 points, or 10.6 % of the old rule's recoveries.
- **Search:** validation passed, and the 200-star search found no flagged run.

**Provenance:**
- **Pre-registration:** `PLAN.md`, `scripts/w3_kepler_hard.py`, the local rule in
  `scripts/w3_kepler.py`, the tests and `sample_hard.json` were committed at 10:53:04 UTC (db00218),
  before any new flux or comparison.
  - A wording-only fix of the §3.1 timing followed 15 s later, before any data.
  - The model check behind §3.1 is now stored (`local_model_check_script.py.txt` →
    `local_model_check.json`, identical numbers).
- **Labels:** injections are simulated; light curves and their flags are observed; outcomes are
  derived.
  - Hard cases give **conditional** false passes, not sky rates.
  - Repeated injections on one star are one source; source counts are given.
  - Numbers taken from `posthoc.json` are post hoc and marked "(post hoc)".
- **Deviations:** none in method, thresholds or seeds.
  - **Phase A run twice:** the second run, at 11:09, came after the val stars were read (11:03).
    It only moved the per-window runs into `phaseA_runs.json.gz` (the plain file exceeded 1 MB)
    and rewrote `phaseA.json` and `decision.json` with identical contents.
  - **Ledger:** tracked as a gzip snapshot.
  - **Usage log:** the `search` entry in `source_usage.jsonl` lists old_bg, but the search's
    common-mode comparison used search + val, as planned.

## 1. Phase A: hard tests on seen backgrounds and real EBs (`phaseA.json`)
Both rules were evaluated on the same injections, 100 per family on the 23 old_bg stars.

| family (category) | old rule passed / 100 (Wilson 95 %) | local rule | sources with a pass (old) |
|---|---|---|---|
| real brightening + single eclipse (real background + deliberately placed event; 4 seen sources × 25) | **56** (0.46–0.65) | 0 | 4 / 4 (per source 10, 23, 22, 1 of 25) |
| single-point anomalies (one +20 % cadence per side) | **72** (0.63–0.80) | 0 | 22 / 23 |
| synthetic 10–20 % brightening + eclipse | **86** (0.78–0.92) | 0 | 23 / 23 |
| long dips 0.1–40 d | 0 | 0 | — |
| quarter-boundary dips | 0 | 0 | — |
| the six #500 controls | 0 each | 0 each | — |

- **Real deep EBs** (real ordinary events; `eb_deep`, the 14 stars with a DR25 TCE depth ≥ 0.6):
  - 10 stars have flagged windows: 2,892 windows.
  - Of these, 2,826 pass common mode with both edges "fast" (post hoc).
  - Under the old rule, no window has a shoulder on either side. Under the local rule, no window
    has an "ok" side; 176 windows have an unevaluated ("insufficient") side
    (`phaseA_runs.json.gz`; 137 of the 2,826 bounded windows are local-insufficient, post hoc).
  - So real deep eclipses do not produce shoulders here. This does not favour either rule.
  - The full chain passes 0 stars.
  - One window (KIC 5262664) passes the one-off test (post hoc). It is 1,543 cadences (31.5 d,
    Q5) long, while the TCE duration is 15 h, so it is probably not an eclipse; it is unverified
    and possibly instrumental. Its edges are not bounded and it has no shoulder.
- **Real single events** (Kirk et al. 2016 table 9, 32 stars): none is flagged.
  - Their deepest observed relative flux is 0.607, i.e. 0.54 mag (post hoc, `eb_single_min_rel`),
    below the 1-mag threshold.
  - This was knowable before the run: the table's Depth column, read as flux at minimum, has a
    minimum of 0.60. These events therefore did not test the shoulder rule.
- **Decision:** the old rule fails 3 families, so the local rule is used (`decision.json`).
- **W3 on seen backgrounds:** the old rule recovers 77–90 / 100 per cell, the local rule 66–84.

## 2. What the local rule's 0 / 100 does and does not show
Most local-rule rejections are guaranteed by the family definitions:
- **single_point** injects exactly one bright cadence per side; the rule needs ≥ 2.
- **synth_bright_eclipse** has a 2–10 d plateau and an eclipse of ≤ 2 d, so the baseline window
  lies inside the brightening on at least one side.
- **long_dip and quarter_edge** contain no brightening.
- **real_bright_eclipse** uses the 4 already seen brightenings that motivated the rule (17–931
  shoulder-like cadences per star, `../posthoc2.json`). The star with the shortest stretch
  (3834001, 17 cadences) gave the old rule 1 pass in 25.

The case the local rule would pass is **≥ 2 cadences that are bright relative to the local
baseline, confined to the shoulder window on both sides of a dip**. One example is two flares
bracketing a dip. No family tested this. The gate therefore has little power against it, and it
remains open (§6).

## 3. Validation: the local rule on 25 unread stars (`validation.json`, `gate.json`)
- **W3:** all 20 cells recover **55–86 / 100**; the total is 1,494 / 2,000. Each cell is ≥ 0.50, so
  that part of the gate passes.
  - One cell's Wilson lower bound is below 0.5: ρ 0.1, blend 0.7, t_E 10 d at 55 / 100
    (0.45–0.64).
- **Controls:** all 10 families pass **0 / 100** on val, from 24–25 sources each (Wilson upper
  bound 0.037). Most of these are by construction (§2).
  - The required conditional controls also hold: real brightening 0 / 100 on 4 sources; real EB
    windows 0 / 2,892 and stars 0 / 10 (≥ 5 required).
- **Gate:** passes.
- **The old rule on the same val injections** (not gated): W3 76–91 per cell, 1,671 / 2,000;
  single-point 66 / 100 and synthetic brightening 84 / 100 pass.
- **Uninjected val stars:** 0 flagged runs.

## 4. Cost of the change (paired, same injections; post hoc, `posthoc.json`)

| | val (unread) | phase A (seen) |
|---|---|---|
| W3 recovered, old → local | 1,671 → **1,494** / 2,000 (−177 injections, −10.6 % of old recoveries; gained 0) | 1,656 → 1,491 (−165) |
| lost per cell | 2–21 | 1–20 |
| lost at t_E = 10 d (4 cells) | 16–21 per cell, 72 of 177 in total | 13–20 per cell |
| lost where the local test was "insufficient" | 141 of 177 | 126 of 165 |
| false passes removed (hard families) | single-point 66 → 0, synthetic brightening 84 → 0 | real brightening 56 → 0, single-point 72 → 0, synthetic 86 → 0 |

- **"Insufficient":** on at least one side, the baseline window held fewer than 20 good cadences in
  the same quarter. The other 36 losses fail the locally-bright criterion.
- **Not counted:** where these windows end (quarter edge or gap inside a quarter) was not counted.

## 5. Search: 200 new stars (`search.json`)
- **Result:** 200 stars, 0 missing, **0 flagged runs**, so 0 survivors.
- **Vetting counts:** 0 checked real dimmings, 0 artefacts and 0 undetermined. No vetting was
  needed, and the pixel-file reserve was unused.
- **Rules not exercised:** with no flagged run, neither the old nor the local shoulder rule acted
  on the search stars.

## 6. Post hoc observations and what is not covered
- **Depth (post hoc):** none of the 225 new stars (val + search) has any good cadence ≥ 1 mag below
  its quarter median. The lowest relative flux is 0.546.
- **Real brightenings (post hoc):** 9 of the 225 have real shoulder-like cadences (rel ≥ 1.10 at
  ≥ 5σ), similar to 4 / 100 in #500.
- **Figure (`cases.png`, looked at):**
  - the old rule's false passes: an eclipse inside a real or a synthetic brightening, and single
    bright cadences;
  - a W3 event lost by the local rule (t_E = 10 d);
  - one W3 event that both rules keep.
- **Not covered:**
  - ≥ 2-cadence brightenings confined to both shoulder windows (e.g. flares on both sides of a
    dip);
  - spotted-star dips, heartbeat stars, dippers, microlensing of blends, and pixel-level
    artefacts (rolling band, SPSD);
  - real events shallower than 1 mag;
  - real deep single eclipses (none in the sample).
- **Single band:** the achromatic test is "not tested". Efficiencies apply from the SAP product
  onward.
- **What a pass means:** it does not exclude every ordinary explanation. It shows that the tested
  families do not pass, under the stated conditions.
