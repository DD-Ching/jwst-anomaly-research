# W3 / Kepler real-background false-pass test and small search: frozen plan

Owner authorization 2026-10-25, following #500. Start 10:42 UTC. Limits: 4 h (to 14:42 UTC), 3 GB
through `scripts/net_ledger.py` (`net_ledger.jsonl` here), 4 cores.

**Committed before any comparison or new flux:** this plan, `scripts/w3_kepler_hard.py`, the local
rule in `scripts/w3_kepler.py` (`shoulder_local`, `screen_star(local=True)`; the default path is
unchanged), the tests and `sample_hard.json`. Before this commit, the ledger holds only VizieR /
NExScI / Villanova metadata requests.

**Hypothesis:** the existing finite-source negative-mass point lens (W3) only.

**Unchanged from #500:** the W3 model, the 20 cells, the exposure average (21 sub-samples per
29.4-min cadence), the quality mask, the quarter join and the run finding. The #500 seeds are not
reused. The paused lines (ZTF, ATLAS, J2244, X4) are untouched.

## 0. Orientation (10:42–10:50 UTC)
- **Repository:** `main` at b397de4. No open PR or branch claims Kepler. #500 is merged
  (`results/w3_kepler/`, D-107).
- **Seen data:** all 100 #500 stars and the interface star are seen. They are used here only as
  real backgrounds for hard tests, never as validation or search stars.

## 1. Groups (`sample_hard.json`; metadata only)
- **old_bg (23, seen):** the 20 #500 dev stars, plus the 3 search stars with real in-quarter
  brightenings (2854411, 2856758, 3834001; 9137992 is a dev star; `../posthoc2.json`).
- **eb_deep (14):** every star with a DR25 TCE of depth ≥ 600,000 ppm (NExScI
  `q1_q17_dr25_tce`; `inputs/tce_depth_ge_0.4.csv`). These are real deep eclipses.
- **eb_single (32):** every system in Kirk et al. 2016 table 9, "no repeating events (long)"
  (VizieR J/AJ/151/68/table9; `inputs/kirk2016_table9.csv`). These are real single events.
- **val (25) and search (200):** the #500 pool (DR25 FGK dwarfs at Kp 12–14, all 17 quarters) minus
  the #500 101 and every EB star, sorted by sha256(`w3-kepler-hard-2026-10-25:<kepid>`). The first
  25 are val and the next 200 search. The groups are disjoint (asserted in `sample`).
- **Read order:**
  - old_bg and the EB stars first, for phase A;
  - val only after `decision.json` exists;
  - search only if `gate.json` passes.
  - `source_usage.jsonl` logs each read. A star whose download fails is recorded in
    `files_<group>.json`, stays in the counts as "missing", and is not replaced.
- **Budget (≈ 7.7 MB per star):** about 2.1 GB of light curves. 0.3 GB is reserved for vetting
  (pixel files).

## 2. Hard tests (fixed now)
Every injection uses the #500 rule of F × raw SAP before the quality mask, join and screen, with
added noise. Every model control is exposure-averaged.

- **Hit:** a flagged run that is not already flagged in the uninjected carrier.
- **Pass:** a hit passes the rule under test.
- **Seeds:** `default_rng([seed, j])`. Phase A uses W3 5100 + c, family i 5200 + i, and real
  brightening star i 5290 + i. Validation uses W3 5300 + c and family i 5400 + i. The family index i
  follows `FAMILIES` = the six #500 controls, then `HARD`.

**Families, with their categories:**

| family | category | definition |
|---|---|---|
| real_bright_eclipse | real background + simulated event (deliberately placed) | On each of the 4 stars with real in-quarter brightenings, 25 single box eclipses (1.5 mag). Centre: a uniformly drawn real "shoulder-like" cadence (rel ≥ 1.10 at ≥ 5σ, uninjected) ± U(−0.5, 0.5) d. Duration log-uniform 0.1–10 d. 4 sources; results reported per source. |
| long_dip | real background + simulated event | Box 1.5 mag, duration log-uniform 0.1–40 d (covers the W3 umbra range up to 4 t_E = 40 d), start U(good range). It can lower the quarter median. |
| single_point | real background + simulated event | Box 1.5 mag, duration U(0.1, 2) d, plus one unflagged +20 % cadence drawn uniformly inside each side window (0.25 × max(duration, 1 d)). |
| quarter_edge | real background + simulated event | Box 1.5 mag, duration U(0.1, 2) d, at a uniformly drawn boundary between consecutive quarters with data. Three variants with equal probability: ends 0–5 cadences before a quarter end; starts 0–5 cadences after a quarter start; spans the gap (0.05–1 d on each side). |
| synth_bright_eclipse | pure model brightening + model eclipse, on real background | A trapezoid brightening of amplitude U(0.10, 0.20), ramps U(0.5, 3) d and plateau U(2, 10) d, with a box eclipse (1.5 mag, duration log-uniform 0.1–2 d) inside the plateau. |
| the six #500 controls | pure model on real background | eclipse, single_eclipse, fade, pspl, common, quarter_step (`w3_kepler.control_flux`, unchanged). |
| eb_deep, eb_single | **real ordinary events** | Uninjected screen of the real light curves. For each flagged window: the single-window outcome without the one-off test (common mode, bounded fast edges, shoulder), and the full chain. Reported per window and per star. |

- **Counts:** n = 100 per family and phase, except real_bright_eclipse (4 × 25, phase A only, on
  seen data). Carriers are drawn from the phase's group.
- **Interpretation:** deliberately chosen hard cases give conditional false passes, not sky false
  positive rates. Repeated injections on one star are one source; source counts are reported.

## 3. The only allowed method change: one local shoulder rule (frozen now, used only if needed)
`shoulder_local`, with parameters chosen from the model only (§3.1) and never from injection truth.
- **Per side, in the quarter of that run edge:**
  - shoulder window S = (0, 0.25 × span] from the edge;
  - baseline window B = (0.25 × span, 0.25 × span + max(1 d, span)] beyond S;
  - span = max(run length, 1 d), as in #500.
- **Baseline and noise:** b = median of rel in B; s_B = 1.4826 × MAD of B.
- **Locally bright cadence in S:** rel − b ≥ max(0.10 b, 5 max(σ, s_B)).
- **Side status:**
  - "ok" if ≥ 2 cadences are locally bright;
  - "insufficient" if B has < 20 good cadences or S has none;
  - otherwise "no".
- **Run:** passes the shoulder test only if both sides are "ok". The full rule is common mode +
  bounded fast edges (the #500 per-quarter rule) + one-off + local shoulder.

### 3.1 Model check (no data, computed 10:50–11:00 UTC; not stored)
- **Method:** noise-free, exposure-averaged W3 flux; 160 u0 values in [0, 1.995]; same windows.
- **Result:** the share of u0 values whose two sides do not both reach 2 locally bright cadences is
  0–0.006 for ρ = 0.01 and 0–0.053 for ρ = 0.1 (grazing u0 ≳ 1.9).
- **Status:** this is a model expectation, not a recovery measurement.

## 4. Decision rule, validation gate, search
- **Phase A** (old_bg carriers, EB light curves): every family is evaluated under both rules on the
  same injections.
  - Comparison stars for common mode: old_bg.
  - **The old rule passes the hard tests only if all of these hold:**
    - every family has passed / n ≤ 0.10;
    - eb windows passing the single-window test / flagged windows ≤ 0.10;
    - eb stars with any passing window / eb stars with flagged windows ≤ 0.10;
    - ≥ 3 real-brightening sources;
    - ≥ 5 eb stars with flagged windows.
  - If all hold, the old rule is kept. Otherwise the local rule is used (`decision.json`).
  - W3 recovery on old_bg (seen data) is reported, not gated.
- **Validation** (25 unread val stars, the chosen rule only, gated). Common-mode comparison stars:
  val + old_bg. **The gate passes only if all of these hold:**
  - every one of the 20 W3 cells has recovered / 100 ≥ 0.50;
  - every family on val has passed / 100 ≤ 0.10;
  - real_bright_eclipse (phase A rows under the chosen rule; seen backgrounds; conditional) ≤ 0.10
    with ≥ 3 sources;
  - eb windows ≤ 0.10 and eb stars ≤ 0.10, with ≥ 5 eb stars flagged.
  - If a required control lacks its sources, the gate fails. It is not dropped.
  - Also reported: the uninjected val screen, and the other rule's numbers (not gated).
- **Search:** only if the gate passes, on the 200 search stars, with the chosen rule.
  - Comparison stars for common mode: search + val.
  - Survivors are taken in search-rank order, then time. **At most 3 are vetted**, using:
    - SAP_QUALITY (all bits) within ±1 d;
    - flux-weighted centroids (MOM_CENTR1/2) in the run against the run's quarter;
    - same-exposure controls (other val / search stars on the same channel and quarter at the same
      cadences);
    - PDCSAP, which is never an independent confirmation;
    - the official target pixel file only if the above leave the case open.
  - Each vetted survivor is classified as a real dimming, an artefact or undetermined. A W3-like
    shape is never called wormhole evidence.
- **Gate fails:** search is not read. The failure is attributed from the stored rows.
- **No second try:** no seed, threshold or rule change, and no second rule. An implementation error
  may be fixed and the affected parts re-run, with both results kept.

## 5. Not covered (stated in advance)
- Ordinary phenomena not modelled here:
  - flares next to dips;
  - spotted-star quasi-periodic dips;
  - heartbeat stars;
  - occultations by circumstellar dust (dippers);
  - microlensing of blends;
  - instrumental features only visible in pixels (rolling band, SPSD).
- The real-event tests cover only the 46 EB stars chosen above.
- **Single band:** the achromatic test is "not tested".
- Efficiencies apply from the SAP product onward only.
- **Labels:** injections are simulated; light curves and flags observed; outcomes derived.
