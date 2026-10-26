# W3 TESS 120 s product-level injections: real light curves keep 0.477 and 0.459 of the ideally identifiable simulated events at the two #561 best masses; both intervals cross one half (not a robust pass); the low masses are insufficient

- **Plan:** `PLAN.md`, sample, scenarios, detector and tests were committed at 14:04:59 UTC (9f7de1e), before any
  flux.
- **Downloads:** 362,291,257 B of the 1 GB cap (`net_ledger.jsonl.gz`).
  - 184 SPOC light curves, with their sha256 in `files.json`. Each was deleted after processing.
  - Plus the metadata queries and the two-file re-computation check (`recheck.json`).
- **Cores:** 2. **Wall clock:** from 13:59 UTC.
- **Labels:** every dimming and every recovery here is a simulated injection. The TESS light curves are observed
  data, used only as a background. The efficiency is product level: after SPOC processing, not from sky to result.
- **#561, kept as is:** its frozen investment gate passed only at |M| ≲ 2 × 10⁻⁵ M☉, where the Kepler long-cadence
  baseline has no reach. That is a new mass range, not a demonstrated tenfold gain in discovery chance.
- **Deviations:**
  - The download route changed before any flux: MAST's `Download/file` answers with a 302 to MAST's public AWS
    copy, and the ledger does not follow redirects. The same files are read from that copy (commit after
    9f7de1e).
  - Before the run finished, `summarize` was dry-run on the first 95 star-sectors to check the code. Its output was
    deleted, and no rule, seed or sample was changed.
  - **Post hoc:**
    - `diagnostics.json` (where the last step loses events; sample medians);
    - the derived fields `kept_per_step`, `density_for_1_x_retention_p5_p95` and
      `ratio_to_baseline_same_mass_x_retention` in `summary.json` (bootstrap unchanged);
    - `recheck.json` (stage `recheck`): 2 star-sectors re-downloaded; sha256 and the re-processed results are identical.

## 1. Sample and method (`sample.json`, `PLAN.md`)
- **Pool:** 14,301 single-sector 120 s records of the #561 TIC sample with R_* and D_S.
- **Draws:** 200 with replacement, ∝ span × D_S^1.5 (seed 564). That is each star-sector's rate of umbra crossings,
  the Level-0 term of #561's integral. They cover 184 star-sectors of 149 stars.
- **The draws are dominated by distant giants** (`sample_per_draw`): median D_S 1,265.1 pc, median R_* 10.178 R☉,
  68.5 % with R_* > 2 R☉. That is what the D_S^1.5 weight of #561's reach selects.
- **Events:** per draw and scenario, 30 drawn from the crossing population.
  - x ∝ √(x(1 − x)); u0 uniform in 0–2.65; t0 uniform over the light curve's full span, gaps included.
  - t_E and ρ follow from |M|, x, v, R_* and D_S.
- **Identifiability is evaluated once,** in the paired ideal, never in the draw.
  - Check: the ideal share of drawn events should equal #561's Level-1/Level-0 × 2/2.65.
  - It is 0.1387 against 0.1279 and 0.1078 against 0.09991 for the two best masses.
  - Events cluster by star-sector, so no formal test is made.
- **Paired levels for each event** (same x, u0, t0 and exposure phase):
  - L0 ideal: the #561 noise-free criterion on a gap-free grid;
  - L1a: noise-free on the file's real cadence times;
  - L1b: as L1a, plus our quality screen (PDCSAP finite and QUALITY = 0);
  - L2: injected into the published PDCSAP flux before the screen, baseline and detection.
- **The detector is blind:**
  - a 3.0 d running-median baseline;
  - a dark run of ≥ 3 contiguous cadences ≥ 1 mag and ≥ 5σ;
  - ≥ 2 shoulder cadences per side ≥ 10 % and ≥ 5σ, within max(run duration, 4 cadences).
  - The truth is used only afterwards, to match a detection to the umbra interval.
- **PDCSAP:** crowding-corrected (blend 1) and cotrended. A post-PDC injection cannot test whether SPOC would
  remove a real event (onboard cosmic-ray mitigation, PDC cotrending, flagging).
- **Cadences that pass our screen:** 5 / 50 / 95 % of each light curve's cadences: 0.656 / 0.899 / 0.999.

## 2. How much is kept (`summary.json`; 6,000 simulated events per scenario)

| scenario | ideal identifiable (share; #561 expectation) | L1a | L1b | L2 | retention L2 / ideal (5–95 %) | frozen verdict |
|---|---|---|---|---|---|---|
| 220 km/s, 5.623 × 10⁻⁶ M☉ | 832 (0.1387; 0.1279) | 708 | 594 | 397 | 0.4772 (0.4075–0.5499) | interval crosses one half (not a robust pass) |
| 220 km/s, 10⁻⁷ M☉ | 26 (0.004333; 0.00492) | 14 | 13 | 3 | 0.1154 (0–0.3) | insufficient (< 30 ideal) |
| 30 km/s, 3.162 × 10⁻⁶ M☉ | 647 (0.1078; 0.09991) | 536 | 412 | 297 | 0.459 (0.3762–0.5407) | interval crosses one half (not a robust pass) |
| 30 km/s, 10⁻⁸ M☉ | 4 (0.0006667; 0.001935) | 2 | 2 | 0 | 0 | insufficient (< 30 ideal) |

- **Kept per step** (`kept_per_step`, best masses, 220 / 30 km/s):

| step | 220 km/s | 30 km/s |
|---|---|---|
| gaps, sector edges and the detector's contiguity and window (L1a/ideal) | 0.851 | 0.8284 |
| quality flags (L1b/L1a) | 0.839 | 0.7687 |
| real background: noise, variability, baseline, 5σ (L2/L1b) | 0.6684 | 0.7209 |

- **Where the last step loses** (`diagnostics.json`, post hoc):
  - 220 km/s: 199 events lost from L1b to L2, of which 87 still had the dark run (their shoulders failed). The
    other 112 lost the dark run.
  - 30 km/s: 117 lost, of which 31 still had the dark run.
  - By ρ bin, L2 / ideal ranges 0.45–0.50 (220 km/s) and 0.39–0.54 (30 km/s). No cause is inferred.
- **By t_E** (cadences; L2 / ideal):

| t_E bin (cadences) | 220 km/s, 5.623 × 10⁻⁶ M☉ | 30 km/s, 3.162 × 10⁻⁶ M☉ |
|---|---|---|
| 2.5–5 | 28 / 44 = 0.6364 | — |
| 5–15 | 191 / 315 = 0.6063 | 1 / 1 |
| 15–50 | 177 / 472 = 0.375 | 82 / 138 = 0.5942 |
| ≥ 50 | no ideal event | 214 / 508 = 0.4213 |

  - The low-mass scenarios have almost no ideally identifiable events in these rate-weighted draws: 26 and 4 of
    6,000. Their t_E is mostly < 5 cadences at 220 km/s and < 15 at 30 km/s (`by_tE_cadences`).
- **Background** (the same detector on the 184 un-injected light curves): 0 dark runs and 0 dark-plus-shoulder flags.
  No extra flag burden was seen at these thresholds. No source was examined.
- **Simulated recoveries without an ideal event** (`L2_not_ideal`): 1, 1, 5 and 0 per scenario. None overlaps a
  background detection (`bg_overlap_in_ideal` = 0).

## 3. What this does to #561's reach (`summary.json` → `reach_561`; only these masses and speeds)
- The same TIC sample and scaling as #561. The retention is the rate-weighted mean over this sample's draws, not a
  per-star efficiency.

| scenario | #561 density for 1 event, ideal (M☉ pc⁻³) | × retention (5–95 %) | ratio to the Kepler LC baseline at the same mass |
|---|---|---|---|
| 220 km/s, 5.623 × 10⁻⁶ M☉ | 0.08734 | 0.183 (0.1588–0.2143) | 2,761 (the baseline barely reaches this mass) |
| 30 km/s, 3.162 × 10⁻⁶ M☉ | 0.6146 | 1.339 (1.137–1.634) | 1.914 |
| 220 km/s, 10⁻⁷ M☉ | 0.3028 | 2.624 (insufficient events) | baseline 0 |
| 30 km/s, 10⁻⁸ M☉ | 1.785 | no recovery (insufficient events) | baseline 0 |

- **The baseline column is still #561's noise-free ideal.** Kepler LC product-level recovery was measured only on
  its own (t_E, ρ) grid (#546, #554) and is not carried over.
- **Scale references only (Read 2014; not negative-mass densities):** local baryons 0.0914; extrapolated dark
  matter 0.005–0.015 M☉ pc⁻³.
  - The best corrected case, 0.183, is about twice the baryon reference (not stored).
  - A density for one expected event is not a limit or an exclusion.

## 4. What a real search would need (`summary.json` → `search_requirement`)
- **Data:** every TESS 120 s light curve is 1,737,708 files × 1,920,741 B (mean here) ≈ 3.338 TB. The top 10 % of
  star-sectors by crossing weight hold 67.91 % of it, about 334 GB (not stored).
  - That weight is Level 0; the identifiable share per star-sector is not in it.
- **Time** (estimate, not stored): this run used about 2 × 22 CPU-minutes for 184 light curves with 24,000
  injections. A background-only search costs well under 0.1 s per light curve. The download dominates.
- **Physical sensitivity:**
  - At the most favourable tested population (220 km/s, 5.623 × 10⁻⁶ M☉, uniform density), a full null search would
    bear on densities of about 0.18 M☉ pc⁻³ and above.
  - At 30 km/s, about 1.3. At 10⁻⁷ or 10⁻⁸ M☉ the measured recovery is not usable.
  - These are product-level numbers and miss the upstream SPOC steps.

## 5. Judgement: a search is not worth starting now
- **Retention:** retaining at least half is not robustly shown (both intervals cross 0.5). The point values are
  just under half.
- **Corrected reach:** about twice the baryon reference even at the best tested mass. A null result over 3.3 TB
  would constrain only densities at that scale.
- **Low masses:** the masses where short cadence opens new ground (10⁻⁸–10⁻⁷ M☉) give too few identifiable events
  in the rate-weighted sample to measure. The draws are dominated by giants, whose ρ is large.
- **What the round adds:** a reproducible product-level injection chain for TESS 120 s, and the measured losses
  per step. It is not physical evidence.

## 6. Not measured
- Whether SPOC would keep a real event: onboard cosmic-ray mitigation, PDC cotrending, and flags set because of the
  event itself.
- Crowding and blends below 1, limb darkening, other masses and speeds, a population model.
- TESS 20 s and Kepler SC.
- **Three separate questions:**
  - (i) The model produces this shape: yes, noise-free (#561).
  - (ii) Data support the model: not addressed (no search).
  - (iii) A controllable signal or channel: not addressed.
