# W3 / Kepler effective search volume of the #542 sample — frozen measurement plan

Owner authorization 2026-10-25 (23:41 UTC). Limits: 4 h, 6 GB through `scripts/net_ledger.py` (this round's
ledger `net_ledger.jsonl` here, cap 6,000,000,000 B), 4 cores. This file is committed before any data are
restored and before any injection. Code: `scripts/w3_kepler_effsearch.py` and its tests, committed before the
injection run.

**What this is:** an injection-recovery measurement of search efficiency on the 500 stars that #542 already
searched (`results/w3_kepler/expand/`). It is not a new blind test: these light curves have been looked at.
No screening rule, threshold or test changes; no new star is read; no shoulder repair; no second method.
Labels: injections are simulated; light curves observed; recoveries, efficiencies and exposures derived.

## 1. Sample and data
- **Stars:** the 500 stars of `expand/sample_expand.json`, all kept. That includes the 5 stars outside the #542
  applicability proxy (`expand/applicability.json`) and the 10 stars in the Kirk et al. 2016 EB catalogue
  (`expand/inputs/kirk2016_catalog.csv`; among them KIC 3954798). No star is dropped or replaced because of a
  result.
- **Restore by list and hash:** each star's Q1–Q17 `_llc.fits` files are fetched by the names recorded in
  `expand/files_expand.json.gz`, and each file's sha256 must equal the recorded one. A mismatch or a failed
  fetch is recorded; the star is not replaced.
  - **Download change (efficiency rule 4.3):** the 245 Q0 files (19.1 MB in #542) are not fetched. #542
    never used them (`w3_kepler._fetch_star` drops Q0), so the column cache is unchanged.
  - No directory listing is fetched (the file names are known).
  - Column extraction and the cache layout are those of `w3_kepler._fetch_star`.
- **Restore check (before any injection):** the uninjected #542 screen is run again on the restored cache
  (`w3_kepler_expand` common-mode counts over the 500 stars, `frac_loo`, `w3_kepler.screen_star(local=True)`).
  - It must reproduce, for every star, `good_cadences` of `expand/screen_expand.json.gz` and the flagged runs
    of `expand/runs_expand.json.gz`: `t_start`, `t_end`, `n` and every pass flag.
  - On any mismatch no injection runs; the round reports the mismatch.

## 2. Screening chain (unchanged, #529 / #542)
- Quality mask `QMASK`; per-quarter median normalisation (a quarter with < 40 good cadences is left out);
  1-mag dim threshold; runs of ≥ 3 consecutive dim good cadences with ≥ 20 non-dim cadences on each side.
- Then:
  - common mode < 0.2 (leave-one-out over the uninjected 500 stars);
  - bounded fast edges (v23 bounds rule per quarter, EDGE_FRAC 0.25, 1-d span floor);
  - one-off test;
  - local shoulder rule.
- Survivor = `pass_local`.
- The 20 W3 cells of `w3_kepler.CELLS`:
  - ρ ∈ {0.01, 0.1};
  - blend ∈ {1.0, 0.7};
  - t_E ∈ {0.1, 0.3, 1, 3, 10} d.
- The model is a finite-source n = 1 negative-mass point lens (`exotic_sim`, via `w3_kepler.w3_flux_avg`),
  averaged over the 29.4-min exposure with 21 sub-samples.

## 3. Injections (fixed before the run)
- **Count:** 2 independent injections per star per cell: 500 × 20 × 2 = 20,000. Each starts from the star's
  original, uninjected SAP series; injections are never stacked.
- **Event distribution:**
  - Event centre t0 ~ uniform on the star's sampling window W_i = [T_lo, T_hi]:
    - T_lo = the first finite TIME among the star's Q1–Q17 cadences (any quality flag);
    - T_hi = the last such TIME.
    - W_i includes the inter-quarter and monthly downlink gaps, safe modes and flagged cadences.
  - Events whose centre falls in a gap, at a quarter boundary or in sparse data stay in the denominator.
  - Impact parameter u0 ~ uniform on [0, 2): umbra-crossing events only. Events with u0 ≥ 2 never dim, so the
    dimming screen cannot find them by construction; they are outside the event definition.
- **Noise:** `w3_kepler.inject_sap`. SAP × F plus Gaussian noise, so that the variance is err² × max(F, 0.09);
  cadences with F = 1 are unchanged.
- **Order:** the injection is applied to the raw SAP flux before the quality mask, the per-quarter
  normalisation and every baseline estimate. Quality flags come from the archive file and are not changed by
  the injection (see section 7).
- **Common mode:** the injected star's own cadences are left out (`frac_loo`); the other 499 stars are
  uninjected.
- **Seeds:** `numpy.random.default_rng([54200, rank, cell, rep])`:
  - rank is the star's index in `sample_expand.json`;
  - cell is 0–19 in `w3_kepler.CELLS` order;
  - rep is 0 or 1.
  - Draw order: t0, u0, then the noise.
- **Run order:** stars in `sample_expand.json` order. Each job is one star with all 40 injections, so a run cut
  short leaves whole stars.

## 4. Recovery (fixed)
- **Matched run:** a flagged run of the injected screen that meets all three conditions:
  - (a) its (t_start, t_end) is not among the star's uninjected runs;
  - (b) it overlaps the umbra interval [t0 − h, t0 + h], h = t_E √(4 − u0²) (#529);
  - (c) more than half of its cadences lie in [t0 − h − Δ/2, t0 + h + Δ/2], with Δ = 29.4 min.
    - This stops an existing eclipse run that the injection lengthens from counting.
- **Recovered:** some matched run has `pass_local`.
- **Spurious:** a new `pass_local` run that is not matched. It never counts as a recovery and is reported
  separately.
- **Loss stage of each unrecovered event,** tested in this order (F = the noiseless exposure-averaged model flux
  factor, dim = F ≤ 10^(−0.4)):
  1. `short`: fewer than 3 dim cadences on an ideal gap-free 29.4-min grid. The grid is phased to the star's
     nearest finite cadence time. Even full coverage cannot make a 3-cadence 1-mag run: a by-construction loss
     of the existing rule.
  2. `gap`: ≥ 3 dim cadences on the ideal grid, but fewer than 3 dim cadences with finite TIME and SAP_FLUX.
  3. `quality`: ≥ 3 such cadences exist, but fewer than 3 of them are good (QMASK, error ≤ 0, or a quarter left
     out).
  4. `not_flagged`: ≥ 3 dim good cadences (covered), but no matched run. Causes: noise, the 20-cadence bracket,
     or a run broken by brighter cadences.
  5. Otherwise, the first failed test of the deepest matched run, in the order common mode → edges (stored as
     `slow` or `insufficient`) → one-off → local shoulder (`no` or `insufficient`).
- **Per event, stored in `events.json.gz`:** rank, kepid, cell, rep, t0, u0, h, the three dim-cadence counts,
  new / matched runs, recovered, spurious, loss stage and edge / shoulder sub-classes.

## 5. Effective search volume and its uncertainty (fixed)
- **Per star i:** window T_i = (T_hi − T_lo) / 365.25 yr; k_ic ∈ {0, 1, 2} recoveries in cell c; p̂_ic = k_ic / 2.
- **Stored separately per cell:**
  - observing time T_obs (good cadences × 29.4 min, the #542 1,732 star-years);
  - sampling time ΣT_i;
  - recovery rate ε̄_c = Σ T_i p̂_ic / Σ T_i;
  - effective search volume E_c = Σ T_i p̂_ic, in star-years;
  - the loss-stage counts.
- **No double counting:** E_c uses the same window as the t0 draw. It is never T_obs × ε̄_c, because ε̄_c
  already contains the gap losses.
- **Uncertainty for these 500 stars** (finite simulations):
  - Given the stars, each cell has 2 independent Bernoulli trials per star (1,000 per cell).
  - E_c's 95 % interval = ΣT_i × the Wilson score interval of ε̄_c, with Kish's effective trial number
    n_eff = 2 (ΣT_i)² / ΣT_i².
  - This is conservative for equal T_i, because Σ p_i(1 − p_i) ≤ n p̄(1 − p̄).
- **Uncertainty per new star** (for extrapolation; star = cluster):
  - e_i = T_i p̂_ic; mean m_c and cluster standard error sd(e_i) / √500. This includes the between-star spread
    and the 2-draw noise, and treats the 20,000 injections as 500 clusters, not 20,000 independent stars.
  - 95 % interval: m_c ± 1.96 SE, widened where needed to the MC interval divided by 500.
- **Per-star efficiencies are not estimated:** 2 draws per star and cell cannot measure one star's efficiency.
  Stored k_ic is only for re-computation; subgroup sums (the 5 out-of-range stars, the 10 EB stars) are
  descriptive.

## 6. What a zero supports (fixed)
- #542 found 0 survivors in these 500 stars.
- **Per-cell conditional upper limit, product level only:** Γ_95,c = 2.996 / E_c,low, the classical Poisson
  95 % limit for zero events (−ln 0.05), at the lower end of E_c's MC interval.
  - Units: W3 signatures of that cell's shape per star-year, with t0 uniform in time and u0 uniform on [0, 2),
    as they would appear in the SAP product.
- **Conditions,** stated with every limit:
  - events are independent and Poisson, at a constant rate over stars and time;
  - every event has that one cell's parameters;
  - the efficiency is the measured one;
  - a star rarely hosts two events. At the limit, λ_i = Γ_95,c T_i, and the share of event-hosting stars with
    ≥ 2 events, 1 − λ e^(−λ) / (1 − e^(−λ)), is reported. Repeating or periodic events are outside the
    model, because the one-off test rejects them.
- **The 20 limits are separate statements,** one per hypothetical single-cell population. They are not a joint
  limit, and no mixture limit is made.
- **Cells:** a cell with E_c,low = 0 has no limit.
- **Not done:** no conversion to a physical lens rate, wormhole number, cosmic abundance or anything else. The
  SAP response to a real event is not validated (section 7). The limits are labelled product-level, and the
  missing pieces are listed.

## 7. Known limits of the measurement
- The efficiency starts at the SAP product. Not covered: aperture photometry of a real event, crowding beyond
  the blend parameter, saturation, background, and whether the archive's quality flags would have flagged
  real dim cadences. For example, the impulsive-outlier bit 512 is set from the data; an injection cannot set
  it.
- The event distribution (t0 and u0 uniform, one cell at a time) is a choice, not a physical model.
- The sample is 500 FGK pool stars chosen by hash after earlier exclusions. The extrapolation in section 8 is
  valid for the same pool under the same rule.

## 8. Expansion estimate (fixed formulas)
- **Gain:** for N new stars (N = 500 and 5,000) from the same pool under the same hash rule,
  ΔE_c = N × m_c, with the 95 % interval from section 5.
  - Limit factor if the result stays zero: E_c / (E_c + ΔE_c).
  - The pool has 18,240 − 500 = 17,740 unused stars (`expand/sample_expand.json`).
- **Cost per 500 stars:**
  - download bytes (#542 ledger, minus Q0);
  - fetch and screen time (#542 `cost.json`);
  - this round's measured injection time.
- **Not given:**
  - no probability of finding an event, which needs an event-rate model;
  - no claim that fainter or more crowded stars would do better.

## 9. Stop rules
- Stop at 4 h, at the 6 GB ledger cap, or if the restore check fails.
- The seeds are never changed for a better result, and no low-efficiency cell is dropped.
- Implementation errors may be fixed, and the affected part re-run, with both records kept.
- If the run cannot finish, the completed stars are reported as a subset, never as the batch.
- The column cache is deleted after the review re-check.
