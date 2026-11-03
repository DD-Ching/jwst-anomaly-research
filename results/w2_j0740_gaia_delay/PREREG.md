# J0740+2926: can Gaia DR3 G epoch photometry identify a delayed replay? Pre-registration

Frozen 2026-11-03 ~06:15 UTC, before any light curve of the pair was fetched. Owner-requested single-pair test,
at most 60 minutes. Any later change is listed as a deviation in RESULT.md.

**Pair:** A = Gaia DR3 878723302007625344 (G 18.24), B = 878723302009111040 (G 19.57); repo-verified identities
(`results/w2_desi_pairs/evidence/identity.json`). SDSS J074013.42+292645.8, separation 2.64″, DESI z 0.9800 / 0.9779.

**Seen before freezing:** `gaia_source`, `vari_summary` and `vari_agn` rows only. Both are VARIABLE with epoch
photometry; 22 (A) and 25 (B) selected G FoV epochs over 971.6 d; std 0.092 / 0.078 mag; std / rms error 5.5 / 3.0;
Abbe 0.21 / 0.29. No light-curve values and no delay result.

## 1. Data and quality cuts
- **Q1:** G per-transit rows from the Gaia DataLink `EPOCH_PHOTOMETRY` service (`DATAMODEL_STANDARD`). Rows are
  dropped if flagged `variability_flag_g_reject` or `rejected_by_photometry`, or if time, magnitude or
  flux-over-error is missing. This is the filter of `gaia_mulens.parse_epoch_csv`.
- **Q2:** error = 1.0857 / flux-over-error, multiplied by s = max(1, √χ²_ν). χ²_ν is pooled over transits that lie in
  groups ≤ 1.0 d long, about each group's mean, for both sources. No other clipping.
- **Q3:** times are `g_transit_time` (BJD(TCB) − 2455197.5); no binning.

## 2. Delay convention and range
- Δt > 0 means B repeats A's light curve Δt later: m_B(t) = m_A(t − Δt) + const, plus slow microlensing.
- Search range |Δt| ≤ 400 d, which leaves ≥ 570 d of overlap in the 971.6 d span.
- |Δt| < 5 d is not interpreted: same-transit common errors and A→B crosstalk act at zero lag.
- Ranges named for the exclusion verdict:
  - R_lens = +20 … +400 d: the fainter B lags, as for an ordinary galaxy lens (SIS value given in RESULT.md);
  - R_all = −400 … +400 d.

## 3. Estimator (PyCS3 3.0.4, spline technique)
- `pycs3.spl.topopt.opt_fine`, knot step 100 d.
- B carries a linear microlensing polynomial (`pycs3.gen.polyml.addtolc`, nparams = 2, one season); A carries none.
- Starts from 17 initial delays (−400 … +400, step 50). The lowest-r² solution inside the range is Δt̂.
- Statistic D(Δt) = r²_joint(Δt) − r²_A − r²_B:
  - r²_joint is the PyCS3 r² of the joint spline with the delay held at Δt;
  - r²_A and r²_B come from separate spline fits with the same knot step;
  - D is the cost of forcing one intrinsic curve on both images.
- Cross-check on the target only: PyCS3 regdiff (`pycs3.regdiff.multiopt.opt_ts`, default kernel).

## 4. Controls (actual epochs and Q2 errors of A and B)
**Simulation inputs:**
- Intrinsic curve: damped random walk, τ = 300 d (observed frame) and σ = 0.10 mag [assumption; the scale roughly
  matches the summary std].
- Common errors [assumption]: transits of A and B within 0.01 d share a Gaussian term of 0.01 mag. B also receives
  1 % of A's flux (crosstalk).

**Controls:**
- **C1, recovery:** same-source pairs with Δt_true ∈ {−200, −100, +50, +100, +200, +300}, N = 40 each.
- **C2, false positives:** independent pairs with the same parameters, N = 100. They include the sampling gaps and
  the common errors. The distribution of Δt̂ is reported, to show aliases from the gaps.
- **C3, real independent sources:** the two `vari_agn` sources nearest to J0740 within 1° with 18.0 ≤ G ≤ 19.8,
  excluding the pair. They are run through the same pipeline as a pair.
- N may be reduced to fit the 60 minutes; any reduction is reported.

## 5. Thresholds and verdict
**Thresholds:** T95 and T99 are the 95th and 99th percentiles of D(Δt_true) over all C1 simulations.

**Pass conditions:**
- C1 passes at a delay if ≥ 68 % of its simulations give |Δt̂ − Δt_true| ≤ 20 d.
- C2 passes if ≤ 5 % of independent simulations give D(Δt̂) ≤ T95 with 5 ≤ |Δt̂| ≤ 400.
- C3 passes if the real independent pair gives D(Δt̂) > T95.

**Verdict, decided in this order:**
1. **Supports a same-source delay** only if all of these hold:
   - C1 passes at some delay, and C2 and C3 pass;
   - the target gives D(Δt̂) ≤ T95 with 5 ≤ |Δt̂| ≤ 400;
   - regdiff agrees with Δt̂ to within 20 d.
2. **Excludes:** the delays on a 5-day grid with D_target(Δt) > T99.
   - Excluded intervals are reported. R_lens or R_all counts as excluded only if fully covered.
   - The power is reported as the fraction of the range excluded for the C2 independent pairs.
3. **Otherwise: data insufficient.**

A "supports" verdict would still only make J0740+2926 a lens candidate.
