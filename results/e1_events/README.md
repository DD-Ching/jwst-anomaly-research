# E1 causal event network (GBM × ICECAT-1 × GWTC × CHIME/FRB Cat 2)

This tests hypothesis E1 (owner idea 4): are time-stamped events from *different* directions statistically
dependent at time lags no ordinary path explains? The statistic is a pair-count cross-correlation in lag (0–10 s,
10–100 s, 100 s–1 h, 1 h–1 d, 1–7 d), split into same-direction (≤ 3σ combined) and wide (> 3σ and > 0.1°) pairs.
GW events have no position in the GWTC CSV, so their channels are lag-only. A dependence would be an anomaly, never
evidence of anything exotic. Decision: DECISIONS.md "E1 causal event network" (D-074). Sources: SOURCES.md "E1
causal event network".

Reproduce (about 20 min single-process; the null ensembles are cached untracked next to the data):

    python scripts/e1_events.py --cpu 1                 # 3 nulls x 1e4 scrambles, injections, limits
    python scripts/e1_events.py --cpu 1 --chime-vet 2000

| File | Provenance | Content |
|---|---|---|
| `counts.ecsv` | derived | per channel × lag window × class: observed pairs; `perm`, `jit` and `jitday` null mean, sd, z and empirical p; 95 % upper limit on added pairs. Rows with `tested = False` are the vetting windows (k × 95.6 min orbit, k × sidereal day, k × solar day) |
| `limits.ecsv` | derived | per wide (GW: any-separation) cell: limit on dependent pairs, injection efficiency, smallest injected n detected in ≥ 50 % of trials, limit per anchor event |
| `injections.ecsv` | simulated | synthetic wide-separation lagged pairs injected into the real catalogues |
| `chime_vet.json` | derived | CHIME–CHIME lag excess vetting: excluded_flag, the 2023-08-25 same-position episode, the in-day null, and an injection calibration of each null (the in-day null absorbs most of an injected signal, so it does not discriminate) |
| `chime_flag_tests.json` | derived | Ordinary-explanation tests on the CHIME–CHIME 1 h–1 d wide flag (`scripts/e1_chime_flag.py`). (a) Busy-day concentration and a 30-day jackknife. "Pairs touching" counts cell pairs with a member on one of the observed top days, so its "share of excess" can exceed 1; the drop-days z is the direct test. (b) Years, the Catalog 1 period and seasons. (c) DM, fluence, S/N, Dec and sidereal-phase differences of cell pairs against scrambled pairs. (d) A rate-modulated null (7-day and 3-day running means of daily counts) with injection calibration. Result: the calibrated 7-day null keeps about 87 % of an injected signal and gives z = −0.9 on the data |
| `chime_rate_null_100s_1h.json` | derived | The same rate-modulated null on the CHIME–CHIME 100 s–1 h wide cell (`--cell 100s-1h`, 2,000 scrambles), plus a variant keeping each event's time of day. 7-day null: z = 2.5 (family-wise p = 0.26), keeps ≥ 97 % of 300 injected pairs (D-074 addendum) |
| `antipodal.json` | derived | Near-antipodal pairs (180° − sep ≤ max(10°, 3σ_comb)) per lag bin in the six localized channels against the `jit` null (2,000 scrambles, `scripts/e1_antipodal.py`), with per-cell 95 % upper limits on extra pairs and `n50_injected` (smallest injected antipodal count detected at family-wise 3σ in ≥ 10 of 20 trials; null = not reached at 100). Global p = 0.56 (D-074 addendum 2) |
| `signed_lag.json` | derived | Signed-lag asymmetry D = N(B after A) − N(B before A) per cross channel × lag bin × class (same/wide; all for GW) against the `jit` null (2,000 scrambles, `scripts/e1_signed_lag.py`); Gaussian tail, or Skellam where the null sd < 1; `n50_injected` from one-sided injections. Global p = 0.68 (D-074 addendum 3) |
| `chime_exposure_dec_profile.ecsv` | derived | CHIME Cat 2 exposure maps (216 MB, sha256 in `data/manifests/e1_chime_exposure.ecsv`) reduced to a 0.1° Dec profile per transit. The file is time-integrated, with no time axis, so it cannot model per-day uptime. Downloaded at the PR #112 review's request, reduced by `scripts/e1_chime_exposure.py` (needs `h5py`), and the raw file deleted |
| `same_dir_pairs_*.ecsv` | derived | same-direction pairs within 7 d (re-triggers, unflagged repeaters, duplicates) |
| `summary.json` | derived | event counts, global (trials-corrected) p per null and family, analytic tails, reachability, positive controls |

Nulls (both keep each event's declination and hour angle, and permute times within catalogue and calendar year):
`perm` keeps the time multiset, so it only tests whether separation depends on lag. `jit` also shifts each event
(GBM: k × orbit ± 5 min within ±3 d; the others: uniform ±3 d) and is the null for lag clustering. `jitday` shifts
the ground instruments and GW by whole days (keeping the time of day). The empirical p is floored at 1/(n+1) = 10⁻⁴.
Claims beyond that use the analytic tail (`analytic_p`, model_prediction) with a Bonferroni factor.
