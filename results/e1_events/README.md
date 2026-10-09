# E1 causal event network (GBM × ICECAT-1 × GWTC × CHIME/FRB Cat 2)

This tests hypothesis E1 (owner idea 4): are time-stamped events from *different* directions statistically
dependent at time lags no ordinary path explains? The statistic is a pair-count cross-correlation in lag (0–10 s,
10–100 s, 100 s–1 h, 1 h–1 d, 1–7 d), split into same-direction (≤ 3σ combined) and wide (> 3σ and > 0.1°) pairs.
GW events have no position in the GWTC CSV, so their channels are lag-only. A dependence would be an anomaly, never
evidence of anything exotic. Decision: DECISIONS.md "E1 causal event network" (D-TBD). Sources: SOURCES.md "E1
causal event network".

Reproduce (about 20 min single-process; the null ensembles are cached untracked next to the data):

    python scripts/e1_events.py --cpu 1                 # 3 nulls x 1e4 scrambles, injections, limits
    python scripts/e1_events.py --cpu 1 --chime-vet 2000

| File | Provenance | Content |
|---|---|---|
| `counts.ecsv` | derived | per channel × lag window × class: observed pairs; `perm`, `jit` and `jitday` null mean, sd, z and empirical p; 95 % upper limit on added pairs. Rows with `tested = False` are the vetting windows (k × 95.6 min orbit, k × sidereal day, k × solar day) |
| `limits.ecsv` | derived | per wide (GW: any-separation) cell: limit on dependent pairs, injection efficiency, smallest injected n detected in ≥ 50 % of trials, limit per anchor event |
| `injections.ecsv` | simulated | synthetic wide-separation lagged pairs injected into the real catalogues |
| `chime_vet.json` | derived | CHIME–CHIME lag excess vetting: excluded_flag, the 2023-08-25 same-position episode, the in-day null |
| `same_dir_pairs_*.ecsv` | derived | same-direction pairs within 7 d (re-triggers, unflagged repeaters, duplicates) |
| `summary.json` | derived | event counts, global (trials-corrected) p per null and family, analytic tails, reachability, positive controls |

Nulls (both keep each event's declination and hour angle, and permute times within catalogue and calendar year):
`perm` keeps the time multiset, so it only tests whether separation depends on lag. `jit` also shifts each event
(GBM: k × orbit ± 5 min within ±3 d; the others: uniform ±3 d) and is the null for lag clustering. `jitday` shifts
the ground instruments and GW by whole days (keeping the time of day). The empirical p is floored at 1/(n+1) = 10⁻⁴.
Claims beyond that use the analytic tail (`analytic_p`, model_prediction) with a Bonferroni factor.
