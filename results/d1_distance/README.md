# D1 distance self-consistency (hypothesis round 1, owner idea 1; D-TBD)

Question: does any single sightline carry two independent distance measures that no ordinary model reconciles?
Result: **no**. Every per-lens and per-FRB pull is far below the 5σ trials-corrected threshold.

Reproduce: `python scripts/d1_distance.py lenses --h0licow <clone> --tdcosmo <clone>` and
`python scripts/d1_distance.py frb --frb <clone> --n-pred 100000` (inputs pinned in
`data/manifests/d1_distance.ecsv`; seed 20261009).

| File | Provenance | Content |
|---|---|---|
| `lenses.ecsv` | model_prediction | per lens: R = D_dt/((1+z_d) D_d) quantiles, pulls A (flat ΛCDM prior), B (wCDM prior), C (LOO R, common ln-scale free), C0 (LOO R, no scale), D (LOO D_dt, H0 + Ωm free); trials-corrected σ |
| `lens_summary.json` | model_prediction | trials, threshold, shuffled-redshift null, minimum detectable factor per lens |
| `injections.ecsv` | simulated | one lens's D_dt × f (0.7–2.0): recovered Δln and pulls |
| `lenses.png` | model_prediction | R vs flat ΛCDM (Ωm = 0.3), and all pulls with the 5σ trials-corrected lines |
| `frb.ecsv` | model_prediction | per localized FRB: predictive median DM, one-sided low/high tails, P(below the MW-only floor) |

Labels: the posteriors are *derived* (published lens models); the predictions are *model_prediction*; the priors,
scatter models and thresholds are *assumptions* (`distance_consistency.Params`, `FRBParams`).
