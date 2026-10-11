# D1 distance self-consistency (hypothesis round 1, owner idea 1; D-073)

Question: does any single sightline carry two independent distance measures that no ordinary model reconciles?
Result: **no**. Every per-lens and per-FRB pull is far below the 5σ trials-corrected threshold.

Reproduce: `python scripts/d1_distance.py lenses --h0licow <clone> --tdcosmo <clone>` and
`python scripts/d1_distance.py tdcosmo --tdcosmo <clone>`,
`python scripts/d1_distance.py kinematic --tdcosmo <clone>` (seed 20261010),
`python scripts/d1_distance.py frb --frb <clone> [--ism ymw16]` (inputs pinned in
`data/manifests/d1_distance.ecsv`; seed 20261009).

| File | Provenance | Content |
|---|---|---|
| `lenses.ecsv` | model_prediction | per lens: R = D_dt/((1+z_d) D_d) quantiles, pulls A (flat ΛCDM prior), B (wCDM prior), C (LOO R, common ln-scale free), C0 (LOO R, no scale), D (LOO D_dt, H0 + Ωm free); trials-corrected σ |
| `lens_summary.json` | model_prediction | trials, threshold, shuffled-redshift null, minimum detectable factor per lens |
| `injections.ecsv` | simulated | one lens's D_dt × f (0.15–8, f = 1 is the baseline): pulls and the whole-chain flag; the detectable factors in `lens_summary.json` come from where these cross the threshold |
| `lenses.png` | model_prediction | R vs flat ΛCDM (Ωm = 0.3), and all pulls with the 5σ trials-corrected lines |
| `tdcosmo_lenses.ecsv` | model_prediction | statistic D on the TDCOSMO 2025 power-law D_dt chains, 8 lenses × (`kext`: own κ_ext PDF; `nokext`: κ_ext = 0); D-073 addendum 2 |
| `tdcosmo_summary.json`, `tdcosmo_injections.ecsv` | model_prediction / simulated | trials, threshold, shuffled-z null, detectable factors; D_dt × f injections |
| `tdcosmo_*_j1206_final_composite.*`, `tdcosmo_*_j1206_final_power_law.*` | as above | `tdcosmo --j1206 final_composite` / `final_power_law`: SDSS1206 replaced by TDCOSMO's final chain (κ_ext included), `kext` only; composite = model-choice check, power law = validation; D-073 addendum 3 |
| `kinematic_lenses.ecsv`, `kinematic_summary.json`, `kinematic_lenses.png` | model_prediction | per lens (76: TDCOSMO 8, SLACS KCWI 13, SLACS SDSS 41, SL2S 14): kinematic ln D_s/D_ds (λ = 1) from the hierArc likelihood pickles, y = minus ln D_s/D_ds(ΛCDM), LOO pull vs the others' offset + intrinsic scatter (Ωm 0.3/0.1/0.5), λ = 1 pull (diagnostic), D_s/D_ds × f injection factors, shuffled-z null; D-073 addendum 4 |
| `kinematic_*_rmax400.*` | model_prediction | the same with the D_s/D_ds prior bound at 400 (`--ratio-max 400 --om 0.3`): robustness of the power-law upper tail |
| `frb.ecsv` | model_prediction | per localized FRB: predictive median DM, one-sided low/high tails (grid convolution), P(below the MW-only floor), the DM at which each tail reaches the 5.81σ flag |
| `frb_injections.ecsv` | simulated | DM_obs = 0.9 × low limit and 1.1 × high limit per FRB, through the same chain |
| `frb_summary.json` | model_prediction | FRB flags, extremes, sensitivity and injection recovery |
| `frb_ymw16.ecsv`, `frb_injections_ymw16.ecsv`, `frb_summary_ymw16.json` | model_prediction / simulated | the same with DM_ISM from YMW16 (`--ism ymw16`; D-073 addendum) |
| `frb_ism_compare.ecsv` | model_prediction | NE2001 vs YMW16 per burst: DM_ISM ratio, pulls, `flag_either` |

Sensitivity (from injections): a lens is flagged only if D_dt/D_d is off by ×2–8 up or ×0.2–0.4 down (R tests), or
D_dt by ×1.3–1.85 / ×0.54–0.78 against the other lenses. Every FRB injection beyond its limit is flagged (94/94 per
side). The low-side FRB limit depends on the sharp lower cutoff of Macquart's p(Δ), which is an assumption.

Labels: the posteriors are *derived* (published lens models); the predictions are *model_prediction*; the priors,
scatter models and thresholds are *assumptions* (`distance_consistency.Params`, `FRBParams`).

YMW16 (`--ism ymw16`) needs pygedm's compiled `ymw16` extension. Where `pip install pygedm` fails (no libf2c for
its NE2001 part): `pip download --no-deps --no-binary :all: pygedm==3.3.0`, delete the `ne21c` Extension from
`setup.py`, wrap `import ne21c` in `pygedm/ne2001_wrapper.py` in try/except, then
`pip install pybind11 "setuptools<81" && pip install --no-build-isolation .`. The script imports only `ymw16`.
