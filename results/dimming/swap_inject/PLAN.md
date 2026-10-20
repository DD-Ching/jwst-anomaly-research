# W3 dimming: the D-039 persistence swap in `inject` and the re-run `combine`

Frozen 2026-10-20 04:12 UTC, before any swapped efficiency was computed (handoff 1 of `forced20/RESULT.md`,
`expcount_eff/RESULT.md` "Next" 1). Screen, thresholds, `forced` and the D-039 test unchanged from `main` (06917c1).

- **Question.** How much do the three-field W3 limits (`wht_coverage/limits_combined.ecsv`, the current headline)
  gain when an injected copy that the catalogue detects in one epoch only is no longer vetoed outright but goes
  through the D-039 persistence test, as real flags now do (#379, #390)?
- **Chain emulated (production order).** A real compact single-epoch flag is cleared only if it is a `vanish` flag
  (the persistence-run target selection), the D-039 test calls it `on_sky` in its detection epoch, and the
  `forced` + depth / footprint step then runs (already emulated by the `forced` inflation proxy in `_recovered`).
  So per injected copy:
  - detected in ≥ 2 epochs and flagged (with the inflation proxy): weight 1 (unchanged);
  - detected in exactly one epoch k, `vanish` flag at both noise levels: weight ε(i, k);
  - otherwise 0. Efficiency per bin = mean weight (expected recoveries; `n_recovered` becomes a float).
- **ε(i, k) (ASSUMPTION, frozen in `expcount_eff/PLAN.md`):** 0.412 if source i has ≤ 4 covering NIRCam `_cal`
  exposures in epoch k, 0.961 if > 4. `n_exposures` per (source, epoch) is read from the committed
  `persistence_run/<field>/verdicts.ecsv` (monitored rows; uid = `<field>-d<row>`, the light-curve row index,
  which the weight-map term does not change). Pairs without a committed count (NEXUS epoch 2, never run for
  controls; or any pair absent from the table) get **ε = 0** (the veto stays: conservative).
- **Variants.** `swap` (above, primary); `swap_pooled` (≤ 4 class at the pooled Wilson lower bound 0.374;
  sensitivity); `swap_unknown_low` (unknown pairs at 0.374 instead of 0; sensitivity, not headline).
- **Gates (computable; `run.py`).**
  - G0: with the swap off, `efficiency_table` reproduces `wht_coverage/runs/<field>/limits_veto_on_on.ecsv`
    `exposure_source_yr` (rtol 1e-9) per field. Failing G0 stops the run (no swapped number is reported).
  - G1: per field and (t_E, ρ) cell, veto-on ≤ swap ≤ veto-off (`limits_veto_off_on.ecsv`) within rtol 1e-9.
    A failing cell is a bug; no headline is quoted until it is explained.
- **Output.** `runs/<field>/limits_<variant>.ecsv`, `limits_combined_<variant>.ecsv` (calibrated fields only, as
  `combine`), `summary.json` (gates, gain per t_E at ρ = 0.1, monitored detected (source, epoch) pairs with and without a committed count).
- **What would change the plan.** A gain < 10 % at every t_E would mean the swap is not worth carrying in
  production `inject` beyond this record; if `swap_unknown_low` raises the
  joint exposure by > 10 % over `swap` at some t_E, the pairs without a count (NEXUS e2) matter and NEXUS e2
  controls become the next step.
