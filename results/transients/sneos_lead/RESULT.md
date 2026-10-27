# SN Eos leading pair: no public epoch can test it (2026-10-27)

Plan `PLAN.md` (committed before the run). Tool `scripts/sneos_lead_epochs.py`, rule functions tested in
`tests/test_sneos_lead_epochs.py`. Numbers from `epochs.json`. Epoch rows are MAST metadata (`observed`); the leading
images' positions, μ and Δt are model_predictions of the VENUS AstroLensPy v2.1 model (`../sneos/summary.json`).

## Result (frozen rule T)
- **Outcome U.** 138 epoch × filter × image rows (HST and JWST, both leading images, 2009–2026); 0 pass rule T,
  and 0 pass the extended-window variant.
- Between the CLASH epochs (2012-04 to 2012-06) and HST PID 17504 (2024-03-11) MAST lists two epochs over both
  positions: HST ACS FR505N (2019-10-01) and JWST MIRI F2100W (2023-09-30).
  - FR505N fails T1: it lies blueward of Lyα at z 5.133 (`lya_obs_A` 7457.7 Å).
  - F2100W fails T2 in the deciding window and T3 (101.1 has no F2100W measurement); it is the only filter that
    enters the extended window.
- Every epoch from 2024-03-11 on maps to 101.1-equivalent dates after the window (T2 fails). Every 2012 epoch maps
  to dates before it.

## Reading
- With this model's delays (1,238 / 1,321 d, ±20 % allowance, an ASSUMPTION), the leading pair's bright phase,
  as 101.1 shows it, fell in years without usable imaging. So the leading images can neither be confirmed nor
  found missing in public data. This is not a lens-model residual and tests nothing about the model.
- No future epoch can help either: in this model the leading images have already passed.

## Limits and deviations
- **Deviation:** the first run started after the commit chain stopped on a docstring line-length lint error, so it
  ran on the uncommitted tool. Rule functions and tests were unchanged; the tool was then committed and re-run. Only
  the committed run's `epochs.json` is reported, and its counts equal the first run's (not stored).
- One model; no Δt uncertainty from the model itself. Ground-based or non-MAST imaging was not searched.
- The rows include HAP skycell and HLA mosaic products, so the 138 rows are not independent epochs.
