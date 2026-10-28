# wh_twoway — finite-energy return pulses through a collapsing Ellis–Bronnikov wormhole

The canonical record (question, fixed problem, deviations, checks, results, limits) is
[`docs/results/wh_twoway.md`](../../docs/results/wh_twoway.md). All results are simulations [sim] of one model and
one configuration; nothing here is an observation.

| File | Content | Produced by (`python scripts/wh_twoway.py …`) |
|---|---|---|
| `checks.json` | static-background, Kain26 Fig. 9 benchmark, GGS09 growth-rate and resolution-noise checks | `checks` |
| `scan_n400.json`, `scan_n800.json` | 6 return energies × 9 emission times, per-run station energies, times, diagnostics | `scan --n 400` / `scan --n 800` |
| `converge_n1600.json` | the 12 boundary runs at n = 1600 | `scan --n 1600 --e-ret 1e-8 0.1 0.3 0.5 --u1 9.5 9.75 10.0 --out converge_n1600.json` |
| `validity_recheck.json` | one n = 400 run re-scored with the deadline-limited station check | `scan --n 400 --e-ret 1e-3 --u1 9.75 --workers 1 --out validity_recheck.json` |
| `boundary.json` | 50 % boundaries per n and energy, Richardson, deadline sensitivity, pass/fail case table | `analyze` |
| `static_return.json` | energy balance of a return pulse in the unperturbed wormhole | `scan --n 400 800 --e-out 0 --e-ret 1e-8 --u1 3.6 --workers 2 --out static_return.json` |
| `imprint.json` | ghost-field imprint of the pre-loaded return pulse at B | `imprint` |
| `fig_*.png` | spacetime diagrams, received fraction vs. emission time, benchmark | `figures` |
