# Grid G1 failures on MACS J0308 / J2129: recovery by epoch-1 exposure count (result, 2026-10-27)

Plan `PLAN.md` (pushed before the changed tool read any pixel). Tool `scripts/transient_grid.py` with the new
`injections.ecsv` output and `--coverage` flag. Numbers: `m2129/summary.json`, `m0308/summary.json` (frozen split
`by_e1_n_exp_min`), `m*/injections.ecsv` (one row per injection), post hoc `posthoc.json`
(`scripts/grid_g1_coverage_posthoc.py`, tested in `tests/test_grid_g1_coverage_posthoc.py`). Inputs are the eight
level-3 mosaics and catalogues of #576 / #577 (observed). We streamed them from S3 into the gitignored data root
and deleted them after the run. Every number below is derived. Thresholds are the D-027 ASSUMPTIONs.

## Reproduction check (passes)
Same seed, positions and thresholds as the frozen runs: overall −0.75 mag recovery 0.518 (J0308) and 0.78 (J2129),
27 / 16 hits in 27 / 13 groups. These are identical to #576 / #577.

## Frozen result: outcome U (a subset has < 100 injections)
The plan gives U when any subset has fewer than 100 injections. The single-exposure subsets hold 67 (J2129) and 50
(J0308), so the frozen outcome is **U**, although the reproduction check passes. The S / R criterion reads only the
≥ 2-exposure subset (433 / 450 injections), and there both grids fall short of 0.85. Read alone, that subset would
give R. No injection falls on zero epoch-1 exposures (`by_e1_n_exp_min` "0": n 0). Recovery at −0.75 mag (frac 1.0),
split by the smallest epoch-1 exposure count over the two bands:

| Grid | e1 n_exp = 1 | e1 n_exp ≥ 2 | Frozen bar |
|---|---|---|---|
| MACS J2129 (F277W / F444W) | 0.612 (67) | 0.806 (433) | 0.85 |
| MACS J0308 (F300M / F410M) | 0.320 (50) | 0.540 (450) | 0.85 |

Single-exposure coverage lowers recovery, but restricting to ≥ 2 exposures does not reach 0.85 in either grid.
Single-exposure strips alone therefore do not account for either G1 failure. There is no rule change: the plan
allows one only under S.

## Post hoc (not in the plan; `posthoc.json`)
- Splitting further by what set the injected flux (`nominal`: the fraction; `floor`: 10 × the band's median
  epoch-2 aperture error, which takes over at faint positions) separates the misses:

  | Grid | nominal, n ≥ 2 | nominal, n = 1 | floor, n ≥ 2 | floor, n = 1 |
  |---|---|---|---|---|
  | J2129 | 0.947 (304) | 0.800 (45) | 0.473 (129) | 0.227 (22) |
  | J0308 | 0.785 (297) | 0.516 (31) | 0.065 (153) | 0.000 (19) |

- The run's floor flux (the smallest injected flux per band, `floor_flux`) at the median aperture errors (epoch 1
  on the epoch-2 zero point) reaches a calibrated significance of 5.53 / 5.61 (J2129 F277W / F444W) and 3.68 / 3.74
  (J0308 F300M / F410M), with no systematic floor. The threshold is 5. On J2129 the floor therefore sits just above
  the threshold at a median position, so positions with larger errors fall below it. On J0308 it sits below the
  threshold even at the median. This matches the low floor-subset recovery (0.473 J2129, 0.065 J0308 on ≥ 2
  exposures). The floor subset is 151 / 500 (J2129) and 172 / 500 (J0308) injections. On J2129 the nominal
  injections on ≥ 2 exposures pass 0.85. On J0308 they do not (0.785), so that pair falls short even where the
  fraction sets the flux.
- Reading: G1 at "frac 1.0" mixes two quantities. One is a −0.75 mag brightening of a source. The other is a
  fixed-flux floor injection, whose significance depends on the noise scale (1.25–1.57 over these four bands) and
  on both epochs' errors. This is a reading of these two grids only. It is not tested on the grids that passed G1,
  and per-position errors are not compared.

## Limits
- Two grids only. The floor-significance numbers use median errors, not each position's error.
- Changing G1 (e.g. scoring only nominal injections, or setting the floor in calibrated σ) would be a D-027 rule
  change. This PR does not make it.
