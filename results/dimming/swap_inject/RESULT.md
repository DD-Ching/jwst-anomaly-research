# W3 dimming: the D-039 persistence swap in the injection chain (2026-10-20)

Run 2026-10-20 04:12–04:22 UTC (re-run once for strict JSON, identical numbers) under the frozen `PLAN.md` (one wording amendment at 04:13, before any output:
the unknown-pair criterion in "What would change the plan" was made computable). Screens rebuilt with the
production chain on `main` (weight-map coverage term on; catalogues from S3 by sha256, Gaia from the ESA TAP), the
`forced` inflation from each field's committed `inject_summary.json`. Numbers: `runs/<field>/summary.json`,
`runs/summary.json`, `runs/limits_combined_<variant>.ecsv`, `runs/<field>/targets_check.json`. Injections are
simulated; limits are derived; ε by exposure count is an ASSUMPTION (`expcount_eff/PLAN.md`).

Commands (repo root; `$S` a scratch directory):
```
results/dimming/persistence_run/fetch_cats.py macs0416 abell2744 nexus
results/dimming/persistence_scope/gaia_esa.py $S macs0416 abell2744 nexus
results/dimming/swap_inject/run.py $S results/dimming/swap_inject/runs <field>     # per field
results/dimming/swap_inject/run.py combine results/dimming/swap_inject/runs
results/dimming/swap_inject/check_targets.py $S $S/cache <field>
```

## Gates
- **G0** (swap off reproduces `wht_coverage/runs/<field>/limits_veto_on_on.ecsv`, rtol 1e-9): passes in all three
  fields and for the joint table (`gate_G0`, `gate_G0_combined`).
- **G1** (veto-on ≤ swap ≤ veto-off per cell): passes for every variant in every field (`gate_G1_lower`,
  `gate_G1_upper`).

## Coverage of ε
Monitored detected (source, epoch) pairs with a committed `n_exposures`: MACS0416 94 / 94, Abell 2744 106 / 106,
NEXUS 234 / 350; the 116 NEXUS pairs without a count are all in epoch 2 (`o014`, controls never run), and take ε = 0
in `swap`. Every MACS0416 and Abell 2744 pair has > 4 exposures, so `swap_pooled` changes NEXUS only.

## Joint limits (calibrated fields MACS0416, Abell 2744, NEXUS; ρ = 0.1; `limits_combined_<variant>.ecsv`)

| t_E (yr) | 0.01 | 0.03 | 0.1 | 0.3 | 1 | 3 |
|---|---|---|---|---|---|---|
| exposure, veto on (= current headline), source-yr | 8.36 | 21.4 | 53.1 | 87.9 | 122 | 144 |
| exposure, `swap` | 12.9 | 38.0 | 106 | 221 | 280 | 294 |
| gain `swap` / veto on | ×1.54 | ×1.78 | ×1.99 | ×2.51 | ×2.29 | ×2.05 |
| gain `swap_pooled` | ×1.54 | ×1.77 | ×1.99 | ×2.50 | ×2.28 | ×2.04 |
| gain `swap_unknown_low` | ×1.62 | ×1.87 | ×2.13 | ×2.71 | ×2.44 | ×2.16 |
| gain veto off (upper bound) | ×1.80 | ×2.12 | ×2.54 | ×3.37 | ×3.24 | ×2.85 |
| 95 % limit per source per yr, `swap` | 0.233 | 0.079 | 0.0284 | 0.0136 | 0.0107 | 0.0102 |

Per field (ρ = 0.1, `runs/<field>/summary.json`): MACS0416 ×1.04–1.47, NEXUS ×1.19–1.84; Abell 2744 goes from
0–1.37 to 2.83–93.3 source-yr (its veto-on recoveries are almost all removed by `single_epoch`, as in
`persistence_scope`).

## Zero surviving events still holds
The real flags the swap would clear are the compact single-epoch `vanish` flags of the production screen. With the
coverage term on there are 11 (MACS0416), 5 (Abell 2744) and 11 (NEXUS); every one was tested by the D-039 run
(`not_tested` empty in `targets_check.json`). Their `on_sky` members (MACS0416 11, Abell 2744 5, NEXUS 8) were all
explained by `forced` + depth / footprint tests (`forced20/RESULT.md`, `nexus_o014/RESULT.md`: 0 survivors); the
NEXUS `inconclusive` (2) and `undetected` (1) flags stay vetoed, as the swap's ε only credits `on_sky`.

## Reading (frozen rules)
- Gain ≥ 10 % at every t_E: the swap is worth carrying. `inject --persistence-verdicts <verdicts.ecsv>` now runs
  it in production; `swap` is the corrected headline, with `swap_pooled` as the sensitivity case (≤ 1 % lower).
- `swap_unknown_low` over `swap`: ×1.05–1.08 at ρ = 0.1 (`runs/summary.json` `unknown_low_over_swap_rho0.1`), below the
  frozen 10 %, so NEXUS e2 controls are not the next step on this criterion.

## Limitations
- ε is a class average by exposure count (two classes), not per copy: an injected copy's surviving epoch may be
  fainter than the S/N ≥ 10 controls the class was measured on (post hoc; not quantified).
- The > 4 class (0.961) is not tested on held-out NEXUS data (`expcount_eff/RESULT.md`).
- The `forced` confirmation stays an error-inflation proxy, as before; the NEXUS `o014` forced run had controls
  8 / 300 flagged, above its frozen 1 % (`nexus_o014/RESULT.md`).

## Next
1. Look at the 11 held-out NEXUS `undetected` control stamps (`expcount_eff/RESULT.md` Next 2).
2. NEXUS e2 controls in a run of their own (tests the > 4 class on NEXUS; would credit the 116 pairs).
3. New fields or epochs (JADES, new NEXUS epochs) now pay double through the swap.
