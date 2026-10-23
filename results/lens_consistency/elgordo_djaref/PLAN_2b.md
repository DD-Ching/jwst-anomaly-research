# El Gordo reference-catalogue extra-image test: stage 2b (plan, frozen)

Frozen 2026-10-23 after stage 2 stopped at rule 2a (D-034 addendum p), before any DJA F200W flux row, mosaic pixel,
frame table or stage-2 prediction was read. Script: `scripts/elgordo_djaref.py --stage 2b`.

**Disclosed before freezing.** The outcome of the new rule 2a is already known: the stage-2 post hoc check
(`summary.json` `stage2.posthoc_rule_2a`) gives 23.9 mag in every row. The rule decides only the flux unit, not
any class or flag.

**Deviation from `PLAN_2.md`, only this:** rule 2a uses u = median of `mag_auto + 2.5 log10(flux_auto × tot_corr)`
over rows with finite positive `flux_auto` and `tot_corr` and finite `mag_auto` (`flux_unit_jy(t,
with_tot_corr=True)`, unit-tested); same thresholds (23.9 → µJy, 28.9 → 10 nJy, ± 0.05 mag, otherwise stop).
Every other rule, input, threshold and outcome of `PLAN_2.md` applies unchanged (rules 1, 2b, 2–7, data volume,
budget). Outputs: `summary.json` `stage2b`, the frame tables in `frame/`, the control / predicted tables in this
directory and `without_contain/`.
