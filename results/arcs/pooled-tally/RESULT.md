# `arcs`: pooled anti count over the S1063 replication fields (exploratory, post hoc)

Plan: `PLAN.md` (written before the computation, after the per-field counts had been read; hence exploratory).
Script: `scripts/arcs_pooled.py`; every number below is in `summary.json`. No new data: each field's stored
deciding-run `arcs.ecsv` is re-read with the cuts from its `arcs.json`; the script stops unless the recomputed
n, anti and expected (to the stored 3 decimals) match each run's stored shape-noise summary (`stored_check`).

## Frozen outcome: **pooled excess** (derived)
| Variant | Fields | Anti | Expected (σ_ε 0.26) | P(≥ anti) |
|---|---|---|---|---|
| A | MACS0417, MACS1423, WHL0137, RXC J0600, MACS0308, RXC J0018 | 10 | 3.222 | 0.0019 |
| B | A without MACS1423 | 6 | 2.273 | 0.029 |

Both variants are below 0.05 (ASSUMPTION threshold), so the frozen label is "pooled excess". Each field alone
"does not replicate" (except MACS1423 at screen level, explained by its clump-pair census); together they hold
2.6–3.1 times the shape-noise expectation.

## Secondary checks (pre-registered; they do not change the outcome)
- **σ_ε scan** (ASSUMPTION grid 0.26–0.40): variant B reaches P ≥ 0.05 at σ_ε 0.28 (P 0.064), variant A at
  0.32 (P 0.051). The pooled excess is about the size of a modest error in the null's intrinsic ellipticity
  dispersion. (σ_ε is not measured on these catalogues; that check is not done.)
- **Abell S1063 at each σ_ε** (from its stored `arcs.ecsv`): 4 anti vs 0.52–1.14 expected, P 0.0021–0.029. Inflating
  σ_ε alone does not remove the S1063 excess anywhere on the grid.
- **g_min split of the pooled anti:** 7 at g_min < 1 and 3 at g_min ≥ 1 (MACS0417 2, MACS1423 1); variant B: 4 and 2.
  Most pooled anti rows lie at g_min < 1, outside the critical curves.

## Post hoc: S1063 against the pooled fields' empirical rate (not in the plan)
`s1063_vs_pooled_posthoc`: scaling S1063's expectation by the pooled observed / expected ratio (B 2.64, A 3.10; 95 %
Poisson ranges 0.97–5.75 and 1.49–5.71) gives 1.38 / 1.63 expected and P(≥ 4) 0.052 / 0.083 at the point value;
over the ratio's 95 % range P runs from 0.002 to 0.36, so the scaled test cannot decide. The conditional test (given
all anti in S1063 plus the pool, the chance that S1063 holds ≥ 4 at its expected share) gives P 0.100 (B) / 0.121
(A) at σ_ε 0.26 and 0.078–0.121 over the σ_ε grid (`sigma_scan`: σ_ε does not scale every field's λ alike).

## Reading (narrowest statement)
- The shape-noise null at σ_ε 0.26 under-predicts `anti` counts across the independent fields by a factor of about
  2.6–3.1 (derived; 6 and 10 counts, so the factor is poorly known). The likelier explanations are ordinary
  (null calibration: σ_ε, blends and clump pairs, photo-z contamination, model error); which one is not tested here.
- Measured against that empirical rate, S1063's 4 anti are not below 0.05 in the conditional test (P 0.078–0.121,
  post hoc). The S1063 anti excess is therefore weaker than its frozen P 0.0021 suggested. It is not explained by this check either:
  1316 and 2184 stay open as in `results/arcs/abells1063-cats/RESULT.md`.

## Limitations
- Exploratory: the plan was written after the per-field numbers were seen.
- Counts are small (6–10); the ratio's 95 % range is wide (about 1–5.7).
- MACS0417's deciding model differs from its HST-era model, and the fields mix GLAFIC, Lenstool and CANUCS maps,
  photo-z sources and depths; the pool treats the ratio as common to all.
- The tuning fields (Abell 2744, El Gordo, MACS1149, MACS0416, Abell 370) are left out by the plan; they could
  calibrate the null further, but they were used to choose `arcs` rules.
