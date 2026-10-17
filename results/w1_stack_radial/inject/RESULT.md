# W1 lenses injected at the `radial` centres through the stack: result (2026-10-17)

Plan: `PLAN.md` (committed before any run). Code: `scripts/w1_stack_inject.py` (painting through
`inject_shear.ShearInjector.lensed_sample`, the D-050/D-053 chain). Every number below is in `summary.json`
(`stack`, `stack.cells`, `fields[*].injections`). Labels: injected shapes `simulated`; S, stacks, shares `derived`;
thresholds ASSUMPTIONs.

## Gates (as planned)
- **Z passed:** θ_E = 0 reproduces S_obs at every valid centre (max |ΔS| 0 in both fields).
- **R passed:** at 2 × 10¹³ M☉, f = 1, the stack responds (frozen median 3.72 > null 95th percentile 0.43).
- Sample check: 19 valid centres (Abell S1063 5, MACS0417 14), observed mean S 0.5095 as in
  `../vet_no_radial_arcs/summary.json`; rotation null redrawn (seed 316): p 0.025, 95th percentile 0.43.

## Frozen rule: failed approach (degenerate)
The planned rule adds each injection to the *observed* shapes and asks whether ≥ 95 % of subset draws exceed the
observed mean. The noise realisation is then fixed, so any positive S shift "excludes": 2 × 10¹⁰ M☉ (the grid's smallest mass;
0–2 lensed sources per centre, `fields[*].injections`) is "excluded" at f ≥ 0.5 with a stack change of only 0.51 → 0.53. The rule measures
whether a lens raises S at these centres, not whether the data could tell; its cells are in `summary.json`
(`excluded`, `share_above_observed`) but **no limit is quoted from them**.

## Noise-aware version (**post hoc**, deviation from PLAN)
Each centre's injected change ΔS = S_inj − S_obs is added to one rotation-null realisation per draw (2000 draws at
f = 1, 1000 subsets otherwise; `stack.cells[*].noise_aware`). ASSUMPTION: ΔS is additive on the null (the painted
images keep their sources' ε plus the lens term). Same exclusion level (≥ 95 % of draws above the observed 0.51).
Centres whose injected S is undefined (too few sources left in the aperture: 2 centres at 2 × 10¹², 1 at 8 × 10¹²,
all MACS0417) drop out of the mean, as in the real stack (post hoc detail; also applied to the frozen cells).

| M (M☉) | f = 0.1 | 0.25 | 0.5 | 1.0 |
|---|---|---|---|---|
| 2 × 10¹⁰ | 0.025 | 0.026 | 0.027 | 0.028 |
| 2 × 10¹¹ | 0.037 | 0.062 | 0.121 | 0.282 |
| 2 × 10¹² | 0.092 | 0.320 | 0.781 | **1.000** |
| 8 × 10¹² | 0.195 | 0.632 | **0.977** | **1.000** |
| 2 × 10¹³ | 0.296 | 0.820 | **0.999** | **1.000** |

Share of draws above the observed mean; **bold** = excluded at the 95 % level (ASSUMPTION).

## Conclusion
Under the post hoc noise-aware rule, W1 lenses of ≥ 2 × 10¹² M☉ at all 19 centres, or of ≥ 8 × 10¹² M☉ at half of
them, would give a larger stack than observed in ≥ 95 % of draws; those cases are excluded at these centres
(conditional on the D-050/D-053 injection model and the additivity ASSUMPTION). At ≤ 2 × 10¹¹ M☉ the stack cannot
tell (efficiency against the null 95th percentile ≤ 0.41 even at f = 1), so the observed excess (p 0.025) is neither
explained nor excluded by such lenses. Not a candidate; no new flag, so no sheet was needed.
The exclusion at f = 1 between 2 × 10¹¹ and 2 × 10¹² M☉ is not resolved by the mass grid.

## Limits of this result
- One lens changes only its own centre's S (PLAN ASSUMPTION, conservative for MACS0417's overlapping apertures).
- Exploratory sample: the centres were chosen by the `radial` screen on the same catalogues; no held-out data.
- The cluster-model-error explanation stays open outside the RELICS area of MACS0417 (`../RESULT.md` addendum).
