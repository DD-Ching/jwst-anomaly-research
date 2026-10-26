# D115b Lenstool stage F2-B: budget probe of Lenstool's optimiser (plan, frozen 2026-10-26 13:16Z, before any run; the first "13:27Z" stamp was a typo, corrected in the tool commit)

## Why a budget probe first
Stage F2 asks whether a re-fit with the cluster-scale potentials free can keep every catalogued image while member
8515's σ is low enough to remove the D115b extra image (RESULT_D14A.md, "Reading"). The CANUCS optimisation used
`inverse 3 0.015 10000` (`input.par`). One fixed-parameter χ² evaluation takes about 0.5 s wall in this session
(`result.json` `runtime_s` 0.5, Lenstool start-up included), so that budget is far beyond one routine run. Before
any F2 Δχ² is read, this stage measures what a short `inverse 3` run costs and returns, with 8515 **fixed** at its
`best.par` σ (f 1.00), so no D115b number is produced here.

## Hypothesis and what would change the plan
- H: a short BayeSys run (`inverse 3 0.5 20`, ASSUMPTION: the smallest budget that still writes `bayes.dat`)
  finishes within the wall cap of 900 s (ASSUMPTION) on 4 OpenMP threads.
- Ordinary expectation: BayeSys starts from prior samples, not from `best.par`, so a short run returns a χ² far
  above 344.42.
- Plan change: if the best χ² of the probe is within 1 % of 344.42 (gate V2, the same tolerance as gate V), a short
  `inverse 3` budget is a usable optimiser and F2 proceeds with it on the f grid. Otherwise F2 needs a local
  optimiser started from `best.par` (to be planned separately); the probe's runtime per sample sizes either route.

## Run file (stage V′ file plus)
- `runmode`: `inverse 3 0.5 20` in place of `inverse 0`; `chi2 1` kept.
- `grille`: `nlens_opt 6` (the first six `potentiel` blocks of `best.par`, O1–O6, are the six `input.par`
  potentials with `limit` blocks).
- `input.par`'s `limit 1`–`limit 6` blocks appended unchanged. The `potfile` stays fixed at `best.par`'s values
  (deviation from the CANUCS free set, stated: the probe measures cost, not the published posterior).
- 8515 fixed at σ 120.600977 km/s.

## Outputs and rule
- Recorded: exit status, wall time, number of samples in `bayes.dat`, the lowest χ² in `bayes.dat` (column read from
  its header), Lenstool's final `chires.dat` `chimul` and `chipenalty` if written.
- Gate V2 (function `gate_v2`, unit-tested before the run): passes iff the run finished within the cap and the
  lowest χ² ≤ 344.42 × 1.01. A stopped run is `undetermined`.
- All χ² are `model_prediction` of re-fitted variants of the pinned `macs0416-canucs` model (D-044); the image list
  is `observed`.
