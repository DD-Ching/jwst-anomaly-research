# D115b vetting: σ(8515) scan with every family's source fitted (frozen before the run)

Vetting record: `docs/candidates/macs0416-canucsref-D115b.md`. Previous step: `../d115_refit/` (#458): at σ × 0.7
(84.4 km/s) the flagged image vanishes with D115a / D115b matched, but Δχ² +12.6 > 6, of which +10.6 (post hoc) sits
in families that kept `validate`'s mean back-traced source. Written 2026-10-22 09:25Z, before any scan run.

## Hypothesis and ordinary explanations
- H0 (`hypothesis`, ordinary): the D115b extra image is a lens-model systematic of member 8515's galaxy-scale mass; a
  σ inside the assumed ±30 % scatter removes it at a χ² cost ≤ the system-51 precedent once every source is fitted.
- H1: every σ that removes the image costs more than that; the extra-image prediction then stays a model tension
  (still not evidence of anything unusual: a full re-fit with other parameters free is untested).
- Result that changes the plan: H0 → D115b closed as ordinary, MACS0416 flag list empty; H1 → the record stays
  `inconclusive` and the Lenstool re-fit (owner machine) becomes the only next step.

## Method (`scripts/m0416_d115_sigma_scan.py`)
- Models: 8515 on its z_spec plane 0.3907 (D-046 two-plane), σ = f × σ_fit with f ∈ {1.00, 0.85, 0.80, 0.75, 0.70}
  (σ_fit 120.6 km/s, so 102.5 / 96.5 / 90.5 / 84.4 km/s). 0.1″ deflection grid, half-width as `validate`.
- **Every** family (not only D115a / D115b) gets a source fitted to the image plane exactly as `../d115_refit`
  (Nelder–Mead on χ²_fam, one-to-one Hungarian assignment, 1e4 penalty when there are fewer predicted than
  catalogued images; starts: the mean and each image's own back-traced source; best kept). Redshifts: `z_used` of
  the back-trace (as `validate`).
- χ²_all(f) = Σ over families of the fitted-source χ²_fam. Δχ²(f) = χ²_all(f) − χ²_all(1.00), summed over the
  families that are matched one-to-one in both models.
- Flag image: a predicted D115b image within 0.7″ of (64.048551, −24.064825) with |μ| ≥ 5 (as before).

## Gates (checked first; a failure stops the scan's interpretation)
- G1 reproduction: f = 0.70 must reproduce `../d115_refit` for D115a / D115b (family χ² within 0.05 of 2.765 / 0.538; no flag image).
- G2 construction: for f = 1.00, every family's fitted χ² ≤ its `validate`-source χ² + 1e-6 (the mean is a start).

## Decision rule (function `verdict`, tested)
- f* = the largest f < 1 at which the flag image is absent and D115a / D115b are matched one-to-one.
- `ordinary` if f* exists, Δχ²(f*) ≤ 6 (ASSUMPTION, system-51 precedent 330.81 → 336.44), and no family that is
  matched one-to-one at f = 1.00 loses that match at f*.
- `inconclusive` otherwise. Reported for every f: Δχ², flag outcome, families losing their match, top-5 family
  changes. The other Lenstool parameters stay fixed (a limit of this test, not of the rule).
