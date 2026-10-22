# D115b vetting follow-up: χ² cost of a lower 8515 σ with a fitted source (frozen before the run)

Vetting record: `docs/candidates/macs0416-canucsref-D115b.md` (PR #458). Written 2026-10-22 07:30Z, before the
refit was run.

## What is already measured (stored in `validate.json` runs, summarised in `validate_runs.json`)
`lens_consistency.py --model macs0416-canucs validate --grid-step 0.1`: χ²_pos 330.81 (fiducial), 330.73 (8515 on
z 0.3907, σ as fitted), 5262.01 (8515 σ 84.4 km/s at z 0.396), 5261.44 (σ 84.4 at z 0.3907). D115a.1–3 share one
predicted image in both σ 84.4 runs. So the 0.25″ grid and the plane redshift do not cause the break.

## Hypothesis and ordinary explanation
- H (`hypothesis`): the break comes from `validate`'s source definition (mean back-traced position). Under σ × 0.7
  that point lies where D115a has a single image; a source position fitted to the image plane (as Lenstool's
  image-plane optimisation does) keeps three matched images.
- Alternative: no source position gives three D115a images near the catalogued ones at σ × 0.7, i.e. the lower
  mass really costs the fit.

## Method (`scripts/m0416_d115_refit.py`)
- Models: fiducial (single plane); 8515 at z 0.3907 with σ as fitted, and with σ 84.4 km/s (two-plane, D-046).
  0.1″ deflection grid, half-width as `validate`.
- Families D115a, D115b (z 2.2999 as in the catalogue): for each model, (1) the `validate` source (mean
  back-traced), (2) a fitted source: Nelder–Mead on χ²_fam(β) = Σ (Δθ_i / σ_i)², σ_i the model's position errors,
  Δθ_i from a one-to-one assignment (Hungarian) of catalogued to predicted images, falling back to the nearest
  predicted image when there are fewer predicted than catalogued. Starts: the mean and each image's own
  back-traced position; the best end point is kept.
- Total χ²_fit(model) = the model's `validate` χ²_pos minus the D115a and D115b terms of that run plus the two
  fitted-source terms. The other families keep `validate`'s source (not refitted, in every model alike).
- Reported per family: number of predicted images, Δθ of each catalogued image, and whether a predicted image lies
  within 0.7″ of the flag position (64.048551, −24.064825) with |μ| ≥ 5.

## Decision rule (computable, function `verdict` in the script, tested)
- ASSUMPTION: Δχ² ≤ 6 counts as an acceptable cost (the MACS0416 system-51 precedent, 330.81 → 336.44).
- Δχ²_fit = χ²_fit(σ 84.4, z 0.3907) − χ²_fit(σ fitted, z 0.3907).
- `ordinary` if Δχ²_fit ≤ 6 **and** the σ 84.4 model has no predicted D115b image within 0.7″ of the flag
  **and** every D115a and D115b catalogued image is matched one-to-one (no shared match).
- `inconclusive` otherwise (then the next step is a σ × 0.75–0.85 scan or a full re-fit, outside this plan).
