# Sunrise arc difference imaging — result

Box: 24″ × 24″ centred on Earendel (Scofield+2025 1.1e, 24.346822, −8.464511), which holds the Sunrise arc's
catalogued images. `arc_diff.py --kernel-px 4 --kernel-stamps --inject 100 --min-bands 2` (full command in
`o010_o052/summary.json`), run with `OMP_NUM_THREADS=1`. Every number is `derived`; thresholds are ASSUMPTIONs.
Sheets: `scripts/arc_diff_sheet.py` (epoch-1 image and S at ±10, 2″ stamps), looked at before writing this.

## Deviations from PLAN.md
- **Deviation:** validation used `--inject 100` (per-band point injections on pixels above the box's 90th
  percentile), not `--inject-grid`; the PLAN's "≥ 0.8 at 5× the box σ" bar was not evaluated as written. The
  injection figures below are per band and before the morphology step.
- **Deviation (postponed):** the o010 → o120 pair did not finish in this session (see below).
- **Post hoc:** "on the arc" is made computable as the distance to the nearest Scofield+2025 multiple image
  (`o010_o052/image_distances.json`); the PLAN did not define it.

## o010 (2022-07-30, jwst 2.0.1) → VENUS o052 (2025-07-08, jwst 3.0.0), F150W + F444W
- Noise scales 1.20 / 1.29; flux scales 1.024 / 1.026 (`summary.json`).
- Injection recovery per band (`summary.json` `injection_recovered`): F150W 0.86 / 0.14 / 0.01 at
  27.5 / 28.0 / 28.5 AB; F444W 0.83 / 0.28 / 0.04.
- **9 hits** (6 `near_bright_core`; `hits.ecsv`). `arc_diff_morph.py` (defaults; `morph.json`): 2 dipole,
  2 sign_flip, 1 edge, 4 pass.
- **The 4 passes, by eye** (`o010_o052/sheet.jpg`):
  - h001, h002: mixed-sign residuals on the cores of a compact galaxy and a bright face-on galaxy (`near_bright_core`);
  - h005: negative residual on a compact galaxy's core in both bands (S −10.0 F150W, −19.0 F444W), 0.50–0.78″ from
    the map edge (`edge_arcsec`); in F444W a positive lobe of +19.9 sits beside it (dipole ratio 1.05; 0.36 in
    F150W, so the two-band dipole rule does not fire);
  - h007: elongated negative residual along an edge-on disk (F150W extent 1.59″).
- All 9 hits lie 3.98″ or more from every catalogued multiple image (`image_distances.json`; nearest h001 to 1.7a).
- **Outcome:** 0 hits within 3.98″ of a catalogued Sunrise-arc image, so 0 candidates on the arc's images in this
  pair. This is not a limit on the whole box below the injection depth.
- **Earendel (control row, not a hit):** S = −5.2 (F150W), −1.5 (F444W) at the reference position (`summary.json`
  `known`); the 2-band rule does not fire and the stamp shows no compact residual at ±10. The catalog forced
  photometry (recentred) gave Δm = +0.05 mag in F150W, 0.3σ calibrated (docs/fields/sunrise.md). The −5.2 is not
  explained here.

## o010 → o120 (2023-01-10, jwst 2.0.1), F277W + F356W (2026-10-16 run)
- **Rules.** They were fixed before the run (PLAN.md addenda 1–2).
- **Deviation:** a four-band attempt gave no output in 90 min, so this pair was run with the long-wavelength (LW)
  bands only. One worker held ~98 % CPU throughout, most likely on a short-wavelength (SW) band; the timing is not
  stored. SW-only changes on this pair are not covered. Files: `o010_o120/` (`summary.json`, `hits.ecsv`,
  `morph.json`, `image_distances.json`, `sheet.jpg`).
- **Run.** Command in `summary.json`, with every other setting as for o052.
  - Shifts: −0.29 / −0.30 px in x.
  - Flux scales: 1.013 / 1.010.
  - Noise scales: 1.03 / 1.14.
- **Validation (`injection_recovered`).**
  - F277W: 1.00 / 0.96 / 0.72 at 27.5 / 28.0 / 28.5 AB.
  - F356W: 0.97 / 0.95 / 0.74.
  - Both bands pass the ≥ 0.8 rule at 27.5 AB.
- **Hits.** There are 12 (`hits.ecsv`). `arc_diff_morph.py --all-hits` gives 10 rosette, 1 dipole, 1 edge and
  **0 pass**.
- **Distances.** All 12 lie 3.87″ or more from every Scofield+2025 image (`image_distances.json`). Under the frozen
  rule (≤ 1.0″), none is on the arc.
- **Looked at (`sheet.jpg`).** The rosettes are the expected ~180° PSF-rotation residuals around compact sources (h001,
  h007–h011) and along one edge-on disk (h003–h006). h000 is a compact-core dipole, and h002 sits on the map edge.
  - Earendel's control row has S = +2.2 (F277W) and +1.5 (F356W), with no compact residual on the stamp.
- **Outcome.** Under the frozen rules there are 0 candidates on the arc's images in this pair.
  - This covers changes brighter than about 27.5 AB in F277W and F356W, at the injection-recovery level above. It is
    not a hard completeness limit and says nothing about SW-only changes.
  - The run was exploratory: the field was examined before (catalogue grid, #257), but this pair was not differenced
    on the arc before.
