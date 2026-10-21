# Image-plane source fit: route passes on MACS1149; Abell S1063 stays map only (2026-10-21)

Frozen plan: `PLAN.md` (21:10Z, before the script ran). Script: `scripts/imageplane_fit.py`; numbers: `summary.json`,
per-image tables `macs1149_fit.ecsv`, `abells1063_fit.ecsv` (provenance `model_prediction`).

| Gate | Rule (ASSUMPTION) | Measured | Outcome |
|---|---|---|---|
| V1 synthetic, S1063 map | rms ≤ 0.05″ in ≥ 9 / 10 | 10 / 10 (rms 3–5 × 10⁻⁴″) | pass |
| V2 MACS1149 CATS `arcs.txt` | all 145 solved, rms ≤ 0.945″ | 145 / 145, rms 0.657″ (quoted 0.63″) | pass |
| R1 S1063 z_spec (S) | all solved, rms ≤ 0.72″ | 36 / 36, rms 0.978″ | fail |
| R1 S1063 all (A) | as S | 108 / 117 paired, rms 1.677″ | fail |
| R2 S1063 z_spec per system | rms ≤ 0.5″ | 0.633–1.460″ over the 12 systems | none open |

**Frozen outcome: map only.** S1063's counter-image test stays off.

**What this shows.**
- The fit removes the route failure: the five z_spec systems whose mean traced source gave 14.2–36.7″ (one-to-one
  pairing; systems 1, 2, 10, 11, 52) now have every image paired at 0.633–1.321″ (`summary.json` `r.per_system`).
  The 18.5″ of `../abells1063_d035/` was a route artefact.
- With the route removed, the CATS v4.1 map still misfits the Chan+2020 z_spec images at 0.98″, about 2 × the
  quoted 0.48″, while the same route reproduces MACS1149 to 1.04 × quoted. Possible causes, not separated: the
  0.2″-pixel deflection map near critical curves, Lenstool's per-image σ weighting (quoted rms is an optimiser
  output), or differences between the CATS and Chan identifications (positions agree to a median 0.061″,
  `../abells1063_d035/summary.json`).
- Post hoc, not a gate: the one-to-one mean-source rms on MACS1149 is 0.920″ (validate's nearest-match measure was
  0.67″); D-035's MACS1149 decision does not change (the fit gives 0.657″).

**Limits.** Point-source images only; no per-image positional errors; Nelder–Mead from the mean and each traced
source (best of up to n_images + 1 starts, 200 evaluations each); 9 of the 117 S1063 images (all outside S) have no fitted partner.

**Next.** None on S1063 image positions from this map. A route with the Lenstool parameter file (not public for
S1063 v4.1: D-035) or another public S1063 model with its own image list would be needed.
