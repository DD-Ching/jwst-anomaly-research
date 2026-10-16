# Status: paused, unfinished (2026-10-16 06:40 UTC)

The owner paused this unit for a process-only round. **Nothing here is a vetting conclusion yet.**

## Done
- **Plan.** `PLAN.md` was committed before measuring.
- **Measurements.** `scripts/orphan_vet.py` measured E1–E4 and E6 detection on the LW stack and 7 bands of the
  11 pairs; the values are in `measurements.json`.
- **Post-run tool fixes** (same plan; deviation):
  - The Monte Carlo error patches are now restricted to source-free sky (|pixel| < 4 × sky std). The first run's
    patches contained galaxies and inflated the E2 errors.
  - The size estimator was checked on synthetic Gaussians. Its deconvolved FWHM is 0.16″ for a true 0.20″ and 0.27″
    for a true 0.30″, and ≤ 0.03″ for a PSF-only source.

## Open
- **E2 validation.** The surface-brightness error model must be checked on known equal-surface-brightness pairs
  (e.g. published multiple images) before any E2 verdict. The current E2 errors (7–72 %) are not validated.
- **E5 spectroscopy.** VizieR cone searches not run: J/ApJ/982/51 (Price+25), J/ApJS/270/12 (Wang+24).
- **E6 neighbours.** Not run: J/ApJS/270/7 (UNCOVER DR3).
- **Records.** The classification table and the CHANGELOG, TASKS and DECISIONS entries are not written.
- **#7.** Member a is below 3σ in every band of our photometry; this needs checking against the DJA detection.
