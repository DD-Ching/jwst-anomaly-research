# Sunrise arc difference imaging (WHL0137−08), o010 → o120 and o010 → VENUS o052 — plan (written before the run)

**Why.** The Sunrise catalog grids (D-027 addenda 2026-10-15, -15i) are null, but a catalog grid misses a transient
blended into an arc (PLATE-23a on Abell 2744, D-027 addendum 2026-10-15j). `scripts/arc_diff.py` has never been run
on the Sunrise arc, the most magnified arc in the sample (Earendel, μ ≳ 1000 near the critical curve).

**Hypothesis (to test).** No point-like flux change on the Sunrise arc beyond what instrument, processing and
known populations explain. Ordinary explanations to test first: PSF/kernel residuals on bright arc knots,
detector persistence (D-039; `persistence_check.py`), spikes, pipeline-version differences (o052 is jwst 3.0.0),
registration dipoles, and microlensing of a caustic-straddling star (an ordinary known population near critical
curves: a single point brightening near the critical curve is expected and is *not* exotic).

**Validation first.** `arc_diff.py --inject-grid` on the same box and pairs; adopt the box only if injected
point brightenings at the arc's surface brightness are recovered as in MACS0416 (≥ 0.8 at 5× the box σ;
ASSUMPTION).

**What changes the plan.** A hit passing the morphology veto (`arc_diff_morph.py`) on the arc in ≥ 1 pair gets
`/vet-candidate`; a hit at a non-Earendel position on the critical curve is a caustic-transient candidate to vet
against microlensing rates, not a lensing-model anomaly. A clean null gives a depth limit on the arc.

## Addendum 2026-10-16: o010 → o120 run, rules frozen before the run
- **Run.** The settings are the o052 run's, unchanged: box `--ra 24.346822 --dec -8.464511 --half-arcsec 12`,
  `--kernel-px 4 --kernel-stamps --inject 100 --min-bands 2`, and the same `known_controls.csv`.
  - Bands: F090W, F115W, F277W and F356W, i.e. every band o120 has.
  - Execution: `--procs 4` with `OMP_NUM_THREADS=1`.
  - The same `--inject 100` validation as o052 (the deviation recorded in RESULT.md).
- **Validation (per band, as reported by `injection_recovered`).** A band's hits are read only if recovery at 27.5 AB
  is ≥ 0.8. Bands below that are reported as shallow and their hits are listed but not interpreted.
- **On the arc (fixed now).** A hit is on the arc if it lies within 1.0″ of a Scofield+2025 multiple image (the same
  image list as `o010_o052/image_distances.json`).
- **Candidate rule.** A hit becomes an ordinary-vetting case if all three hold:
  - it passes `arc_diff_morph.py` (defaults);
  - it is on the arc;
  - it is detected in ≥ 2 bands (`--min-bands 2`).
- **Vetting order.** The first checks are persistence (`persistence_check.py`), the PSF rotation (~180° between the
  epochs, so rosettes are expected), registration dipoles, and caustic-crossing-star microlensing. A single point
  brightening near the critical curve is ordinary until those are excluded.
- **Status of the data.** o010 was used in the o052 difference; o120 was used only in the catalogue grid (#257). This
  run is the first arc difference of this pair, but the field has been examined, so it is not a blind sample.

## Addendum 2 2026-10-16 10:25 UTC (**deviation**, written before any o120 output existed)
- **What happened.** The four-band o010 → o120 run (`--procs 4`, `OMP_NUM_THREADS=1`) gave no output in 90 min. One
  worker held ~98 % CPU throughout, most likely on a short-wavelength (SW) band. This is the second unfinished SW
  attempt on this pair (the first is in RESULT.md); it was stopped, and the timing is not stored.
- **What changes.** The pair is run with **F277W + F356W only**, with every other setting unchanged.
  - `--min-bands 2` therefore requires both long-wavelength (LW) bands.
  - The validation, on-arc and candidate rules of addendum 1 are unchanged.
- **What is lost.** SW-only (blue) changes on this pair are not covered.
- **Not done.** No change to `arc_diff.py`: tool work is out of scope (owner, 2026-10-16).
