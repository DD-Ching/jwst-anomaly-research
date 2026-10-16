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
