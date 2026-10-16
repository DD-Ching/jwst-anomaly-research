# No-injection control at the positions of vetoed SMACS injections (plan, written before any measurement)

**Why.** `results/transients/review_list/RESULT.md`: off the arc mask, morphology vetoes most recovered 27.5 AB
injections in the SMACS ERO → VENUS box (8 rosette, 5 dipole, 3 extended of 29 hits). The cause was not tested.
This is a measurement only; no rule, threshold or veto is changed (owner, 2026-10-16).

**Hypotheses (each a hypothesis).**
- H1 (neighbour residuals): the opposite-sign lobes that veto an injected star already exist in the
  no-injection difference map at that position (residuals of neighbouring sources).
- H2 (injection-induced): they are absent without the injection, so the injection itself makes them (PSF-matching
  kernel applied to the injected star in the convolved epoch, or an injected-PSF shape mismatch).

**Measurement.** Re-run the stored no-injection command of `results/transients/arc_chain/runs/smacs/summary.json`
with `--maps`, check that `hits.ecsv` is byte-identical to the stored one, then run `arc_diff_morph.py --known` on that
map at every injected position of the 3 stored seeds that was a hit and not `near_bright_core`, default parameters.
Per position, record whether an opposite-sign (negative) region with |S| ≥ 3 and ≥ 4 px lies within 0.8″ (the lobe
definition of the morphology code), and `extent_arcsec`.

**Decision rule (fixed now).** Compare vetoed (rosette + dipole) with passed injected hits:
- H1 supported if ≥ 70 % of the vetoed positions have a pre-existing negative lobe and the vetoed fraction exceeds
  the passed fraction by ≥ 30 percentage points.
- H2 supported if ≤ 30 % of vetoed positions have one.
- Otherwise mixed / undetermined. Thresholds are ASSUMPTIONs; small samples (≈ 13 vetoed per magnitude) are stated.

Outcome changes the plan as follows: H1 → the off-arc efficiency loss is a crowding property of the box (any
fix is a morphology-rule change, owner's call); H2 → an injection-procedure defect, so the R2 number itself is
suspect and the injection code is the next thing to inspect.
