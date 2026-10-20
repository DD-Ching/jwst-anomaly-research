# Abell 2744 system 700.1 in the UNCOVER v2.0 maps (result)

Plan: `PLAN.md` (frozen and committed before the run). Script: `scripts/a2744_s700.py --model uncover2`. Numbers:
`summary.json` (`known_cases` per radius with `failed_detail`, `outcome`, `uncover_counterparts`). Computed positions
are `model_prediction`; catalogued positions and redshifts are `observed` (catalogue values).

**Frozen outcome: "tool not validated"; no z scan was run.** On UNCOVER's own list (63 systems with z > 0 and ≥ 2
images, each at its own z) the one-to-one rule reproduces 31 / 63 (0.492) at 1.0″ and 44 / 63 (0.698) at 1.5″, both
below the 0.8 gate (ASSUMPTION). At 1.5″ no failure lacks images (every failing system gets ≥ as many predicted
images as catalogued). 9 of the 19 have every catalogued image within 1.5″ of *a* predicted image but two share one
(`n_failed_shared`). In the other 10, the nearest predicted image is 1.77–2.31″ away (`failed_detail`). So this
rule is too strict for a map model with a 0.84″ own-list rms (D-041): a shared predicted image, the #412 signature
of 700.1 in Bergamini+23, also occurs in known UNCOVER systems.

**Post hoc finding (catalogue look-up, not part of the frozen rule).** 700.1a/b/c are UNCOVER's images 78.21 /
78.22 / 78.23 (0.02″ each; `uncover_counterparts`). UNCOVER lists system 78 (three knots, 78.1–78.3) at z 2.3151,
`Redshift_type` spec, reference "Price+24", JWST/NIRSpec, as given in the v2.0 image table (the paper itself was
not checked here). In the known-case run, 78.2 at that z gets 3 predicted images, every catalogued image within
1.34″ of one, best one-to-one assignment 1.72″ (`failed_detail`; 78.1: 1.08″ / 1.63″; 78.3: 0.27″ / 1.74″).
UNCOVER uses system 78 as a constraint, so this is a fit, not an independent prediction.

**Reading.** 700.1 has a spectroscopic redshift (2.3151 per the UNCOVER table), so the redshift is no longer free.
At z 2.30–2.325 Bergamini+23 fails as at every z (rms 10.86–10.90″; `../a2744_s700_zscan/summary.json`). The
newer UNCOVER model, constrained by this system, places separate images on a and b within ≈ 1.3–1.7″. The most
economical reading: Bergamini+23 lacks this constraint in the infall region (ordinary model limitation). The
one-image reading of the a/b feature is not required by UNCOVER's model. Not a candidate; nothing here bears on
exotic lensing.

**Deviations (labelled).** Post hoc: the per-system `failed_detail` / `n_failed_shared` fields and the
`uncover_counterparts` look-up were added after the first run; the gate numbers are identical across the three runs.
The plan's timestamp line was corrected to the commit time (21:12Z) before the run's results were written.

**Limitations.** One rule (one-to-one, fixed radii). UNCOVER's own-list fit to 78 is not independent. The
Price+24 redshift was taken from the UNCOVER table, not from the paper or the spectrum.
