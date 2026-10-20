# Abell 2744 system 700.1: fine redshift scan in Bergamini+23 (plan, frozen 2026-10-20 20:10Z before any scan)

**Question.** D-036 / docs/fields/abell2744.md: Bergamini+23 (`abell2744-bergamini23`, z free, best-fit 1.217) never
reproduces 700.1a/b/c as three separate images at the sampled z = 1, 1.5, 2, 3, 5 ("a finer z scan is open"). Is the
non-reproduction a fitted-redshift problem (some z reproduces the configuration), or does the model fail at every z?

**Ordinary explanations.** (1) wrong fitted z; (2) model resolution in the northern infall region (160″ from the
core, few constraints); (3) solver grid (fold pair 700.1a/b, 0.25″ grid; D-038).

**Method (model_prediction; image positions observed).** z = 1.00–4.00 in steps of 0.025. At each z: one source
(magnification-weighted least squares of the three back-traced positions, `a2744_uncover`'s `wls_source`),
`find_images` on a 0.25″ grid refined to 0.02″ next to critical curves, then a one-to-one assignment of 700.1a/b/c to predicted images (minimum total
separation). Recorded per z: number of images, the three separations, rms, parities of the assigned images.

**Frozen rule (ASSUMPTION thresholds).**
- *reproduced at z*: all three catalogued images get distinct predicted images within 1.0″ each (the D-037 addendum
  pair criterion; ≈ 2.3 × the model's 0.427″ rms).
- Outcome A "fitted-z artefact": at least one z is reproduced. Report the z range.
- Outcome B "not reproduced at any z": record the minimum rms and where the failure lies (a/b merged, c missing).
  The system stays "model-limited"; next step an independent model (UNCOVER v2.0 maps, D-041) at its own z.
- Neither outcome is a candidate on its own: B is a model residual in a weakly constrained region.

**What would change the plan.** Outcome B together with UNCOVER also failing at its own z would make 700.1 a
two-model lens-model residual worth `/vet-candidate` (are 700.1a/b/c one source? colours, photo-z).
