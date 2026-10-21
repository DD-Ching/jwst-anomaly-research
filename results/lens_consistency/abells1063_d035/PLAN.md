# Abell S1063 D-035 revisit: CATS v4.1 image-plane rms with the Chan+2020 redshifts (frozen 2026-10-21 20:25Z)

**Question.** D-035 gated the `abells1063-cats` image list off (rms 11.8″ vs 0.48″ quoted) because `arcs.txt` carries
placeholder redshifts. With published system redshifts (Chan+2020, `../abells1063_images/chan20_images.ecsv`), does
the CATS v4.1 deflection map reproduce the images well enough to open the list for `images` (counter-image test)?

**Ordinary explanations of a failure:** wrong redshifts (Chan Geoz is another model's value), identification
differences between Chan and CATS, frame offset between the lists, map resolution near critical lines.

**Route** (`scripts/s1063_d035.py`): `lens_consistency.load_model("abells1063-cats")`; images = Chan list without
`dup_position` rows; `backtrace_images` + `imageplane_residuals` with the same solver grid step (0.25″) as `validate`.
Pre-check: nearest CATS `arcs.txt` image to each Chan image (median offset; a frame offset > 0.1″ is reported).

**Frozen rules** (functions in the script, unit-tested before the run):
- R1 (whole-list gate, D-035, per subset): rms(dtheta) over solved images ≤ 1.5 × 0.48″ = 0.72″ AND every image of
  the subset solved. Evaluated on subset S (systems with z_spec) and subset A (all systems).
- R2 (per-system opening, `systems_passing_rms`, 0.5″ ASSUMPTION): applied to S only (Geoz systems inherit
  another model's z, so A-only systems are not opened).
- Outcome: R1(S) passes → list opens for S; else the S systems passing R2 open; none → S1063 stays map only.
**What changes the plan:** an opened subset makes `images` on S1063 (the `arcs` anti-excess cluster) the next step.
