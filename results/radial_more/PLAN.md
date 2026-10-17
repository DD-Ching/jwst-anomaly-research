# Radial (W1 exotic-lens) discovery screen on the shear-only fields: plan (2026-10-17, frozen before any run)

**Question.** Do MACS0417, MACS1423, WHL0137, RXJ2129 and MACS0647 (W1 shear-screened, never `radial`-screened)
show anti-tangential arcs converging on a dark centre (docs/exotic_lensing.md W1; D-031/D-049)?

**Hypothesis (exotic, not assumed).** A compact negative-mass lens produces radially stretched, demagnified image
pairs; `exotic_screens.py radial` flags convergence centres of ≥ 3 `anti` arcs with p_random < 0.05.

**Ordinary explanations to test first for any flag.** Diffraction spikes and streaks (Gaia `--spike-stars`; the
W1 streak masks of `results/w1_streak_veto`), cluster radial arcs (model-radial veto), shape noise (null draws),
lens-model errors (a model that fails the `arcs` alignment control gives no result), edge and blend artefacts.

**Inputs.** `inject_shear.field_spec(<field>)` (models, F200W catalogues, photo-z, Gaia veto centres) as used for W1.
Screen thresholds: the CLI defaults (`radial_defaults()`); `--max-radius` 100″ (ASSUMPTION, as for the shear runs).
RXJ2129 (Lenstool maps 8 / 12) and MACS0647 (CLASH 67 % / 50 %) failed the `arcs` alignment control: they are run
but reported as "no result".

**Known-case check first.** Re-run the screen on Abell S1063 with its stored inputs and reproduce the counts in
docs/fields/abells1063.md (51 `anti`, 0 flags).

**What changes the plan.** Any flag on a control-passing field → `/vet-candidate` before any reading. A null → the
radial screen covers eleven fields; the TASKS item moves on.
