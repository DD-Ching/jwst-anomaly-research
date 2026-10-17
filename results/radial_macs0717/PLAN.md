# Radial (W1 exotic-lens) screen on MACS0717 with ASTRODEEP-FF photo-z: plan (2026-10-17, frozen before any run)

**Question.** TASKS "re-run MACS0717 radial when a photo-z catalogue exists": the 2026-10-08 run
(docs/fields/macs0717.md step 5: 74 `anti`, 68 not model-radial, max 5 lines, p_random >= 0.70) had no photo-z, so
cluster members and foreground galaxies entered. ASTRODEEP-FF photo-z (`scripts/astrodeep_zout.py`, #319) now exist
over the HST WFC3 footprint.

**Hypothesis (exotic, not assumed).** As in `results/radial_more/PLAN.md`: a compact negative-mass lens gives
radially stretched, demagnified image pairs; `exotic_screens.py radial` flags convergence centres of >= 3 `anti` arcs
with p_random < 0.05.

**Ordinary explanations to test first for any flag.** Spikes (Gaia `--spike-stars`), cluster radial arcs
(model-radial veto), shape noise (null draws), lens-model errors (CATS v4.1), streaks, edges and blends.

**Run (frozen).** `scripts/radial_more.py --fields macs0717 --out results/radial_macs0717` (unchanged code; inputs
`inject_shear.field_spec("macs0717")`: VENUS o029 F200W catalogue, ASTRODEEP zout, CATS v4.1, Gaia veto,
`--max-radius` 100", CLI defaults). Thresholds are ASSUMPTIONs (CLI defaults).

**Known limitation stated in advance.** Photo-z cover only ~17 % of the shape rows (855 of 4,953, not stored), so
the photo-z cut removes few members; outside the WFC3 footprint rows enter as in the unchanged code.

**What changes the plan.** Any flag -> `/vet-candidate` before any reading. A null -> the TASKS item is closed.
