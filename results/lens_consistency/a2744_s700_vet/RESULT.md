# Abell 2744 system 700.1 (UNCOVER 78): vetting (result)

Plan: `PLAN.md` (frozen and committed before any measurement, with a computable specification added before the run).
Script: `scripts/a2744_s700_vet.py`; numbers: `summary.json`; sheet: `sheet.png` (looked at). Fluxes and noise are
`derived` from the observed `jw02561-o001_t003` mosaics; redshifts are catalogue values (`observed` spec-z, eazy
photo-z `model_prediction`); lens-model components are `model_prediction` inputs.

**Frozen outcome (from `decide`): "undetermined (colour test not validated; member perturber present)".**

1. **Colour test not validated.** Of 56 spec-z UNCOVER systems, the first 20 pairs with F200W S/N ≥ 10 in both
   images were tested; only 10 / 20 (0.5) read "consistent" at p ≥ 0.01, below the 0.8 gate (ASSUMPTION). With
   r = 0.15″ high-pass apertures and random-aperture noise, true multiple images often differ by more than the
   noise (χ² up to 254 for 7 bands; `known_pairs`). So "consistent" colours for 78 decide nothing.
   For the record (not a decision): the three 78 a/b pairs give p 0.96, 0.93, 0.68 (`target_pairs`).
2. **Member perturbers present.** Within 5″ of the 78.2 a/b midpoint, the DJA v7.2 table has two spec-z cluster
   members brighter than F200W 22 (crude r = 0.5″ aperture): DJA 45749 (z_spec 0.305, 2.11″, ≈ 20.4) and DJA
   45736 (z_spec 0.301, 4.43″, ≈ 18.6) (`perturbers`).
3. **Sheet** (`sheet.png`, by eye): the six 78 knots lie on one thin, straight arc about 2.5″ long, elongated
   roughly perpendicular to the direction of DJA 45749, which lies about 2″ south.

**Post hoc (after the frozen outcome; `post_hoc_bergamini23_components`).** Both members are already mass
components of Bergamini+23 (`best.par`): DJA 45749 is member potential 20036 (0.11″; σ 119 km/s, cut 21.9 kpc),
and DJA 45736 is the centre of the group halo O4 (σ 839 km/s) plus its galaxy-scale component O9 (σ 388 km/s)
(0.05″). So 700.1 lies 4.4″ from the centre of a group-scale halo of the northern substructure and 2.1″ from a modelled
member. A critical curve there depends on these components' parameters (and on the UNCOVER model's treatment of the
same galaxies, not checked). That leaves an ordinary explanation (mis-scaled local mass) open. It is not
established; no re-fit was made.

**Reading.** 700.1 stays undetermined and is not a candidate. The ordinary explanations are not excluded: one image of
a single galaxy (untested, because the colour test is not validated), or a fold pair whose critical curve the models
misplace near a modelled group core and member (plausible, untested). Only a local re-fit of the northern group
(O4/O9/20036 free) with 78 as a constraint, or resolved kinematics or spectroscopy of both arc ends, would decide it.

**Failed approach.** High-pass r = 0.15″ aperture colours with random-aperture noise do not reproduce colour equality
of known UNCOVER multiple images (10 / 20). Do not use pair colour χ² for one-image-or-two questions without a
systematic floor calibrated on known pairs first.

**Deviations (labelled).** Post hoc: the `--model-check` block and the sheet's stretch (vmax 10 → 5 σ, display only).
The run was made on a script version whose only later changes are the sheet, the model check and the
docstring (no change to measurement or decision code).

**Limitations.** One aperture size, one noise model. The crude magnitude is a 0.5″ aperture, not a total magnitude.
Photo-z membership was not used for any perturber that mattered (both are spec-z). The UNCOVER model's member list
was not checked.
