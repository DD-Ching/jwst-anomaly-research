# MACS0417 and MACS1423: published multiple images against the pinned CANUCS v1 maps (frozen plan, 2026-10-21 03:27Z)

**Question.** Do the published multiple images of MACS J0417.5-1154 and MACS J1423.8+2404 map back to one source per
system and re-solve to the catalogued positions in the pinned CANUCS DR1 v1 Lenstool maps (`macs0417-canucs-maps`,
`macs1423-canucs-maps`; second model for MACS0417: `macs0417-relics-lenstool`)? Mission item: image positions and
multiple images where none are predicted (or predicted where none are catalogued), on two clusters not yet tested.

**Inputs (observed; transcribed from arXiv LaTeX sources 2026-10-21).**
- MACS0417: Jauzac et al. 2019 (arXiv:1811.02505, MNRAS 483, 3082) Table 4, the gold + silver set of Mahler et al. 2019:
  56 images, 20 families (knots 1a/1b/1c/1m, 2a/2b/2c count as families). Spectroscopic: systems 1 (z 0.8718), 2 and 3
  (z 1.046), 25 images. The other systems carry the Jauzac model z.
- MACS1423: KALEIDOSCOPE (arXiv:2405.04577, MNRAS 533, 4500) Table "tab:mulimg": 17 images in 4 systems; 1 and 2
  (z 2.84, one galaxy, two components) and 3 (z 1.779) spectroscopic; 4 model z 1.797.
- CANUCS v1 image catalogues are not public (readme: with v2); the CANUCS models used 30 / 5 spectroscopic systems,
  so the spectroscopic images here are probably among their constraints (ASSUMPTION; not verifiable until v2).

**Already seen before freezing (stated, not hidden).** A loading check at 03:22Z ran `image_plane` on the 25 MACS0417
spectroscopic images: CANUCS reproduced 24 within 0.8″, 1b.4 (a fourth image of knot 1b near the BCG) at 5.3″;
RELICS v1 missed most by 3–13″. The gates below were written after that look; their thresholds are the ones used in
the MACS1149 / MACS0717 pins (G2), unchanged.

**Frozen rules (ASSUMPTIONs; `scripts/m0417_images.py`).**
- Gate K, per model and cluster, on the spectroscopic images: image-plane rms (back-trace to the family mean source,
  re-solve, nearest predicted image) ≤ 1.0″ and ≥ 0.8 of the images matched within 3″. A model that fails K gives no
  reading.
- Non-spectroscopic systems: z scanned on a grid from z_lens + 0.2 to 10 (60 log-spaced steps); the z with the
  smallest family rms among those matching every image within 3″ is the fitted z; if none, the z with the smallest
  rms. Per image: `reproduced` (≤ 3″ at that z) or `not_reproduced`.
- Spectroscopic images are read at their z_spec with the same 3″ rule.
- Reading per image across the gate-passing models: `not_reproduced` in all of them → `lensing_residual_input` (goes
  to `/vet-candidate`, no interpretation); otherwise `no_flag`. Families whose model predicts more images than are
  catalogued are listed (`extra_predicted`) for a later forced-photometry step; they are not flags here.

**Ordinary explanations to test first for any flag.** Catalogue identification (wrong image or system), model core /
BCG profile (central images depend on the core radius and the BCG truncation), frame offset (HST vs JWST, < 0.3″),
transcription errors, fitted-redshift degeneracy.

**What changes the plan.** CANUCS failing K → loading or frame problem, fix before any reading. RELICS failing K →
record and run its κ map check (G1 of the earlier pins) to separate a loading error from model quality; MACS0417
readings then rest on CANUCS alone (single-model, labelled so). A flag → `/vet-candidate` in the next cycle.
