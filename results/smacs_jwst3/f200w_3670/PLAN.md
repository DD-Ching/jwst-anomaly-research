# Plan: view and test SMACS f200w_3670 (TASKS Then 4), frozen 2026-11-07 before any image is read

**Object.** `jw02736-o001_t001_nircam_f200w_3670` (RA 110.836848, Dec −73.458635), stellar-locus "star", rank 24 in
the star stratum, 0.369″ from NED `SMACS J0723-73:[MJR2023] 002.3` (G_Lens, z 1.3715; Mahler+2023 system 2).

**Question.** Is it (A) a Galactic star projected next to arc image 2.3, or (B) a compact part of the lensed galaxy
(system 2) that the stellar locus mislabels? (B) would be a lens-consistency test: a compact lensed knot should have
counterparts in images 2.1 / 2.2 with the same colours.

**Ordinary explanations considered:** foreground star (A), compact lensed knot (B), cluster-member compact galaxy,
deblend fragment of the arc.

**Rules (ASSUMPTIONs):**
1. Look first: F200W / F444W / colour cutouts 6″ around the source and around images 2.1, 2.2.
2. (A) if a PSF diffraction pattern (six spikes) is visible in F200W at the source **and** its F090W−F444W / F150W−F277W
   colours differ from the 2.3 arc light by > 0.3 mag.
3. (B) if no spike pattern and its colours agree with arc 2.3 within 0.3 mag (aper_total, pipeline catalog).
4. Otherwise undetermined; recorded as such.
5. Any (B) outcome adds a counterpart check in 2.1 / 2.2; an unmatched compact knot goes to `/vet-candidate`.
