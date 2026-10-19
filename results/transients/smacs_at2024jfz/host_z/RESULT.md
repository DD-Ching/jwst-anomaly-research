# Result: host redshift of AT 2024jfz (2026-10-19)

Plan: `PLAN.md` (frozen first). Script: `scripts/smacs_at2024jfz_hostz.py`; rule functions `classify` and
`pick_host` unit-tested (`tests/test_smacs_at2024jfz_hostz.py`). Numbers: `result.json`.

**Frozen outcome: `untested`.** No DJA v7.4 source lies within R_MAX = 1.0″ of the TNS position (`candidates`
empty), and no Caminha+22 MUSE source lies within 3″ (`muse_within_3arcsec` empty). The host redshift stays
unknown, so this step does not decide which counter-image window of #360 applies.

**Post hoc (not part of the frozen rule; `posthoc_nearest5`).**
- The nearest DJA source is id 9171 at 1.632″ (r_half 0.209″, d / r_half 7.79), eazy z_phot 1.626 (z160–z840
  1.597–1.634, derived; 20 filters). It sits 0.12″ from catalogued image 6.1 of Mahler+22 ICLv2 (not stored),
  i.e. it is that image. Its interval lies inside the undetermined window 1.55–1.70; under the frozen rule it is
  not the host (beyond R_MAX).
- The next DJA sources are ≥ 3.05″ away.
- The ERO F150W template stamp (`../niriss4368/stamps_tns.png`) shows faint flux at the position (31.3 ± 3.2 nJy,
  #360) on a bright vertical band. A faint host blended into that band is not excluded by an empty catalogue
  match; "no catalogued host" is the narrowest statement.

**Labels.** z_spec: observed (none found); photo-z: derived; R_MAX, MUSE_MATCH, R_FLOOR: ASSUMPTION; image 6.1
identity: catalogue match (not stored).

**What would decide it.** A spectrum at the transient position (NIRSpec MSA or MUSE re-reduction), or a
deblended photo-z of the 31 nJy template flux after modelling the vertical band. If the transient belongs to the
system 6 galaxy (post hoc reading, 7.8 r_half from its catalogued centre), z_s ≈ 1.60–1.63 falls in the
undetermined window, and the next step is a PSF-matched difference image at the z_s 1.55–1.70 counter-image track.

**Deviations.** The rule functions and their tests passed before the first run, but their commit landed after it
(a lint failure stopped the commit). The `posthoc_nearest5` block and the removal of a local path from
`sources` were added after the first run; the re-run gave the same outcome.
