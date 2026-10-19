# AT 2024jfz: NIRISS 4368 F150W check of the z_s 1.35–1.40 counter-image (plan, frozen before any pixel read)

## Question
#358 (`../RESULT.md`, post hoc) found one archival epoch whose timing matches a repeat of AT 2024jfz at a
counter-image predicted by Mahler+22 ICLv2: NIRISS imaging, program 4368, CLEAR;F150W, MJD 60557.96 (MAST
`jw04368-o001_t001_niriss_clear-f150w`, 859 s, centre 110.82741, −73.45474: the predicted positions are ≈ 9″
from it, inside the 2.2′ field), if z_s ≈ 1.35–1.40. There the counter-image lies 6.4–6.8″ from the TNS position,
μ_k/μ_0 ≈ −0.84, delay 94–114 d (`../posthoc.json`).

- **Hypothesis (H1):** AT 2024jfz is a lensed transient at z_s 1.35–1.40 and the model is right to ≲ 0.5″; then a
  point source of ≈ 2579 × 0.84 ≈ 2170 nJy (ASSUMPTION: F150W ≈ the F090W flux; the rule allows a factor 8 below
  this, see the threshold) appears at the predicted track in the 4368 epoch.
- **Ordinary explanations of a non-detection:** z_s outside 1.35–1.40 (most of the z_s range predicts no epoch);
  model position error > 0.5″; the event faded faster than the ≤ 30 d plateau assumed in #358; it is not lensed
  by the cluster in the way modelled (e.g. a foreground or cluster-member transient).
- **Ordinary explanations of a detection:** a static source not in the template (resolution or PSF differences
  between NIRISS and NIRCam), a detector artefact or cosmic ray, misregistration of a bright host.
- **What changes the plan:** "present" → vet (other 4368 data, any later epoch, morphology) before any
  interpretation; "absent" → the z_s 1.35–1.40 branch of the lensed-transient hypothesis is disfavoured under this
  model and plateau assumption, and the AT 2024jfz line closes until a host redshift exists; "untestable" →
  record why and close the same way.

## Data (observed)
- Epoch E: level-3 `_i2d.fits` of `jw04368-o001_t001_niriss_clear-f150w` (S3 byte-range box; nothing stored).
- Template T: the ERO 2736 NIRCam F150W level-3 `_i2d.fits` (2022), the pre-event reference in the same filter
  name (bandpasses differ slightly; ASSUMPTION: negligible for a ≥ 270 nJy aperture excess).

## Positions (model_prediction, then ASSUMPTION grid)
Track: the z_s 1.35 and 1.40 nearby counter-image positions of `../posthoc.json` and 4 points evenly between them
(linear; ASSUMPTION). Grid: every point of a 0.1″ square grid (centred on the track midpoint) within
**0.5″** (ASSUMPTION: model position error) of the track segment.
Also measured, not part of the rule: the TNS position itself (t_E − t_0 = +116 d).

## Measurement (derived)
`hit_lightcurve.aperture_photometry` (aperture r = 0.15″, annulus 0.35–0.6″) at each grid point in E and in T;
σ_E and σ_T = `empirical_error` (300 random apertures in each box). Δ = F_E − F_T; σ_Δ = hypot(σ_E, σ_T).

## Frozen rule (function `outcome` in `scripts/smacs_at2024jfz_niriss.py`, unit-tested before the run)
Combined over the grid (one outcome):
- T_ap = 0.125 × 2170 nJy = 271 nJy (ASSUMPTION: aperture fraction ≥ 0.5 and a factor ≤ 4 for fading or μ error).
- **untestable** if more than 20 % of grid points are NaN in E or in T, **or** the injection control fails: a
  Gaussian point source (FWHM 0.08″, ASSUMPTION) of 2170 × 0.25 = 543 nJy total injected into E at each track
  point must give Δ ≥ T_ap **and** Δ ≥ 5 σ_Δ at ≥ 5 of the 6 track points;
- else **present** if any grid point has Δ ≥ T_ap **and** Δ ≥ 5 σ_Δ;
- else **absent**.
A "present" is a flag for vetting, not a detection. Stamps of E, T and E − T (by eye) are made and looked at
before anything is written.
