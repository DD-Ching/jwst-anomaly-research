# `arcs` anti row MACS0308 2157: MUSE redshift attempt — gate M failed as frozen; no redshift

Plan and the pre-data amendment: `PLAN.md`. Script: `scripts/arcs_muse_2157.py` (stream reader, aperture, noise,
z scan, `classify`; unit-tested in `tests/test_arcs_muse_2157.py` before the run). Numbers: `summary.json`
(`streamed_utc` 2026-11-08T20:12), spectra and scores `spectra.json.gz`, post hoc `posthoc.json`.
Cube: ESO MUSE `ADP.2020-01-04T23:56:39.902` DATA extension (1.55 GB) streamed once, two 80 × 80-spaxel boxes kept in
memory, nothing stored (Cloud disk decision; the portal ignores HTTP Range).

## Gate M (validation): **failed as frozen**
- Control (central galaxy, NED WISEA J030855.93+264537.9 position): z scan best 0.3568, but 0 features with
  |S/N| ≥ 5 (Mg b 4.5, Na D 4.1, G band 3.7, Ca H 2.6, Ca K 2.3 at the best z; `summary.json`), so gate M (≥ 2
  features, |z − 0.356| < 0.01) fails. Per the plan the target reading is not used for a membership decision.
- Contact sheet (`contact_sheet.png`, looked at): the 0.6″ aperture at the NED/WISE position lies about 1″ from the
  galaxy's white-light peak (by eye, not stored), so it samples the outer light; the 4000 Å break and absorption
  troughs are visible by eye.

## Target (observed / derived; not a decision, gate failed)
- 2157 is visible as a faint source at the aperture position in the white light; median continuum S/N 2.0 per
  1.25 Å pixel; 94 sky-masked pixels.
- Frozen function output: z_best 0.4174, 0 features ≥ 5 → `undetermined`. At z 0.356 and at the photo-z 0.756 no
  feature exceeds |S/N| 1.9 (`features_at_z_cluster`, `features_at_z_phot`).

## Post hoc (not in the plan; `posthoc.json`)
- Peak contrast (best score / highest score > 0.01 away): control 1.79 (0.3568 vs 0.1552), target 1.07 (0.4174 vs
  0.317). The control's scan identifies the cluster redshift; the target's scan has no preferred redshift.
- Narrowest reading: this MUSE exposure (2629 s, `../spec-coverage/summary.json`) gives no emission line at |S/N| ≥ 5 in a 0.6″ aperture at 2157 at
  any z in 0–1.6 under these windows; absorption-line redshifts at its continuum S/N are not reachable with this
  tool. Membership stays undecided.

## Limitations
- Noise is the robust scatter of random apertures in the box, not the STAT extension.
- The control position is a WISE position (offset from the MUSE peak); no recentring was in the plan.
- Line windows and the S/N 5 cut are ASSUMPTIONs; the per-line cut is too strict for absorption-line galaxies at this
  depth (control), so the frozen gate tests the cut as much as the cube.
- No injection-recovery of an [O II] line at 2157's position was run, so no line-flux limit is quoted.

Not a lensing anomaly claim and not a candidate.
