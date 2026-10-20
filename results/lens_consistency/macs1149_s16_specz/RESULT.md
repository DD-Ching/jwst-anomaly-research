# MACS1149 system 16: published spectroscopic redshift check: result

Plan `PLAN.md` (frozen before any catalogue was read). No images used; catalogue positions only. Labels: observed
(catalogue values), derived (separations), assumption (0.5″ match radius).

## Frozen outcome: **not decided** (no spectroscopic redshift for 16.1 or 16.2); D-037 unchanged

| Source | Match to CATS v4.1 16.1 (177.40971, 22.387662) / 16.2 (177.40989, 22.387828) | Redshift |
|---|---|---|
| CATS v4.1 `arcs.txt` | itself | z column 0.0 (no spec-z; the model fits z 4.419) |
| Treu et al. 2016 (VizieR J/ApJ/817/60, table 3; compiles HST-GO-13459/14041 grism and MUSE 294.A-5032 redshifts) | images 27.1 / 27.2 at 0.0″ / 0.0″ (derived), set "Silver" | `z` and `zspec` empty |
| DJA NIRSpec emission-line table v4.4 (`dja_msaexp_emission_lines_v4.4.csv.gz`, CANUCS 1208 PRISM) | 16.1: none within 0.5″; 16.2: `macs1149-v4_prism-clear_1208_5114276` at 0.04″ | `z_best` −1, no visual grade; automated `zline` 0.441 at median S/N 1.3 (not a redshift) |

- Not system 16: `1208_5106431` (Treu 26.2) 1.6″ from 16.1 has grade 2, z 3.189; it is a different catalogued
  system and is not used.
- CANUCS 5106117 (the possible far image, D-037) has no row in the DJA table.
- The 16.2 spectrum file could not be read (S3 path 403), so it was not inspected by eye.

## Post hoc (exploratory; does not change the frozen outcome)
- DJA's EAZY photo-z for 16.2 (grizli v7.0 photometry, `id_phot` 14171) is **2.79** (point estimate, no interval in
  the table). With CANUCS DR1 (2.25, 95 % 0.23–2.33), two independent photometric catalogues place 16.2 inside the
  redshift ranges where both models merge the flagged extra image onto 16.2 (CATS 2.1–2.9, Sharon 2.5–3.3;
  `../macs1149-cats/sys16_track.json`), and neither near the fitted 4.419. This supports D-037's reading
  (model error from a fitted redshift) without deciding it.

## What would decide it
A graded NIRSpec redshift of 16.1 or 16.2 (a DJA re-grade of 5114276, or a deeper spectrum), as D-037 already says.

## Cost
One 133 MB catalogue download (scratch, deleted after use; reason: the only public JWST spectroscopic table covering
CANUCS MACS1149), one VizieR query, one 3 kB `arcs.txt`.
