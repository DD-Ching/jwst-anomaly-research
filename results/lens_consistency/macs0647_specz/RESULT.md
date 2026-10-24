# MACS0647: CLASH LTM-Gauss maps at DJA NIRSpec redshifts (2026-10-24)

Plan `PLAN.md` frozen in f537735; rule functions S1–S4 and their tests (`scripts/macs0647_specz.py`,
`tests/test_macs0647_specz.py`) in b91478c, both before the table was read. Numbers from `summary.json`. Labels:
observed (DJA catalogue values, Meena+2023 positions), derived (separations), model_prediction (image-plane offsets),
assumption (thresholds). No imaging and no spectra were read: this is a position test of a published model at
catalogued redshifts.

## Frozen outcome: **not decided** (gate S3 fails)

- DJA v4.6 (99,478 rows read, `rows_read`; sha256 `ee64bb46…c6ddf3`, same bytes as in `../macs1149_s16_specz/dja_v46/`),
  streamed and matched within 0.5″. The first run read a copy cut to a box around the cluster; the stored
  `summary.json` is from the full stream and is byte-identical to that run's before the review fixes (not stored).
- 19 Meena+2023 images have a DJA row within 0.5″ (`S1`); 10 have a grade-3 z (S1). Families with an S2 z: 1, 2, 3,
  6-1, 6-2, 10, 11 (pre-JWST controls) and 13, 15 (held out); 16.4's grade-3 z 6.13 sits on an uncertain row and is
  not used. Each S2 z rests on one image (`S2.images`); no family is `inconsistent` and no image `ambiguous`.
  6-1.2 and 6-2.2 take their z from the same DJA source 3349 (`shared_sources`; two components of MACS0647-JD in
  Meena+2023), so the 9 families rest on 8 spectra.
- **S3 (controls at their spec-z):** 11 / 19 certain images within 3″ (0.58 < 0.80, ASSUMPTION) → fails. The misses
  are 1.1 (4.76″), 6-1.2 / 6-2.2 (10.8″ / 10.9″), 6-1.3 / 6-2.3 (6.5″), 11.1–11.3 (4.8–11.9″) (`per_image`).
- S4 therefore reads `not_decided` for the six held-out images of families 13 and 15 (offsets 1.1–3.9″, not readings).
- Post hoc (review): with each DJA source counted once (6-2 left out), S3 still fails, 10 / 16 (0.625;
  `post_hoc_S3_sources_once`).

## Reading
- With each family's z fixed at its catalogued spectroscopic value, the HST-era LTM-Gauss maps place 8 of 19
  pre-JWST control images more than 3″ from any predicted image (MACS0647-JD at z 10.17 and family 11 at z 2.59
  miss by up to 10.9″ and 11.9″). At fitted z the same maps matched 37 / 38 images of all pre-JWST families 1–12 (addendum w).
  So under this model the free z hides these offsets, and it cannot test the redshift structure of the JWST-first
  families.
- Not separated here: model error, a wrong spec-z match (each family z comes from one image within 0.5″, e.g. 11.1,
  which itself misses by 11.9″) or a misidentified family member. None is tested; all are ordinary.
- No lensing-residual input: with the control failed, the held-out images are not read.
- Post hoc, not rules: fitted z (addendum w) vs spec-z: family 1 1.81 / 2.14, 2 4.77 / 4.75, 3 2.18 / 3.05, JD
  5.00 / 10.17, 10 6.29 / 7.47, 11 1.81 / 2.59, 13 3.15 / 3.60, 15 3.46 / 3.24 (`summary.json` here and
  `../macs0647_clash_heldout/summary.json`).

## What would decide it
A JWST-era MACS0647 model with public maps (fitted with these spectroscopic redshifts), run with S3 as its control.

## Reproduce
```
curl -sS https://s3.amazonaws.com/msaexp-nirspec/extractions/dja_msaexp_emission_lines_v4.6.csv.gz | gunzip -c \
  | python scripts/macs0647_specz.py /dev/stdin
```
(The run read a copy cut to `box_prefilter_deg` first, which holds every 0.5″ match.)
In cloud sessions `archive.stsci.edu` needs a CA bundle with Mozilla's roots (`SSL_CERT_FILE`; CHANGELOG 2026-10-23).

## Cost
One 172 MB table streamed (not stored; cloud disk rule; reason: the only public graded JWST spectroscopic table
covering MACS0647), two CLASH map files (cache, not committed); < 1 min of compute.
