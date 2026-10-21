# MACS1149 system 16: DJA NIRSpec table v4.6 re-check: result

Plan `PLAN.md` (frozen before the table was read). Catalogue rows only; no spectra or images. Labels: observed
(catalogue values), derived (separations), assumption (0.5″ match radius, grade 3 = spectroscopic z).

## Frozen outcome: **not decided** (no grade-3 redshift for 16.1 or 16.2 in v4.6); D-037 unchanged

| Row (v4.6, 99,478 rows, sha256 `ee64bb46…6ddf3`) | Match | v4.6 values | Change vs v4.4 (`../summary.json`) |
|---|---|---|---|
| `macs1149-v4_prism-clear_1208_5114276` | 16.2 at 0.045″ (derived) | `grade` empty, `z_best` −1, automated `zline` 0.441 at `sn50` 1.34, EAZY `z_phot` 2.79 | none in the values v4.4 recorded (v4.6 msaexp version 0.9.5.dev2+g80b81b4; v4.4's was not recorded) |
| 16.1 | no row within 0.5″ | — | none |
| CANUCS 5106117 (D-037's possible far image) | no row | — | none |
| `…_1208_5106431` (Treu 26.2, not system 16; context only) | 1.6″ from 16.1 | grade 2, `zline` 3.189 | none |

- v4.5 (2026-02-23) was not read: it is older than v4.6, which supersedes it.
- The 16.2 row is still ungraded with identical automated values; the 2026-10-20 recheck condition ("when DJA
  grades `1208_5114276`") still holds.

## What would decide it
Unchanged: a graded NIRSpec redshift of 16.1 or 16.2 (a DJA re-grade of 5114276 or a deeper spectrum). Next check
only when a DJA table newer than v4.6 appears (HEAD `extractions/dja_msaexp_emission_lines_v<ver>.csv.gz`).

## Reproduce
`curl -sS https://s3.amazonaws.com/msaexp-nirspec/extractions/dja_msaexp_emission_lines_v4.6.csv.gz | gunzip -c |
python -I <filter>`, where the filter reads the CSV from stdin and keeps rows with `file`, `ra`, `dec` within 0.5″
of the two CATS v4.1 positions, or whose `file` contains `5114276`, `5106117` or `5106431` (`summary.json` lists
the columns kept).

## Cost
One 172 MB table streamed through a filter (not stored; cloud disk rule; reason: the only public graded JWST
spectroscopic table covering CANUCS MACS1149), < 1 min.
