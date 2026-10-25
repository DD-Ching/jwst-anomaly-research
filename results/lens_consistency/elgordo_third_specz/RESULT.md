# El Gordo third-image candidates at DJA NIRSpec redshifts (result)

Plan `PLAN.md` frozen in 1a9e05e1 with its rule functions and tests (`scripts/elgordo_third_specz.py`,
`tests/test_elgordo_third_specz.py`), before the table was read. Numbers: `summary.json` (frozen rules),
`posthoc.json` (post hoc). Labels: observed (DJA rows, Diego+2023 positions, MUSE z), assumption (0.5″, 0.01 (1 + z),
gate Z thresholds).

## Frozen outcome: **not decided** (no input)
- DJA v4.6 read in full (99,478 rows, `rows_read`). No row lies within 0.5″ of any of the 218 Diego+2023 images (`diego23_images.ecsv` of addendum ay)
  (`matches` empty), so gate Z has 0 control images and fails (`gate_Z`), and all seven targets read
  `not_decided`; each also has no spectrum (`status` none).

## Post hoc (not in the plan)
- The table has no El Gordo spectra at all: 0 rows within 10′ of the cluster, nearest 207′ (`posthoc.json`).
- MAST (2026-10-25): the 12 JWST observations centred in a 0.10° × 0.08° box on El Gordo are all NIRCam imaging
  (programs 1176 PEARS, 6882 VENUS); 0 NIRSpec (`posthoc.json` `mast_query`). The box holds observation centres,
  so a NIRSpec pointing centred further out is not excluded.

## Reading
- Whether 6c / 13c / 20c are misidentified or Caminha's counter-image positions are off by 3–6″ stays open (addendum
  ay). No public JWST spectrum can decide it today.

## What would decide it
A NIRSpec or MUSE redshift at 6c, 13c or 20c (none public in DJA v4.6 or MAST as queried), or a JWST-era El Gordo
model with public maps.

## Reproduce
```
curl -sS https://s3.amazonaws.com/msaexp-nirspec/extractions/dja_msaexp_emission_lines_v4.6.csv.gz | gunzip -c \
  | python scripts/elgordo_third_specz.py /dev/stdin
```
