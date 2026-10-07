# Sources

Enough information to recover each external source later: URL, version or DOI where relevant, and
the date it was checked. Sections marked with a unit number are owned by that bootstrap unit;
edit only your own section in parallel work.

## Data (checked 2026-10-07)

- **MAST** (Mikulski Archive for Space Telescopes), CAOM API: https://mast.stsci.edu/api/v0/_services.html ;
  Python examples: https://mast.stsci.edu/api/v0/pyex.html
- **JWST program 2736**, ERO SMACS J0723.3-7327: https://www.stsci.edu/jwst/science-execution/program-information?id=2736
  Level-3 imaging: NIRCam F090W/F150W/F200W/F277W/F356W/F444W (obs `o001`), MIRI F770W/F1000W/F1500W/F1800W
  (obs `o002`), all PUBLIC (release MJD 59773.625). NIRCam F200W: `_cat.ecsv` 3,444,394 B (3,145 rows,
  59 columns, made with jwst 2.0.1 / photutils 2.3.0), `_i2d.fits` 1,762,162,560 B. MIRI F770W `_i2d.fits` 36,083,520 B.
- **JWST program 1345** (CEERS): https://www.stsci.edu/jwst/science-execution/program-information?id=1345
  Control pointing `jw01345-o001_t021_nircam_clear-f200w`: `_cat.ecsv` 4,017,583 B (3,695 rows, jwst 3.0.0),
  `_i2d.fits` 1,248,422,400 B.
- **STScI public S3 bucket** `stpubdata` (anonymous): https://stpubdata.s3.amazonaws.com/ — level-3 key pattern
  `jwst/public/jw{PPPPP}/L3/t/o{OOO}/{obs_id}_{suffix}` (verified for 2736 and 1345).

## Core software

- **Astropy**: https://www.astropy.org — cloud FITS byte-range access:
  https://docs.astropy.org/en/stable/io/fits/usage/cloud.html
- **astroquery.mast**: https://astroquery.readthedocs.io/en/latest/mast/mast_obsquery.html
- **JWST pipeline source catalog step** (defines `_cat.ecsv` columns):
  https://jwst-pipeline.readthedocs.io/en/latest/jwst/source_catalog/main.html

## Agent tooling (unit 8 extends)

- Claude Code skills: https://code.claude.com/docs/en/skills ; routines (`/schedule`):
  https://code.claude.com/docs/en/routines
- `anthropics/claude-code-action`: https://github.com/anthropics/claude-code-action

## Archive access (unit 1)

_Pending._

## Catalogs and cross-band matching (unit 2)

_Pending._

## Features and anomaly-detection methods (unit 3)

Checked 2026-10-07. Library versions are the ones resolved for the project environment that day.

- **scikit-learn** 1.9.1. `IsolationForest`:
  https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.IsolationForest.html ;
  `LocalOutlierFactor`: https://scikit-learn.org/stable/modules/generated/sklearn.neighbors.LocalOutlierFactor.html ;
  overview: https://scikit-learn.org/stable/modules/outlier_detection.html
- **SciPy** 1.18.1. `median_abs_deviation`:
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.median_abs_deviation.html ;
  `rankdata`: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.rankdata.html
- Liu, Ting & Zhou 2008, "Isolation Forest", ICDM, 413–422, doi:10.1109/ICDM.2008.17
- Breunig, Kriegel, Ng & Sander 2000, "LOF", SIGMOD, 93–104, doi:10.1145/342009.335388
- JWST `source_catalog` step: aperture EE defaults of 30/50/70 and the CI star thresholds,
  https://jwst-pipeline.readthedocs.io/en/latest/jwst/source_catalog/arguments.html ; APCORR reference
  file (NIRCam rows: filter, pupil, eefraction, radius),
  https://jwst-pipeline.readthedocs.io/en/latest/jwst/source_catalog/reference_files.html
- **PyOD** 3.6.6 (PyPI, released 2026-09-17; evaluated and not adopted, see D-004):
  https://github.com/yzhao062/pyod . ECOD: Li et al. 2022, arXiv:2201.00382,
  doi:10.1109/TKDE.2022.3159580
- **ADBench**, Han et al. 2022, arXiv:2206.09426 (30 algorithms on 57 datasets)
- **Astronomaly**, Lochner & Bassett 2021, arXiv:2010.11202, doi:10.1016/j.ascom.2021.100481 ;
  https://github.com/MichelleLochner/astronomaly (last push 2026-05-12)

## Imaging, cutouts and visualization (unit 4)

_Pending._

## External catalog services (unit 5)

_Pending._

## Candidate store and run provenance (unit 6)

_Pending._

## Open-source project tooling (unit 7)

_Pending._

## Future milestones: models, frameworks, lensing, HLSPs, tracking (unit 9)

_Pending._
