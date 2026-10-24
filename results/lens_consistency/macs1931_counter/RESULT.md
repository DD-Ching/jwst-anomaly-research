# MACS J1931: VENUS imaging footprint for the counter-image route (2026-10-24)

Script `scripts/macs1931_counter.py footprint`; tests `tests/test_macs1931_counter.py`; numbers from `summary.json`.
Only MAST metadata and catalogued (observed) image positions were read; no map value, predicted position or pixel.

## Result
- GO 6882 (VENUS) has one level-3 NIRCam tile on the field, `jw06882-o038_t041`, in **10 bands** (F090W, F115W,
  F150W, F200W, F210M, F277W, F300M, F356W, F410M, F444W). Each band's `s_region` contains **all 51 used
  images**, the 26 spec-z images included.
- F200W level-3 products: `_i2d.fits` 1,291,757,760 B (over the 200 MB rule, so S3 byte-range cutouts only),
  `_cat.ecsv` 3,927,313 B, `_segm.fits` 184,645,440 B; all on `s3://stpubdata/jwst/public/jw06882/L3/t/o038/`.
- The F200W SCI header reads by byte range from this session: 10077 × 4577 px, MJy/sr, 0.0312″ / px.

## Reading
Imaging is not the blocker for MACS J1931: the field is covered in every band and the mosaic can be cut by byte
range. The blocker is the reference catalogue (`PLAN.md`): no public deblended catalogue, and pipeline segments are
rejected in cluster cores. `PLAN.md` proposes an own photutils catalogue on a background-filtered cutout, validated
first against the MACS0417 gate C2. It is not frozen, and nothing has been run on it.
