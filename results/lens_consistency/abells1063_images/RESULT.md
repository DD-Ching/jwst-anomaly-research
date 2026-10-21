# Abell S1063: published multiple-image list for the D-035 revisit (data step, no lens-model reading)

Script: `scripts/s1063_imagelist.py`. Table: `chan20_images.ecsv` (observed positions; z per system). Numbers:
`summary.json`.

- **Source:** Chan et al. 2020 (ApJ 888, 35), VizieR `J/ApJ/888/35` table1, streamed 2026-10-21 (already used
  for the `arcs` anti positions on 2026-10-16, `results/arcs/abells1063-cats/RESULT.md`; the 2026-10-21 18:08Z
  keyword search for Caminha+2016 / Karman+2017 missed it).
- **Content:** 121 images in 42 systems; 12 systems (36 images) carry a MUSE z_spec, the rest get Chan's model
  geometric z (`Geoz`, z_spec False). Chan confidence classes: 72 / 20 / 29 images in classes 1 / 2 / 3.
  No system has images with different z_spec. Median |z_spec − Geoz| over the images of z_spec systems 0.129
  (derived; Chan's model may be fitted to these z, so this is not an independent check).
- **Source defect:** 27.1, 27.2 and 27.3 carry the same coordinates as 26.1 in the VizieR table
  (`dup_position_ids` in `summary.json`; found by `/code-review`). All four rows are flagged `dup_position` and are to be
  left out until checked against the paper's table; systems 26 and 27 are then untestable by position.
- **Not done here:** `validate` / `images` against `abells1063-cats`. The pinned CATS v4.1 files are on
  `archive.stsci.edu`, which still fails TLS from cloud sessions (2026-10-21 19:15Z, curl exit 000), and the MAST
  Download API returns 404 for the same `mast:HLSP/frontier/...` URIs (also for the MACS0416 CATS `arcs.txt`).

## Next
- Local session: `lens_consistency validate abells1063-cats` with this list as the image list (z_spec systems
  first, then all), then the D-035 rms gate and `images`; the Chan class ≥ 2 subset as a robustness variant.
