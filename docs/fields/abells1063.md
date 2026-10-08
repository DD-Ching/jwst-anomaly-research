# Abell S1063 (Hubble Frontier Fields)

Field note for the exotic-lens screens on the HFF CATS map model `abells1063-cats` (D-035). Abell S1063 is at
z = 0.348. Every number below is `observed`, `derived` or `model_prediction` as labelled, and every threshold is an
ASSUMPTION (D-029, D-031, D-034, D-035). No physical interpretation is made.

## Data (accessed 2026-10-08)

- **JWST (`observed`):** GLIMPSE, program 3293 (PI Atek), NIRCam observation 1, association `t001`, all 9 PUBLIC
  level-3 bands: `jw03293-o001_t001_nircam_clear-{f090w,f115w,f150w,f200w,f277w,f356w,f410m,f444w,f480m}`
  ([`configs/abells1063.yaml`](../../configs/abells1063.yaml), `abells1063_v1`; `stages:` identical to
  `configs/reference_sample.yaml`).
  - One MAST position query (342.1832, −44.5309, r = 0.02°; 4.7 s) returned 19 NIRCam level-3 image observations:
    GLIMPSE 3293 o001 (9 bands, F200W 70 ks, F444W 141 ks), 1840 o010 (6 bands, ~1 ks each) and VENUS 6882 o059.
    3293 was chosen for depth and band coverage. Every band's `s_region` contains the cluster position.
  - jwst 3.0.0 for every band. `_cat.ecsv` 1.3–5.3 MB; `_i2d.fits` are 5.8–9.2 GB (short wavelength) and 0.5–0.8 GB
    (long wavelength). They are never downloaded; any cutout must read S3 byte ranges. Manifests:
    `data/manifests/abells1063.ecsv` (9 catalogs with sha256) and `data/manifests/abells1063_products.ecsv`
    (18 level-3 products with S3 URIs).
- **DJA v7.5 photo-z:** the DJA v7 index lists an `abells1063` mosaic. `abells1063-grizli-v7.5-fix.photoz.tar.gz`
  (74,481,914 bytes, sha256 `088a1954e357d96e9e3b0ab4af1f474eb9d35676774c1ad990f40140e81b1319`, last-modified
  2025-01-09) holds `abells1063-grizli-v7.5-fix.eazypy.zout.fits`: 12,597 rows, 246 `z_spec` > 0, sha256
  `03bc4c289982dc50a4540305b8858b279aca7c7cdb64c809d69a855e631b6428`. It was used directly as `--photoz`
  (`attach_photoz` matches within 0.3″). **Caveat:** 13 of the 46 `arcs.txt` images with a DJA match within 0.3″
  get z_phot < 0.6 (blends with members and intracluster light), so the background cut drops many real arcs here.
  The config has no `matched_photometry` block (DJA photometry is not wired in).
- **Lens model:** CATS v4.1 HLSP maps and `arcs.txt`, pinned in `lensmodel.HFF_CATS` (D-035). No `params.txt`
  (404): placeholder redshifts, image list gated off, so only `radial` applies.

## Steps and results

| Step | Command | Wall time |
|---|---|---|
| 1 | MAST position query + footprint/product check | 5 s + ~1 min |
| 1 | `fetch_reference_sample.py --config configs/abells1063.yaml --catalogs-only` | 66 s |
| 2 | `lens_consistency.py --model abells1063-cats validate` | 5 s |
| 3 | frame offset (scratch: `arcs.txt` against the F200W/F150W/F277W catalogues) | 2 s |
| 4 | DJA index check + photo-z tarball download (75 MB) and extraction | ~10 s |
| 5 | `exotic_screens.py --model abells1063-cats radial --catalog <F200W cat> --max-radius 100` (without / with `--photoz`) | 18 s / 14 s |

**2. Validate (`model_prediction` against the published products).** κ median relative difference 1.8e-2 (p95
1.8e-2, 1,522 points; median ratio 0.982); μ(z = 2) 5.8e-3 (p95 4.1e-2, 1,639 points). Image plane: 72 of 72
images solved, rms 11.81″, max 46.9″ against the quoted 0.48″. The gate fails (`image_list_gate.passes: false`,
`agrees_with_models_setting: true`), as in D-035. Shared predicted images: 1.1–1.3, 2.1–2.3, 11.1–11.3.

**3. Frame offset (`derived`).** Median offset (JWST − `arcs.txt`) to the nearest catalogue source within 0.5″:
F200W dRA cos δ −0.059″, dDec −0.046″ (n = 42 of 72; bootstrap σ 0.009″/0.010″); F150W −0.057″, −0.044″ (n = 43);
F277W −0.068″, −0.065″ (n = 16). |offset| = 0.075″ < 0.1″ (ASSUMPTION threshold): **no `frame_offset_arcsec`**,
no code change.

**5. Radial screen (`derived`).** ASSUMPTIONs are the script defaults (as in docs/fields/macs0416.md).

| Run | elongated (spike segs dropped) | not behind lens | `anti` (not model-radial) | centres | null mean (p95) | max lines | p_random_max |
|---|---|---|---|---|---|---|---|
| no photo-z | 436 (26) | 0 | 51 (46) | 10 | 8.2 (13) | 4 | 0.435 |
| DJA v7.5 photo-z | 436 (26) | 82 | 36 (32) | 6 | 3.2 (7) | 3 | 0.95 |

No centre has p_random < 0.05, and the strongest has 4 lines (< 6). Under the vetting rule, no cutout vetting is
needed. **Null.**

The two 4-line centres in the run without photo-z are at (53.0, 7.5)″ (342.16257, −44.52881; dark, 1.2″ from
source 1142) and (−26.0, 16.0)″ (342.19335, −44.52644; 0.18″ from source 3818, isophotal S/N 10,595, a bright
cluster galaxy). Neither is significant (p 0.435).

**Counts.** Screened: 51 `anti` arcs (36 with photo-z), of which 46 (32) are not model-radial. Flags: 0. **Result:** no exotic-lens candidate in Abell
S1063 with the CATS v4.1 map model; nothing for `/vet-candidate`.

## Limits

- Radial screen only (image list gated off, D-035).
- The DJA photo-z put 13 of 46 matched `arcs.txt` images in the foreground (z < 0.6), so the photo-z run removes many
  genuine background arcs; the run without photo-z is the more complete one.
- Gaia-seeded veto (`--spike-stars`, 11 Gaia DR3 stars G < 17 within 4′; D-043): the result is unchanged. 26 segments
  dropped, 51 `anti`, max 4 lines, p 0.435; with photo-z max 3, p 0.95. No saturated or off-mosaic star adds segments
  here.
