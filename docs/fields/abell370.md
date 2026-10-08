# Abell 370 (Hubble Frontier Fields)

Field note for the exotic-lens screens on the HFF CATS map model `abell370-cats` (D-035). Abell 370 is at z = 0.375.
Every number below is `observed`, `derived` or `model_prediction` as labelled, and every threshold is an ASSUMPTION
(D-029, D-031, D-034, D-035). No physical interpretation is made.

## Data (accessed 2026-10-08)

- **JWST (`observed`):** CANUCS GTO program 1208, NIRCam observation 2, association `t001`, all 8 PUBLIC level-3
  bands: `jw01208-o002_t001_nircam_clear-{f090w,f115w,f150w,f200w,f277w,f356w,f410m,f444w}`
  ([`configs/abell370.yaml`](../../configs/abell370.yaml), `abell370_v1`; `stages:` identical to
  `configs/reference_sample.yaml`).
  - One MAST position query (39.9714, −1.5822, r = 0.02°; 1.9 s) returned 76 NIRCam level-3 image observations:
    1208 o002 (8 bands, 6,399 s each), MAGNIF 2883 and 3538 (medium bands only), 5324 (three shallow 4–6-band
    visits), JUMPS 5890 (F150W + medium bands), VENUS 6882, and exclusive-access 7345 (Dragon) and 9371. 1208 was
    chosen: the only deep broad-band set with 8 bands. Every band's `s_region` contains the cluster position.
  - jwst 3.0.0 for every band. `_cat.ecsv` 2.4–5.7 MB; `_i2d.fits` 1.65 GB (short wavelength), never downloaded:
    cutouts read S3 byte ranges. Manifests: `data/manifests/abell370.ecsv` (8 catalogs with sha256) and
    `data/manifests/abell370_products.ecsv` (16 level-3 products with S3 URIs).
- **DJA:** none (no Abell 370 mosaic in the DJA v7 imaging index, checked 2026-10-08).
- **CANUCS DR1 photo-z:** `hlsp_canucs_jwst-hst_multi_a370-clu_multi_v1_photometry-cat.fits.gz` (32,153,881 bytes,
  sha256 `f5622f2867aa3094df6861981ee67f7f1879b8381b348224748c7381436ae17b`; 13,567 sources, 509 `Z_SPEC`). Converted
  to zout columns as in docs/fields/macs1149.md (`Z_SPEC` where > 0 for `z_phot`/`z160`/`z840`, else
  `Z_ML`/`Z160`/`Z840`; `ra`/`dec`; 0.3″ match). Blend caveat: 3 of the 68 `arcs.txt` images with a CANUCS match
  within 0.3″ get z < 0.6, so the background cut drops a few real lensed images.
- **Lens model:** CATS v4 HLSP maps and `arcs.txt`, pinned in `lensmodel.HFF_CATS` (D-035). No `params.txt`
  (404), so `arcs.txt` redshifts are placeholders and the image list stays gated off: only `radial` applies.

## Steps and results

Wall times are for this machine (4 cores, S3 from the session's network).

| Step | Command | Wall time |
|---|---|---|
| 1 | MAST position query + footprint/product check | 2 s + ~1 min |
| 1 | `fetch_reference_sample.py --config configs/abell370.yaml --catalogs-only` | 15 s |
| 2 | `lens_consistency.py --model abell370-cats validate` | 8.5 s |
| 3 | frame offset (scratch: `arcs.txt` against the F200W/F150W/F277W catalogues) | 2 s |
| 5 | `exotic_screens.py --model abell370-cats radial --catalog <F200W cat> --max-radius 100` (without / with `--photoz`) | 51 s / 43 s |
| 6 | vetting: line lists, cutouts (S3, 30″, 4 bands, two positions), Gaia DR3 query, sensitivity re-runs | ~6 min |

The first `validate` attempt failed the sha256 check of the z = 2 magnification map (a truncated download through
the proxy); the immediate retry matched the pinned sha256.

**2. Validate (`model_prediction` against the published products).** κ median relative difference 3.0e-2 (p95
3.0e-2, 2,471 points; median ratio 0.970); μ(z = 2) 2.1e-3 (p95 2.8e-2, 6,487 points). Image plane: 101 of 101
images solved, rms 10.66″, max 73.4″ (placeholder redshifts, no quoted rms): image list stays off, as in D-035.
Shared predicted images: 2.2–2.5, 5.2–5.4, 7.1/7.2, 10.1/10.2, 14.1/14.2, 16.1/16.3, 21.1–21.3.

**3. Frame offset (`derived`).** Median offset (JWST − `arcs.txt`) to the nearest catalogue source within 0.5″:
F200W dRA cos δ **−0.121″**, dDec −0.015″ (n = 41 of 101; bootstrap σ 0.009″/0.006″); F150W −0.107″, −0.014″
(n = 32); F277W −0.114″, −0.006″ (n = 20). Re-matching after shifting by the median gives the same values (F200W
−0.121″, −0.017″, n = 42). |offset| > 0.1″ (ASSUMPTION threshold), so it is pinned on this branch:
`_HFF_FRAME_OFFSET["abell370"] = (-0.121, -0.015)` in `scripts/lens_consistency.py`, and every radial run below
uses it (D-040 applies the offset to map models).

**5. Radial screen (`derived`).** ASSUMPTIONs are the script defaults (as in docs/fields/macs0416.md).

| Run | elongated (spike segs dropped) | not behind lens | `anti` (not model-radial) | centres | null mean (p95) | max lines | p_random_max |
|---|---|---|---|---|---|---|---|
| no photo-z | 949 (13) | 0 | 164 (153) | 71 | 94.8 (118) | **15** | **0.0** |
| CANUCS photo-z | 949 (13) | 121 | 135 (124) | 51 | 64.3 (82) | **15** | **0.0** |
| no photo-z, aper50 S/N ≥ 5 (3,806 of 5,266 rows) | 524 (9) | 0 | 66 (64) | 12 | 20.4 (31) | 8 | **0.01** |
| no photo-z, Gaia-spike veto (166 rows dropped) | 882 (13) | 0 | 140 (129) | 73 | 77.6 (95) | 5 | 0.945 |
| CANUCS photo-z, Gaia-spike veto | 882 (13) | 121 | 111 (100) | 53 | 47.7 (61) | 5 | 0.495 |

Two centres were flagged (p_random < 0.05) and both were vetted:
- the 15-line centre at (38.0, −8.0)″ in the default runs;
- the 8-line centre at (−67.0, 45.5)″ in the S/N ≥ 5 run. It is the same structure as the 10-line centre at
  (−72.5, 48.5)″ (p 0.66) in the default run.

## Flag vetting (step 6)

Cutouts F150W/F200W/F277W/F444W (30″, SCI and WHT, S3 byte ranges) were inspected for both flags. "aper50 S/N" is
the catalogue's 0.5″-aperture flux over its error.

| Flag | Centre (RA, Dec) | Line sources | Ordinary explanation tested | Verdict |
|---|---|---|---|---|
| 15 lines, p 0.0 (both default runs) | 39.960747, −1.584486 | Labels 308, 268, 3688, 357, 218, 3685–3687, 3677–3683. Isophotal S/N 14–165, but aper50 S/N 1.0–14.5 (aper50 AB 28.1–31.1). e 0.51–0.94; semimajor up to 45 px | **Chain:** all 15 lie on one line (perpendicular rms 0.42″, PA 62.3°). Their `pa_obs` is 44–73°, along the chain: a chain, not convergence. **Mosaic axis:** the i2d column axis is at PA 63.4°. **Seam:** in the cutouts the chain is a faint straight streak along the columns in F150W/F200W. It sits at a SW weight step (relative WHT 0.3 → 0.7; 10 of 15 sources at WHT 0.3, below the pipeline's 0.5 gate, D-011), 3–10″ from the mosaic's southern edge (cutout flags `edge,nan,low_weight`). It is weaker in LW (F444W catalogue: 0 sources on the chain). **Spike:** the line passes 1.1″ from Gaia DR3 G = 13.71 (39.955714, −1.586726), 13″ beyond the chain end and off the F200W mosaic (no catalogue source within 10″). | **ordinary: instrumental.** Diffraction spike of an off-mosaic G = 13.7 star along the column/V3 spike axis, at an SW low-weight seam. With S/N ≥ 5 the centre falls to 6 lines (p 0.65); with the Gaia-spike veto it disappears |
| 8 lines, p 0.01 (S/N ≥ 5 run); 10 lines, p 0.66 (default) | 39.989925, −1.569625 | Labels 3203, 4716, 3036, 3014, 2973, 2963, 2926, 4719 (4719: semimajor 114 px, e 0.98); aper50 S/N 5.7–53 | **Chain:** one line (perpendicular rms 0.29″, PA 66.1°), `pa_obs` 59–64°. **Cutouts:** a bright straight stripe along the mosaic columns in all four bands, with parallel ridges in F444W. It points at Gaia DR3 G = 12.71 (39.996246, −1.566737), 1.0″ off the line, 12–37″ beyond the line sources. The star has no catalogue counterpart within 4.8″ (saturated core), so `spike_segments` (D-034), which needs a catalogued point source brighter than AB 20 and caps the spike length at 20″, cannot see it. Both flags lie on one straight line at PA ≈ 63°, about 125″ apart. | **ordinary: instrumental (diffraction spike of a G = 12.7 star)**. With the Gaia-spike veto the centre disappears |

**Gaia-spike veto (scratch sensitivity test, not in the code).** Gaia DR3 stars with G < 17 within 4′ of the
cluster (7 stars, VizieR I/355/gaiadr3). A catalogue row was dropped when such a star lay 0.5–60″ away, the row's
major axis was within 7° of the direction to it, and that direction was within 7° of a spike axis (63.4° + 0/60/90/120°).
This dropped 166 of 5,266 rows. Max 5 lines, p 0.945 (no photo-z) and 0.495 (CANUCS photo-z): **null**.

**Counts.** Screened: 153 `anti` arcs without photo-z, 124 with. Flags: 2 (radial); both are diffraction spikes
of bright Gaia stars. Surviving every ordinary test: **0**.

**Result:** no exotic-lens candidate in Abell 370 with the CATS v4 map model. Nothing for `/vet-candidate`.

## Limits

- Radial screen only: no `params.txt`, so the image list (and `images`, `fluxratio`) is gated off (D-035).
- `spike_segments` misses spikes of stars that are saturated (not catalogued) or off the mosaic, and its 20″
  length cap is too short for G ≈ 13 stars (segments 12–37″ out). The Gaia-veto test above is scratch code;
  folding a Gaia-seeded spike veto into `exotic_screens.py` is proposed (TASKS).
- The null distribution of the default runs is inflated by the same spike segments (random max 7–14 lines), so
  the default-run p values for other centres are too high (conservative); the vetoed runs are the cleaner numbers.
- CANUCS photo-z blends (3 of 68 matched images at z < 0.6) drop a few real arcs from the photo-z runs.
