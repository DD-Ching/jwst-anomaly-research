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
| no photo-z, Gaia-spike veto (`--spike-stars`, committed code) | 880 (82) | 0 | 140 (129) | 73 | 77.6 (95) | 5 | 0.945 |
| CANUCS photo-z, Gaia-spike veto (committed code) | 880 (82) | 121 | 111 (100) | 53 | 47.7 (61) | 5 | 0.495 |

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

**Gaia-spike veto (in the code since D-043).** It was first a scratch test, which dropped 166 of 5,266 rows. It is
now `exotic_screens.py radial --spike-stars`: Gaia DR3 stars join `spike_segments` with spikes up to 60″, and the
field's hexagonal spike axes are estimated as before. The reproducing commands:

    python scripts/gaia_stars.py 39.97134 -1.58226 --radius-arcmin 4 --out stars.ecsv   # 7 stars, G < 17
    python scripts/exotic_screens.py --model abell370-cats radial --catalog <F200W cat> --max-radius 100 \
        --spike-stars stars.ecsv [--photoz <CANUCS zout>]

The committed code reproduces the scratch result: 82 spike segments dropped, 140 `anti`, 73 centres against 77.6,
max 5 lines, **p 0.945**. With CANUCS photo-z: 111 `anti`, 53 against 47.7, **p 0.495**. **Null.**

The "aper50 S/N ≥ 5" run in the table above was a scratch catalogue cut (rows with aper50 S/N < 5 removed before
`radial`). The aper50 floor is not in the code yet (TASKS).

**Counts.** Screened, after the Gaia veto: 140 `anti` arcs without photo-z, 111 with. Flags: 2 (radial); both are
diffraction spikes of bright Gaia stars. Surviving every ordinary test: **0**.

**Result:** no exotic-lens candidate in Abell 370 with the CATS v4 map model. Nothing for `/vet-candidate`.

## Limits

- Radial screen only: no `params.txt`, so the image list (and `images`, `fluxratio`) is gated off (D-035).
- Without `--spike-stars`, `spike_segments` misses spikes of stars that are saturated (not catalogued) or off the
  mosaic, and its 20″ cap is too short for G ≈ 13 stars (segments 12–37″ out). Use the Gaia-seeded veto (D-043).
- The null distribution of the default runs is inflated by the same spike segments (random max 7–14 lines), so
  the default-run p values for other centres are too high (conservative); the vetoed runs are the cleaner numbers.
- CANUCS photo-z blends (3 of 68 matched images at z < 0.6) drop a few real arcs from the photo-z runs.

## Arc orientations against the predicted shear (`arcs`, 2026-10-13, D-024 addendum 3)

Hypothesis and thresholds as in `macs1149.md` (D-024 defaults, ASSUMPTIONS). Model: CANUCS DR1 Lenstool
(`abell370-canucs`, frame offset as in `MODELS`). Inputs: F200W pipeline catalogue `jw01208-o002_t001` (sha256
`f4774b39…`, manifest) and CANUCS DR1 photo-z (sha256 `f5622f28…`, SOURCES; zout columns as in `macs1149.md`).
Outputs: `results/arcs/abell370-canucs/`.
- **Known cases:** 59 elongated catalogued images: median offset 8°, 53 within 30°. The two beyond 60° are knots
  6b.1 and 6b.3 (77°) of system 6 (z 1.0633), a resolved spiral.
- **Background, strong shear:** 47 of 56 aligned (p = 8e-15), **4 anti, 0.94 expected (P(≥ 4) = 0.016)**. All four
  looked at (6″ F200W cutouts):
  - 1565 (z 2.94, g_min 1.2, μ(z = 2) 98) and 1788 (z 5.2, g_min 4.2; aligned at z = 2, 89° off across its photo-z
    range): on critical curves, where the predicted direction flips and depends on z and on the exact curve
    position. The shape-noise null does not model this.
  - 4540 (m 20.7, z_spec 1.063, g 0.23) is catalogued image **6.1** (0.8″ from knot 6a.1): a resolved spiral whose
    disk shape (e 0.56) exceeds the shear; the model fits its knot positions. 340 (m 25.3, z 1.11, g 0.26): a
    compact source, e 0.55. Both are intrinsic shapes in weak shear.
- **Sensitivity (derived; `anti_expected` on the stored tables):** without the g_min ≥ 1 sources, 2 anti vs 0.61
  (P = 0.13); with σ_ε = 0.30 / 0.35, 0.98 / 1.44 expected. Over all five fields run so far (Abell 2744, El Gordo,
  MACS1149, MACS0416, Abell 370): 7 anti vs 2.6 (P = 0.017) with every source; 4 vs 2.13 (P = 0.17) for g_min < 1;
  4 vs 3.5–5.1 with σ_ε 0.30–0.35. The excess sits entirely at critical curves.
- **Members/foreground:** 5 of 13 anti (P = 1/3 each if unlensed). **No CANUCS match (77):** 15 anti, near random.
- **Result (superseded for 1565 by the scan below): null.** No anti-tangential arc that a critical-curve flip or an intrinsic shape does not explain.

### Redshift scan of the anti sources (2026-10-13, D-024 addendum 4)
The predicted axis turns by 90° only where κ crosses 1. `results/arcs/anti_flip_scan.json` gives κ, g, the lensing
axis ratio q_lens and the shape-noise P(anti) at each source's z16/z/z84 in CANUCS and CATS (`model_prediction`
and `derived`). CANUCS DR1 95 % ranges (Z025–Z975, `observed` catalogue values): 1565 0.40–3.47, 1788 0.98–5.88.
- **1788:** κ = 1 at z 3.16 (CANUCS) / 2.95 (CATS); at its photo-z κ 1.04–1.06, q_lens 0.62–0.71, P(anti) 0.05–0.09.
  A weak, model-sensitive radial prediction: explained.
- **1565 (39.977168, −1.583366):** κ 0.81–0.82 (CANUCS) / 0.76–0.77 (CATS) at z 2.8–3.2, below 0.9 up to z = 8:
  tangential at every redshift, predicted PA 58°, q_lens 0.02–0.14, P(anti) ≤ 2e-4. Weighted moments (Gaussian
  σ 0.08–0.12″, local plane background; noise from ERR): F150W 150–158° (S/N 2–4), F200W 150–155° (±7–11°, S/N 8),
  F277W 163–168° (±7°, S/N 10). In the F150W–F356W cutouts 1565 is a knot on a faint ~1″ structure at PA ≈ 160°
  that reaches neighbour 1562 (0.73″ at PA 162°; CANUCS z_ML 0.27, 0.09–2.85); a cluster member (CANUCS z 0.36,
  F200W 1.5 μJy) lies 1.2″ to the SE. The model's counter-image is 12″ NNW (39.979289, −1.580843; μ ≈ +14, so ≈
  28.9 mag): CANUCS sources there at 0.52″ and 0.85″ have F200W 8 ± 5 and 21 ± 4 nJy and colours unlike 1565's
  (F150W/F200W 0.47, F277W/F200W 1.1). Inconclusive.
- **Three more published models** (HFF v4 κ/γ maps at one pixel, scaled by D_LS/D_S; SOURCES): Sharon κ 0.79–0.81,
  GLAFIC 0.76–0.78, Keeton 0.75–0.77 at z 2.8–3.2 (≤ 0.87 at z = 8); q_lens 0.00–0.08. Shear axis PA 56° / 58° / 63°
  if γ1 is along the map's pixel x axis (that convention reproduces CANUCS's 58°). All five models predict a strong
  tangential stretch; model choice does not explain 1565.
- **CANUCS deblending (`observed` catalogue and segmentation values):** 1565's CANUCS counterpart has
  `FLAG_DEBLEND` = True. In the DR1 segmentation map (`…a370-clu-40mas_multi_v1_segmentation.fits.gz`, 1,239,107 B,
  sha256 `9baa4407…`) its segment 25146 has 31 pixels of 0.04″ (0.36″ long), major axis PA 165°, axis ratio 0.46, and
  touches the member's segment 922 (375 px in the 2″ box) at its southern end; 1562 is a separate segment (20618).
  Its colours differ from the member's (F150W/F200W 0.47 vs 0.90; F277W/F200W 1.1 vs 0.84), so it is not member
  light alone. At z 2.94, [O III] 5007 Å falls in F200W and Hα in F277W, which fits its F200W/F277W excess
  (`hypothesis`). A deblended fragment whose shape is set by the cut against the member's halo was the leading
  ordinary explanation; the next item rules it out for the image moments.
- **Member modelled out (`derived`):** a Sérsic + constant fit to the member on a 4.4″ stamp (n 0.9–1.1, r_eff
  0.22–0.23″; 1565 and 1562 masked) leaves the axis unchanged: F150W 151–156°, F200W 150–155°, F277W 163–168°
  (±7–12°). The member's smooth light does not set the axis; the deblending only affects the catalogue segment.
- **Noise axis (`derived`):** the same 0.08″ window at 150 blank F277W positions within 15″ (no catalogue source within
  1.2″) gives no preferred axis: 17 % within 15° of the mosaic row axis (153.4°), 14 % of the column axis, 17 %
  expected. The F200W run timed out (not run). Correlated noise does not obviously set the F277W axis.
- **Status: unvetted flag, not a candidate.** Ordinary explanations still open: a member tidal feature or a blend
  (the structure through 1562), a low-z photo-z solution (Z025 0.40; P(anti) 0.10 at z 0.6 and 0.24–0.33 below z ≈ 0.45), and an unmodelled
  local perturber (five models agree, so only a perturber none of them contains). Next: `/vet-candidate`.
