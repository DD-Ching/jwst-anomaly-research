# W3 dimming: cost and value of the D-039 persistence test in place of `single_epoch` (scoping)

Run 2026-10-19 15:14–15:20 UTC under the frozen `PLAN.md`. Code and thresholds of the screen unchanged.

## Data route (deviation)
- Catalogues: S3 `cloud_uri` from the committed `data/manifests/dimming_*_products.ecsv` (MAST API skipped).
- Gaia star mask: VizieR timed out (60 s read timeout, then a 20 s probe), so the same cone (6′) was taken from the
  ESA Gaia archive TAP (`gaiadr3.gaia_source`, J2016.0; `gaia_esa.py`). Reproduction check below shows the mask is
  equivalent for these fields.

## Reproduction (derived)
With the veto on, `whatif.py` reproduces every committed `exposure_source_yr` (both ρ, all t_E) of
`macs0416_regen`, `abell2744_calib` and `nexus_calib` exactly (`compare.py`, `np.allclose` rtol 1e-9).

## Workload (derived; `<field>/workload.json`)
Real screen flags removed by `single_epoch` and by no other catalogue veto:

| Field | flags | single-epoch only | of which vanish | compact (point-like) | compact vanish |
|---|---|---|---|---|---|
| MACS0416 | 854 | 203 | 201 | 12 | 12 |
| Abell 2744 | 530 | 267 | 267 | 8 | 8 |
| NEXUS | 3,793 | 875 | 873 | 17 | 17 |

**37 compact single-epoch vanish flags** over the three fields would go to the persistence test (plan limit 50).

## Value (simulated; what-if with the per-copy `single_epoch` veto off, ρ = 0.1)
Exposure (source-yr), joint over the three fields:

| t_E (yr) | 0.01 | 0.03 | 0.1 | 0.3 | 1 | 3 |
|---|---|---|---|---|---|---|
| veto on (committed) | 6.85 | 26.8 | 69.1 | 103 | 146 | 191 |
| veto off | 17.6 | 53.3 | 159 | 347 | 457 | 505 |
| gain | ×2.57 | ×1.99 | ×2.31 | ×3.36 | ×3.14 | ×2.64 |

- Per field: MACS0416 ×1.6–2.3, NEXUS ×1.7–2.7, Abell 2744 ×17–53 (its W3 recoveries are almost all full vanishes
  of two-epoch sources, which the veto removes; veto-on efficiency 0.000–0.013).
- 100 % dimming recovered: 240 → 680 / 680 (MACS0416), 20 → 1,238 / 1,240 (Abell 2744), 1,140 → 2,779 / 2,780
  (NEXUS); 20 % and 50 % unchanged (the veto only bites when the source disappears from the catalogue).

## Reading
- The plan's rule (≥ 20 % gain, ≤ 50 compact flags) is met: **build the swap.**
- The gain is an **upper bound**: it assumes the persistence test calls every real full vanish `on_sky`. D-039's
  `on_sky` needs two clean exposure-level detections at S/N ≥ 5 (ASSUMPTION); a source at mosaic S/N ≈ 10 in four
  dithers sits near S/N 5 per exposure, so the faint bins will lose part of the gain, and the first exposure of a
  visit cannot be cleared. That efficiency is not measured yet.

## Next (in order)
1. Measure the persistence test's `on_sky` efficiency against mosaic S/N on monitored real sources (known on-sky)
   in their epochs; it multiplies the what-if recoveries per magnitude bin.
2. Run `persistence_check.py` on the 37 compact single-epoch vanish flags; `on_sky` ones go through `forced` and the
   exposure-level check, then `/vet-candidate`.
3. Only then swap the veto in `efficiency_table` (per-copy efficiency from step 1) and re-run `combine`.
