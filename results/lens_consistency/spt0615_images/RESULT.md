# SPT-CL J0615−5746 published images vs the RELICS Lenstool v1 and GLAFIC v2 maps (2026-10-25)

Frozen plan `PLAN.md` (written before any map was read; pushed in a6ba478c, see Deviation). Numbers come from
`summary.json` (frozen run) and `posthoc_3b_zspec.json` (post hoc), written by `scripts/spt0615_images.py` and
(post hoc file) by session code built on the same `m0417_images.Solver` / `lensmodel` calls, not committed as a script. Image list: Paterno-Mahler et al. 2018 Table 2 (observed), `spt_images.ecsv`.
Maps are observed products; fitted z, offsets, factors and classes are derived (model_prediction).

## Frozen result
| | Lenstool v1 | GLAFIC v2 |
|---|---|---|
| Gate P (κ, μ(z 2) median ratio) | 1.000003, 1.0000004 → pass | pass |
| Gate K0 (families 1–3, fitted z; 8 images) | rms 0.447″, matched 1.0 → pass | rms 0.283″, matched 1.0 → pass |
| Frame F (median ΔRA cos δ, ΔDec) | −0.224″, +0.088″ → pass | −0.074″, +0.048″ → pass |
| Reading B (families 1, 3 at z_spec; 6 images) | rms 0.472″, 0 misses | rms 0.197″, 0 misses |
| Reading A (3.3, 3.4, 3.6 in family 3b) | 3 misses (3b fitted at the grid edge z 1.172, 1 image predicted) | 0 misses |
| Rule 7 factors (families 1 / 2 / 3) | 1.026 / 0.987 / 0.984 → median 0.987, `near_1` | 1.026 / 0.932 / 0.997 → median 0.997, `near_1` |
| Extra predicted images (catalogued, predicted) | family 2 (2, 3), 3 (3, 5) | family 2 (2, 3), 3 (3, 7), 3b (6, 7) |

**Q1: `H_none`.** Both maps are usable and reproduce the spectroscopic families at their z_spec with no amplitude
factor. **Q2: no joint flag** in reading B or A (GLAFIC v2 reproduces 3.3, 3.4, 3.6; one usable model suffices).

## Post hoc (not in the plan)
- Family 3b at z_spec 4.013 without a z fit: GLAFIC all six images ≤ 0.27″ (7 predicted); Lenstool v1 ≥ 28″ (1
  predicted: the mean source is pulled by 3.3, 60″ from the others).
- From the secure family-3 source at each model's fitted z: Lenstool v1 predicts 5 images, nearest to 3.3 / 3.4 / 3.6
  at 36.3″ / 7.4″ / 2.7″; GLAFIC 7 images at 0.27″ / 0.28″ / 0.47″. So the two maps differ in system 3's image
  configuration, consistent with the paper's model 1 omitting 3.3, 3.4, 3.6 (Table 2 note) while GLAFIC v2 matches
  them; whether GLAFIC v2 used them as constraints is not documented in the HLSP (not tested).
- Position sheet looked at (observed images and 3b misses per model; not committed, regenerable from the ECSVs).

## Reading
- Addendum bk's SPT0615 "intermediate" (R 0.772 over |μ_ref| ≥ 10) is not an amplitude offset of either map at the
  image redshifts: both need factor ≈ 1. The difference over S goes with a different system 3 configuration
  (5 vs 7 images), i.e. a shape difference between the models (not measured beyond that).
- Not a lensing anomaly: 3.3, 3.4, 3.6 are "less secure" in the paper and one usable model reproduces them.
- Limits: families 1–3 are the models' own constraints (reproduction, not prediction); single plane (the z ≈ 0.4
  foreground structure is in neither map); HLSP Lenstool v1 = paper model 1 is an ASSUMPTION.

## Deviation
- The GLAFIC `validate` (gate P input, κ and μ maps only, no image position) ran about a minute before the plan
  commit was pushed; the plan text was written before it and not changed after.
- Hosts: archive.stsci.edu with the proxy + certifi bundle in `SSL_CERT_FILE` (one connection reset, retried);
  arXiv e-print served. Map files cached under the gitignored data root; no data committed.
