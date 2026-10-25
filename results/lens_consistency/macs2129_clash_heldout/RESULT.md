# MACS J2129: CLASH HST-era maps on the SLICE JWST-era image list (2026-10-25)

Frozen plan `PLAN.md` (pushed in 834486b before any map was read); the pinned maps, `scripts/macs2129_clash_heldout.py`
and its rule tests (`tests/test_macs2129_clash_heldout.py`; the rule functions are those of
`macs0647_clash_heldout`, tested there) were pushed in 3bf8bec before the run. Numbers come from `summary.json`;
per-image residuals are in `imageplane_<model>.ecsv` (provenance `model_prediction`; positions observed, SLICE
Table 8). No imaging was read: this is a position test of published maps against a published image list.

## Frozen result
| | CLASH PIEMD + eNFW | CLASH LTM-Gauss |
|---|---|---|
| Gate K0, families 1–14 (44 images, fitted z) | **fails**: rms 2.05″ (> 2.0″), 91 % within 3″ | **fails**: rms 2.03″, 89 % |
| Frame check F (median observed − predicted, RA·cos δ / Dec) | **fails**: +0.03″ / +0.59″ | **fails**: −0.14″ / +0.74″ |
| Reading A, families 16–18 (9 images, fitted z), reported | rms 3.01″, 7 / 9 within 3″ | rms 3.16″, 7 / 9 |
| Reading B, 13 spectroscopic families at z_spec (41 images), reported | rms 3.20″, 59 % within 3″ | rms 7.96″, 61 % |
| Usable (K0 and F) | no | no |

Rule 4 / 5 outcome: **no usable model, so no image gets a reading** (`flags_A` and `flags_B` are empty). Both
models miss gate K0 by less than 0.06″ in rms and miss check F on the Dec axis.

## Post hoc (not in the plan; no reading)
- Reading A's misses are the same two images in both models: 18.1 (6.45″ / 6.72″) and 18.2 (5.23″ / 4.88″). Source
  18 has SLICE model z 5.53 (no spectrum).
- Fitted z against z_spec (`z_fit`, `z_spec`): PIEMD + eNFW fits a lower z than the spectroscopic one for 12 of 13
  families (all but 12), so at z_spec its images are over-deflected; LTM-Gauss scatters both ways and its z_spec
  rms is dominated by families 3 and 9 (20.9–26.5″, z_spec 2.24).
- The Dec offsets (+0.59″, +0.74″) are of the size of an HST-era vs JWST-era astrometric frame difference; that
  cause is not tested here.

## Reading
- The frozen test gives no lensing-residual input: neither 2015 CLASH model is a usable reference on the SLICE
  image list under the K0 / F gates frozen for MACS0647.
- **Limit:** a frame correction or a looser K0 would be a post hoc change and is not applied. The question needs
  the SLICE Lenstool model (rms 0.79″), which the LAM atlas holds; the atlas was closed on 2026-10-25 (security
  suspension, internal network only).
