# SMACS 0723 vs CEERS t021 ranked-galaxy features: colours and morphology differ strongly, carried by z_phot cluster members

Plan: `PLAN.md` (pushed before any fetch; rule implemented as `verdict()` with a unit test before the run). Run
`20261106T220924Z-df2b1117` of `configs/reference_sample.yaml` at main `e3eef1a`; script
`scripts/science_control_features.py`; numbers in `summary.json`; plot `feature_cdfs.png` (looked at).

## Samples (derived)
- Ranked stratum (the `scores` table): SMACS 2938 (galaxy stratum; stars split by D-012), CEERS 3156.
- **Deviation:** CEERS `classify_sources` failed (Gaia XMatch table name refused, Gaia TAP 503 / read timeouts),
  so CEERS was ranked as one mixed stratum (stars not split). SMACS classification succeeded.
- Shared features: `ref_mag`, three colours (F150W−F200W, F200W−F277W, F277W−F356W), the F200W morphology set and the
  dropout / gap counts. The plan's bins (`ref_mag` ≤ 24, 24–26, 26–28) hold 1119 of the 2938 SMACS rows; the rest
  is NaN (451) or > 28 (1368); that bin is reported post hoc (`post_hoc_faint_bin`), not judged.

## Frozen rule (ASSUMPTION thresholds: max binned KS D ≤ 0.1 comparable, > 0.2 differ strongly)
| Feature | max binned D | Verdict |
|---|---|---|
| `ref_mag` | 0.118 | intermediate |
| F150W−F200W / F200W−F277W / F277W−F356W | 0.306 / 0.540 / 0.589 | differ strongly |
| `ref_log_isophotal_area` (size) | 0.298 | differ strongly |
| `ref_ellipticity`, `ref_log_ci_50_30`, `ref_log_ci_70_50` | 0.233, 0.372, 0.363 | differ strongly |
| `ref_sharpness` (26–28 only; point-like sources) | 0.485 | differ strongly |
| `ref_roundness` | 0.193 | intermediate |
| `ref_log_nn_dist` | 0.244 | differ strongly (intermediate without the bright bin) |
| dropouts, `n_gaps` | ≤ 0.055 | comparable |

The bright-bin H2 exemption does not change any colour or morphology verdict (`verdict_excl_bright`). By the plan,
the fields "differ strongly". Caveat: at the bright end the bins are small (SMACS 104 / CEERS 55), where a KS D of
about 0.23 is reached by chance at the 5 % level (asymptotic formula; not stored), so bright-bin D values alone are
weak; the AD p values (`ad_p`) are in `summary.json`.

## Post hoc: which population carries the difference
SMACS ranked rows matched to the DJA v7.4 eazy catalogue (0.2″; SOURCES.md, sha256 pinned): 2340. Members
(0.25 ≤ z_phot ≤ 0.55; ASSUMPTION): 614; background (z_phot ≥ 0.8): 1274. Each group vs all CEERS ranked rows, by bin
(`post_hoc_zphot_split`):
- **Members** differ in every colour at every magnitude: F277W−F356W D 0.79–0.86 (median −0.451 to −0.511 vs CEERS
  −0.07 to +0.08), F200W−F277W D 0.52–0.81. This is the expected colour of an old stellar population at z ≈ 0.39
  (the 1.6 µm bump falls in F200W–F277W); in the plot, the SMACS F277W−F356W curve piles up near −0.5. Members are also
  larger at ≤ 26 (area D 0.26–0.31) and rounder at 26–28 and > 28 (ellipticity D 0.28 / 0.30).
- **Background** galaxies match CEERS much more closely: colours D 0.09–0.12 at 24–28 except F277W−F356W at 24–26
  (D 0.28, n 59 vs 130) and > 28 (D 0.19); area D 0.06–0.17; CI_50_30 D 0.08–0.14; ellipticity D 0.05–0.13. The
  ≤ 24 bin has only 19–21 background rows and is not judged.
- Hence a processing difference at fixed magnitude (plan H3) is not supported by the background sample beyond these
  residuals; the strong differences track the cluster-member population (H2), at all magnitudes, not only the
  bright end. Caveat: z_phot is itself colour-based, so selecting members by z_phot partly selects by colour; the
  morphology differences (size, ellipticity) are not subject to that circularity.

## What this changes
- D-022's "Revisit if the two fields' size distributions differ strongly" holds by the frozen rule: no
  science-vs-control top-k statement without matching. A fair control comparison should use SMACS background
  galaxies (z_phot ≥ 0.8) or member-matched weighting, in `ref_mag` bins (TASKS item; D-022 addendum).
- SMACS rankings themselves are unchanged (each sample is ranked on its own).

## Limitations
- CEERS stars not split (deviation above); CEERS has few catalogued stars (TASKS Then 8), so the effect is expected to
  be small but is not measured.
- Pixel scales of the two F200W level-3 catalogs were not checked here (size compares `isophotal_area` in pixels).
- CEERS photo-z were not used (the CEERS DJA photo-z tarball is 396.8 MB); the CEERS redshift mix differs from the
  SMACS background by construction.
- Data: SMACS photo-z member (13 MB) and the DJA catalogues were streamed to scratch and deleted after the run.
