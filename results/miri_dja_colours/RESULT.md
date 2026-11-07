# SMACS MIRI colours from DJA fixed apertures: result (2026-11-07)

Plan: `PLAN.md` (frozen before measuring). Script: `scripts/miri_dja_colours.py` (rules V, C, B and `decide` are
unit-tested in `tests/test_miri_dja_colours.py`). Numbers: `summary.json`. Contact sheet: `contact_sheet.png`.

**Runs (scratch outputs, not stored).** `20261107T151133Z-675b4a49` (`configs/smacs_jwst3_miri_dja.yaml`) and the
baseline `20261107T151133Z-0490f84e` (`configs/smacs_jwst3_miri_confirm.yaml`), the same MAST L3 jwst 3.0.0 MIRI
catalogues and DJA v7.4 `fix_phot.fits` (sha256 361a4eb9…). CDS (SIMBAD, XMatch) was unreachable in both runs; the
Gaia TAP route classified stars, as in the earlier runs today. 160 sources ranked in each run.

**Frozen rules (ASSUMPTIONs) and outcome (derived).**
| Rule | Measure | Value | Threshold | Outcome |
|---|---|---|---|---|
| V | MAD of (`dja07` − `aper50`) F770W−F1000W, lowest area quartile (32 rows, area ≤ 42 px) | 0.058 (median −0.031) | ≤ 0.15 | pass |
| C | ranked sources with a finite `dja07` F770W−F1000W | 134 / 160 = 0.8375 | ≥ 0.80 | pass |
| B | Spearman(log F770W isophotal area, F770W−F1000W), 124 common rows | `aper50` +0.067 (p 0.46); `dja07` +0.053 (p 0.56); isophotal +0.152 (p 0.093) | `aper50` \|ρ\| ≥ 0.10 needed for a bias | no bias |

**Verdict (frozen function): not adopted.** The MIRI pipeline `aper50` colours show no size dependence above the
0.10 floor in F770W−F1000W, so there is no D-013-type bias for fixed apertures to remove; the ranking keeps the
pipeline colours.

**Reported either way.** Top-20 Jaccard vs the baseline 0.60: new `f770w_193`, `_109`, `_177`, `_231`, `_161`;
dropped `f770w_78`, `_136`, `_246`, `_156`, `_38` (`summary.json`).

**Post hoc (not in the plan's rules).**
- F1000W−F1500W (47 rows): ρ −0.335 (p 0.021) `aper50`, −0.367 (p 0.011) `dja07`, −0.254 isophotal. A fixed aperture
  does not reduce it. The compact-quartile offset (`dja07` − `aper50`) is −0.292 here (MAD 0.041), i.e. the 0.7″
  aperture is not PSF-matched between F1000W and F1500W. The size trend's cause is not tested (hypothesis: a real
  population trend; a PSF-size effect would make compact sources bluer, i.e. ρ > 0, the opposite sign).
- Contact sheet of the `dja07` run's top 20 looked at: the largest `color_f770w_f1000w` z-scores (+7.6 to +11.5) are
  F770W−F1000W ≈ 0.9–1.55 mag against a ranked median −0.37 (MAD 0.11; not stored), and they agree with `aper50` to
  ≤ 0.13 mag where both exist (11 rows) (not stored): red compact sources and galaxies, the same population as the
  baseline. No new entry is vetted; none is an anomaly claim.

**Limitations.** One field; 124 rows for the main test; DJA MIRI apertures are not PSF-homogenised.
