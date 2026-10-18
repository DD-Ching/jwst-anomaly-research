# W1 shear screen in CEERS (EGS, program 1345, no lens model): result (2026-10-18)

Plan: `PLAN.md` (frozen before any run). Code: the `CEERS_FIELDS` entries in `scripts/inject_shear.py`, the DJA
v7.4 zout and catalogue pins in `scripts/w1_stage_inputs.py`, `scripts/w1_footprints.py` (pointing choice by
`drop_overlapping`, `footprints.json`; the function was written after the choice by hand and agrees with it, post
hoc). Labels: shapes and pixels `observed`; S, p, efficiencies and limits `derived`; injected lenses
`simulated`; R and every threshold an ASSUMPTION.

## Gate P (`gate_p.json`, decided by `passes_gate_p`)
All ten pointings pass (photo-z share 0.695–0.813). **Deviation:** as in GOODS-N, the counts come from the Z = 0.1
runs that also hold the screen; the function uses only row counts, never S.

## Screens (`summaries/z0.1_<field>.json`, `screen`; the screen does not depend on Z)
| Pointing | S_max | p_rot | max \|S_×\| | E/B rule | area (arcmin²) | streak-vetoed rows |
|---|---|---|---|---|---|---|
| ceers-o001 | 3.99 | 0.05 | 4.13 | no | 10.49 | 232 |
| ceers-o002 | 3.38 | 0.51 | 3.89 | no | 10.39 | 528 |
| ceers-o003 | 3.57 | 0.295 | 3.94 | no | 10.14 | 28 |
| ceers-o004 | 4.96 | < 0.005 (0 of 200) | 4.17 | **passes** | 10.12 | 56 |
| ceers-o061 | 3.32 | 0.555 | 3.42 | no | 9.66 | 124 |
| ceers-o064 | 3.31 | 0.64 | 3.78 | no | 9.82 | 24 |
| ceers-o066 | 3.40 | 0.62 | 3.35 | no | 10.06 | 15 |
| ceers-o068 | 3.58 | 0.17 | 3.77 | no | 9.74 | 99 |
| ceers-o069 | 4.07 | 0.02 | 3.62 | **passes** | 10.01 | 0 |
| ceers-o071 | 3.21 | 0.725 | 3.21 | no | 9.79 | 200 |

Streak masks (`masks/*.json`, `*.png`): status-2 trails in nine of ten pointings (none in o069). The o002 sheet
shows two diffuse scattered-light patches inside detectors that the mask does not cover (seen by eye; not measured).

## The two flags (`sheet_smax_flags.png`, `vet/vet_<field>.json`)
The 40″ F200W sheet was looked at before the vetting numbers were read:
- in o004 the peak sits on the edge of a detector gap lined with catalogue detections, with a bright star's
  diffraction spikes inside the aperture (the GOODS-N o004 pattern);
- in o069 a horizontal spike-like line with detections along it runs through the peak.

Under the frozen rule (`w1_shear_more` addendum A):
- **o004 is explained:** without V1 ∪ V2, S = 1.54 and p_rot = 1.0;
- **o069 is explained:** without V1 ∪ V2, S = 1.85 and p_rot = 1.0 (V2 carries it: without V1 alone S = 4.14).

Neither is a candidate. Both stay out of the headline (D-098). **Deviations:** the plan recorded the tarball sha256 at the first fetch; only the member sha256 was
recorded (the archive was streamed unpinned). The plan ran injections on every pointing that passes gate P; the flagged two were not injected, because D-098 keeps them out of every limit here.

## Limits (`w1_volume.py`; zero detections; Planck18)
The sixteen stored headline summaries reproduce `results/w1_shear_goodsn/limits_16fields_unflagged.json` exactly
(checked in the session, not stored separately).

| Set | file | area (arcmin²) | n₉₅ at 2 × 10¹² / 8 × 10¹² / 2 × 10¹³ M☉ (h³ Mpc⁻³) | × Takahashi & Asada |
|---|---|---|---|---|
| sixteen fields (before) | `results/w1_shear_goodsn/limits_16fields_unflagged.json` | 119.8 | 2.46 × 10⁻³ / 4.29 × 10⁻⁴ / 1.82 × 10⁻⁴ | 24.6 / 4.29 / 1.82 |
| **twenty-four: + eight unflagged CEERS (headline, D-098)** | `limits_24fields_unflagged.json` | 199.9 | 1.50 × 10⁻³ / 2.53 × 10⁻⁴ / 1.09 × 10⁻⁴ | 15.0 / 2.53 / **1.09** |
| CEERS only (eight) | `limits_ceers_only.json` | 80.1 | 3.83 × 10⁻³ / 6.18 × 10⁻⁴ / 2.69 × 10⁻⁴ | 38.3 / 6.18 / 2.69 |

The headline double-counts the pairwise overlaps among the eight: at least 1.33 arcmin², the sum of their pairs
listed in `footprints.json` (pairs ≤ 0.05 arcmin² are not listed). The area is the screened grid area, the same
definition as every earlier field; for the eight it is 80.1 arcmin², 7 % more than their MAST footprints (74.9
arcmin², `footprints.json`; sum not stored), so the CEERS V_eff may be overstated by up to that much.

Plan check ("what would change the plan"):
- the headline does not reach 10⁻⁴ h³ Mpc⁻³ at 2 × 10¹³ M☉ (1.09 × 10⁻⁴);
- V_eff at 2 × 10¹³ M☉ grows by 68 % (53,057 → 88,969 Mpc³), above 20 %, so blank-field additions continue (T4).

## Responsivity sensitivity (`r_sensitivity.json`, Z = 0.3, eight pointings)
| Mass (M☉) | R = 0.414 | R = 0.231 | ratio |
|---|---|---|---|
| 2 × 10¹² | 68 | 19 | 0.279 |
| 8 × 10¹² | 270 | 51 | 0.189 |
| 2 × 10¹³ | 435 | 101 | 0.232 |

Per pointing and mass the ratio spans 0.0–0.5 (24 cells with a non-zero adopted count). **The CEERS contribution,
like GOODS-N's, depends strongly on the borrowed R.**

## Limitations
- R is borrowed from another program (CANUCS Abell 370 F200W fit); CEERS F200W is shallower (2.8 ks) than JADES 1181.
- Overlaps among the eight (≥ 1.33 arcmin²) are double-counted; the screened area exceeds the MAST footprints by 7 %.
- The headline is 1.09× the Takahashi & Asada bound at 2 × 10¹³ M☉ under that R; with R = 0.231 the recoveries at Z = 0.3 fall to
  0.19–0.28 (CEERS, eight-pointing totals) and 0.21–0.26 (GOODS-N, four-pointing totals) of the adopted ones.
