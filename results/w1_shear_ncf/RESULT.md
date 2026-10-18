# W1 shear screen in the CANUCS NIRCam flanking fields (no cluster model): result (2026-10-18)

Plan: `PLAN.md` (frozen before any run). Code: `inject_radial.NullLensModel`, the NCF entries in
`scripts/inject_shear.py`, `scripts/w1_ncf_gate_r.py`, `scripts/w1_peak_vet.py` (post hoc V1r and its control
added). Labels: shapes and pixels `observed`; R, S, p, efficiencies and limits `derived`; injected lenses
`simulated`; every threshold an ASSUMPTION.

## Gate R: passed (`gate_r.json`)
F150W isophotal moments on the CANUCS cluster pointings, with the Gaia spike veto as in the F200W fits: R = 0.425 ±
0.035 (MACS0417) and 0.618 ± 0.037 (MACS1423), against the stored F200W 0.488 and 0.604 (|ΔR| ≤ 0.15). Adopted
R = 0.425 (`adopted_r`). **Deviation (corrected before the final grid):** the first gate run left the spike veto out
(R 0.421 / 0.614; found in code review); every number below comes from the re-run grid with R = 0.425.

## Screens (`summaries/z0.3_<field>.json`, `screen`; the screen does not depend on Z)
| Field | rows used | S_max | p_rot | max \|S_×\| | E/B rule | area (arcmin²) |
|---|---|---|---|---|---|---|
| macs0416-ncf | 3,472 | 4.64 | 0.005 | 4.33 | **passes** | 10.50 |
| macs1149-ncf | 3,379 | 3.45 | 0.64 | 3.81 | no | 10.48 |
| abell370-ncf | 3,155 | 3.94 | 0.065 | 4.81 | no | 10.33 |
| macs0417-ncf | 4,040 | 6.02 | < 0.005 (0 of 200) | 4.01 | **passes** | 10.62 |
| macs1423-ncf | 3,643 | 4.03 | 0.06 | 3.82 | no | 10.36 |

Streak masks (`masks/*.json`, `*.png`; 0–4 status-2 trails per field; the sheet titles say F200W because the title
is hard-coded in `w1_streak_mask.py`, the image is F150W). The macs1423-ncf mask misses a faint diagonal streak in
its upper-left module (seen on the sheet).

## The two flags (`sheet_smax_flags.png`, `vet/vet_<field>.json`)
Looked at on the 40″ F150W sheet before any test: the macs0416-ncf aperture is crossed through its centre by a
diagonal streak that the trail mask missed; the macs0417-ncf aperture is crossed by horizontal stripes and a chip
gap.
- **macs0416-ncf: explained under the frozen rule** (`w1_shear_more` addendum A): without V1 ∪ V2 (V2: 12 sources
  pointing at the brightest catalogue source, 9.3″ from the peak) S = 3.53, p_rot 0.575.
- **macs0417-ncf: not explained under the frozen rule** (S 5.04 without V1 ∪ V2, p_rot < 0.005): V1 looks for
  stripe *columns*, and these stripes run along rows. **Post hoc** V1r (the same test on rows; 86 stripe rows in the
  30″ cutout, 48 used sources within 1″): S = 1.35, p_rot 1.0. Control (`control_v1r_rotated`): the same row set
  rotated about the peak by 200 random angles gives median S 5.84, and S ≤ 1.35 in 9 % of the angles (angles near
  the stripes' own orientation are among them by construction). The low S needs the band along the stripes, which
  supports the stripe explanation only weakly; the rows next to the chip gap may be flagged by partial-row medians
  (not tested). Kept out of the headline; it goes to `/vet-candidate` only if a frozen row/diagonal streak veto
  leaves it standing.

## Injections (`summaries/z<Z>_<field>.json`; 100 lenses per field, mass and Z; seed 50)
Recovered of 100, at 2 × 10¹² / 8 × 10¹² / 2 × 10¹³ M☉ (fields in the order of the table above):

| Z | macs0416-ncf | macs1149-ncf | abell370-ncf | macs0417-ncf | macs1423-ncf |
|---|---|---|---|---|---|
| 0.1 | 41 / 75 / 66 | 49 / 70 / 69 | 19 / 55 / 53 | 42 / 80 / 68 | 30 / 72 / 52 |
| 0.3 | 9 / 26 / 64 | 19 / 48 / 74 | 9 / 20 / 54 | 11 / 46 / 72 | 14 / 46 / 74 |
| 0.6 | 2 / 18 / 33 | 5 / 23 / 45 | 2 / 7 / 15 | 5 / 30 / 46 | 3 / 15 / 39 |
| 1.0 | 0 / 3 / 13 | 2 / 7 / 19 | 0 / 4 / 5 | 1 / 12 / 19 | 1 / 6 / 18 |
| 1.5 | 0 / 2 / 3 | 0 / 2 / 8 | 0 / 0 / 2 | 0 / 4 / 5 | 0 / 1 / 3 |

## Volume density limits (95 %, zero detections, Planck18; `limits_*.json`)
The eleven stored cluster summaries reproduce `results/w1_volume/limits.json` exactly (checked in the session; not
stored). n₉₅ in h³ Mpc⁻³ and the ratio to Takahashi & Asada's 10⁻⁴:

| Set | Area (arcmin²) | 2 × 10¹² M☉ | 8 × 10¹² M☉ | 2 × 10¹³ M☉ |
|---|---|---|---|---|
| eleven cluster fields (`results/w1_volume`) | 66.3 | 6.5 × 10⁻³ (65×) | 8.5 × 10⁻⁴ (8.5×) | 3.65 × 10⁻⁴ (3.6×) |
| **14: + the three NCF without a flag (headline)** | 97.5 | 3.03 × 10⁻³ (30×) | 5.35 × 10⁻⁴ (5.4×) | 2.35 × 10⁻⁴ (2.35×) |
| 15: + macs0416-ncf (flag explained, frozen rule) | 108.0 | 2.75 × 10⁻³ (27×) | 4.81 × 10⁻⁴ (4.8×) | 2.11 × 10⁻⁴ (2.1×) |
| 16: + macs0417-ncf (flag explained post hoc, weakly) | 118.6 | 2.32 × 10⁻³ (23×) | 3.97 × 10⁻⁴ (4.0×) | 1.86 × 10⁻⁴ (1.9×) |

The headline keeps both flagged fields out, as `results/w1_shear_more` kept its explained MACS1423 hit out.
NCF-only V_eff (`limits_ncf_only.json`) is 1.82 / 1.15 / 0.97 times the eleven cluster fields' at the three masses,
from 52.3 against 66.3 arcmin² (the cause of the higher efficiency per area is not measured).

## Against the plan's criteria
- "Within 2× at 2 × 10¹³ M☉ → blank deep fields next": **not met by the headline** (2.35×); met only by the
  16-field set (1.9×), which rests on a weak post hoc explanation.
- "NCF add < 20 % to V_eff at every mass → stop": not met (headline V_eff +116 % / +59 % / +55 %). Blank fields are the
  cheapest W1 area so far; continuing is allowed.

## Limitations
- R for the NCF is borrowed from the cluster pointings (ASSUMPTION; gate R bounds the band change only).
- The NCF residual cluster shear (3–6′ from the cores) is not removed; assumed smooth over a 10″ aperture.
- macs0417-ncf: horizontal-stripe rows are not in any frozen veto; V1r is post hoc and its 86 rows cover ~9 % of
  the cutout's rows; its rotated-band control is weak (9 %).
- Single-plane, Z-grid interpolation and 100-lens binomial errors as in `results/w1_volume`.
