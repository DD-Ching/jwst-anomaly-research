# W1 shear screen: aperture-coverage cut on grid centres: result (2026-10-17)

Plan: `PLAN.md` (committed before the runs). Labels: coverage, S maps, b_max `derived`; injections `simulated`;
limits `derived` from `simulated`; c = 0.5 and every D-053 threshold are ASSUMPTIONs. Inputs: the five fields'
committed inputs (`inject_shear.py` DEFAULT_FIELDS; every file sha256-checked against SOURCES.md / manifests).

## Diagnostic (`diagnose.json`, cut off)
| Field | S_max coverage | max \|S_×\| | its coverage | valid centres with coverage < 0.5 |
|---|---|---|---|---|
| Abell 2744 | 1.00 | 4.21 | 1.00 | 1.9 % |
| MACS0416 | 0.497 | 4.04 | 0.47 | 5.4 % |
| MACS1149 | 1.00 | 3.76 | 1.00 | 3.6 % |
| Abell 370 | 0.59 | 3.46 | 0.51 | 6.2 % |
| Abell S1063 | 1.00 | 3.33 | 0.68 | 8.8 % |

Only MACS0416's S_max (3.79) and B-mode maximum lie at centres below 0.5. The highest floor (Abell 2744, 4.21) is fully covered,
so the cut cannot lower it.

## Runs (`decide.json`; `combined_c0.json`, `combined_c05.json`, `summary_c05_<field>.json`)
- **Reproduction check:** c = 0 gives the committed recoveries in all five fields exactly (`reproduces_committed`).
- **Real fields at c = 0.5 (rule 1):** no field passes the D-053 rule (MACS0416 S_max 3.79 → 3.76, p_rot 0.05, b_max 4.01;
  the others keep their S_max and b_max). Null; nothing for `/vet-candidate`.
- **b_max changes:** MACS0416 4.04 → 4.01; the other four are unchanged.
- **Recovered (of 200), c = 0 → 0.5:** Abell 2744 8 × 10¹² 72 → 71; MACS0416 7/50/112 → 8/52/111 (2 × 10¹² / 8 × 10¹² /
  2 × 10¹³); Abell S1063 8 × 10¹² 69 → 67, 2 × 10¹³ 91 → 89; MACS1149 and Abell 370 unchanged.
- **Five-field 95 % limits (deg⁻²), c = 0 → 0.5:** 2 × 10¹² 5,765 → 5,672; 8 × 10¹² 1,082 → 1,088; 2 × 10¹³ 567 → 570.
- **Rule 2 (adopt only if lower at each of the three masses): not met.** c = 0 stays the default; the cut is a
  failed approach. The limits change by ≤ 1.6 %; each field's recoveries move by at most 2, below the binomial
  standard deviation of those cells (about 2.6–6.8 of 200; not stored).

## Limits of this result
- One cut value (0.5) tested; coverage uses the footprint definition of the limits (a catalogue row within 4″), not
  the image weight map.
- The B-mode floor is set elsewhere (Abell 2744's 4.21 at full coverage): PSF anisotropy, blends or streaks remain the
  candidates; none is tested here.
