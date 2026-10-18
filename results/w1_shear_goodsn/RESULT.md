# W1 shear screen in GOODS-N (JADES 1181, no lens model): result (2026-10-18)

Plan: `PLAN.md` (frozen before any run). Code: the `GOODSN_FIELDS` entries, `passes_gate_p` and `--fixed-r` in
`scripts/inject_shear.py`, the shared-zout staging in `scripts/w1_stage_inputs.py`, `scripts/w1_peak_sheet.py`
(display only). Labels: shapes and pixels `observed`; S, p, efficiencies and limits `derived`; injected lenses
`simulated`; R and every threshold an ASSUMPTION.

## Gate P (`gate_p.json`, decided by `passes_gate_p`)
| Pointing | rows used | with photo-z | share | gate P |
|---|---|---|---|---|
| goodsn-o001 | 2,259 | 1,954 | 0.865 | pass |
| goodsn-o003 | 2,329 | 1,817 | 0.780 | pass |
| goodsn-o004 | 3,187 | 2,265 | 0.711 | pass |
| goodsn-o005 | 3,234 | 2,308 | 0.714 | pass |
| goodsn-o009 | 2,945 | 646 | 0.219 | **fail, dropped** |
| goodsn-o010 | 2,748 | 49 | 0.018 | **fail, dropped** |

The DJA v7.3 GOODS-N zout matches few rows in the two eastern pointings (o009, o010). Why it does not cover them
was not checked. **Deviation:** the counts come from the same runs as the screen, not from an earlier step. The
function uses only row counts, never S.

## Screens (`summaries/z0.1_<field>.json`, `screen`; the screen does not depend on Z)
| Pointing | S_max | p_rot | max \|S_×\| | E/B rule | area (arcmin²) | streak-vetoed rows |
|---|---|---|---|---|---|---|
| goodsn-o001 | 3.29 | 0.735 | 3.43 | no | 11.55 | 0 |
| goodsn-o003 | 4.77 | < 0.005 (0 of 200) | 3.70 | **passes** | 11.47 | 0 |
| goodsn-o004 | 4.35 | < 0.005 (0 of 200) | 3.72 | **passes** | 10.76 | 0 |
| goodsn-o005 | 3.38 | 0.745 | 3.39 | no | 10.77 | 0 |

Streak masks (`masks/*.json`, `*.png`): no status-2 trail in o001–o005. The o004 sheet shows faint diagonal streaks
and row striping that the mask does not cover (seen by eye; not measured).

## The two flags (`sheet_smax_flags.png`, `vet/vet_<field>.json`)
The 40″ F200W sheet was looked at before any test. Both peaks sit on a mosaic edge:
- in o003 a bright star's vertical diffraction spike carries a line of catalogue detections into the aperture;
- in o004 a column of catalogue detections runs along the mosaic edge, next to a bright galaxy.

Under the frozen rule (`w1_shear_more` addendum A):
- **o003 is explained:** without V1 ∪ V2, S = 3.80 and p_rot = 0.075;
- **o004 is explained:** without V1 ∪ V2, S = −0.18 and p_rot = 1.0. Without V1 alone, S = 0.10.

Neither is a candidate.

## Limits (`w1_volume.py`; zero detections; Planck18)
The fourteen stored headline summaries reproduce `results/w1_shear_ncf/limits_14fields.json` exactly (not stored
separately).

| Set | file | area (arcmin²) | n₉₅ at 2 × 10¹² / 8 × 10¹² / 2 × 10¹³ M☉ (h³ Mpc⁻³) | × Takahashi & Asada |
|---|---|---|---|---|
| fourteen fields (before) | `results/w1_shear_ncf/limits_14fields.json` | 97.5 | 3.03 × 10⁻³ / 5.35 × 10⁻⁴ / 2.35 × 10⁻⁴ | 30 / 5.4 / 2.35 |
| **sixteen: + o001, o005 (headline, D-098)** | `limits_16fields_unflagged.json` | 119.8 | 2.46 × 10⁻³ / 4.29 × 10⁻⁴ / 1.82 × 10⁻⁴ | 24.6 / 4.29 / **1.82** |
| eighteen: + o003, o004 (frozen PLAN rule) | `limits_18fields.json` | 142.0 | 1.95 × 10⁻³ / 3.39 × 10⁻⁴ / 1.54 × 10⁻⁴ | 19.5 / 3.39 / 1.54 |
| GOODS-N only (four) | `limits_goodsn_only.json` | 44.6 | 5.45 × 10⁻³ / 9.27 × 10⁻⁴ / 4.48 × 10⁻⁴ | 54.5 / 9.27 / 4.48 |

**Deviation (headline):** the plan let a flag that the frozen rule explains enter the headline. D-098 admits only
blank fields whose screen does not pass the E/B rule. D-098 wins, so the headline is sixteen fields. The
eighteen-field value is reported, not adopted.

Plan check ("what would change the plan"):
- the headline does not reach 10⁻⁴ h³ Mpc⁻³ at 2 × 10¹³ M☉;
- V_eff at 2 × 10¹³ M☉ grows by 29 % with o001 and o005 (41,089 → 53,057 Mpc³), above 20 %, so blank-field additions
  continue (T4).

## Responsivity sensitivity (`r_sensitivity.json`, Z = 0.3, four pointings)
The same injections with R = 0.231 (the MACS0647 value) recover:

| Mass (M☉) | R = 0.414 | R = 0.231 |
|---|---|---|
| 2 × 10¹² | 33 | 7 |
| 8 × 10¹² | 152 | 32 |
| 2 × 10¹³ | 254 | 67 |

That is 0.21–0.26 of the recoveries in the four-pointing totals, and 0.13–0.33 per pointing and mass (11 cells with a non-zero adopted count; `ratio` in the file). **The GOODS-N contribution depends strongly on the assumed R.** R cannot be
fitted in a blank field, and program 1181 has no cluster pointing to fit it on (a deviation from D-098's
same-program rule; D-099).

## Limitations
- R is borrowed from another program (CANUCS Abell 370 F200W fit).
- The o004/o005 overlap (≈ 0.5 arcmin² of MAST footprint) is double-counted in the headline: less than 1 % of the area.
- The two eastern pointings are lost to photo-z coverage.
- The faint unmasked o004 streaks were found by eye only.
- The summaries store absolute `catalog` / `photoz` paths of the session machine; resolve them through `inject_shear.field_spec` and the data root (as `w1_peak_sheet.py` does). The zout is the pinned DJA member.
- Gate P counts rows with a valid z160 (stricter than the screen's finite z_phot; `inject_shear.py`).
