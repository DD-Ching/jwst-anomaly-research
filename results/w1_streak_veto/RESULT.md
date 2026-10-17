# W1 calibrated streak veto: result (2026-10-17)

Plan: `PLAN.md` (committed before the runs; addendum A, a revised mask normalisation, committed after looking at
the first two masks and before any screen or injection run: **deviation** from the first frozen rule). Labels: pixels
and shapes `observed`; masks, S maps and limits `derived`; injections `simulated`; every threshold an ASSUMPTION.

## Masks (`masks/<field>_streak_mask.{fits.gz,json,png}`; `scripts/w1_streak_mask.py`)
STScI `acstools.findsat_mrt` (median Radon transform) on the streamed F200W mosaic (block 4, clip ±10σ, footprint
eroded 1.5″, rebin 2; status-2 trails only). Persistent trails per field: Abell 2744 0, MACS0416 2, MACS1149 0,
Abell 370 2, Abell S1063 3, SMACS 0723 3, El Gordo 0, MACS0417 4, WHL0137 1, RXJ2129 1, MACS0647 9, MACS1423 3.
Sheets looked at: the trails are vertical stripes and diagonal streaks or spikes. Missed by eye: a faint diagonal
streak in MACS1423, the fainter crossing streak at the MACS0647 S_max, and faint parallel diagonal streaks in the
S1063 S_max aperture (`sheet_abells1063_veto_smax.png`). False trails: the MACS0647 mosaic gap and halo (PLAN addendum A).

## Known cases and adoption (`masks/adoption.json`, computed by `w1_streak_mask.evaluate`)
- K1: the MACS1423 mask covers the 2026-10-17 hit position: **yes**.
- K2: MACS1423 with the veto: S_max 3.41, p_rot 0.215, b_max 3.63, D-053 rule **not passed** (was 3.88, 0.01, passed).
- **Adopted** (K1 and K2 hold).

## Screen and injections with the veto (`inject_summary_<field>.json`; seed 50, 200 per mass)
| Field | rows vetoed | S_max | p_rot | b_max | rule | recovered at 2e12 / 8e12 / 2e13 M☉ |
|---|---|---|---|---|---|---|
| Abell 2744 | 0 | 3.39 | 0.58 | 4.21 | no | 14 / 72 / 146 |
| MACS0416 | 224 | 3.83 | 0.025 | 4.04 | no | 5 / 48 / 105 |
| MACS1149 | 0 | 3.40 | 0.19 | 3.76 | no | 3 / 30 / 74 |
| Abell 370 | 599 | 3.30 | 0.375 | 3.72 | no | 7 / 20 / 50 |
| **Abell S1063** | 187 | 3.61 | 0.045 | 3.33 | **yes** | 15 / 65 / 103 |
| SMACS 0723 | 382 | 2.88 | 0.79 | 3.19 | no | 5 / 51 / 114 |
| El Gordo | 0 | 2.95 | 0.62 | 3.25 | no | 0 / 14 / 48 |
| MACS0417 | 401 | 3.50 | 0.095 | 3.42 | no | 10 / 46 / 97 |
| WHL0137 | 78 | 2.75 | 0.68 | 3.16 | no | 1 / 25 / 41 |
| RXJ2129 | 12 | 2.91 | 0.375 | 3.00 | no | 8 / 54 / 92 |
| MACS0647 | 381 | 2.89 | 0.315 | 3.27 | no | 1 / 5 / 13 |
| MACS1423 | 312 | 3.41 | 0.215 | 3.63 | no | 7 / 66 / 119 |

No recoveries at 2 × 10¹⁰ or 2 × 10¹¹ M☉ in any field. "rows vetoed" counts all catalogue rows dropped, not only
rows the screen uses (the used-row share was not computed; **deviation** from PLAN "Reported regardless"); no field
exceeds the 20 % rule on all rows (largest: Abell 370, 599 of 5,266 rows staged, 11 %; row count from the session, not stored).
Against the earlier summaries of the same fields, the veto changes the recoveries in both directions (e.g. MACS1423
9 / 72 / 140 → 7 / 66 / 119; MACS0647 0 / 1 / 1 → 1 / 5 / 13; WHL0137 1 / 12 / 44 → 1 / 25 / 41). The changes also
include the different rotation-null draws, so the veto's own cost is not separated.

**Abell S1063 passes the D-053 rule with the veto.** S_max (3.61), its position and b_max (3.33) are the same as
without the veto (`results/w1_shear_s1063/inject_summary_abells1063.json`); only p_rot moved, 0.07 → 0.045, at the
0.05 threshold. Its aperture holds faint parallel diagonal streaks that the mask misses (sheet, by eye). Per PLAN it
is left out of the headline and goes to vetting (next run). It is not a candidate.

## Limits (`combined_11fields_veto.json`; `combined_12fields_veto.json` with S1063 for reference only)
- **Headline, eleven fields with the veto (S1063 out, MACS1423 in), 62.5 arcmin²:** Σ(W1) < 5.3 × 10³ / 8.0 × 10² /
  3.9 × 10² deg⁻² at 2 × 10¹² / 8 × 10¹² / 2 × 10¹³ M☉ (border-corrected 6.0 × 10³ / 9.2 × 10² / 4.4 × 10²).
- Earlier eleven-field limits without a veto (`results/w1_shear_rxj_m0647/combined_11fields.json`, S1063 in,
  MACS1423 out): 4.5 × 10³ / 7.6 × 10² / 3.8 × 10². The new headline is weaker at 2 × 10¹² M☉ because S1063 left;
  it is the first W1 limit that runs a streak test through the injections.
- Still blind at ≤ 2 × 10¹¹ M☉.

## Limits of this result
- The mask finds bright and long trails. It misses faint parallel streaks, such as those at the S1063 peak.
- The veto rule was revised once after the first two masks were looked at (addendum A).
- p_rot rests on 200 rotation draws, so a peak at the threshold can flip with the draws (S1063).
- Inputs were staged by `scripts/w1_stage_inputs.py` (sha256 checked where recorded; the RXJ2129 photo-z tarball
  has no recorded hash, its member does).
- The mask FITS files are not tracked (`*.fits.gz` is gitignored): re-make them with
  `scripts/w1_streak_mask.py --fields <field> --out results/w1_streak_veto/masks` (about 30–60 s per field) and
  check `mask_pixels` against the committed JSON before `inject_shear.py --streak-masks`.

## Open code-review findings (not fixed in this run; they need re-runs)
- **Spike-veto coupling:** streak rows are dropped before `shear_screen`, so the catalogue spike-segment veto and
  `psf_sigma_px` are estimated from the reduced table. Some of the S_max / b_max / p_rot changes above may come
  from a changed spike veto, not from the streak rows. Next run: apply the streak veto as a row flag after the spike
  estimate, then re-run the twelve fields.
- The used-row share for the 20 % rule is not computed: the summaries need n_used without and with the veto.
- The veto distance is measured from the nearest mask-pixel centre, so it is exact only to about half a mask
  pixel (≈ 0.18″) (documented in `veto_function`).
- `w1_stage_inputs.py` fetches the files without a recorded sha256 through a sentinel that compares equal to any
  hash. The fetch helpers should get an explicit no-verify path instead.
- `w1_streak_mask.py` streams each mosaic one strip at a time, with no concurrent range reads and no
  download/compute overlap (scripts/CLAUDE.md parallel topology).
