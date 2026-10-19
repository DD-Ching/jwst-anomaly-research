# E1 + Swift BAT: result (2026-10-18)

Plan: `PLAN.md` (frozen first). Code: `scripts/e1_swift.py` (frozen rules `control_passes`, `outcome`, unit-tested
in `tests/test_e1_swift.py`), post hoc diagnostics `scripts/e1_swift_diag.py`. Inputs pinned in
`data/manifests/e1_swift.ecsv` (all five re-pinned 2026-10-18; the live `fermigbrst` no longer matches D-074's pin).

**Events (observed, `summary.json`):** BAT 1,635 (on-board Swift triggers, MJD 53359.96–61318.90), GBM 4,398,
ICECAT 340, GW 391, CHIME 3,641 (one node per repeater).

**Positive control (frozen rule passes):** BAT–GBM same-direction pairs at 0–100 s: 466 vs a jit-null mean of 0.10
(sd 0.31); analytic p underflows to 0.

**Tested family (35 cells, 10⁴ jit + 10⁴ jitday scrambles; `counts.ecsv`): outcome "null".**
- Smallest cell: BAT–CHIME 10–100 s wide, 4 pairs vs 1.33 ± 1.16; analytic p 0.047; Bonferroni p 1.0.
- Pooled empirical global p 0.54 (jit), from `summary.json`.
- Same-direction pairs within 7 d (`same_dir_pairs_*.ecsv`): BAT–BAT 0; BAT–ICECAT 2 (both with IC190629A,
  σ_comb 19.5°); BAT–CHIME 1 (GRB 220412A / FRB20220415D, 2.8 d, 0.52° at σ_comb 0.25°). Nothing to vet: no
  cell passes the rule.

**Limits (`limits.ecsv`, D-074 method; all 20 rows valid):** 95 % upper limits on dependent wide (GW: any-separation)
pairs per eligible anchor:

| Channel | 0–10 s | 10–100 s | 100 s–1 h | 1 h–1 d | 1–7 d |
|---|---|---|---|---|---|
| BAT–BAT (per BAT burst) | 1.8 × 10⁻³ | 1.8 × 10⁻³ | 1.9 × 10⁻² | 1.7 × 10⁻² | 6.3 × 10⁻² |
| BAT–ICECAT (915 anchors) | 3.3 × 10⁻³ | 3.3 × 10⁻³ | 8.5 × 10⁻³ | 3.1 × 10⁻² | 3.2 × 10⁻² |
| BAT–CHIME (343 anchors) | 8.7 × 10⁻³ | 2.7 × 10⁻² | 6.6 × 10⁻² | 0.16 | 0.35 |
| BAT–GW (255 anchors) | 1.2 × 10⁻² | 1.2 × 10⁻² | 3.7 × 10⁻² | 6.3 × 10⁻² | 0.10 |

**Control-channel observations (outside the family; post hoc):**
- 8 BAT–GBM pairs at 0–10 s fall in the "wide" class (`counts.ecsv`; 0.22 expected). By name and time they are the
  same GRBs (e.g. GRB 090423 / bn090423330); their GBM rows carry an external position (error_radius < 0.5°, σ
  floored at 0.05°), and the separations are 0.17–1.02°. This is a position-convention effect of the D-074 σ
  floor, not a dependence.
- BAT–GBM all-separation pairs show a deficit at 1 h–4 d (1 h–1 d: 1,603 vs 1,784.3, z −4.9; `counts.ecsv`). It
  is in all 19 BAT years (`control_deficit_by_year.json`), symmetric in sign and gone beyond 8 d
  (`control_lag_curve.json`; ratio 0.75 at 1–3 h). A null that keeps each satellite event's position relative to
  the Sun leaves it (1 h–1 d wide z −5.0, `control_sun_frame_null.json`).
  With the 466 joint bursts removed from both catalogues, the 1 h–4 d sum is 3,844 vs 3,863.7 ± 37.9
  (`control_lag_curve_no_counterparts.json`): the deficit is carried by the joint bursts. Why joint bursts carry
  it is not determined (not tested).
  Update 2026-10-19 (post hoc, `selfpair/RESULT.md`, D-074 addendum 10): the long-lag deficit behaves as the
  jit null displacing the joint bursts predicts; frozen rule "undetermined" (100 s–1 h cell). No tested channel includes GBM.

**Limitations:** the BAT 1σ comes from the BAT 90 % radius only (XRT positions not used). The Swift-orbit quantum
(5,700 s, ASSUMPTION) keeps orbit phase only approximately over ±3 d. GW channels have no sky positions (D-074).
The tested cells are one-sided upper tails; deficits are not tested.
