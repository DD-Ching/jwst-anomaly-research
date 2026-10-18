# E1 + IceTracks-DR2: result (2026-10-18)

Plan: `PLAN.md` (frozen first). Code: `scripts/e1_icetracks.py` (frozen rules `control_passes`, `outcome`,
unit-tested in `tests/test_e1_icetracks.py`). Inputs: `data/manifests/nf_icetracks_dr2.ecsv` (D-079) and
`data/manifests/e1_swift.ecsv`, both unchanged.

**Events (observed, `summary.json`):** DR2 7,583 northern tracks (Dec > −5°, log10 E ≥ 4.0, one per readout;
MJD 54562.7–59722.4, 2008-04 to 2022-05), GBM 4,398, BAT 1,635, ICECAT 340, GW 391, CHIME 3,641.

**Positive control (frozen rule passes):** DR2–ICECAT same-direction pairs at 0–10 s: 270 vs a jit-null mean of
0.0013 (analytic p underflows to 0). 270 distinct ICECAT alerts and 270 distinct DR2 tracks
(`same_dir_pairs_DR2-ICECAT.ecsv`): the alerts are recovered one-to-one at the same time.

**Tested family (35 cells, 10⁴ jit + 10⁴ jitday scrambles; `counts.ecsv`): outcome "null".**
- Smallest cell: DR2–CHIME 1–7 d same, 35 pairs vs 26.51 ± 5.11; analytic p 0.065; Bonferroni p 1.0.
- Pooled empirical global p 0.79 (jit), 0.78 (jitday). Largest deficit among tested cells: DR2–BAT 100 s–1 h
  wide, 113 vs 133.0 (z −1.7).
- Short lags (observed): DR2–GBM same 0 pairs at 0–100 s (0.47 expected), DR2–BAT 0 (0.00), DR2–CHIME 0 (0.01),
  DR2–GW 0 at 0–100 s (0.13). Same-direction pairs within 1 h: DR2–GBM 19 (16.8 expected in the 0–1 h same cells),
  DR2–BAT 0, DR2–CHIME 0 (`same_dir_pairs_*.ecsv`). No cell passes the rule, so nothing is vetted.

**Limits (`limits.ecsv`, D-074 method; 19 of 20 rows valid, DR2–GW 1–7 d not reached):** 95 % upper limits on
dependent wide (GW: any-separation) pairs per eligible anchor:

| Channel | 0–10 s | 10–100 s | 100 s–1 h | 1 h–1 d | 1–7 d |
|---|---|---|---|---|---|
| DR2–GBM (7,321 anchors) | 4.1 × 10⁻⁴ | 2.7 × 10⁻³ | 4.4 × 10⁻³ | 2.0 × 10⁻² | 4.8 × 10⁻² |
| DR2–BAT (7,583) | 4.0 × 10⁻⁴ | 1.2 × 10⁻³ | 2.5 × 10⁻³ | 1.2 × 10⁻² | 3.6 × 10⁻² |
| DR2–CHIME (769) | 3.9 × 10⁻³ | 8.2 × 10⁻³ | 2.3 × 10⁻² | 0.10 | 0.19 |
| DR2–GW (435) | 6.9 × 10⁻³ | 6.9 × 10⁻³ | 1.8 × 10⁻² | 6.0 × 10⁻² | — |

**Control-channel observations (outside the family; post hoc; found in review):**
- 6 DR2–ICECAT pairs at 0–10 s fall in "wide" (0.05 expected): the same neutrinos with DR2 and ICECAT positions
  more than 3σ_comb apart (different reconstructions and error conventions). Not a dependence.
- DR2–ICECAT wide pairs show a **deficit** at long lags: 1 h–1 d 334 vs 415.5 ± 18.1 (z −4.5), 1–7 d 2,057 vs
  2,241.7 ± 28.3 (z −6.5) (`counts.ecsv`). Six distinct null variants (`scripts/e1_icetracks_diag.py`,
  `control_null_variants.json`, 10³ scrambles each; its `both` and `in_uptime` select the same 303 alerts) leave it at z −4.4 to −5.3 and −6.1 to −7.1: dropping alerts
  after the last DR2 track, keeping only alerts inside DR2 good runs (303 of 340), scrambling ICECAT with the DR2
  uptime jitter, and moving either side alone. So it is not downtime, the DR2 end date or one side's scramble;
  around an ICECAT alert time there are fewer DR2 tracks within 1 h–7 d than the null predicts. Cause not
  determined (not tested further). Like addendum 6's BAT–GBM deficit, it appears in a channel whose two catalogues
  share events. No tested channel shows a deficit beyond z −1.7.

**Limitations:** `jitday` scrambles DR2 exactly as `jit` (only the partners change), so DR2 time-of-day
structure is not separately tested; the null is not validated against a same-detector channel (the control deficit above); the cut (log10 E ≥ 4, north) keeps mostly atmospheric neutrinos (D-079: f_astro not computed), so
the limits are per track above the cut, not per astrophysical neutrino. AngErr is used as 1σ without a
systematic floor. CHIME overlaps DR2 only 2018-08 to 2022-05 (769 anchors). The ±3 d uptime jitter keeps each
track's hour angle and good-run membership, not seasonal rate structure finer than ±3 d. One-sided upper tails
only.
