# SMACS `f356w_1243`: two-epoch F444W check (2736 o001 → VENUS 6882 o057)

Plan: `PLAN.md` (frozen in the claim commit; the two nearest-row magnitudes were seen during the coverage check
and are disclosed there). Script `scripts/smacs_epochs_1243.py`; `decide` unit-tested before the run
(`tests/test_smacs_epochs_1243.py`). Numbers: `summary.json`. Contact sheet looked at (`contact_sheet.png`).

## Data (observed)
- Level-3 F444W CATs read from `s3://stpubdata` (sha256 in `summary.json`); 2022-06 and 2026-06-05, baseline 4.0 yr.
- F444W is the only band with a target row in both epochs; 12774 o002 does not cover the target (PLAN.md).

## Result (derived; F444W only)
| Quantity | Value |
|---|---|
| Frame tie | 241 references (S/N ≥ 30, ≤ 60″), offset (4.32, 20.54) mas |
| Target S/N | 121.7 (2736), 48.4 (6882) |
| Target Δm (minus reference median) | −0.0329 mag |
| Peers (|Δmag| ≤ 0.5, size ±30 %, ≤ 60″) | 15; Δm σ_peer 0.0187 mag |
| Target residual | 14.13 mas; peer median 8.63, p95 20.40, p99 21.45 mas; 2 / 15 peers larger |

**Decision (frozen rules):** variability `constant` (|Δm| 0.0329 ≤ 2 σ_peer = 0.0374; fade bound 3 σ_peer = 0.056 mag);
motion `no motion detected` (≤ peer p95; bound 20.40 mas / 4.0 yr = 5.10 mas/yr).

## Reading
- H-SN: an ordinary supernova in 2022 would have faded by 2026; the F444W flux is within 0.056 mag (3 σ_peer) of its
  2022 value, so a 2022 transient dominating the F444W light is disfavoured. Long-lived events and slowly varying
  AGN at the < 0.06 mag level are not tested.
- H-BD: a brown dwarf moving ≳ 5.1 mas/yr is disfavoured; slower (more distant or low v_t) ones are not excluded.
- The record's verdict (`known population: catalogued red, point-like high-redshift candidate`) stands; the
  transient open question is answered for F444W.
- The D-017 tool significance (6.71) is not evidence of motion: for faint compact targets the tool error is
  under-estimated (D-017 addendum 1); the frozen rule uses the peer null.
- **Post hoc (centring):** the frozen rule centres Δm on the bright-reference median. Centred on the peer median
  instead, |Δm| is 0.0365 (`post_hoc.dmag_minus_peer_median`) against the same 0.0374 threshold: still `constant`, but
  by under 0.001 mag, so the `constant` verdict is marginal. Against the bright-reference σ (0.0338,
  `compare.reference_dmag_sigma`) |Δm| is 0.97 σ (`post_hoc`). Either way any change is ≤ 0.06 mag.
- The contact sheet marks the record position only; a 14 mas residual (about a quarter of an F444W pixel) is not
  visible by eye, so the sheet checks presence, morphology and the neighbour, not the motion number.
- **Not tested:** F356W (no 2026 epoch), colour change, sub-5 mas/yr motion.
