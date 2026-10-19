# DR2–ICECAT long-lag deficit: result (2026-10-19)

Plan: `PLAN.md` (frozen first). Code: `scripts/e1_icetracks_selfpair.py`, frozen rule `outcome` unit-tested
in `tests/test_e1_icetracks_selfpair.py` before the run. Numbers: `selfpair.json` (derived; 1,000 scrambles
per variant, the run's main null, seeds from 8,000,000). Post hoc verification of a control-channel
observation; not part of the E1 tested family.

**Outcome of the frozen rule: "explained".**

| Variant (DR2 tracks) | wide 100 s–1 h | wide 1 h–1 d | wide 1–7 d |
|---|---|---|---|
| `all` (7,583, as the run) | 12 vs 17.82 (z −1.38) | 334 vs 415.22 (z −4.53) | 2,057 vs 2,242.99 (z −6.48) |
| `rm_shared` (7,307: 276 tracks within 10 s of an alert dropped) | 12 vs 12.33 (z −0.09) | 274 vs 292.84 (z −1.24) | 1,817 vs 1,792.77 (z +1.01) |
| `rm_random` (7,307: 276 random non-shared tracks dropped) | 12 vs 17.29 (z −1.27) | 324 vs 406.24 (z −4.63) | 1,993 vs 2,179.73 (z −6.98) |

`all` re-runs the main null with 10³ scrambles so that all three variants share seeds; its z values differ
slightly from `counts.ecsv` (10⁴ scrambles: −4.5, −6.53).

**Reading (derived).** The 276 shared tracks (270 same + 6 wide pairs at 0–10 s) are the same neutrinos as
ICECAT alerts. In the data each sits at lag ≈ 0; the null moves the two members apart by up to a few days with
RA following sidereal time, so it places one track per shared event in a long-lag "wide" cell. Dropping the
shared tracks lowers the null mean by 122.38 (1 h–1 d) and 450.22 (1–7 d) pairs and the data by 60 and 240
(`rm_shared_wide_drop`); the null-only parts, 62.38 and 210.22, are within 22 % and 11 % of the plan's
estimate of displaced pairs (≈ 80, ≈ 190; ASSUMPTION of uniform ±3 d shifts on both sides). Dropping as many
random tracks leaves the deficit. Under the frozen rule the deficit is attributed to the null's treatment of
shared events, not to the sky.

**Consequences.**
- With the shared events removed, the jit null matches this same-detector channel in all five wide lag cells
  (|z| ≤ 1.24, `wide_z.rm_shared`).
- The tested E1 channels have no shared events, so their results are unchanged.
- For long-lag wide cells under the jit null, shared events raise the null mean (here by 62–210 pairs), which
  lowers the sensitivity to an excess in those cells; the 0–10 s cells, where the shared pairs sit, are the
  control and are read with them included.
- The BAT–GBM deficit (D-074 addendum 6; it vanished when the 466 joint bursts were removed) is consistent with
  the same mechanism; not recomputed here.

**Limitations.** One null type (`jit`); 1,000 scrambles per variant; "shared" is a 10 s time match
(ASSUMPTION), so a shared event with a larger time offset would stay in. The frozen rule treats only negative
residuals as "not explained"; a large positive `rm_shared` residual would read "undetermined" (observed
maximum +1.01). The fix is not in the shared null code (`event_network`, `e1_icetracks`): a later channel with
shared events must drop them for its long-lag cells (D-074 addendum 9). Review fixes (z for zero-spread
cells, chunked seeds) were made before this re-run; the first run gave the same outcome (not stored).
