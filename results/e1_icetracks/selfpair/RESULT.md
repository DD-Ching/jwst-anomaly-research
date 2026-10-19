# DR2–ICECAT long-lag deficit: result (2026-10-19)

Plan: `PLAN.md` (frozen first). Code: `scripts/e1_icetracks_selfpair.py`, frozen rule `outcome` unit-tested
in `tests/test_e1_icetracks_selfpair.py` before the run. Numbers: `selfpair.json` (derived; 1,000 scrambles
per variant, the run's main null, seed 8,000,000). Post hoc verification of a control-channel observation;
not part of the E1 tested family.

**Outcome of the frozen rule: "explained".**

| Variant (DR2 tracks) | wide 100 s–1 h | wide 1 h–1 d | wide 1–7 d |
|---|---|---|---|
| `all` (7,583, as the run) | 12 vs 18.02 (z −1.47) | 334 vs 414.57 (z −4.41) | 2,057 vs 2,242.85 (z −6.65) |
| `rm_shared` (7,307: 276 tracks within 10 s of an alert dropped) | 12 vs 12.23 (z −0.07) | 274 vs 292.90 (z −1.23) | 1,817 vs 1,793.63 (z +0.97) |
| `rm_random` (7,307: 276 random non-shared tracks dropped) | 12 vs 17.21 (z −1.27) | 324 vs 406.62 (z −4.66) | 1,993 vs 2,180.58 (z −6.76) |

**Reading (derived).** The 276 shared tracks (270 same + 6 wide pairs at 0–10 s) are the same neutrinos as
ICECAT alerts. In the data each sits at lag ≈ 0; the null moves the two members apart by up to a few days with
RA following sidereal time, so it places one track per shared event in a long-lag "wide" cell that the data do
not have. Dropping the shared tracks lowers the null by 122 (1 h–1 d) and 449 (1–7 d) but the data by only 60
and 240; the differences (62, 209) match the plan's estimate of the displaced pairs (≈ 80, ≈ 190; ASSUMPTION
of uniform ±3 d shifts). Dropping as many random tracks lowers data and null alike, and the deficit stays.
So the deficit is a property of the null on a channel whose catalogues share events, not of the sky.

**Consequences.**
- With the shared events removed, the E1 null reproduces a same-detector channel at every lag (|z| ≤ 1.23 in
  the five wide cells): the validation that `RESULT.md`'s Limitations said was missing.
- The tested E1 channels have no shared events, so their results are unchanged. The bias lowers the null's
  long-lag mean wherever shared events remain, so it can hide an excess but never create one.
- The BAT–GBM deficit (D-074 addendum 6; it vanished when the 466 joint bursts were removed) fits the same
  mechanism; not recomputed here.

**Limitations.** One null type (`jit`); 1,000 scrambles per variant; "shared" is a 10 s time match
(ASSUMPTION), so a shared event with a larger time offset would stay in.
