# BAT–GBM long-lag deficit: result (2026-10-19)

Plan: `PLAN.md` (frozen first). Code: `scripts/e1_swift_selfpair.py`, frozen rule `outcome` unit-tested in
`tests/test_e1_swift_selfpair.py` before the run. Numbers: `selfpair.json` (derived; 1,000 jit scrambles per
variant, seeds from 8,000,000, both catalogues scrambled as in the run). Post hoc verification of a
control-channel observation; not part of the E1 tested family. All counts are all-separation (same + wide).

**Outcome of the frozen rule: "undetermined"** (the 100 s–1 h cell of `rm_shared` is z +4.05, outside |z| < 2).
**Post hoc outcome (deviation, the same rule on `rm_joint_1h` / `rm_random_1h`): "explained".**

| Variant (BAT bursts) | 100 s–1 h | 1 h–1 d | 1–7 d |
|---|---|---|---|
| `all` (1,635, as the run) | 94 vs 74.56 (z +2.22) | 1,603 vs 1,785.18 (z −4.91) | 10,321 vs 10,629.12 (z −5.26) |
| `rm_shared` (1,159: 476 bursts within 100 s of a GBM burst dropped) | 72 vs 43.95 (z +4.05) | 1,031 vs 1,038.74 (z −0.27) | 6,445 vs 6,476.67 (z −0.65) |
| `rm_random` (1,159: 476 random others dropped) | 61 vs 56.11 (z +0.67) | 1,177 vs 1,356.28 (z −5.61) | 7,729 vs 8,010.63 (z −5.30) |
| post hoc `rm_joint_1h` (1,135: also same-direction GBM within 1 h) | 44 vs 42.54 (z +0.22) | 999 vs 1,003.36 (z −0.16) | 6,270 vs 6,283.75 (z −0.29) |
| post hoc `rm_random_1h` (1,135: 500 random others dropped) | 80 vs 56.09 (z +3.18) | 1,163 vs 1,353.92 (z −5.81) | 7,626 vs 7,936.96 (z −5.74) |

**Reading (derived).**
- Removing the BAT bursts that have a GBM burst within 100 s removes the long-lag deficit (1 h–1 d z −0.27,
  1–7 d z −0.65); removing as many random bursts keeps it (−5.61, −5.30). The null mean falls by 746.44 and
  4,152.45 pairs, the data by 572 and 3,876; the null-only parts are 174.44 and 276.45 pairs
  (`rm_shared_drop`). As for DR2–ICECAT (D-074 addendum 9), the jit null moves each joint burst's two members
  apart and so over-predicts long-lag pairs.
- The frozen rule fails only on the 100 s–1 h cell. Of those 72 `rm_shared` pairs, 25 are same-direction
  (session listing; the 24 BAT bursts are stored as `joint_1h_only_bat_names`), most at separation ≤ 0.1° with
  σ_comb 0.1° (the GBM row carries the Swift position, as in `results/e1_swift/RESULT.md`), e.g. GRB 221009A (GBM trigger 3,197 s before BAT's) and GRB 130925A (two GBM triggers).
  These are joint GRBs whose BAT and GBM triggers are 100 s–1 h apart, which the frozen 100 s window misses
  (the plan's own "Limitations" case). Dropping them too (post hoc, these 24 bursts) gives |z| ≤ 0.29 in all
  three cells while the random control keeps the deficit (1 h–1 d z −5.81).

**Consequences.** The BAT–GBM deficit is attributed, post hoc, to the same null artefact as DR2–ICECAT;
D-074 addendum 10. A joint-event rule needs a position match over ≥ 1 h, not a 100 s time window alone.
No tested E1 channel includes GBM–BAT, so no E1 conclusion changes.

**Limitations.** One null type (`jit`); 1,000 scrambles per variant; the 1 h joint window and its use are
post hoc (chosen after the frozen rule failed), so the "explained" reading is exploratory. The `all` row
re-runs the null with 10³ scrambles; its z differ slightly from `counts.ecsv` (10⁴).
