# BAT–GBM long-lag deficit: result (2026-10-19)

Plan: `PLAN.md` (frozen first). Code: `scripts/e1_swift_selfpair.py`, frozen rule `outcome` unit-tested in
`tests/test_e1_swift_selfpair.py` before the run. Numbers: `selfpair.json` (derived; 1,000 jit scrambles per
variant, seeds from 8,000,000; GBM is scrambled first from each seed, so its scrambles are identical across
variants). Post hoc verification of a control-channel observation; not part of the E1 tested family. All counts
are all-separation (same + wide).

**Outcome of the frozen rule: "undetermined"** (`rm_shared` 100 s–1 h z +4.2, outside |z| < 2).
Post hoc (deviation; the same rule applied to `rm_joint_1h` / `rm_random_1h`): "explained", but see the
circularity note below.

| Variant (BAT bursts) | 100 s–1 h | 1 h–1 d | 1–7 d |
|---|---|---|---|
| `all` (1,635, as the run) | 94 vs 74.32 (z +2.31) | 1,603 vs 1,783.77 (z −4.78) | 10,321 vs 10,631.22 (z −5.04) |
| `rm_shared` (1,159: 476 bursts within 100 s of a GBM burst dropped) | 72 vs 43.33 (z +4.2) | 1,031 vs 1,040.29 (z −0.33) | 6,445 vs 6,476.41 (z −0.66) |
| `rm_random` (1,159: 476 random others dropped) | 61 vs 55.78 (z +0.69) | 1,177 vs 1,358.2 (z −5.83) | 7,729 vs 8,007.44 (z −5.31) |
| post hoc `rm_joint_1h` (1,135: also a same-direction GBM burst within 1 h) | 44 vs 42.84 (z +0.17) | 999 vs 1,003.48 (z −0.16) | 6,270 vs 6,282.03 (z −0.26) |
| post hoc `rm_random_1h` (1,135: 500 random others dropped) | 80 vs 56.17 (z +3.21) | 1,163 vs 1,353.94 (z −5.81) | 7,626 vs 7,937.7 (z −6.01) |

**Reading (derived).**
- The 476 shared bursts are 466 with a same-direction GBM pair within 100 s (the joint bursts of
  `results/e1_swift/RESULT.md`) plus 10 with only a wide-class pair (`n_shared_bat_with_same_direction_pair`).
- Long lags: dropping the shared bursts removes the deficit (1 h–1 d z −0.33, 1–7 d z −0.66); dropping as
  many random bursts keeps it (−5.83, −5.31). The null mean falls by 743.48 and 4,154.81 pairs, the data by
  572 and 3,876; the null-only parts are 171.48 and 278.81 pairs (`rm_shared_drop`). This is what the
  displacement mechanism of D-074 addendum 9 predicts, and these cells are not circular: the dropped bursts
  were chosen by a 0–100 s match, not by their long-lag pairs.
- 100 s–1 h: the frozen rule fails here. 25 of the 72 `rm_shared` pairs are same-direction pairs from 24 BAT
  bursts (`joint_1h_only_pairs`); 23 of them lie within 0.1° (σ_comb ≈ 0.05°, the GBM row carrying the Swift
  position), e.g. GRB 221009A (GBM trigger 3,197.0 s before BAT's) and GRB 130925A (two GBM triggers, −900.7 s
  and −117.3 s). Such positional identity reads as the same GRB with BAT and GBM triggers 100 s–1 h apart,
  which a 100 s window does not catch (a gap in the frozen plan, found after the run).
- **Circularity:** `rm_joint_1h` drops exactly the bursts with same-direction pairs at ≤ 1 h, so its 100 s–1 h
  z near 0 is close to built in and is not evidence. Only its long-lag cells (−0.16, −0.26) add anything, and
  they repeat `rm_shared`. `rm_random_1h` keeps the 24 bursts (they are excluded from the random draw), so its
  100 s–1 h z +3.21 is the same 25 pairs.

**Consequences.** The long-lag BAT–GBM deficit behaves as the shared-event displacement predicts; the frozen
rule stays "undetermined" because of the 100 s–1 h cell, whose excess consists of position-identical joint
GRBs (identity by position, not tested by an independent rule). No tested E1 channel includes BAT–GBM, so no
E1 conclusion changes. A future joint-event rule should be pre-registered with a position match over a
longer window (hypothesis for a frozen plan, D-074 addendum 10).

**Limitations.** One null type (`jit`); 1,000 scrambles per variant; one random-control seed per variant; the
1 h window and its variants are post hoc. The `all` row re-runs the null with 10³ scrambles; its z differ
slightly from `counts.ecsv` (10⁴).
