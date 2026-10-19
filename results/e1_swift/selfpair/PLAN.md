# BAT–GBM long-lag deficit: shared-event displacement test (plan, frozen 2026-10-19 before the run)

Post hoc verification of the BAT–GBM control-channel deficit in `results/e1_swift/RESULT.md` (all-separation
1 h–1 d 1,603 vs 1,784.3, z −4.9; carried by the 466 joint bursts, cause "not determined"). Same mechanism test as
`results/e1_icetracks/selfpair/PLAN.md` (D-074 addendum 9). Not part of the E1 tested family.

**Hypothesis H_self (ordinary, processing).** A joint burst (one GRB in both BAT and GBM) sits at lag ≈ 0 in the
data. The jit null shifts the two members independently, so in the null each joint burst leaves a BAT–GBM pair
at a random direction and a lag of up to days, inflating the null mean at long lags. Removing the joint bursts'
BAT rows should then remove the deficit; removing as many random non-joint BAT rows should not.

**Alternatives:** a real anti-correlation (hypothesis); orbit or Sun-constraint structure finer than the null
keeps (the Sun-frame null already leaves the deficit, `control_sun_frame_null.json`).

**Test (`scripts/e1_swift_selfpair.py`, the run's jit null, 1,000 scrambles, seed 8,000,000):**
- V_rm: drop every BAT burst within 100 s of a GBM burst (the positive-control window; any class), recount, rerun.
- V_ctl: drop the same number of random BAT bursts not in such a pair (seed 5), recount, rerun.

**Frozen rule (`outcome`, unit-tested before the run):** on the all-separation (same + wide) counts:
"explained" when V_rm has |z| < 2 in each of 100 s–1 h, 1 h–1 d, 1–7 d **and** V_ctl keeps z ≤ −3 in 1 h–1 d;
"not explained" when V_rm has z ≤ −3 in 1 h–1 d or 1–7 d; otherwise "undetermined". Thresholds are ASSUMPTIONs.

**What changes the plan:** "explained" closes the BAT–GBM deficit as the same null artefact as DR2–ICECAT;
"not explained" keeps it open (orbit-phase structure).
