# DR2–ICECAT long-lag deficit: shared-event displacement test (plan, frozen 2026-10-19 before the run)

Post hoc verification of the control-channel deficit in `results/e1_icetracks/RESULT.md` (wide 1 h–1 d
334 vs 415.5, z −4.5; 1–7 d 2,057 vs 2,241.7, z −6.5). Not part of the E1 tested family.

**Hypothesis H_self (ordinary, processing).** About 276 DR2 tracks are the same neutrinos as ICECAT alerts
(270 same + 6 wide pairs at 0–10 s). In the data each such pair sits at lag ≈ 0. Every null scramble shifts
the two members independently (ICECAT per-year permutation + ±3 d jitter; DR2 ±3 d uptime jitter), with RA
following sidereal time, so in the null each shared event leaves a DR2 track at a random direction ("wide")
at a lag of up to a few days next to an alert time. The null therefore over-predicts long-lag wide pairs by
about N_shared × P(lag bin) (ASSUMPTION: uniform ±3 d shifts on both sides give a triangular lag on ±6 d:
P(100 s–1 h) ≈ 0.013, P(1 h–1 d) ≈ 0.29, P(1–7 d) ≈ 0.69, i.e. ≈ 4 / 80 / 190 of 276; observed deficits
5.9 / 81.5 / 184.7). The deficit is then not a property of the sky.

**Alternatives:** a real anti-correlation (hypothesis), seasonal or run-level rate structure finer than
±3 d, alert-selection effects. The earlier diagnostic already excluded downtime, the DR2 end date and
either side's scramble alone (`control_null_variants.json`).

**Test (`scripts/e1_icetracks_selfpair.py`, main null as the run: 1,000 scrambles, seed 8,000,000):**
- V_rm: drop from DR2 every track that is in a DR2–ICECAT pair at 0–10 s (any class), recount, rerun null.
- V_ctl: drop the same number of random DR2 tracks that are not in such a pair (seed 5), recount, rerun null.

**Frozen rule (`outcome`, unit-tested before the run):** "explained" when V_rm has |z| < 2 in each of the
three wide cells 100 s–1 h, 1 h–1 d, 1–7 d **and** V_ctl keeps z ≤ −3 in the 1–7 d wide cell (removing tracks
alone does not remove the deficit); "not explained" when V_rm has z ≤ −3 in the 1 h–1 d or 1–7 d wide cell;
otherwise "undetermined". Thresholds are ASSUMPTIONs.

**What changes the plan:** "explained" closes the deficit (and is a same-detector validation of the E1
null once shared events are removed); "not explained" keeps it open as a property to vet (alert lag curve,
realtime vs archival alerts).
