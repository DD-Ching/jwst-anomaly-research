# W3 dimming: NEXUS-Center-Tile-6 (o015) as an unseen ≤ 4-exposure sample for the k = 1 aperture (frozen 2026-10-20 08:30 UTC)

Handoff 1 of `results/dimming/nan_k1/RESULT.md` ("a k = 1 adoption needs new wide-epoch data"). Written before any
o015 pixel was read.

**New input (observed, MAST query 2026-10-20 08:15 UTC):** program 5105 has one public NIRCam observation not in
`configs/dimming_screen.yaml`: `jw05105-o015_t005` (target NEXUS-Center-Tile-6, s_ra 268.7782, s_dec 65.1786,
F200W t_exptime 1,245 s, t_min MJD 60878.76, 80 `_nrc` `_cal` files across bands). Its centre is 3.4′ from the
NEXUS-Center pointing, so it partly overlaps the Center mosaics. No NEXUS-Center epoch after MJD 61127.37 (o012) is
public yet. Its short exposure makes ≤ 4 covering exposures likely; it has never been read at any k.

**Hypothesis.** (H_k1) as in `nan_k1/PLAN.md`: the k = 1 fill raises the ≤ 4-exposure `on_sky` efficiency to about
0.58–0.72 on unseen data. Ordinary alternatives: the gain is specific to the Center wide epochs (o015 has other
dithers and detector positions), or the overlap holds < 10 ≤ 4 pairs.

**Sample (deviation from the monitored-pair definition, stated before the run).** o015 has no catalogue in the
screen, so the known on-sky sources are NEXUS monitored sources detected in the deep e2 (`jw05105-o014`) mosaic,
measured on the o015 F200W `_cal` exposures that cover them. Depth proxy: e2 mosaic S/N × √(1,245 / 31,137) ≥ 10,
i.e. e2 S/N ≥ 50 (ASSUMPTION; o015 S/N ≈ 0.2 × e2 S/N). The 30 d gap means a real source can fade; at this S/N
that is assumed negligible (ASSUMPTION; any `undetected` pair is looked at on a stamp before scoring is read as
efficiency only if the stamp shows the source).

**Measurement.** `run_o015.py <gaia> <out> max_nan=K` for K = 0 and K = 1 (in parallel), `persistence_check.run`
unchanged. Budget: both stopped at 09:00 UTC if unfinished (an unfinished run writes nothing and is reported as
not done). `score.py`: pairs with n_exposures ≥ 1, `on_sky` / n by class (≤ 4 / > 4), Wilson 95 %, k = 0 beside
k = 1; `persistence` counts; G0 = finite k = 0 fluxes identical at k = 1.

**Decision rules (ASSUMPTION; same values as `nan_k1/PLAN.md`).**
- R1 (adopt k = 1 for the next PR): G0 passes; n(≤ 4) ≥ 10; the k = 1 ≤ 4 Wilson interval contains 0.584 or lies
  above it, and its lower bound is > 0.412; `persistence` Wilson lower bound ≤ 2 %. Then the next PR re-runs exact
  k = 1 on every NEXUS epoch and the 37 targets, sets ε(≤ 4) = 0.584 in the swap and re-runs `inject` / `combine`.
- R2 (not adopted): any condition of R1 fails. The swap keeps 0.412.
- R3 (not tested): a run unfinished, or n(≤ 4) < 10. The swap keeps 0.412.
- > 4 class (second NEXUS held-out test of 0.961): holds if the k = 0 Wilson interval contains 0.961 or lies above.

**What would change the plan:** G0 failing; a filled exposure making a pair `persistence`; the k = 1 ≤ 4 rate
below 0.412.
