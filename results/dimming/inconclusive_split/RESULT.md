# W3 dimming: every measurable held-out NEXUS `inconclusive` pair fails D-039 only for lack of an earlier exposure (2026-10-20)

Plan frozen first (`PLAN.md`). `split.py` reads the committed `persistence_run/nexus` tables (no pixels, no
network; `split.json`). Provenance: derived. Thresholds are the frozen D-039 ones (S/N 5, ratio 20, 3 h lookback;
ASSUMPTIONs); 60 s as "same exposure" is an ASSUMPTION (SW and LW starts differ by 0.06 s here).

## Frozen split (held out, e5–e7, ≤ 4 exposures, mosaic S/N ≥ 10, ≥ 2 finite fluxes)
| cause | pairs |
|---|---|
| C1 `n_prior` = 0 (no earlier exposure of the run covers the pixel; no failed prior read) | **12 / 12** |
| C2 a finite exposure below S/N 5 | 0 |
| C3 other | 0 |

- Every one of the 12 is detected at S/N 20–538 in each finite exposure; none is `suspect`.
- **11 of 12 are a single exposure**, the first dither (`_00001`) of its visit, read in its SW and LW detector
  (two "exposures" 0.06 s apart). The twelfth (d47796) has one finite flux in each of two dithers; the second is
  clean, the first has no prior.
- **Frozen decision (PLAN.md): C1 ≥ 6**, so the remaining held-out ≤ 4 deficit is a test-design effect of visit
  structure, not depth. With measurability (`undetected_stamps`) it accounts for the whole held-out ≤ 4 deficit:
  33 `on_sky` + 12 C1 = 45 / 45 measurable pairs.

## Post hoc (train epochs were used to set the D-039 rules; not held out)
- Train e0, e1, e3, e4: 13 such pairs, 12 C1 (11 single-exposure), 1 C3 (d44759 e0: one clean, one `suspect`
  detection at S/N 242–259).
- `posthoc_onsky.py` (`posthoc_onsky.json`): 67 of 176 NEXUS monitored `on_sky` pairs rest on clean detections at
  one exposure time (SW + LW); MACS0416 and Abell 2744 have none. All 32 target (`flag`) `on_sky` verdicts have
  clean detections at ≥ 2 times. So D-039 "two clean detections" does not mean two times; targets and controls
  get the same rule, so ε stays consistent with how targets are judged.

## What this changes
- No limit changes: `swap_inject` already counts these pairs as misses (ε(≤ 4) = 0.412, conservative).
- The pairs are untestable by D-039 as built: its lookback sees only the exposures given to the run (the same
  program epoch), and the first dither of a visit has none. Crediting them would need an earlier exposure on the
  same detector from the preceding visit, whichever program took it. That changes the D-039 test (frozen plan,
  unit test, controls and targets re-run in all three fields), so it is not done here.

## Next
1. Optional, method change: a cross-visit prior for first-dither exposures (preceding exposure on the same
   detector from MAST, any program), validated on the train controls first. Upper bound on gain: the 12 + 12 C1
   pairs, i.e. held-out ≤ 4 measurable 33 / 45 → at most 45 / 45 (Wilson 95 % lower 0.921); joint gain not
   estimated.
2. k = 1 aperture still waits for a new NEXUS-Center wide epoch or JADES (unchanged).
