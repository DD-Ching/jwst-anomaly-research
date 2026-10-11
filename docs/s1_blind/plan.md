# S1 blind spots: single-pulse, faint and short bursts (TTE)

Written 2026-10-11 ~22:25 UTC, before any single-pulse score or TTE light curve was computed.

**Why 2,975 of 4,389 reduced bursts are excluded** [derived from `results/s1_twins/bursts.ecsv.gz`]:
- 846 have no resolved pulse (no peak at 5σ prominence). They are faint; TTE adds no photons.
- 2,129 have a single resolved pulse. 1,197 of these have T90 > 10 s.
- **Short bursts (T90 < 2 s):** 699 excluded against only 42 eligible. At 64 ms bins they span a few bins, and
  TTE's microsecond tags can resolve their structure. 147 of them have peak S/N > 20.

**Step A: discrimination test, bcat only, no download.** Can single-pulse bursts be searched at all once stretch,
lag and a per-band scale are free?
- **Sample:** npulse = 1 with peak S/N ≥ 10. Split in halves with seed 2. Thresholds come from the calibration
  half's real-pair background at a per-pair FAR of 10⁻⁴, using the frozen M2 (grid step 1.05, |Δk| ≤ 2). The
  consistency test is the refined per-band χ² with c = 10⁻³.
- **Injections:** the frozen deformation set (s0, E_peak, ratio) on the other half, with min_pulses = 1.
- **Unrelated pass rate:** among the 2,000 most similar calibration-half pairs.
- **Decision rule.** Single-pulse bursts are worth a full screen only if both hold:
  - (a) the expected false candidates over all single-pulse pairs (pairs above threshold × unrelated pass rate)
    are ≤ 3;
  - (b) the total recovery of deformed injections is ≥ 10 %.
- Otherwise a single smooth pulse does not carry enough shape information for this test, with bcat or with TTE
  (the 64–1,024 ms bins already resolve pulses lasting ≥ 10 s). The blind spot is then recorded as
  **structural**, not as a data-access problem.

**Step B: TTE pilot for short bursts** (only if Step A's answer, or the short-burst case, justifies it):
- Stream TTE for the bcat detectors of a small set of short bursts (about 30, peak S/N > 20), one burst at a
  time, deleting after reduction. No archive is stored, and the reason and volume are logged.
- Build two-band light curves (50–300 keV; 10–50 + 300–1,000 keV) at 2–16 ms.
- Count resolved pulses, and measure seconds and MB per burst.
- **Continue to all 699 short bursts** only if ≥ 25 % gain ≥ 2 pulses and the extrapolated volume streams in
  ≤ 2 h. Otherwise record the cost and stop.
