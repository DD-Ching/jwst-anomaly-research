# W3 TESS 120 s product-level injection trial — plan

Owner authorization 2026-10-26 (about 13:59 UTC), following #561. Limits: 2 h wall clock, 1 GB through
`scripts/net_ledger.py` (ledger `net_ledger.jsonl` here), 2 cores; no extension. Committed before any flux is
read. Injected events and their recoveries are simulations; no real-event search; no other paused branch.

## 0. What #561 showed (kept as is)
The frozen investment gate passed only at |M| ≲ 2 × 10⁻⁵ M☉, where the Kepler long-cadence baseline has no reach.
That is a new mass range, not a demonstrated tenfold gain in discovery chance. All of it was noise-free.

## 1. Sample (`sample.json`, metadata only)
- **Pool:** the #561 TESS 120 s TIC sample (TIC ID ending in 37) with TIC R_* and D_S. One CAOM record per star and
  sector: single-sector SPOC `timeseries`, t_exptime 120. Excluded: segment 1751, and any span > 40 d. 14,301
  records.
- **Weight:** span × D_S^1.5. This is each star-sector's rate of umbra crossings (u0 < 2.65) for any one |M| and v,
  the Level-0 term of #561's integral.
- **Draws:** 200, with replacement, with probability ∝ weight (`numpy.random.default_rng(564)`). This gives 184
  unique star-sectors (149 stars). A record drawn k times gets k independent sets of injections. Download order:
  first draw.
- No star or sector is added, removed or swapped after any flux is read.

## 2. Scenarios and events (frozen; from #561's non-zero TESS 120 s intervals)

| scenario | \|M\| (M☉) | v (km/s) | why |
|---|---|---|---|
| v220_M5.6e-6 | 5.623 × 10⁻⁶ | 220 | #561 best mass, 220 km/s |
| v220_M1e-7 | 10⁻⁷ | 220 | inside the 220 km/s interval (10⁻⁸ – 1.778 × 10⁻⁵), shorter t_E |
| v30_M3.2e-6 | 3.162 × 10⁻⁶ | 30 | #561 best mass, 30 km/s |
| v30_M1e-8 | 10⁻⁸ | 30 | inside the 30 km/s interval (10⁻⁹ – 1.778 × 10⁻⁶) |

- **Per draw and scenario, 30 events from the crossing population** (`draw_events`; seed [564, draw, scenario]):
  - x ∝ √(x(1 − x)) (the crossing rate ∝ R_E(x));
  - u0 uniform in 0–2.65;
  - t0 uniform over the light curve's full time span, gaps included.
  - t_E = R_E/v and ρ = x R_*/R_E follow from |M|, x, v, R_* and D_S.
- **Identifiability is not in the draw.** It is evaluated once, in the paired ideal, so it is not counted twice.
  Draws, x and u0 weight every event equally, as in #561's reach integral.

## 3. Paired levels for each injected event (same x, u0, t0 and exposure phase)
- **L0 ideal:** the #561 criterion (D and S, noise-free, no window) on a gap-free 120 s grid with the light
  curve's cadence phase.
- **L1a:** the noise-free model on the file's own cadence times (orbit gaps and sector edges are real), through the
  detector below with σ = 0.
- **L1b:** as L1a, after our quality screen: PDCSAP_FLUX finite and QUALITY = 0. Every set SPOC flag masks the
  cadence.
- **L2:** the model multiplied into the published PDCSAP flux at every cadence, before our quality screen, baseline
  and detection; then the full detector.
- **Exposure integration:** the mean of A(t) over each contiguous 120 s cadence (TIMEDEL), from
  `exotic_sim.finite_source_magnification` (uniform disk, blend 1).
- **PDCSAP:** crowding-corrected (consistent with blend 1) and cotrended. Injecting into the published PDCSAP
  measures a product-level efficiency only. Whether SPOC (cosmic-ray mitigation, PDC cotrending, flagging) would
  remove a real event upstream is not tested. This is not a sky-to-result efficiency.

## 4. One fixed detector (`detect`), blind to the truth
- **Baseline:** centred running median over 3.0 d (2,160 cadences) of the screened PDCSAP flux; rel = flux /
  baseline.
  - Set from the time scales before any flux: the longest frozen-scenario dark phase is ≤ 4 t_E ≈ 1.1 d (#561
    `events_at_best_mass`, 30 km/s, 95 % t_E 24,200 s), below half the window.
- **σ:** 1.4826 × MAD of rel over the sector.
- **D:** ≥ 3 contiguous cadences (spacing ≤ 1.5 cadences) with rel ≤ 10^−0.4 and 1 − rel ≥ 5σ.
- **S:** on each side of the run, within max(run duration, 4 cadences), ≥ 2 cadences with rel ≥ 1.10 and
  rel − 1 ≥ 5σ.
- **Recovery (post hoc, truth used only here):** a D-and-S detection overlapping the umbra interval
  t0 ± √(4 − u0²) t_E.
- **Background:** the same detector on each un-injected light curve. Its D-runs and D-and-S flags are counted as
  the extra flag burden only. They are not candidates, and no source is examined.

## 5. Reported, and the frozen verdict (`summarize`)
- **Per scenario (speed and |M|) and per t_E bin** (in cadences: < 2.5, 2.5–5, 5–15, 15–50, ≥ 50):
  - events, ideal-identifiable events, and those kept at L1a, L1b and L2;
  - retention = L2 / ideal over the same events;
  - a 5–95 % bootstrap over draws (2,000 resamples).
- **Losses per step:** L0 → L1a (gaps, edges, contiguity, window), L1a → L1b (quality flags), L1b → L2 (noise,
  variability, baseline, 5σ).
- **Verdict per scenario** (function `verdict`, unit-tested; "retain at least half"):
  - fewer than 30 ideal events → insufficient;
  - 5 % bound ≥ 0.5 → retains at least half;
  - 95 % bound < 0.5 → retains less than half;
  - otherwise the interval crosses one half, and this is not written as a robust pass.
- **Effect on #561:**
  - #561's TESS 120 s Level-1 reach at that scenario's mass × retention, with the same TIC sample and scaling, and
    only for these four scenarios. No efficiency is carried to other masses or other TESS stars.
  - The data volume a real search would need (from the pool's weight distribution and the file size).

## 6. Bookkeeping
Each file's sha256 and size are kept in `files.json`, and per-star-sector results in `per_ss/`. The raw files are
deleted after processing. One independent review and one targeted re-check follow.
