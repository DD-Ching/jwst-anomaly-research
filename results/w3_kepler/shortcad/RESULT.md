# W3 short-cadence physical reach: all three products pass the frozen investment gate, but only at masses the long-cadence baseline cannot resolve (|M| ≲ 2 × 10⁻⁵ M☉ at 220 km/s); none reaches 10× where the baseline has reach; one noise-included injection test is proposed, not run

- **Plan:** `PLAN.md` was committed at 11:47 UTC (a985360), before any inventory or model number.
- **Downloads:** 19,938,373 B of the 150 MB cap (`net_ledger.jsonl.gz`).
  - All of it is metadata: MAST CAOM records, TIC and DR25 rows, and 3 product lists.
  - A first inventory run was stopped when the TIC lookups of every 20 s star proved too slow. Its bytes are
    counted.
- **No:** target flux, pixels, injection or search.
- **Labels:**
  - Catalogue R_* and D_S are model values (DR25 isochrone fits; TIC v8.2).
  - Coverage spans are upper bounds.
  - Light curves are noise-free model predictions.
  - The reach numbers are derived scenario requirements. They are not abundances, limits or exclusions.
- **Deviations:**
  - TESS 20 s was sampled (TIC ID ending in 7, about 10 %), not taken whole: the plan's fallback.
  - S2 (30 km/s) is an unsourced assumption, as in #559.
  - Post hoc, not in the plan:
    - the D_L ≤ 1 kpc sensitivity (section 3);
    - the TESS light-curve file sizes for the cost estimate (`lc_file_sizes.json`, product metadata);
    - the per-product distributions (`summary.json`).
  - The plan did not fix the bootstrap size; 100 resamples were used.
- **Review fixes** (independent review, `/code-review`):
  - Multi-sector CAOM records (an obs_id with two sector tokens, such as s0001–s0096 spanning 2,579 d) had entered
    the TESS per-star spans and inflated TESS 120 s coverage by about 2.9×. They are excluded now (`tess_fix`).
  - Below the ρ grid (ρ < 0.001), the lookup returned w = 0 instead of the point-source limit. It now uses the first
    ρ row. This changed only the high-|M| end; no best mass changed.
  - The first version's statement "TESS 120 s most robust" is withdrawn.

## 1. How many usable observations exist? (`inventory.json`, `summary.json`)

| product | cadence | stars | records | stars with R_* and D_S used | coverage per star, 5 / 50 / 95 % (upper bound) | Σ T (star-years, upper bound) | D_S median (5–95 %) pc | R_* median (5–95 %) R☉ |
|---|---|---|---|---|---|---|---|---|
| Kepler LC baseline (#559) | 1,765.5 s | 1,971 of the 2,000 searched | — | 1,971 | 1,459 d (window) | 7,876 | 587.9 (265.6–918.9) | 1.08 (0.7265–1.64) |
| Kepler SC (prime) | 58.85 s | 5,337 | 24,381 SC months | 5,069 (DR25) | 1 / 1 / 24.2 months; 35 / 35 / 875 d | 2,210 | 467.9 (127.1–1,777) | 1.739 (0.7342–7.773) |
| TESS 20 s (sample 6,454) | 20 s | 63,967 | 157,309 | 4,049 of the sample (63 %) | 1 / 2 / 7 sectors; 25.21 / 51.67 / 199.0 d | 7,889 (scaled) | 207.1 (24.36–1,862) | 1.723 (0.18–60.39) |
| TESS 120 s (sample 5,812) | 120 s | 592,890 | 1,737,708 | 4,930 of the sample (85 %) | 1 / 2 / 9 sectors; 24.58 / 28.8 / 229.9 d | 106,600 (scaled) | 196.6 (50.4–1,023) | 1.146 (0.3027–17.5) |

- **Sources:** MAST CAOM (`dbo.ObsPointing`, SPOC `timeseries`), TIC v8.2 (MAST TIC TAP), Kepler DR25 stellar table
  (NExScI).
  - TESS sectors are 1–107 (20 s from 27). One 2026 segment of about 7 d is labelled 1751 in CAOM; it is included.
- **Kepler SC coverage:** CAOM keeps one record per star. The obs_id quarter string gives SC months per quarter.
  - The reading matched the number of `_slc.fits` files for 3 stars (49 / 49, 1 / 1, 1 / 1).
  - Each month is counted as ≤ 35 d (Q0 ≤ 10 d), an upper bound. CAOM's `t_min/t_max` span the whole mission and
    are not used.
- **TESS coverage:** T = Σ over distinct sectors of (t_max − t_min). That includes the orbit gaps, so it is an upper
  bound.
  - Multi-sector records are excluded: 1,836,854 − 1,737,708 = 99,146 at 120 s (difference not stored), none at 20 s.
    After that, no star has two records for one sector.
  - The longest remaining single-sector span is 57.4 d (`max_record_span_d`). Sectors 97 and 98 are 55–57 d long in
    CAOM.
  - Segment 1751 lies inside sector 99 (MJD 61055–61062, within 61045–61073; one CAOM query, not stored). A star in both has up to 7.23 d
    counted twice. At most about 0.04 % (120 s) or 0.15 % (20 s) of Σ T (not stored).
- **Not added together:** a 20 s star's 120 s light curve is the same photons. Products are compared separately.
  The mission length is never used as a star's time.
- **Stars without values contribute 0 (conservative):**
  - 268 Kepler SC stars lack DR25 R_* or D_S, and 13 of the 5,069 used have D = 0.00 in DR25;
  - 37 % of the 20 s sample and 15 % of the 120 s sample lack TIC R_* or D_S; what kind of stars they are was not
    counted.
- **Kepler SC row:** months per star come from all 5,337 rows, days and Σ T from the 5,069 used rows.
- **TESS scaling:** totals are the sample sums × (all stars / sampled stars): 9.911 for 20 s and 102.0 for 120 s.
  The 5–95 % bootstrap over sampled stars is in `reach.json`.

## 2. Can short cadence see the shape? (`shape.json`)
- **Coupled, not free:** a lens of |M| at D_L = x D_S with transverse speed v has:
  - R_E = √(4G|M|/c² · x(1 − x) D_S);
  - t_E = R_E/v;
  - ρ = x R_*/R_E.
  - So for given R_*, D_S and v, choosing |M| and x fixes both t_E and ρ.
- **Light curve:** `exotic_sim.finite_source_magnification` (uniform disk, blend 1), averaged over each
  contiguous exposure. Random phase (16 phases), u0 from 0 to 2.65. Noise-free.
- **Identifiable shape** (ASSUMPTION; the Kepler pre-screen and local-shoulder numbers, noise-free):
  - (D) ≥ 3 consecutive exposures at ≥ 1 mag below baseline;
  - (S) ≥ 2 exposures at ≥ 10 % above baseline on each side.
- **w = phase-averaged u0 half-width with D and S (θ_E units; 2 would mean every umbra crossing):**

| ρ | t_E = 1 cadence | 1.58 | 2.51 | 3.98 | 10 | ≥ 100 |
|---|---|---|---|---|---|---|
| 0.001–0.107 | 0 | 0.125–0.131 | 1.32–1.37 | 1.88–1.95 | 1.93–1.98 | 1.93–1.98 |
| 0.405 | 0 | 0.0375 | 1.067 | 1.703 | 1.725 | 1.725 |
| 0.790 | 0 | 0 | 0.178 | 1.491 | 1.525 | 1.525 |
| 1.103 | 0 | 0 | 0 | 0.0094 | 1.375 | 1.375 |
| 1.539 | 0 | 0 | 0 | 0 | 0 | 0 |

- **Reading:**
  - The shape needs t_E ≳ 2.5 exposures (about 4 for ρ ≈ 0.8, 10 for ρ ≈ 1.1). It needs ρ ≲ 1.1: at
    ρ = 1.54, finite-source smoothing keeps the shoulders below 10 %.
  - The dark run alone (D) survives to larger ρ; `w_D` is in `shape.json`.
  - Between grid points w is interpolated in (log ρ, log t_E); the cut between ρ 1.10 and 1.54 is not resolved
    finer.
- **Thresholds per cadence (t_E ≳ 2.5 exposures):** about 74 min for Kepler LC, 2.5 min for Kepler SC, 5 min for
  TESS 120 s and 50 s for TESS 20 s (not stored).
  - #559's "three minutes" (the t_E that would give one expected event at 0.015 M☉ pc⁻³ under ε = 1) is
    therefore within reach of the 20 s and 58.85 s cadences. Whether real lenses have such t_E is what section 3
    computes; it is not assumed.

## 3. How much does the physical reach improve? (`reach.json`; post hoc `sensitivity_dl1kpc.json`)
- **Population** (the same for every product): lenses of one |M| and one transverse speed v, with a uniform lens
  mass density |ρ_neg| along each sight line.
  - Expected events per unit density: Σ_i T_i (v D_S/|M|) ∫ 2 R_E(x) w(x) dx.
  - Level 0 uses w = 2 (every umbra crossing). Level 1 uses w from section 2 (every crossing with an identifiable
    shape, ε = 1 inside it).
  - No Kepler LC recovery rate is used for any product, the baseline included.
  - |M| runs from 10⁻⁹ to 10 M☉ in 0.25 dex.
- **Best Level-1 reach per product** (density for one expected event, M☉ pc⁻³; TESS 5–95 % bootstrap over sampled
  stars in brackets):

| product | 220 km/s: best \|M\| (M☉) | density for 1 event | 30 km/s: best \|M\| | density for 1 event | same, lenses only at D_L ≤ 1 kpc (post hoc), 220 / 30 km/s |
|---|---|---|---|---|---|
| Kepler LC baseline | 5.623 × 10⁻⁵ | 0.6684 | 5.623 × 10⁻⁶ | 1.958 | 0.6959 / 2.026 |
| Kepler SC | 1.778 × 10⁻⁶ | 0.1511 | 1.778 × 10⁻⁶ | 1.107 | 0.8977 / 6.365 |
| TESS 20 s | 3.162 × 10⁻⁸ | 0.05005 (0.03492–0.07133) | 1.778 × 10⁻⁸ | 0.3579 (0.2628–0.5485) | 0.1177 / 0.7616 |
| TESS 120 s | 5.623 × 10⁻⁶ | 0.08734 (0.07605–0.1) | 3.162 × 10⁻⁶ | 0.6146 (0.5456–0.7067) | 0.1155 / 0.8299 |

- Densities are 1/N (`densities.json`, from `reach.json` and the post hoc `sensitivity_dl1kpc.json`); the
  brackets are 1/(bootstrap N).
- **Against the baseline's best density (`baseline_best_over_this_best` in `densities.json`):**
  - TESS 20 s is 13.35× lower, TESS 120 s 7.653× and Kepler SC 4.424× (220 km/s).
  - With lenses only within 1 kpc (post hoc): 5.912×, 6.025× and 0.7752×.
  - These compare different masses, so they are not the same-population ratio the gate uses.
- **Events at the best mass (220 km/s; rate-weighted 5 / 50 / 95 %, `events_at_best_mass`):**
  - baseline: t_E 4,100 / 5,670 / 7,140 s; ρ 0.128 / 0.288 / 0.581; D_L 105 / 312 / 660 pc;
  - Kepler SC: t_E 968 / 4,210 / 5,420 s; ρ 0.247 / 0.705 / 1.25; D_L 174 / 4,200 / 10,700 pc;
  - TESS 20 s: t_E 133 / 298 / 478 s; ρ 0.221 / 0.649 / 1.22; D_L 165 / 1,220 / 4,830 pc;
  - TESS 120 s: t_E 620 / 1,990 / 4,470 s; ρ 0.0829 / 0.581 / 1.24; D_L 21.0 / 311 / 2,300 pc.
- **Why the post hoc 1 kpc column:** the rate grows as D_S^1.5. Under a uniform density along the sight line, much
  of the short-cadence reach comes from distant stars and lenses at kpc distances (Kepler SC median D_L 4.2 kpc).
  - A density that is uniform for several kpc is a strong assumption, above all for a disk-like population.
  - Restricting lenses to D_L ≤ 1 kpc changes the baseline by 4 % but removes Kepler SC's gain in best density.
- **Same |M|, same v (ratio of Level-1 reach to the baseline's, `ratio_vs_baseline_N1`):**
  - Where the baseline has reach (|M| ≥ 5.623 × 10⁻⁵ M☉ at 220 km/s, ≥ 5.623 × 10⁻⁶ at 30 km/s), no product
    reaches 10×:
    - TESS 120 s 4.69–8.748 (220 km/s) and 3.083–8.748 (30 km/s); with lenses only within 1 kpc (post hoc), 3.696–7.733 from
      3.162 × 10⁻⁵ M☉ up (220 km/s);
    - Kepler SC 1.235–1.606 and 1.158–1.47;
    - TESS 20 s 0.6726–1.11 and 0.6244–1.11.
  - At lower |M| the LC exposures are too long for the shape, the baseline's reach falls to 0, and every product's
    ratio grows without bound. That is where the gate passes (section 4).
- **Scale references only (Read 2014; not negative-mass densities):** local baryons 0.0914 M☉ pc⁻³; extrapolated
  dark matter 0.005–0.015.
  - TESS 20 s's best case, 0.05005, and TESS 120 s's, 0.08734, lie between them. Within 1 kpc, both (0.118 and 0.116)
    are above the baryon reference.
  - One expected event is a scenario requirement. It is not a limit or an exclusion, and no light curve was read.
- **Optimistic by construction:**
  - ε = 1 inside the noise-free shape criterion;
  - coverage spans include gaps;
  - one mass and one speed;
  - uniform density.
  - Each of these can only lower the real reach. The exception is the stars left out for missing values, which
    would raise it.

## 4. Is a next measurement worth it? Gate: passes, only at low masses (`reach.json` → `gate`)
- **Gate result:** by the frozen rule (≥ 3 contiguous grid points with Level-1 reach ≥ 10 × the baseline at the
  same |M| and v), all three products pass at both speeds:

| product | 220 km/s interval (M☉) | 30 km/s interval (M☉) |
|---|---|---|
| Kepler SC | 1.778 × 10⁻⁸ – 1.778 × 10⁻⁵ | 1.778 × 10⁻⁹ – 10⁻⁶ |
| TESS 20 s | 10⁻⁹ (grid edge) – 10⁻⁵ | 10⁻⁹ (grid edge) – 10⁻⁶ |
| TESS 120 s | 10⁻⁸ – 1.778 × 10⁻⁵ | 10⁻⁹ (grid edge) – 1.778 × 10⁻⁶ |

- **What the pass means:**
  - Every interval lies where the baseline resolves the shape rarely or never: negative masses from 10⁻⁹ to
    about 2 × 10⁻⁵ M☉ (up to about 6 Earth masses, not stored), with t_E of minutes to about an hour. Short cadence
    opens these masses; it does not improve by 10× the masses the baseline already reaches.
- **Which product:**
  - TESS 20 s gives the lowest best density (0.050 at 220 km/s). Within 1 kpc it matches TESS 120 s (0.118 against
    0.116).
  - TESS 120 s has the largest same-mass gain where the baseline has reach (up to 8.748×, below the gate factor).
    Its files are about 1.8 MB per star and sector, against about 11 MB at 20 s (`lc_file_sizes.json`).
  - Kepler SC adds little: about 1.2–1.6× at the baseline's masses, and nothing within 1 kpc.
- **What the gate does not show:**
  - The shape criterion is noise-free. Real photometry (noise, momentum dumps, scattered light, SPOC quality flags)
    may remove the 10 % shoulders or the 3-exposure dark run.
  - This is the main unknown and the reason for the test below.

### Minimal measurement plan (not run; needs a new authorization)
- **Question:** with real noise and flags, what fraction of the noise-free Level-1 events stays identifiable in
  TESS 120 s SPOC light curves?
  - 120 s is chosen on cost: 20 s reaches lower masses at about 6× the bytes per star and sector. The same plan can
    run on 20 s at that cost.
- **Sample:** 200 star-sectors drawn at random, with a frozen seed, from the 4,930 sampled 120 s stars that have R_*
  and D_S. There is no real-event search; the flux is used only as a carrier for injections.
- **Injections:** 20 per star-sector, multiplied into PDCSAP flux. Each is drawn from the best-mass population
  (5.623 × 10⁻⁶ M☉ at 220 km/s, 3.162 × 10⁻⁶ at 30 km/s):
  - x from that star's rate weight 2 R_E(x) w(x);
  - then t_E and ρ from x;
  - u0 uniform in 0–2.65 and a random phase.
- **Criterion:**
  - D: ≥ 3 consecutive exposures that are both ≥ 1 mag below baseline and ≥ 5σ below the local scatter
    (`LOCAL_NSIG`);
  - S: ≥ 2 exposures per side that are both ≥ 10 % above baseline and ≥ 5σ above it;
  - quality-flagged exposures are masked.
- **Gate for any later search (frozen before any flux):** the injection-recovered share is at least half of the
  noise-free share for the same draws. Otherwise stop.
- **Cost:**
  - 200 × 1.85 MB ≈ 370 MB, streamed and deleted. That is above this round's 150 MB cap (not stored).
  - Estimates: about 1–2 h wall clock on 2 cores, mostly downloads; a new 120 s adapter of about 100 lines around
    the existing W3 code.
- **Not included:** a real-event search, Kepler SC, or any other source.

## 5. Not done, and unknown
- Photometric noise and systematics (section 4 plan), crowding and blends below 1, limb darkening.
- Mass and speed distributions, clustering, and a Galactic density model; uniform density is the only population
  model.
- The lens's physical size against R_E.
- Whether a negative-mass object of 10⁻⁹–10⁻⁵ M☉ is compact enough to act as a point lens is not addressed.
- The steps from sky to SPOC / Kepler products (onboard cosmic-ray mitigation of 120 s data, flags).
- Stars without catalogue values.
- **Three separate questions:**
  - (i) The model can produce this shape at these cadences: yes, noise-free, section 2.
  - (ii) Data support the model: not addressed; no flux was read.
  - (iii) A controllable signal or channel: not addressed.
