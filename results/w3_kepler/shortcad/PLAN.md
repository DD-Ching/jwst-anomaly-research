# W3 short-cadence physical reach — plan

Owner authorization 2026-10-26 (about 11:40 UTC): one inventory of whether existing short-cadence products could
substantially improve the physical reach for the n = 1 negative-mass point lens. Limits: 90 min, 150 MB
(`scripts/net_ledger.py`, ledger `net_ledger.jsonl` here), 2 cores. Only documents, catalogues, observation
metadata and pure model calculations: no target flux, pixels, injection or search. Not a test of general
wormholes or of any communication channel. Committed before any inventory number or model number is computed.

## 1. Products and inventory (question 1)
- **Baseline:** #559's 1,971 Kepler long-cadence (LC, 1765.5 s) stars with DR25 R_* and D_S, T = 1459.49 d each
  (the window #559 used).
- **Kepler short cadence (SC, 58.85 s):** every prime-mission SC light curve in MAST CAOM (`obs_collection =
  'Kepler'`, `t_exptime` 60).
  - CAOM keeps one record per star. Per-month SC coverage is decoded from the obs_id quarter string, checked
    against the product lists of 3 stars.
  - Months are counted at ≤ 35 d each (Q0: ≤ 10 d), an upper bound. CAOM `t_min/t_max` span the whole mission
    and are not used.
  - R_* and D_S come from DR25 (`keplerstellar`), as in #559.
- **TESS 120 s and 20 s:** SPOC `timeseries` records in CAOM (`t_exptime` 120 or 20). One record per star and
  sector. T = Σ over distinct sectors of (t_max − t_min), an upper bound because the orbit gaps are included.
  - Counts and spans are server-side aggregates.
  - R_* and D_S come from TIC v8.2 (MAST TIC TAP). For 120 s, they come from a sample of the stars whose TIC ID
    ends in a fixed two digits (about 1 %), scaled to the total, with a bootstrap error. For 20 s, all stars if
    the volume allows, else the same sampling.
- **Dedupe:** within a product, one record per star and sector or month. A 20 s star's 120 s light curve is the
  same photons. Products are compared separately and are not added together.
- Stars without R_* or D_S are left out and counted. Missing quantities are reported as unknown or as bounds.

## 2. Shape after finite source and exposure integration (question 2)
- **Model:** `exotic_sim.finite_source_magnification` (n = 1, sign −1, uniform disk, blend 1), with the
  exposure-averaged flux F_k equal to the mean of A(t) over each contiguous exposure (exposure = cadence).
  Noise-free. Random exposure phase.
- **Identifiable shape** (ASSUMPTION, the Kepler pre-screen and local-shoulder numbers applied noise-free):
  - (D) ≥ 3 consecutive exposures with F ≤ 10^−0.4 (≥ 1 mag);
  - (S) on each side of that run, ≥ 2 exposures with F ≥ 1.10.
- **w(ρ, t_E/Δ):** the expected half-width in u0 (θ_E units) over which D and S hold, averaged over phase.
  Tabulated once on a (ρ, t_E/Δ) grid; w for D alone is reported too.
- **Coupling, not free parameters.** For a lens of |M| at D_L = x D_S with transverse speed v:
  - R_E = √(4G|M|/c² · x(1 − x) D_S);
  - t_E = R_E/v;
  - ρ = x R_*/R_E.
  - "Three minutes" (#559) is only a scale to check, not a prediction.

## 3. Reach (question 3)
- **Population** (the same for every product): lenses of one |M| and one transverse speed v, uniform lens mass
  density |ρ_neg| along each sight line. Speeds: S1 220 km/s (Read 2014) and S2 30 km/s (ASSUMPTION, as #559).
- **Expected events per unit density for star i:**
  N_i/|ρ_neg| = T_i (v D_S/|M|) ∫₀¹ 2 R_E(x) w(x) dx.
  - Level 0: w = 2, every umbra crossing.
  - Level 1: w from section 2, every crossing with an identifiable shape (ε = 1 inside it).
  - No Kepler LC recovery rate is used.
- **Grid:** |M| from 10⁻⁹ to 10 M☉ in 0.25 dex. Reported: Σ_i N_i per unit density per product, and the density
  for one expected event.
  - Ordinary and dark-matter densities (Read 2014) are scale references only. One expected event is not a limit
    or an exclusion.

## 4. Investment gate (question 4; a gate for spending, not a discovery standard)
Decided by function `gate()` (unit-tested):
- Passes for a product if, at one of the two speeds, a contiguous |M| interval of ≥ 3 grid points (≥ 0.5 dex) has
  Level-1 reach ≥ 10 × the baseline's Level-1 reach at the same |M| and v. A baseline of 0 counts as ≥ 10×.
- The ratio against the baseline's Level 0 is reported beside it.
- **If it passes:** one minimal test plan with its cost. Nothing is run.
- **If it fails, or the data are insufficient:** stop. No third data source is added.

## 5. Not done
Photometric noise, crowding, blending below 1, limb darkening, mass and speed distributions, clustering, lens
size versus R_E, TESS / Kepler sky-to-product steps, and pipeline flags. Each is listed as unknown.
