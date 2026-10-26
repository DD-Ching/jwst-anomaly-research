# W3 / Kepler physical mapping and value of a larger search — plan

Owner authorization 2026-10-26 (about 10:35 UTC). Limits: 2 h, 150 MB through `scripts/net_ledger.py`, 2 cores. No new
light curve or image, no injection re-run, no threshold change, no new search or photometry framework. Committed
before the mapping and population numbers are computed.

## 1. Model (as coded, `src/jwst_anomaly/exotic_sim.py`)
- **Lens:** n = 1, sign −1, a negative-mass point lens (Cramer et al. 1995; Safonova, Torres & Romero 2002).
  - θ_E = √(4G|M|/c² · D_LS / (D_L D_S)) (`einstein_radius`, Kitamura et al. 2014).
  - The umbra is β < 2 (`caustic_beta(1)` = 2).
  - β(t) = √(u0² + ((t − t0)/t_E)²) (`impact_track`).
  - ρ = θ_*/θ_E, in Einstein radii (`finite_source_magnification`).
  - The blend is the lensed share of the baseline flux.
- **Not Ellis:** the Ellis wormhole (n = 2) and general wormholes are other models and are not tested here.
- **Implied:** t_E = θ_E/μ_rel. So ρ · t_E = θ_*/μ_rel = t_*, the source-radius crossing time.
  - For a lens at D_L = x D_S moving at transverse speed v relative to the line of sight: μ_rel = v/D_L, and t_* =
    x R_*/v.
  - (ρ, t_E) cannot be chosen independently of (R_*, D_S, x, v).
  - |M| = R_E² c² / (4G D_L (1 − x)), with R_E = v t_E.

## 2. Inputs
- **Stars:** the 2,000 stars searched under one rule (#542 / #546 500 + #554 1,500). R_* and D_S with their
  errors come from the official Kepler DR25 stellar table (NExScI `keplerstellar`, `q1_q17_dr25_stellar`). These
  are model values (isochrone fits), not Gaia distances. No light curve or image is read.
- **Efficiencies:** `results/w3_kepler/s2000/combined.json`, `round_plus_546`: the 20 measured cells, 2,000 stars.

## 3. Outputs (`scripts/w3_kepler_physical.py`)
- **`mapping.json`:** for every star and cell, θ_*, θ_E, μ_rel and t_*. For each scenario speed
  v ∈ {30, 220} km/s, the lens-distance fraction x = v t_*/R_*. If x < 1: D_L, R_E and |M|. Per-cell summaries over
  stars: median and 5–95 % ranges, and the share with x < 1. A θ_* range from the catalogue errors (R high and D
  low, and the reverse) is carried along.
- **`population.json`:** the population requirement.
  - For a population whose events have one t_E (a tested value), with ρ inside the tested range and u0 < 2 (umbra
    crossings), the expected number of detected events is N = Σ_i ε · (4/π) τ_i T_i / t_E, with:
    - τ_i = (4πG/c²) |ρ_neg| D_S,i² / 6, for a uniform lens mass density |ρ_neg| along each sight line;
    - T_i, the star's window;
    - ε, the lowest measured efficiency among that t_E's four cells (pessimistic) and the highest (optimistic).
  - Reported: the |ρ_neg| needed for N = 1 and N = 3. It is compared with scale references (not population
    models): Read 2014, the local dark-matter density extrapolated from the rotation curve, 0.005–0.015 M☉ pc⁻³, and
    the local baryonic density, 0.0914 ± 0.009 M☉ pc⁻³ (Flynn et al. 2006 via Read 2014).
  - The same for +500 and +5,000 stars. Assumptions: zero result, similar efficiency, and the same
    D_S distribution.
- **Two scenarios** (only for the mapping to masses and distances; the density requirement depends on t_E, not on
  them):
  - S1 halo-like, v = 220 km/s (local circular speed; Read 2014);
  - S2 disk-like, v = 30 km/s (ASSUMPTION; a typical thin-disk speed, not sourced here).
  - Both use a uniform spatial density along the sight line. That is approximate for a disk population at the
    Kepler field's b ≈ 13°.
- **Not done:**
  - no abundance, no physical event-rate limit, no probability of discovery;
  - no efficiency between or outside the tested grid points.

## 4. Sky-to-SAP steps
Only the official documentation and existing evidence are used. The Kepler Archive Manual (KDMC-10008-006) Table
2-3, quality bits, is checked against the mask `w3_kepler.QMASK`. Unverified steps are listed; nothing is simulated
at pixel level.

## 5. Recommendation
Exactly one of A (expand), B (fill one gap) or C (pause), with its reason and evidence.
