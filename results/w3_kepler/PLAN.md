# W3 / Kepler dense-cadence pilot: frozen plan

Owner authorization 2026-10-24. Start 07:39 UTC. Limits: 3 h (to 10:39 UTC), 2 GB through
`scripts/net_ledger.py` (`net_ledger.jsonl` here), 4 cores. This plan, `scripts/w3_kepler.py` and
`tests/test_w3_kepler.py` are committed **before any target flux is read**.

**Hypothesis tested:** the existing finite-source negative-mass point lens (W3; n = 1, ε < 0;
`exotic_sim`). Not all wormhole models, and not a faster-than-light communication test.

**Unchanged:** paused lines (ZTF / ATLAS, J2244, X4 and the rest).

## 0. Duplicate check (07:39–07:50 UTC)
- **`main` at 180c0ac:** no Kepler light-curve work.
  - "Kepler" appears only as Keplerian orbits (X1, `wh_farside.py`) and a PDC literature citation
    (`docs/hypotheses/round-3/B-on-R3A.md`).
- **Open PRs and `claude/*` branches:** none on Kepler. W3 branches cover OGLE, MOA, Gaia, ZTF,
  ATLAS and the short-t_E study only.
- **Earlier TESS check:** the coverage check of the 40 ZTF stars was not a Kepler search, and it is
  not used here.

## 1. Data
- **Product:** MAST Kepler prime-mission DR25 long-cadence light curves (`kplr<id>-<ts>_llc.fits`,
  `https://archive.stsci.edu/missions/kepler/lightcurves/`).
  - Interface check on one non-target star, KIC 2578941 (rank 100, excluded from every group).
    Its Q5 header reads `DATA_REL` 25 and SOC 9.3.23; the file is 492 kB with 4,634 cadences.
- **Quarters:** Q1–Q17 only; Q0 is not used.
- **Flux:** `SAP_FLUX` with `SAP_FLUX_ERR` is primary. `PDCSAP_FLUX` is a diagnostic on flagged runs
  only (`pdc_in_run`).
  - SAP and PDCSAP come from the same observation, so they are never two independent confirmations.
- **Storage:** only the needed columns are cached under the data root (gitignored), and the cache
  is deleted at the end. Per file, size and sha256 go to `files_<group>.json`.
- **Expected volume:** 100 × 18 files × ~0.49 MB ≈ 0.9 GB.

## 2. Sample (metadata only; `sample.json`)
- **Pool:** NExScI `keplerstellar`, delivery `q1_q17_dr25_stellar`, with logg ≥ 4.0,
  4500 ≤ Teff ≤ 6500 K and 12.0 ≤ Kp ≤ 14.0.
  - The query result's sha256 is pinned in `sample.json`; it has 24,474 rows.
  - Kept: `st_quarters` = 17 ones, i.e. observed in Q1–Q17. That leaves 18,570 stars.
- **Selection:** sort by sha256(`w3-kepler-2026-10-24:<kepid>`).
  - Ranks 0–19 are **dev**, 20–39 **val** and 40–99 **search**.
  - A star belongs to exactly one group. Nothing is swapped after any flux is read.
- **Read order:** dev first; val only after the dev phase; search only if `gate.json` passes.
  `source_usage.jsonl` logs each read.

## 3. Processing (SAP product onward; fixed)
- **Quality:** a cadence is excluded if any of these `SAP_QUALITY` bits is set: 1, 2, 4, 8, 16, 32,
  64, 128, 512, 1024, 16384, 65536 (`QMASK`).
  - Also excluded: non-finite time or flux, or an error ≤ 0.
  - Straylight (2048, 4096) and planet-search exclude (8192) cadences are kept.
  - An excluded cadence is **absent**: never zero flux, never dim, never an edge.
- **Quarter join:** relative flux = SAP / (median SAP of that quarter's good cadences). Same for σ.
  - A quarter with < 40 good cadences is left out.
  - There is no cross-quarter interpolation and no gap filling.
- **Series** (thin interface to `w3_ztf`):
  - dm = −2.5 log10 rel; rel ≤ 0 is a non-detection, which counts as dim.
  - `lim` is set so that w3_ztf's error formula gives the Kepler σ of dm.
- **Run finding** (`flagged_runs_kepler`):
  - ≥ 3 consecutive dim cadences (dm ≥ 1.0, `w3_ztf.DIM_MAG`);
  - ≥ 20 non-dim cadences on each side.
  - The ground-based "≥ 3 nights" condition is dropped (Kepler is continuous).
- **Single band:** there is no second band. The achromatic test is **not applicable** and is
  recorded as "not tested (single band)", never as passed.

## 4. Injection (before any baseline, join or screen of this round)
- **Order:** model F(t) × raw `SAP_FLUX` at every cadence, then the quality mask, quarter
  normalisation and screen. No baseline is estimated before the injection.
- **Exposure average:** F is averaged over each 29.42-min cadence (21 sub-samples, uniform).
- **W3 model:** `exotic_sim.light_curve` (n = 1, sign −1, finite source), tabulated in u and
  interpolated.
  - The test checks this against direct evaluation (max difference < 2 × 10⁻³).
  - F = blend × A + 1 − blend.
- **Noise:** SAP_inj = SAP × F + N(0, err² (max(F, 0.09) − F²)) where positive;
  err_inj = err √max(F, 0.09). This is the sky-limited floor of `w3_microlensing.inject_w3`
  (ASSUMPTION).
- **W3 grid** (20 cells), fixed now from span, gaps and integration:
  - **t_E ∈ {0.1, 0.3, 1, 3, 10} d:**
    - The umbra lasts 2 t_E √(4 − u0²) ≤ 4 t_E.
    - At t_E = 10 d that is ≤ 40 d, below half a ~93-d quarter, so the per-quarter median survives.
    - t_E ≥ 30 d is **excluded in advance**: its umbra of up to 120–400 d spans whole quarters,
      which the per-quarter join cannot keep.
    - t_E = 0.1 d gives an umbra of ≥ 3 cadences for u0 ≤ 1.96.
  - **ρ ∈ {0.01, 0.1}.**
  - **blend ∈ {1.0, 0.7}**, as in ZTF.
  - **Draws:** u0 ~ U[0, 2), t0 ~ U over the star's good-cadence range (gaps included). A carrier
    star is drawn uniformly from the phase's group.
  - **n = 100** per cell and phase.
  - **Excluded from every cell:** t_E 30–300 d, blend < 0.7 and ρ > 0.1.
- **Ordinary controls** (100 each per phase, the same F → SAP rule and exposure average):
  - `eclipse`: `w3_ztf.control_model` unchanged (periodic box, P U(1, 50) d, duration U(0.1, 2) d,
    1.5 mag);
  - `fade`: `w3_ztf.control_model` unchanged;
  - `pspl`: `w3_ztf.control_model` unchanged (a brightening);
  - `common`: 5 consecutive good cadences at 3 mag with the comparison dim fraction set to 0.5. Its
    rejection by the common-mode test holds by construction, as in ZTF.
  - `single_eclipse` (Kepler's 4-yr span shows single eclipses of long-period binaries): one 1.5-mag
    box, duration U(0.1, 2) d, start U over the good range.
  - `quarter_step` (instrumental level change that resets at the quarter roll): 1.5 mag from
    t1 ~ U(good range) to the end of t1's quarter.
- **Seeds:**
  - dev: W3 cell c 3100 + c, control k 3200 + k;
  - val: 3300 + c, 3400 + k;
  - each injection j uses `default_rng([seed, j])`.
- **Hit:**
  - W3: a flagged run overlapping t0 ± t_E √(4 − u0²);
  - controls: any run.
  - Either way, the run must not already be flagged in the uninjected carrier.
  - **Recovered / passed:** some hit passes §5.
  - Every injection stays in its denominator. "Insufficient" (edges not bounded, other tests
    passed) is listed separately and is not a pass.

## 5. Vetting (pass = all of these; computable, `screen_star`)
1. **Common mode** (`w3_ztf.vet`): the mean leave-one-out dim fraction of the comparison stars at
   the run's cadences must be < 0.2 (`COMMON_FRAC`).
   - Comparison stars: dev phase = the other dev stars; val = dev + val; search = all 100.
2. **Bounded fast edges** (v23 `bounds`, `EDGE_FRAC` 0.25, span = max(run length, 1 d)), computed
   per side on that edge's quarter only (`edge_bounds_kepler`).
   - An edge whose nearest normal cadence is in another quarter (separate aperture and
     normalisation) is not bounded.
   - Both sides must be "fast".
3. **One-off** (`w3_ztf.vet`): no other dim run of ≥ 3 cadences anywhere in the light curve.
4. **W3 shoulders** (`shoulder`): on each side, within 0.25 × span of the edge and in the same
   quarter, some cadence with rel ≥ 1.10 and (rel − 1)/σ ≥ 5. Both sides are required.
   - Basis (`model_check.json`, model only): the exposure-averaged W3 excess just outside the edge
     is ≥ 0.393 in every cell (worst case over 12 u0 values and 5 cadence phases).
   - 0.10 is 4× below that. Eclipses, fades and steps have no shoulder.
- **Secondary (reported, not gating):** pass without test 4 (`pass_no_shoulder`).
- **Diagnostic:** PDCSAP in flagged runs.

## 6. Gate, search and stop
- **Dev phase:** the fixed method on dev carriers, plus the uninjected screen of the 20 dev stars.
  - Implementation errors may be fixed, and the affected parts re-run with both results kept.
  - No threshold, seed or method change.
- **Validation phase:** the same on the unread val stars. **Gate (validation only):**
  - **every** one of the 20 cells has recovered / 100 ≥ 0.50, **and**
  - **every** one of the 6 control families has passed / 100 ≤ 0.10.
  - These are pilot thresholds, not discovery significance and not sky false-positive rates.
  - Also reported: the uninjected screen of the val stars (flags and passes still to check).
- **Gate passes:** search the 60 search stars (`search`) without asking.
  - Survivors are taken in sample order, then time. At most 3 are vetted, using:
    - the existing quality flags (all bits within ±1 d);
    - the flux-weighted centroids (`MOM_CENTR1/2`) during the run against its quarter;
    - same-exposure controls (sample stars on the same channel and quarter, and all sample stars
      at the same cadences);
    - PDCSAP;
    - if needed, the official target pixel file.
  - Unresolved cases are "undetermined". A W3-like shape is never called a wormhole candidate.
- **Gate fails:** the search group is not read. The failure is attributed to sampling, product
  processing, ordinary events or implementation, from the stored per-injection rows.
- **No second try:** no seed, threshold or method change after a result is seen.
- **Budget:** if 3 h or 2 GB would be exceeded, completed results and the unfinished scope are
  recorded; denominators are not reduced and the budget is not restarted.

## 7. Limits stated in advance
- **Processing scope:** efficiencies apply to processing from the SAP product onward, not to the
  official pixels-to-light-curve chain. For example, a real event flagged by the pipeline (impulsive
  outlier, PDC exclusions) is not modelled.
- **Grid:** t_E ≤ 10 d, blend ≥ 0.7 and ρ ≤ 0.1. Events covering more than half of a quarter
  (especially short Q1 / Q17) lose their dip in the per-quarter median.
- **Sample:** 100 FGK dwarfs at Kp 12–14, observed in all 17 quarters; not a sky rate.
- **Labels:** injections are simulated; light curves and flags observed; outcomes derived.
