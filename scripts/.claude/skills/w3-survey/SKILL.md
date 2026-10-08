---
name: w3-survey
description: Recipe for the W3 (inverted microlensing - flux vanishes between caustic spikes; negative-mass n = 1 or Ellis lens) search on any light-curve survey - adapter, selection emulation, pre-screen, fit, vet, revet, inject, limit - with the failed-approach rules from OGLE-IV, Gaia DR3 and MOA-II. Use when adding a light-curve survey or changing scripts/w3_microlensing.py, scripts/w3_gaia.py or scripts/w3_moa.py.
---

# W3 survey recipe

Rules: root CLAUDE.md, `scripts/CLAUDE.md`, `src/jwst_anomaly/CLAUDE.md`. Decisions: D-047 (prediction), D-052,
D-054 (signature layer), D-057–D-059 (OGLE-IV), D-061 (Gaia DR3), D-062 (MOA-II). Results:
docs/exotic_limits.md "W3 in the OGLE-IV microlensing samples", "W3 in MOA-II".

## 1. Choose the data: the selection decides everything
- A published, PSPL-selected event list cannot constrain W3: inject W3 through its selection first. OGLE-IV Mróz
  samples pass 0 / 600, Gaia DR3 `vari_microlensing` 0 / 240 (single-bump / constant-outside / PSPL gates).
- Use light curves released **before** any bump or PSPL cut. MOA-II 9-year Cut-0 (sign-agnostic difference-image
  detections) keeps W3; Gaia DR4 epoch photometry is the next candidate when released. KMTNet and OGLE EWS are
  event-finder or terms-limited (EWS needs the owner).

## 2. Adapter
- `src/jwst_anomaly/<survey>.py`: a `signatures.LightCurveSurvey` (`events`, `light_curve`, `efficiency`), inputs
  pinned by sha256 + a manifest; light curves via `signatures.standard_light_curve` (magnitudes) or
  `standard_flux_light_curve` (difference fluxes that can be negative).

## 3. Commands
- OGLE-IV: `python scripts/w3_microlensing.py fit|merge-chunks|vet|revet|sheet|inject|audit|limit|summary|manifest
  [--sample bulge2019|disk2020]` (`fit --chunk K/N` writes tracked tables under `results/w3_ogle/`).
- Gaia DR3: `python scripts/w3_gaia.py fetch|fit|inject|summary|manifest`.
- MOA-II: `python scripts/w3_moa.py prescreen|fit|merge-chunks|vet|sheet|inject|limit|manifest`
  (`fit --chunk K/N` → `results/w3_moa/`).
- Always `OMP_NUM_THREADS=1` with process pools; `JWST_ANOMALY_DATA` set.

## 4. Vetting chain (cheapest first; never drop a test for speed)
refit with every ordinary model (PSPL, FSPL, parallax with |π_E| ≤ 5, D-058) → robust errors / isolated outliers →
variable baseline → free blend per season → season drifts → binary source → binary lens → VSX / Gaia DR3 variability
(CDS XMatch) → arXiv names → `revet`: feature coverage (≥ 3 epochs where the models differ by > 3σ, classified at
the epochs), jackknife (drop ≤ 3 influential epochs, keep ≥ 3 in the feature; skip with exactly 3), two unrelated
PSPL bumps. MOA adds shared-epoch tests (field- and chip-wide Poisson), a neighbour test, an eclipse-dip model and
≥ 3 nights in the feature. Look at the contact sheet before concluding.

## 5. Limit
- Inject through the whole chain on the survey's own carriers (selection emulation, pre-screen, fit, vetting);
  grid t_E ∈ {3, 10, 30, 100, 300} d × ρ ∈ {0.01, 0.1}, u₀ ~ U[0, 2), plus PSPL controls; ≥ 200 per cell when
  luminosity-function weighting reduces n_eff.
- Γ₉₅ = 3 / (N_s T ε) with zero survivors; N_s from published star counts where they exist (Nunota et al. 2024 for
  MOA), else a stated ASSUMPTION; rule-of-three floor only for zero recovered injections.

## Failed approaches (rules)
- Exotic fits started only from the PSPL solution miss the W3 geometry: start with caustic spikes on pairs of maxima
  and on absolute t_E.
- ΔBIC alone is not a candidate test: season blends and feature coverage removed 120 of 127 OGLE flags; in MOA every
  dip flags against ordinary microlensing, so vetting decides.
- A jackknife must keep ≥ 3 epochs in the exotic feature or it kills real short events.
- An exotic fit whose spikes sit on two bumps years apart needs the two-unrelated-events test (bounded, nesting PSPL).
- Feature epochs must be classified with both models evaluated at the epochs, not on a coarse grid.
- The rule-of-three floor applies to "recovered" (selected and flagged), not "selected".
- MOA pre-screen: a plain box deficit against the median passes ~10k light curves (trends, season offsets) — compare
  against the local flanks *and* the median, normalised by each light curve's red noise; a flanks-only notch passes
  PSPL bumps — also require the box below the median; 97 % of shape passes sit at epochs shared across many objects —
  run the shared-epoch test before fitting.
- Unbounded t_E or ρ makes single fits take ~200 s: cap t_E ≤ 1,000 d and ρ ≤ 0.3 for MOA.
- Fitting every shape pass before the shared-epoch cut wasted a chunk (131 / 133 flagged dip-shaped variables).
- The variable-baseline test is the largest MOA efficiency loss (35 % of quiet carriers have χ²/dof > 2 alone):
  calibrate it on the carrier distribution; re-run injections after any pre-screen or vetting change.
- The Gaia Extractor cuts cannot be emulated from the paper (guessed definitions fail 126 / 163 real events);
  a single-id DataLink request returns bare CSV; retry truncated chunked replies.
- `table3.dat` (Mróz 2019) rows overflow their byte ranges: split on whitespace, sexagesimal on colons.
- Pin BLAS threads; two 4-process pools on 4 cores stalled both. Wait on output files, not `pgrep -f`.
- 6 of 212 arXiv name queries hit the rate limit (recorded as −1); re-run before any publication-facing claim.
