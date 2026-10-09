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
- MOA-II: `python scripts/w3_moa.py --field gbF prescreen|merge-prescreen|fit|merge-chunks|vet|sheet|inject|limit|
  summary|combine|manifest|run-field` (D-062, D-068). `prescreen` streams the field tar by HTTP range reads (never
  stored) into one tracked table per 4 GiB in `results/w3_moa/prescreen/`; `run-field` runs every stage, skipping
  finished pre-screen and fit chunks; `combine` sums N_s T ε over the tracked `limits_gb*.ecsv`. Long queues run
  detached (`setsid nohup`, logs under `derived/w3_moa/`): harness background tasks are killed after 30 min.
- Always `OMP_NUM_THREADS=1` with process pools; `JWST_ANOMALY_DATA` set; install the `mulens` extra (fresh cloud venvs
  lack it, and without it the parallax and binary-lens tests were skipped silently; the stages now refuse to run).

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
- A fixed variable-baseline threshold (χ²/dof > 2) lost 73 / 188 vetted MOA injections (35 % of quiet carriers
  exceed it from red noise): use the field's 95th percentile of quiet χ²/dof (gb22 5.3, gb21 6.3; loss 1–2 / 130).
- A far-field exotic fit (u0 ≫ 1) with cancelling giant fs/fb can mimic any smooth dip — require the fit domain
  and a physical source flux (`exotic_in_domain`: u0 < 2, f_s ≤ 3 × DoPHOT or Gaia RP reference) and test a smooth
  Gaussian dip (`smooth_dip`); gb20-R-4-0-49379 (a red giant's ~270-d dimming) passed every other test.
- An exotic feature must be bracketed by baseline on both sides; a one-sided step or secular change is not an umbra
  crossing (`feature_bracketed`, ≥ 20 epochs each side; `step_ramp`, a level change plus ramp with free t_s):
  gb19-R-4-4-31159 (a season-boundary step that never recovers) passed every other test.
- Uniform injection magnitudes re-weighted to the luminosity function leave n_eff ≈ 20–25 of 60: draw them from
  the LF (`--sampling lf`).
- Streaming: a parent that downloads ranges and ships bytes to workers was OOM-killed (8–14 GB); each worker reads
  its own ranges with prefetch threads. Read ranges with `stream=True` and accept only 206: the archive sometimes
  answers 200 with the whole tar (retry it unread). Tracked pre-screen rows must carry every column the fit reads
  (`SCAN_KEYS`). After a killed run, look for orphaned pool workers (`ps`) before benchmarking. `*.log` is
  gitignored (heartbeat file: `results/w3_moa/progress.txt`).
- The Gaia Extractor cuts cannot be emulated from the paper (guessed definitions fail 126 / 163 real events);
  a single-id DataLink request returns bare CSV; retry truncated chunked replies.
- `table3.dat` (Mróz 2019) rows overflow their byte ranges: split on whitespace, sexagesimal on colons.
- Pool sizing, BLAS pinning and how to wait on runs: `scripts/CLAUDE.md` "Parallel topology".
- 6 of 212 arXiv name queries hit the rate limit (recorded as −1); re-run before any publication-facing claim.
