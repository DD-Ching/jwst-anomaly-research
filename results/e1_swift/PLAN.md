# E1 + Swift BAT: frozen plan (2026-10-18, before any count is read)

**Hypothesis (owner idea 4, D-074):** events from different directions are dependent at lags no ordinary path
explains. This run adds Swift BAT GRBs as a fifth node set (TASKS.md Now 0: "IceTracks-DR2, Swift, Einstein Probe").

**Input (observed):** HEASARC TAP `swiftgrbba` (Swift GRB compilation by the burst advocates; 2,043 rows,
2004-12-17 to 2026-10-07 at the first look, which read only row counts, dates and column units). Columns `name`,
`time` (trigger MJD), `trigger_obs`, `ra`, `dec`, `error_radius` (BAT 90 % radius, arcmin). Pinned by sha256 in
`data/manifests/e1_swift.ecsv`. The four D-074 catalogues come through `scripts/e1_events.py`'s `fetch` (its own
pins).

**Selection (ASSUMPTIONs):**
- `trigger_obs == "Swift"` (on-board BAT triggers). Rows triggered by Fermi, INTEGRAL, SVOM, ground or GUANO
  analysis are dropped: their time is another instrument's trigger and their selection differs.
- Rows without finite `time`, `ra`, `dec` or `error_radius` are dropped.
- 1σ = `error_radius` / 2.146 (90 % radius of a 2-D Gaussian, as ICECAT in D-074), in degrees.

**Channels.** Tested family: BAT–BAT (same, wide), BAT–ICECAT (same, wide), BAT–CHIME (same, wide), BAT–GW (all),
over the five D-074 lag bins (0–10 s, 10–100 s, 100 s–1 h, 1 h–1 d, 1–7 d): 35 cells.
BAT–GBM is the **positive control** and is not in the family: the same GRB is seen by both, so its same-direction
cells at 0–100 s must show a large excess, and its wide cells at short lags carry GBM localisation tails of the same
burst (ordinary). The D-074 GBM/ICECAT/GW/CHIME channels are not re-tested.

**Nulls (D-074 definitions, `event_network`):** `jit` (primary; per-year permutation then a ±3 d shift; BAT shifts
are k × 5,700 s ± 300 s, the Swift orbit (≈ 95 min, ASSUMPTION; it decays over the mission, so phase is kept only
approximately)) and `jitday` (vetting), 10⁴ scrambles each. Pooled per-cell and global p as D-074.

**Frozen rules (functions in `scripts/e1_swift.py`, unit-tested before the run):**
- `control_passes`: BAT–GBM same-direction cells 0–10 s and 10–100 s together hold ≥ 50 pairs and the summed
  count's analytic p under `jit` is ≤ 1.8 × 10⁻⁵. If it fails, there is **no result** (time or position
  convention broken).
- `outcome`: "excess" if any tested cell's analytic p (jit) × 35 ≤ 1.35 × 10⁻³ (family-wise 3σ, D-074); otherwise
  "null". An "excess" goes to vetting (jitday null, duplicate rows, Swift retriggers and safe-mode gaps), never to
  interpretation.
- Limits as D-074: injections into each tested cell (40 trials per grid point), the 95 % upper limit on dependent
  pairs and the rate per eligible anchor, valid only where an injected n reached 50 % detection.

**What changes the plan:** a failed control stops the run (record why). An "excess" in a same-catalogue BAT–BAT
cell at ≥ 1 h is expected first to be a rate-modulation effect, as CHIME–CHIME was (D-074 addendum); it would need
the rate-modulated null before any other reading.
