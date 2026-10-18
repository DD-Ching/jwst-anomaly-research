# E1 + Einstein Probe: frozen plan (2026-10-18, before any pair count is read)

**Hypothesis (owner idea 4, D-074):** events from different directions are dependent at lags no ordinary path
explains. This run adds Einstein Probe (EP) WXT transients as a sixth node set (TASKS.md Now 0; handoff of #351).

**Input (observed, then parsed):** no public EP event table exists (2026-10-18 web search: papers carry their own
PDF tables; the EP site indexes GCN Circulars). Source: the GCN Circulars JSON archive
(`https://gcn.nasa.gov/circulars/archive.json.tar.gz`, ~31 MB, rewritten daily), pinned by sha256 in
`data/manifests/e1_ep.ecsv`. `scripts/e1_ep.py build` writes `ep_events.ecsv` (one row per EP name, with the
circular id) and `build_summary.json`; that table is the record the counts use. Parse rules (ASSUMPTIONs, in
`parse_circular` / `ep_events`, unit-tested):
- EP name from the subject (`EPyymmdd[a-z]`); the body mentions WXT before the position.
- Time: the light-curve onset ("began / started ... at <UTC>") if stated, else the first UTC timestamp (often the
  trigger; it can be an observation start up to ~20 min earlier). The onset must lie within 1 d of the name's date.
- 90 % radius: the first "uncertainty [radius] of x arcmin", else the first arcsec one (flagged `r90_unit`);
  σ = r90 / 2.146 as D-074.
- One row per name, from the lowest-numbered circular that parses; names with a circular calling the trigger
  "not a real source" are dropped.
- Build check (before this plan): 230 events, MJD 60384.8–61328.7; 5 random rows read against their circulars,
  5 / 5 correct in time, position and radius (not stored beyond this line).

**Overlap:** EP starts 2024-03; ICECAT-1, CHIME/FRB Catalog 2 and the pinned GWTC end before it, so only GBM and
BAT (both pinned in `data/manifests/e1_swift.ecsv`, reused unchanged) overlap EP. Those three channels are not run.

**Channels.**
- Tested family (18 cells): EP–EP (same, wide) × five D-074 lag bins (10); EP–BAT and EP–GBM (same, wide) at
  1 h–1 d and 1–7 d (8). Short-lag EP–BAT / EP–GBM cells are the same GRBs seen twice (ordinary) and are the
  control, not tested.
- **Positive control:** EP–BAT plus EP–GBM same-direction pairs at 0–10 s, 10–100 s and 100 s–1 h, summed. EP
  time conventions (trigger vs onset vs observation start) can move a joint GRB by minutes, hence 0–1 h.

**Nulls (D-074, `event_network`):** `jit` (primary; per-year permutation, then a shift by k × 5,760 s ± 300 s for
EP, ASSUMPTION: EP's ~96 min orbit at ~600 km; GBM and BAT keep their D-074 / #351 quanta) and `jitday`
(vetting), 10⁴ scrambles each.

**Frozen rules (functions in `scripts/e1_ep.py`, unit-tested before the run):**
- `control_passes`: the summed control count ≥ 5 pairs (ASSUMPTION: EP–GRB joint detections are a minority of EP
  transients) and its analytic p under `jit` ≤ 1.8 × 10⁻⁵. If it fails: **no result** (time or position
  convention broken).
- `outcome`: "excess" if any tested cell's analytic p (jit) × 18 ≤ 1.35 × 10⁻³ (family-wise 3σ, D-074), else
  "null". An "excess" goes to vetting (jitday null, duplicate names, the same transient reported twice, Galactic
  flares and stellar repeats), never to interpretation.
- Limits (D-074 method, 40 trials per grid point) for EP–EP only, valid only where an injected n reached 50 %
  detection.

**What changes the plan:** a failed control stops the run. An EP–EP excess at ≥ 1 h is expected first to be a
repeating Galactic or stellar flarer (same direction) or EP pointing strategy (wide); it would need those checks
before any other reading.
