# E1 + IceTracks-DR2: frozen plan (2026-10-18 23:30 UTC, before any pair count is read)

**Hypothesis (owner idea 4, D-074):** events from different directions are dependent at lags no ordinary path
explains. This run adds IceCube IceTracks-DR2 tracks as the seventh node set (TASKS.md Now 0; handoff of #353).

**Input (observed):** IceTracks-DR2 v3.1 (doi:10.7910/DVN/MMIIZA), the 28 events/uptime files already pinned in
`data/manifests/nf_icetracks_dr2.ecsv` (D-079) and read by `nf_ghost_dr2.fetch` (reduced cache only). The other
five catalogues reuse `data/manifests/e1_swift.ecsv` unchanged (GBM, BAT, ICECAT-1, GWTC, CHIME Cat 2).

**Node set (ASSUMPTIONs, as E-NF1b / D-079):** northern tracks, Dec > −5°, log10 E ≥ 4.0 (the higher-purity cut
D-079 used); one track per (run, event) readout (the first subevent), so a readout cannot pair twice; 1σ = AngErr.
DR2–DR2 pairs are not re-tested here (D-079 covers ≤ 1 h; its forecast shows no reach at longer lags).

**Channels.**
- Tested family (35 cells): DR2–GBM, DR2–BAT, DR2–CHIME (same, wide) × five D-074 lag bins (30); DR2–GW (all,
  no GW positions) × five (5).
- **Positive control:** DR2–ICECAT same-direction pairs at 0–10 s. ICECAT-1 Gold/Bronze tracks are IceCube track
  events, so many are in DR2 with the same event time (ordinary: the same neutrino twice). DR2–ICECAT is never
  tested (neutrino–neutrino pairs belong to E-NF1).

**Nulls:** `jit` (primary): DR2 events get an independent uniform ±3 d shift redrawn until inside a good run
(`neutrino_frontier.jitter_uptime`, the D-079 null; hour angle kept); the other catalogues as in #351
(`scramble_jit`, GBM / BAT orbit quanta). `jitday` (vetting): the other catalogues as #351's `jitday`, DR2 as `jit`.
10⁴ scrambles each.

**Frozen rules (functions in `scripts/e1_icetracks.py`, unit-tested before the run):**
- `control_passes`: ≥ 10 control pairs (ASSUMPTION: at least ~10 % of the ~150 northern ICECAT tracks pass the
  DR2 cuts and match within 3σ_comb) and analytic p (jit) ≤ 1.8 × 10⁻⁵. If it fails: **no result**.
- `outcome`: "excess" if any tested cell's analytic p (jit) × 35 ≤ 1.35 × 10⁻³ (family-wise 3σ, D-074), else
  "null". An "excess" goes to vetting (jitday, detector-configuration and season dependence, the same GRB or FRB
  repeated, atmospheric seasonal rate), never to interpretation.
- Limits (D-074 method, 40 trials per grid point), valid only where an injected n reached 50 % detection.

**What changes the plan:** a failed control stops the run (time or position convention broken). An excess in
DR2–GBM at ≤ 100 s same would first be read as a known GRB–neutrino coincidence search (IceCube GRB papers), and
any 1 h–7 d excess first against seasonal / uptime structure that the ±3 d jitter does not remove.
