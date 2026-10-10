# Neutrino Frontier

Program started by the owner brief of 2026-10-10 (D-077). Neutrinos are probes only. Exotic physics is a
hypothesis; an anomaly score or an excess is never called a wormhole, warp effect or discovery.
Status ladder per framework: A mathematical possibility, B self-consistency, C compatible with observations,
D testable with public data, E evidence, F engineering feasibility.

Files: [round1_invention.md](round1_invention.md) (N-A, literature-blind), [round1_review.md](round1_review.md)
(N-B, adversarial, with references), [data_audit.md](data_audit.md) (sources, sizes, fields, hosts).
Global state: TASKS.md "Neutrino Frontier", CHANGELOG.md, DECISIONS.md (D-077), SOURCES.md.

## Hypothesis registry

| ID | Framework (round 1) | Verdict (N-B) | Ladder | Next action |
|---|---|---|---|---|
| NF-H01 | F1 granular locality: stochastic distance, flare width ∝ L² | KNOWN-REDUNDANT, UNTESTABLE NOW | B/C | Needs ≥ 2 multi-event neutrino flares at different distances; none public. Park. |
| NF-H02 | F2 twin-sheet flavor (pseudo-Dirac partners, δm² 10⁻¹⁹–10⁻¹⁷ eV²) | KNOWN-REDUNDANT | C | Window already excluded by published IceCube data (arXiv:2406.06476). No compute. |
| NF-H03 | F3 second causal cone: ν–γ lag = ε × Shapiro delay | KNOWN-REDUNDANT | B/C | Standard EP test; robust \|ε\| ≲ 10⁻³–5 × 10⁻³ (SN1987A); D-074 GBM × ICECAT null covers it. No rescan. |
| NF-H04 | F4 sparse nonlocal links: delayed "ghost" neutrino at an unrelated position | CONDITIONALLY VIABLE (phenomenology) | B | Experiment E-NF1 below. TXS-flare ghost test rejected (all-sky trials swamp it). |

Assumptions that keep NF-H04 alive: static links (else causality fails), and a coupling that passes neutrinos
but not photons (assumed, not derived). Its only testable content here is a rate R_g of delayed wide-separation
copies of alert tracks.

## Experiments

### E-NF1: ICECAT-1 ghost-pair limit (NF-H04), pre-registered 2026-10-10
- **Hypothesis:** a fraction R_g of astrophysical alert tracks have a delayed copy (ghost) at an unrelated sky
  position after Δt_g ∈ [0, 180 d].
- **Ordinary explanations:** alert-rate and uptime modulation (seasonal, the 2019 realtime-stream change),
  Dec-dependent acceptance, split or re-reconstructed events of one run, atmospheric-muon background pairs.
- **Statistic:** ICECAT-1 v4 (340 tracks, CR_VETO dropped) unordered pairs with separation > 3σ_comb and > 0.1°
  (D-074 "wide"), per lag bin 0–10 s, 10–100 s, 100 s–1 h, 1 h–1 d, 1–7 d, 7–30 d, 30–180 d. Two statistics per
  bin: the pair count and the signalness-weighted count Σ s_i s_j (a ghost pair of an astrophysical parent is
  astrophysical on both sides). 14 cells. ASSUMPTION: the weighting.
- **Nulls:** D-074 `jit` (Dec and hour angle kept, within-year permutation, ±3 d) for bins ≤ 7 d; a cyclic
  ±1 yr jitter inside the catalogue span for 7–30 d and 30–180 d (±3 d is too short there; ASSUMPTION: alert
  rate stationary on 1 yr). Within-year permutation alone keeps every pairwise lag, so it has no power here.
- **Reproduce first:** D-074 wide counts 36 (jit 25.45 ± 4.4) at 1 h–1 d and 138 (149.39) at 1–7 d.
- **Detection:** family-wise 3σ (analytic p × 14 ≤ 1.35 × 10⁻³), as D-074. A hit goes through `jitday`-type
  vetting, run/day split checks and `/vet-candidate` before any interpretation.
- **Sensitivity and limit:** inject ghosts (R_g × signalness of each parent; Δt log-uniform in the bin; RA uniform,
  Dec from the catalogue, error resampled). R_g,50 = smallest R_g detected in ≥ 50 % of trials; 95 % upper limit
  = smallest R_g whose injected statistic exceeds the observed one in ≥ 95 % of trials.
- **Deviations from this pre-registration (made before the final run, recorded in D-077):** per-cell p is
  empirical from 20,000 scrambles (floor 5 × 10⁻⁵ < 9.6 × 10⁻⁵), not analytic, because the weighted statistic is
  not an integer count; the limit is a 95 % CLs limit (the classical one excluded R_g = 0 on a low fluctuation),
  computed as a null draw plus the ghost excess; a ghost's Dec and error come from the same catalogue event.
- **Result (2026-10-10): null.** Reproduces D-074 (36, 138). Pooled global p = 0.058; largest cell 1 h–1 d
  weighted (z = 2.7, Bonferroni p = 0.094; its 36 pairs inspected: all years, 6 same-run, two triplets). 95 % CLs
  limits on R_g: ≤ 0.02 (≤ 10 s), ≤ 0.05 (10 s–1 h), ≤ 0.2 (1 h–1 d), ≤ 0.1 (1–7 d), ≤ 0.5 (7–30 d); 30–180 d
  excludes only R_g = 1. NF-H04 stays at ladder B; the parameter region above these limits is rejected for alert
  tracks. `results/nf/ghost_pairs.json`; D-077.
