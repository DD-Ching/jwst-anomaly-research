# Bullet Cluster CANUCS extra-image test (WLS source): gate C2 passes; no flag (2026-10-24)

Plans: `PLAN.md` (frozen 2a7e851; script and CANUCS sha256 pins 1ae40dd, before the run) and the post-failure
deviation `PLAN_1b.md` (frozen 6d75951; code 72d7b61; both before any CANUCS pixel was read). Outputs:
`summary.json`, `control_bullet.ecsv`, `predicted_bullet.ecsv`, `without_contain/` (rule 8, reported only),
`frame/` (rule 3), `summary_wls.json` (rule 2).

## Gates (derived)
| Gate | Value | Threshold (ASSUMPTION) | Outcome |
|---|---|---|---|
| F (frame, recomputed) | 92 / 135 gold matched ≤ 0.5″; median dx −0.0005″, dy −0.0000″ | ≥ 30; \|median\| ≤ 0.1″ | pass |
| K-WLS (stored, rule 2) | rms 0.391″, 135 / 135 within 3″ | ≤ 1.0″; ≥ 80 % | pass |
| W1 (catalogue vs `W_hdr`) | 20,138 rows; 0-based rms 1e-10″; 1-based rms 0.057″ | rms ≤ 0.01″, p99 ≤ 0.04″ | pass (catalogue is 0-based) |
| W2 (pixel centroids) | 51 / 51 usable; median Δx −0.035 px, Δy +0.205 px | ≥ 50; \|median\| ≤ 0.25 px | pass (Δy close to the limit; 0.008″) |
| C2 shared (rules 5–6) | all 58 / 67 (0.866); arc-referenced 37 / 42 (0.881) | ≥ 8, ≥ 80 %; arc ≥ 3, ≥ 80 % | **pass** |

W1 agrees to 1e-10″, so the catalogue's `X`, `Y` were computed with the same 0.04″ header geometry (the `CD` /
`PIXSCALE` keywords are the inconsistent ones). Only 51 point sources passed the W2 selection (bright, isolated,
unsaturated), one above the minimum.

## Targets (model_prediction vs observed)
22 predicted unmatched counter-image targets (WLS source). Local classes: recovered 1 (N12, 27 σ), confused 1 (K5),
undetectable_local 8, undetectable 5, inconsistent_reference 5, no_reference 2, **absent 0 → no flag** (rule 7).
Controls: recovered 61, confused 7, absent 8, ambiguous 6, undetectable_local 16, inconsistent_reference 19,
no_reference 15, undetectable 2.

## Reported, not rules (rule 8)
- Without the containing-source rule: C2 fails (51 / 67 = 0.761, arc 31 / 42 = 0.738); no flag.
- Plain gate C2: passes (68 / 82, arc 47 / 57).
- Without the local-noise rule (global depth): C2 passes (80 / 98) and 5 targets read `absent`: N2a, N5a ×2, N5b ×2.
  Under the frozen local-noise rule all five are `undetectable_local` (pred / σ_loc 1.0–2.4; σ_loc 36–179 nJy in 0.3″
  against ~6 nJy at N12).

## Contact sheet (looked at, not stored: scratchpad)
The 10 locally tested targets on the F200W bCG-subtracted mosaic (6″ boxes). High σ_loc has visible causes: K7a / K9a
sit in bright-star diffraction spikes; N2a and N11 on bCG-subtraction residuals and masks; N5a / N5b in cluster-galaxy
halos (one N5a position touches a bright galaxy). N12 shows the predicted source. No predicted position shows an
obvious empty patch where a bright image should be.

## Conclusion
On the Bullet Cluster, with the CANUCS DR1 best-fit maps (WLS source) and the CANUCS DR1 reference route, gate C2
passes and no extra-image flag is raised. This is a null with limited reach: 13 of 22 targets are untestable at the
5 σ local depth (undetectable / undetectable_local), mostly near bright stars, bCG residuals and cluster members; the
5 global-depth `absent` readings dissolve under the frozen local-noise rule. No exotic reading is considered.

## Limitations
- `PLAN_1b.md` is a post-failure deviation (header-derived WCS); W2's Δy (0.2 px) is near its threshold.
- System z fixed (ASSUMPTION); WLS source only (the mean source fails gate K on K1b / K1e).
- The local-noise rule decides 5 targets; a deeper or star-/bCG-cleaned image would be needed to test them.

## Failed approaches
- `PLAN.md` rule 5 as written: `WCS(header)` of the `bgsub-sci` mosaic reads 3600″ per pixel (`PC` without `CDELT`),
  and the `CD` keywords say 0.08″; neither is the 0.04″ grid.
- A 10-min background job limit killed the first `PLAN.md` attempt (rules 1–3 take ~16 min); `--plan-1b` reuses the
  stored rule 2–3 outputs and runs in ~9 min under `nohup`.
