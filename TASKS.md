# Tasks

Prioritized queue. Agents pick from **Now** first; in-flight work is visible as open PRs.
Evidence for the current priorities: D-023 (owner, 2026-10-08), D-025/D-026 and the CHANGELOG entries of 2026-10-08.

## Division of labour (owner brief 2026-10-10)
- **Local-only (owner machine):** GitHub releases for derived data; files that exist only locally.
- **Cloud routine:** TXS 0506+056 SkyLLH benchmark on IceTracks-DR2 (done, D-080), E-NF1 on IceTracks-DR2 (short lags done, D-079; rest per Neutrino Frontier item 1), round 2,
  all scans.

## Exploration trial (owner-approved 2026-10-12 → 2026-11-09; D-083)
Rules, weekly split (45 % routine search / 30 % high-risk / 15 % verification / 10 % upkeep), units and ledger:
`docs/state/trial.md`. Provisional conclusions of open needs-human PRs, and what is "done; do not redo":
`docs/state/pending.md`.
- **X1 far-side SV hot spots at Sgr A*** closed 2026-10-13 (gates passed as written, 3.7 / 20 h; the tool is cross-checked, the position comparison is exploratory, the S2 flux-contrast observable is untested; no candidate): conditional results and
  the reopen inputs are in `docs/state/units/X1.md` "Closure". Tracer: `scripts/x1_tracer.py` (angles and image
  weights verified by second routes, `scripts/x1_s4.py`, `scripts/x1_s4b.py`).
- **X2 SV-throat echoes vs published post-merger limits** (authorized 2026-10-13): state file
  `docs/state/units/X2.md`; S1–S5 done. The S5 background is radially unstable, so there is no further echo search on it.
  S6 rotation screen done. S6a done: the 5D merged modes continue as an oscillating growing pair (verified at J/A
  0.121 and 0.142). Writeup drafted (`docs/state/units/X2_s6a_summary.md`); next is a human expert review, the owner's
  call. Paused; no automatic follow-up. Archival brief:
  `docs/results/s6a_5d_rotating_wormhole_modes.md` (+ .pdf).

- **S0-2026-10-16 message-and-reply window in unstable wormholes** (owner check 2026-10-16):
  `docs/state/units/S0-2026-10-16.md`. Stopped. The literature already shows test-ray passages and a possible round trip within a window of
  order r_throat/c in two-universe Ellis models; published maintenance uses active phantom injection. No new unit; reopen conditions are in §7.
- **S0-2026-10-15b Kapitza dynamic stabilisation** (owner screen 2026-10-15): `docs/state/units/S0-2026-10-15b.md`.
  Paused. The proposed Ellis + static radiation-fluid background does not exist (re-check), so no solver. A follow-up
  needs a specified radiation closure (distribution, stress tensor, boundaries, confinement or supply, fixed
  matter theory).
- **X4 wave-optics GW lensing by negative-mass objects** (authorized 2026-10-15 by S0-2026-10-15): state
  file `docs/state/units/X4.md`. S1 done (simulated; 18 / 36 preliminary passes against aligned-spin models).
  S1b done: against precessing, higher-mode binaries (IMRPhenomXPHM), R2 and R3 no longer pass the screening
  threshold (ρ_req > 30). S1c done (R1): residual kept against validated aligned-spin fits (ρ_req 5.8,
  bound-conditional, provisional); undetermined with XPHM (the precessing fitter fails its check; best ρ_req 9.2). S1d done: the reference
  global search (SciPy DE) fails its injection checks (lensed procedure flawed, unlensed budget-limited), so R1
  stays undetermined. S1e done: the joint search makes U reliable; P passes 3/4 runs, so the tool still fails
  and R1 stays undetermined. No pilot proposed. Open, on the owner's go: more P budget (an estimated 4–5 h on 4
  cores) or precession-parametrisation advice; or an R1 U fit testing outcome 1 only (about 0.5 h). S1f done: that U fit finds no counter-example (best FF 0.9135,
  ρ_req 7.4); R1 stays undetermined. S1g done (2026-10-21): with only the budget raised (DE 1,000 generations, 20
  refinement rounds) the P search passes 6 / 6 frozen checks, S1e's failure and 2 new validation injections included;
  on R1 its seeds converge to different optima (0.9476 vs 0.9459), so R1 stays undetermined (P and U search
  reliability). S1h done (2026-10-21): seeds 24–26 per ordinary model reproduce the common optima (P 0.94588, U
  0.89229); the isolated bests (P 0.9476, U 0.9135) stay, so the seeds still disagree (1.75e-3, 0.0212 > 2e-4):
  tool/optimization undetermined, R1 undetermined. Extra-seed route closed. Open, on the owner's go: a search method
  with a global-coverage diagnostic (e.g. nested sampling on the same score) or GW data-analysis advice. S1i done
  (2026-10-21): dynesty nested sampling passes U validation (V7, new V10) and, with the S1e polish, both R1 runs beat the old U best
  (0.9146, 0.9151 vs 0.9135; verified), in a different region, disagreeing by 5.6e-4; U optimum still not
  established, R1 undetermined. Open, on the owner's go: P with the same method needs a validated interpolated lens
  factor first (exact factor about 1 s per call; estimated 12.5–20 CPU-h plus that check). S1j done (2026-10-21):
  spin-bound U validation fails on B2 (near-cap, anti-aligned, q 9.35; ε 6.9e-4 with a converged polish), B1 passes;
  R1 not refined; tool/optimization undetermined. P acceleration not supported until the U search handles R1's spin
  regime; the owner decides (spin parametrisation or GW data-analysis advice).
- **X3 W3 in unselected Gaia epoch photometry** (paused 2026-10-15; reopen when the DR4 epoch-astrometry fields and
  precision are confirmed. 2026-10-16 check: fields confirmed in the draft data model, precision not yet
  published; waiting for the DR4 release documentation) (authorized 2026-10-14 by S0-2026-10-14): state file
  `docs/state/units/X3.md`. G0 passed on dev. G1 done: held-out GAPS screened once (395,764 exposed sources),
  196 flags, every one failing a costed (post hoc) ordinary rule, no candidate; conditional limit Γ₉₅ < 1 × 10⁻⁵ per star per year for
  photocentre-stable events (blended events are mostly vetoed by RUWE, a model estimate). G2 prerequisite done:
  astrometry mostly repeats the ordinary blend effect; only the lens image shift is new, measurable for G ≲ 17,
  θ_E ≳ 3 mas and t_E ≈ 100 d (under the assumed per-transit σ). Next: wait for ESA's DR4 epoch-astrometry documentation (2 Dec 2026); develop only if it carries AL
  coordinates, scan angles and σ for those sources.

## Now (owner direction 2026-10-09: System A/B hypothesis rounds, D-069; docs/hypotheses/)
0. **E1 causal event network** (D-074; owner idea 4): GBM × ICECAT-1 × GWTC × CHIME Cat 2 pair counts by lag and
   separation. The five requested channels are null (family p = 0.27 on the re-pinned GBM snapshot). The
   CHIME–CHIME 1 h–1 d wide excess (z ≈ 3.9) is reproduced by a calibrated rate-modulated null (week-scale detection-rate modulation; its cause is a
   hypothesis). The 100 s–1 h wide cell drops to z = 2.5 under it (family-wise p = 0.26, D-074 addendum). Next:
   - a time-resolved CHIME uptime series (only from the collaboration: owner decision); it would also decide
     the 100 s–1 h residual;
   - signed-lag ("which comes first", addendum 3), antipodal (addendum 2) and GW sky-map channels (addendum 4:
     282 maps, global p = 0.83) and signed-lag GW cells by sky-map class (addendum 5, p = 0.56) are null;
   - ~~Swift~~ done 2026-10-18 (BAT, D-074 addendum 6, `results/e1_swift/RESULT.md`: 35 cells null, control passes);
     ~~Einstein Probe~~ done 2026-10-18 (GCN-circular table, D-074 addendum 7, `results/e1_ep/RESULT.md`: 18 cells
     null, control passes; EP–GBM 1 h–1 d same is the smallest cell, × 18 = 0.039); ~~IceTracks-DR2~~ done
     2026-10-18 (northern tracks, D-074 addendum 8, `results/e1_icetracks/RESULT.md`: 35 cells null, control passes).
     All listed node sets are done; E1 now waits only on the CHIME uptime series.
     ~~DR2–ICECAT control deficit~~ explained 2026-10-19 (D-074 addendum 9, `results/e1_icetracks/selfpair/RESULT.md`:
     the null displaces the 276 shared events; removed, z −1.24 / +1.01).
     BAT–GBM control deficit checked 2026-10-19 (D-074 addendum 10, `results/e1_swift/selfpair/RESULT.md`): long
     lags behave as the displacement predicts; frozen rule undetermined (100 s–1 h joint GRBs; post hoc only).
0. **D1 distance self-consistency** (D-073 + addenda; owner idea 1): H0LiCOW (R and D_dt LOO), 94 localized FRBs
   (NE2001 and YMW16) and TDCOSMO 2025 power-law D_dt LOO on 8 lenses (addendum 2: max 1.12σ with κ_ext): all null.
   Reach: R ×2–8 / ×0.2–0.4; D_dt ×1.36–1.85 / ×0.41–0.78. Composite model (SDSS1206 only): null
   (addendum 3). Kinematic D_s/D_ds from the hierArc likelihood pickles on 76 lenses (TDCOSMO, SLACS KCWI/SDSS,
   SL2S): null, max 2.49σ vs 5.78σ, detectable ×0.21 / ×5.6 (addendum 4). Paused: revisit with λ_int-free
   (spatially resolved JWST/KCWI) kinematics or new time-delay lenses.
1. **S1 burst twins** (D-071): Fermi GBM null, < 7.2 × 10⁻³ twin pairs per eligible burst (ratio 1, 15 % band
   mismatch). **Done on unmerged PRs (provisional; see `docs/state/pending.md`; do not redo):** the s ≠ 1 chain and
   the full M2 re-screen (#167); TTE for short and faint bursts and the catalogue E_peak joint constraint (#169; the
   S1 energy upgrade is paused). CHIME/FRB Catalog 2: route and 307-burst pilot in issue #153 (research note); 30-event
   validation 2026-10-21 (D-071 addendum, `results/s1_chime/RESULT.md`): the criteria do not carry over as frozen
   (recovery 6.9 % at ratio 1; 0 / 133 chance and 0 / 3 same-repeater passes), so no CHIME search. Reopen only with
   real same-signal pairs to calibrate the χ² systematic, or an owner-approved scattering re-convolution in the twin
   model. Remaining: a generative pulse-model null.
2. **S2 flat-kernel SN residuals** (D-070 + addendum): Pantheon+ × LS DR9 null (γ < 0.025), DES-SN5YR × LS DR9 null
   (γ < 0.019 mag per unit T/⟨T⟩, α = 2, one-sided 95 %, dilution- and chain-corrected); α = 1 is kernel-degenerate. The
   lensing positive control is only ~1.7σ (DES) / ~2σ combined. ~~LS DR10 z < 22 column~~ done 2026-10-18, control
   gate failed (γ_L 1.1σ, `results/s2_flat_kernel/dr10/RESULT.md`): **S2 paused**. Reopen only with DES Y3 Gold or a
   spec-z / deeper column, plus the full DES-SN5YR covariance.
3. **S3 hybrid images** (D-072): COSMOGRAIL doubles null, sensitivity only r ≈ 0.4–2.4 (red quasar variability makes
   the copy collinear with the main image). Next, low priority: SN Refsdal / SN H0pe imaging at the model lags.
4. **System A rounds paused** (D-076: rounds 2–3 gave no testable survivor). Restart only with a new owner brief.
   Optional: the B0218+357 VLA polarisation sign test (R2-B).

## Neutrino Frontier (owner brief 2026-10-10, D-077; docs/neutrino_frontier/)
Independent program; takes no resources from MOA-II, S1, E1/D1 or data validation without an allocation decision.
1. **E-NF1 follow-up (NF-H04 ghosts):** ICECAT-1 null (global p = 0.058), CLs limits R_g ≤ 0.02–0.05 (Δt ≤ 1 h),
   ≤ 0.1–0.2 (1 h–7 d), ≤ 0.5 (7–30 d). **IceTracks-DR2 short lags (E-NF1b, D-079): null**, R ≤ 0.002 per northern
   track with log10 E ≥ 4 (≤ 100 s), ≤ 0.01 (100 s–1 h). The DR2 forecast rules out R_g ~ 10⁻³ at every lag and any
   gain beyond 1 h (optimistic floor ≥ 0.055 at 1 h–1 d, ≥ 0.64 at 30–180 d). Next, only if NF-H04 gets a
   quantitative prediction: f_astro per cut from the DR2 effective areas (turns R into R_g per astrophysical
   neutrino); the 30–180 d bin has no public sample that can reach it.
2. **Benchmarks:** TXS 0506+056 2014–15 box with SkyLLH **reproduced** (D-080: n̂s 12.72, γ̂ 2.26 vs 12.7 / 2.3,
   once the box edges sit on the edge events). Next: the time-integrated benchmark (arXiv:2605.19040 Table 8: TXS
   8.8 / 2.0 and NGC 1068 80.1 / 3.2 with SkyLLH, all 14 seasons; needs the IC40/IC59/IC79 IRFs too: 3.4 GB of
   smearing matrices, all on disk at once; cloud only, delete after use).
3. **Round 2** (N-A then N-B): only frameworks with a prediction decidable on IceTracks-DR2, GCN alerts after
   ICECAT-1 (77 event times; clean revisions and retractions first), HESE-12 topology labels or KM3NeT/ANTARES.
4. GW–neutrino ghost corner (N-B: missed by LVK lensing searches): GW sky maps now exist (D-074 addendum 4).
   Define a GW-parent ghost test only if round 2 gives it a quantitative prediction.

## Part 0 speed-ups (owner brief 2026-10-10, D-078)
- **W3 short-event study (D-082):** stage 1 done (frozen t_E = 3 d sample; `eclipse_dip` losses are unsampled
  spikes; the `jackknife_nights` share is undecided). Stage 2: re-vet the dev flags recording the jackknife's own
  feature nights and ΔBIC, then a night-robust likelihood test on dev, an FP check on the real gb12 flags and PSPL
  controls, and one validation evaluation.
- **Injection speed (V2.1 P2):** compiled W3 finite-source integrand (46 % of fit time) with an identical-outcome
  test on fixed injections; then 2-column closed-form solves in place of `lstsq`.
- **Acceleration V2 audit** (`docs/audits/2026-10-10-acceleration-v2.md`): MOA pre-screen reuse done (D-078
  addendum 3). Next, by information per hour: the shared-epoch masking study (dev/validation injection split,
  real flags as the false-positive check), then profile `run_inject` (≈ 60 % of a field run), then the
  short-t_E jackknife study. Same fix for fit chunks (`fit_chunk_problem` still keys on the whole `Params`;
  ~200–290 s of refits per field after a vetting-only change).
- Done: `derived-data-20261010` holds the GW sky maps; a fresh session takes them in ~3 s (D-078 addendum).
- **E1 GBM re-pinned (cloud, D-078 addendum 2):** one snapshot (8cc2…, 4,392 bursts) in both E1 manifests; every
  E1 conclusion unchanged. **Next (local session, soon):** fetch 8cc2… while HEASARC still serves it, publish it to
  the derived-data store and wire `e1_events.fetch` to it (else `--refresh` and publish that). Later: per-catalogue
  RNG streams in the E1 nulls (own PR).
- Profile every null/injection script once before runs > 10 min (S-3); next candidates: `e1_events.py` injection
  loops, `GWMaps` construction (region distances, ~27 s).

## Continuing (owner direction 2026-10-08: evidence of wormholes / negative mass / warp in any public dataset; D-054)
1. **W3 in published microlensing samples** (D-057–D-059, D-061). Mróz 2019 bulge (5,790) and 2020 disk (460) fitted and
   vetted: 127 + 6 flags, **0 survive** (the D-059 chunk survivors BLG667.04.62161 and BLG624.18.69573 fail
   `feature_coverage`). **The published samples cannot limit W3:** 0 of 600 injected W3 events pass the emulated Mróz
   selection (PSPL controls 15–43 %), so the chunk re-fits only re-test a selection that excludes the signal. Next:
   light curves from **before** a PSPL selection. **MOA-II 9-year Cut-0 light curves keep W3 (D-062):** gb22
   null (30 flags, 0 survive). **D-068: eleven fields streamed (gb7, gb11, gb12, gb15–gb22; 712,780 light curves), 582
   flags, 0 candidates (gb17, gb18: one chain survivor each, explained outside the chain); limits withdrawn**
   (injections had a lenient source-flux bound). The four D-068 chain gaps are closed (D-068 addendum
   2026-10-10, `CHAIN_VERSION` 2026-10-10.1). All eleven fields re-run under it (gb18 last): 0 survivors, limits tracked (CHANGELOG 2026-10-10; gb18's t_E = 300 d,
   ρ = 0.1 cell has 0 recoveries). **Combined (21 fields, gb14 added): Γ₉₅ ≈ 1.2 × 10⁻⁸–7.4 × 10⁻⁸ per star per year**
   (`results/w3_moa/limits_combined.ecsv`; no field > 19 % of any cell; gb14 tightened it by 5–10 %). The lenient-reference fields were re-vetted with a working CDS XMatch: 0 survivors (CHANGELOG 2026-10-10).
   Follow-ups (from the gb10 review): the pre-screen per-chunk clock starts at a chunk's first finished segment, so
   chunk windows overlap and `wall_time_s` / `mb_per_s` are overstated (use run totals); ~~manifest rows and
   retrieval dates~~ done (2026-10-17: all 22 streamed tars listed; dates "not recorded" until a re-stream stamps them).
   A fresh clone needs `merge-prescreen`, `merge-chunks` and a re-run of `vet` (~3 min) before `limit`: the tracked
   `results/w3_moa/vetting_gbF.json` is a summary without the chain stamp `limit` checks (gb14 run 5).
   MOA-II is complete except gb13; next: gb13 after the owner decides #145. Each further field now moves the combined
   limit by ≲ 10 % (T4: diminishing; prefer a new survey for W3 coverage). Gaia DR3 `vari_microlensing` is ruled out (D-061: 0 / 240 W3
   injections both selected and flagged; one real flag, not a candidate). Next: Gaia DR3 variables with epoch
   photometry (vari_summary; check by injection whether W3 survives the variability classifier first), KMTNet
   public seasons, OGLE EWS (owner decision, terms). Optional:
   re-fit the bulge under D-058's bounded π_E (`merge-chunks`, then `vet` + `revet`; can only add flags); re-run the
   6 arXiv name queries that errored.
   **W3 in ZTF DR24 light curves** (v13 rule; v22 handoff: 17 qualifying partitions left). Fade-control diagnostic
   2026-10-20 (`results/w3_ztf/fade_diag/RESULT.md`): per-area rates are heterogeneous (χ² 234 / 24) and no inventory
   metric predicts them; post hoc they cluster by HEALPix order-4 parent (between 208 / 11, within 26 / 13 dof).
   v25 (2026-10-21, `results/w3_ztf/v25/RESULT.md`): the frozen sibling prediction holds in 6 / 6 areas, so H-sib
   is the queue rule (untested parents and fade-passing parents' siblings first; siblings of failing parents are
   not fetched). 2346 and 2347 searched (3,926 groups), 0 survivors. v26 (`results/w3_ztf/v26/RESULT.md`): both
   remaining untested parents fail at their heads (2878 fade 11.7 %, 2813 14.3 %), so four siblings were skipped;
   **the v22 inventory is exhausted** (0 searched). Next: a new-area inventory (the v22 round ran under the owner's
   4 h / ≤ 6 GB authorization after #395; it does not cover another: v22, v25 and v26 streamed ≈ 4.1 + 4.5 + 1.1 GB
   per their CHANGELOG cost lines, so a new inventory waits for a new owner authorization). The fade limit stays an
   owner decision.
2. **W1/W2 in wide imaging** (D-056: published lens catalogues null; f_dark < 0.31 (typical, 3/25 after the 2026-10-08
   pair-match and one-entry-per-lens amendment) for quasar/radio-selected lenses, the only selections sensitive to a
   dark lens). Next: HSC PDR3 photometry, or a *validated* PSF-subtracted HST
   image model for the 307 blended or too-close lensed quasars (HST *catalogue* photometry is not decisive: D-060,
   efficiency 0.46). The Moffat N-image chain (2026-10-21, D-060 addendum, `results/w2_hst_blend/RESULT.md`) recovers
   9/10 known lenses but its no-lens negative control gives S/N 13 (≥ 5σ): method validation fails, no target is
   interpreted. It needs a per-image PSF matched to the HAP-MVM drizzle (STScI empirical PSF grids, photutils) or
   lenstronomy lens-light + point-source fits (new dependencies), judged by no-lens negative controls < 5σ and controls
   ≥ 0.9; the CHITAH lens models of the 3 open SuGOHI IX pairs; rejected lensed-quasar pairs: LS DR10 cannot decide them
   (D-064, control efficiency 0/5 at 1.9–2.6″); of the 11 colour-matched ones, J0130+0725 has no lens light in HST to
   F814W ≈ 23 (D-064 addendum: binary vs dark lens needs spectra of both images or two-epoch flux ratios), J0728+2607
   has none either with an empirical PSF (F814W ≈ 23), and 9 have no HST (Euclid DR1 when public); catalogued image separations for the 15
   pair-decided systems the pair check cannot test yet (incl. the 3 CHITAH pairs; HSC lens models); HST lens photometry (CASTLES-type) for the blended quasar systems; ALMA positions for submm systems;
   deeper imaging (HSC PDR, Euclid DR1) for the conservative variant; a lens list with image positions for W1.
3. **W5** (D-063): Legacy Surveys DR10, 340.5 deg², 40 flags, 0 survivors; n₉₅ ≈ 0.012–0.018 deg⁻² at θ_E = 8–32′.
   Euclid Q1 counts would not open θ_E < 6′ (D-065: clustering-limited, ×1.0–1.4 S/N on a fifth of the area).
   **Euclid Q1 radial-shear screen (D-066, D-067):** all three Deep Fields (60 deg²) null; R = 0.56 ± 0.16 from SZ
   clusters (R = 0.5 kept); n₉₅ ≈ 0.049 / 0.051 deg⁻² at θ_E = 2′ / 4′; 1′ not limited (efficiency 0.43–0.58). Next:
   per-tile star-ellipticity gradient test; Euclid DR1 when public; θ_E ≤ 1′ and ≈ 1° are still open.
4. **Warp:** monitor the literature (checked 2026-10-08: no distant-observer EM template — Lentz & Felton 2024 give fluxes only; Clough et al.
   2024 waveform not public; O3 superluminal-burst search already null, Kuwahara & Cannon 2023). Recheck monthly.

## JWST focus (D-047 screens; continues under the direction above)
1. **Injection-recovery, then limits** (Phase 2; `exotic_sim.inject_images` / `inject_light_curve`):
   - W1: the shear screen (D-053, `exotic_screens.py shear`, `inject_shear.py`) is built; thirteen fields null (Abell S1063,
     SMACS 0723 and El Gordo added 2026-10-17; MACS0417, WHL0137, RXJ2129 and MACS0647 too, MACS1423 explained hit, back in under the streak veto; MACS0717 added 2026-10-17 with ASTRODEEP-FF photo-z, D-053 rule not passed, nearly blind (efficiency ≤ 0.045), kept out of the headline: 83 % of its rows lack photo-z (not stored), `results/w1_shear_macs0717/RESULT.md`); twelve-field limits with the streak veto 4.4 × 10³ / 7.2 × 10² / 3.5 × 10² deg⁻² at 2 × 10¹² / 8 × 10¹² /
     2 × 10¹³ M☉, still blind at ≤ 2 × 10¹¹ M☉. Next:
     lower the B-mode floor (Abell 2744 max S_× 4.21 at (88, −42), full coverage: ~~PSF anisotropy from stars~~ done
     2026-10-17, not adopted (`results/w1_psf/RESULT.md`: 3 stars, too few for a map; a constant χ_PSF lowers the
     floor to 3.85, rule ≤ 3.79); ~~model PSF per detector (STPSF)~~ done 2026-10-17, not adopted
     (`results/w1_stpsf/RESULT.md`: axis close to the stars, |χ| depends on the isophote level and over-predicts the
     stars by 0.05–0.07); ~~smoothed STPSF (detection kernel, jitter, drizzle box)~~ done 2026-10-17, rule S3
     fails (`results/w1_psf_smooth/RESULT.md`; PSF-model route stops; next PSF-weighted pixel moments or new fields), blend
     rejection (2026-10-17 diagnostic, `results/w1_floor/RESULT.md`: no streak there; 44 sources near bright
     galaxies carry it, above the rotation null; a bright-neighbour veto is not floor-specific); ~~drop edge apertures with < 50 % coverage~~ done 2026-10-17, not adopted: limits move ≤ 1.6 %,
     `results/w1_edge_cut/RESULT.md`); ~~add SMACS 0723 and El Gordo~~ done 2026-10-17, null
     (`results/w1_shear_smacs_elgordo/RESULT.md`); ~~MACS0417, MACS1423, WHL0137~~ done 2026-10-17 (`results/w1_shear_more/RESULT.md`: MACS1423 hit carried by sources on an image stripe; nine-field limits 4.9 × 10³ / 8.2 × 10² / 4.1 × 10² deg⁻²); ~~RXJ2129 and MACS0647~~ done 2026-10-17, null (`results/w1_shear_rxj_m0647/RESULT.md`; MACS0647 R 0.23, nearly blind, top centres on streaks); ~~MACS0717 with photo-z~~ done 2026-10-17 (ASTRODEEP-FF, VizieR; rule not passed, R 0.23, ≤ 9 / 200 recoveries; not in the headline; the low R is unexplained); ~~a calibrated streak veto~~ done 2026-10-17, adopted (`results/w1_streak_veto/RESULT.md`: MACS1423 back in, no field passes; twelve-field limits 4.4 × 10³ / 7.2 × 10² / 3.5 × 10² deg⁻²; next: faint parallel streaks such as the S1063 peak's, fewer false trails on MACS0647); stack S around `radial` (done 2026-10-17 for four fields, `results/w1_stack_radial/RESULT.md`: excess mostly from the defining arcs; residual p 0.027, exploratory; MACS0417 second model: S unchanged on the 7 common centres; injections at the centres done 2026-10-17, `results/w1_stack_radial/inject/RESULT.md`: post hoc noise-aware rule excludes ≥ 1 × 10¹² M☉ at all centres (mass-point addendum), ≥ 8 × 10¹² at half; ≤ 5 × 10¹¹ not excluded; line done) and orphan-pair
     centres; ~~other lens redshifts → volume density vs Takahashi & Asada~~ done 2026-10-17 (`results/w1_volume/RESULT.md`: eleven fields, Z 0.1–1.5, n₉₅ < 3.65 × 10⁻⁴ h³ Mpc⁻³ at 2 × 10¹³ M☉, 3.6× weaker than T&A; 8.5× / 65× at 8 × 10¹² / 2 × 10¹²; single-plane ASSUMPTION; RXJ2129 out, catalogue 404); ~~blank CANUCS flanking fields~~ done 2026-10-18 (D-098, `results/w1_shear_ncf/RESULT.md`: no lens model, F150W, R 0.425 from gate R; fourteen fields n₉₅ < 2.35 × 10⁻⁴ h³ Mpc⁻³ at 2 × 10¹³ M☉, 2.35× T&A; macs0416-ncf / macs0417-ncf flags on a missed diagonal streak and horizontal stripes, kept out). Next W1 area: ~~a lower findsat_mrt threshold for the NCF streaks~~ failed 2026-10-18 (`results/w1_streak_rowdiag/RESULT.md`: masks identical at 5 / 3 / 2.5; macs0416-ncf is a spike case, macs0417-ncf detector row striping); ~~a field-wide row-segment median veto~~ failed 2026-10-18 on cost (`results/w1_row_veto/RESULT.md`: covers the macs0417-ncf peak but keeps 56–75 % of the rows; a retry needs source footprints masked first), ~~GOODS-N (JADES 1181)~~ done 2026-10-18 (D-099, `results/w1_shear_goodsn/RESULT.md`: four of six pointings pass gate P; o003 / o004 flags on mosaic-edge artefacts, explained by V1 ∪ V2 and kept out; sixteen fields n₉₅ < 1.82 × 10⁻⁴ h³ Mpc⁻³ at 2 × 10¹³ M☉, 1.82× T&A; the gain is R-dependent: R 0.231 keeps 0.21–0.26 of the recoveries in the four-pointing totals). ~~CEERS (program 1345)~~ done 2026-10-18 (D-099 addendum, `results/w1_shear_ceers/RESULT.md`: ten pointings pass gate P; o004 / o069 flags on gap-edge and spike artefacts, explained by V1 ∪ V2 and kept out; twenty-four fields n₉₅ < 1.09 × 10⁻⁴ h³ Mpc⁻³ at 2 × 10¹³ M☉, 1.09× T&A; R 0.231 keeps 0.19–0.28 of the CEERS recoveries). ~~a blank-field R calibration~~ done 2026-10-18 (D-101: borrowed R conservative). ~~PRIMER (program 1837)~~ done 2026-10-18 (D-099 addendum 2, `results/w1_shear_primer/RESULT.md`: twenty 21.8 arcmin² tiles pass gate P; four E/B flags, three explained by V1 ∪ V2, UDS x3y0 explained post hoc by its `_cal` exposures as an optical artefact of a G 10.07 star 58″ away that the 60″ external spike cap misses; forty fields n₉₅ < 6.38 × 10⁻⁵ h³ Mpc⁻³ at 2 × 10¹³ M☉, 0.638× T&A under R = 0.414; R 0.231 keeps 0.15–0.28 of the PRIMER recoveries). Next W1: ~~lower masses: filter pilot~~ done 2026-10-18 (D-103, `results/w1_low_mass/RESULT.md`: the frozen Schirmer 10″ filter beats point-mass and 5″ variants; V0 recovers 3–4 / 100 at 1 × 10¹² M☉, so a forty-field 1 × 10¹² grid would give ≈ 30–45× Takahashi & Asada, low priority), ~~an R measurement on real blank-field images~~ done 2026-10-18 (D-105, `results/w1_r_realmetacal/RESULT.md`: metacal of GOODS-N F200W tiles R = 0.953 ± 0.24, frozen rule "conservative"; too uncertain to adopt R ≈ 1; next for it: an empirical drizzled PSF and a bad-pixel mask before deconvolution, more tiles and fields); ~~v2: mask, PSF σ floor, unique matching, all tiles~~ done 2026-10-18 (D-105 addendum, `results/w1_r_realmetacal/v2/RESULT.md`: 18 tiles, R = 1.056 ± 0.156, R − 2σ 0.744 < 0.8, still "conservative"; PSF_d × 1.15 moves R by −0.18, so next is the empirical drizzled PSF, then a second field); ~~v3: empirical PSF_d from Gaia stars~~ failed 2026-10-18 (D-105 addendum 2, `results/w1_r_realmetacal/v3/RESULT.md`: gate E, 2 usable stars; post hoc ePSF R 0.707 ± 0.169 vs analytic 1.006, simulated gate fails; R = 0.414 kept; route paused, reopen only with a validated drizzled STPSF or a star-rich field);
   - W2 / dark deflectors (D-048 clusters, D-051 deep fields: null, limits; D-055: the flanking-field orphan
     excess was a null (e) artefact, no excess under the fixed (e) or the companion-aware (f); D-048 clusters re-run, null): segmentation-map
     adjacency for the same-galaxy and visible-lens rules (the main efficiency losses); pixel-level injections;
     optionally CEERS / GOODS-S / PRIMER (~600 MB each);
   - W3 (D-052, null; headline limits from MACS0416 + Abell 2744): ~~Abell 2744 `forced` with 300 controls~~ done
     2026-10-19 (calibrated, null, headline 5–11 % tighter; `results/dimming/abell2744_calib/RESULT.md`); ~~NEXUS
     `forced` with 300 controls~~ done 2026-10-19 (calibrated; 5 artefact-test survivors, none a W3 dimming by eye,
     post hoc; headline ×2.6–6.2 tighter if they hold; `results/dimming/nexus_calib/RESULT.md`); next `/vet-candidate`
     on its survivors (d79956 SN-like; three single-band points) done (#371, 0 W3 events); ~~MACS0416 `forced
     --workers` and a joint `combine`~~ done 2026-10-19 (`results/dimming/macs0416_regen/RESULT.md`); ~~fold
     `chromatic_veto` into `forced`~~ done 2026-10-19 (D-052 addendum e; injections unaffected, no re-run needed); replace the single-epoch veto for
     vanish flags with the D-039 persistence test (scoped 2026-10-19, D-052 addendum f,
     `results/dimming/persistence_scope/RESULT.md`: worth it, joint exposure ×2.0–3.4 upper bound, 37 compact flags;
     run 2026-10-19 partial, `results/dimming/persistence_run/RESULT.md`: NEXUS 12 / 17 `on_sky`, controls 61 / 61 `on_sky`; ~~`forced` on the 12~~ done 2026-10-19, `results/dimming/nexus_o014/RESULT.md`: 0 survivors, 11 / 12 at other epochs' mosaic edges, 7 detected in covering epochs; controls 8 / 300 flagged, above the frozen 1 %; all 37 targets tested 2026-10-19: 32 `on_sky` (MACS0416 12, Abell 2744 8 each in one epoch only), 0 `persistence`, controls 194 / 194 at S/N 11–33; next: `forced` + epoch depth / footprint test on the 20, check edge coverage in the injection chain (tried 2026-10-19, `results/dimming/edge_coverage/RESULT.md`: gate V failed as written, mis-specified gate set, no decision; post hoc the proxy counts 9–18 % of injection pairs at zero / low weight; re-gated 2026-10-19: gate V2 on per-epoch stamp flags passed, decision `over_counts`, veto-off ×0.81–0.93 masked; ~~a weight-map coverage term, re-run `inject`, corrected headline~~ done 2026-10-19 (D-052 addendum h, `results/dimming/wht_coverage/RESULT.md`: per-source limits ×1.07–1.33 weaker at t_E ≥ 0.03 yr; swap gain ×1.80–3.52 with the term)), ~~the persistence test's `on_sky` efficiency~~ done 2026-10-20 by exposure count, held out on NEXUS e5–e7 (D-052 addendum j, `results/dimming/expcount_eff/RESULT.md`: ≤ 4 exposures 0.412 holds, > 4 untested); ~~`forced` + depth / footprint test on the 20 MACS0416 / Abell 2744 targets~~ done 2026-10-20 (`results/dimming/forced20/RESULT.md`: 0 survivors, 16 constant-flux catalogue misses, 4 weight edges; all 32 `on_sky` targets explained), ~~fold the swap into `inject` with that efficiency~~ done 2026-10-20 (D-052 addendum k, `results/dimming/swap_inject/RESULT.md`: joint exposure ×1.54–2.51 at ρ = 0.1, new headline); ~~the NEXUS `undetected` control stamps~~ done 2026-10-20 (`results/dimming/undetected_stamps/RESULT.md`: visible sources lost to 1–2 DO_NOT_USE aperture pixels or the edge; measurable held-out ≤ 4 0.59–0.84, decision R2); ~~NaN-tolerant aperture~~ tried 2026-10-20 (`results/dimming/nan_aperture/RESULT.md`: G1 fails at 2 pixels, k stays 0; post hoc k = 1 lifts held-out ≤ 4 to 58 / 89); ~~k = 1 frozen plan on unseen NEXUS e0 / e2 controls~~ done 2026-10-20 (`results/dimming/nan_k1/RESULT.md`: R3, only 8 unseen ≤ 4 pairs, all `on_sky`; e2 > 4 class 143 / 143 holds); next: k = 1 waits for a new wide epoch as held-out data); ~~split the held-out `inconclusive` pairs by cause~~ done 2026-10-20 (`results/dimming/inconclusive_split/RESULT.md`: 12 / 12 lack an earlier exposure, 11 a single first dither; optional next a cross-visit prior in D-039); add JADES and new NEXUS epochs;
     SN/TNS check for any survivor; leave the Sunrise transient track to its owner run;
   - ~~W5: counts N(>S) around `radial` centres~~ done 2026-10-17 (`results/w5_radial_counts/RESULT.md`): no deficit,
     but no power: at JWST depth the model predicts an excess; the post hoc excess (Q ≈ 1.26 at 9″ without arc
     neighbours) is flat in R. An excess test needs a density-matched control (low priority).
2. ~~W3 inside caustic-crossing arcs needs a microlens with macro shear~~ done 2026-10-20 (simulated, `results/w3_macro_shear/RESULT.md`: minimum images keep a zero-flux umbra narrowed to 2√λ1 along the stretch; saddle images have none). Next only with a W3 search on caustic-transient stars: inject through `macro_magnification_map` with each image's (κ, γ) and parity.
3. Warp: recheck only when a paper gives an electromagnetic prediction for a distant observer.

## Supporting (M3 lensing-violation search, D-023; serves the focus above)
1. **Lens-model consistency on SMACS** (D-024: model validated; arc orientations give a null result). Next:
   - counter-images: done for the catalogued systems (D-029: 0 of 11 testable uncatalogued images absent). Left:
     - ~~bright single arcs~~ — done for five fields (D-029 addendum 2026-10-14, `arcs_counter.py`): null; Abell 370
       3881 was CATS system 7 with a wrong photo-z; controls on all fields;
     - ~~a BCG/ICL-subtracted test for system 17~~ — done (D-029 addendum 2026-10-14): `absent` once the halo is
       removed, explained by its fitted redshift (the pair is near a caustic; model-dependent); not a candidate;
     - ~~reference fluxes for systems 11, 16 and 26~~ — done (D-029 addendum 2026-10-14, `smacs_ref_fluxes.py`):
       counter-images ≤ 1–4σ in three models, except Zitrin's system 16 third image (18σ, recovered); the 10σ
       ICLv2 prediction for 26 rests on a neighbour-contaminated reference; no `absent`, not candidates;
     - ~~`--halo-filter-arcsec` on the other fields~~ — done (D-029 addendum 2026-10-14): Abell 2744, MACS1149,
       MACS0416 give no new `absent`; member-buried rows need member subtraction (isophote/Sérsic), not a median filter;
     - ~~member subtraction~~ — done (D-029 addendum 2026-10-14, `member_sub.py`): 1 of 31 rows testable (MACS1149 2,
       `absent`, already explained); unconvolved Sérsic residuals leave the rest untestable;
   - flux ratios and colours: done for SMACS and El Gordo (D-032: no anomaly; most images untested). Left:
     - ~~Abell 2744 `fluxratios`~~ — done (D-032 addendum 2026-10-13): 34 + 35 images tested, 7 flags, all ordinary;
       an ICL-aware total (local background or ICL-subtracted) would remove the halo-light flags;
     - totals for arcs inside cluster-galaxy/ICL light: `scripts/halo_phot.py` explains Abell 2744 602.1c (600.1 is
       a resolved spiral, untestable; D-032 addendum; `_phot_apcorr` rejected). Run on all three fields (2026-10-13: no anomaly; scatter 0.10–0.27 mag on the
       compact subset; El Gordo 1c is a thin arc that passed the concentration gate). ~~Elongation gate, re-run,
       `fluxratios --halo-phot --mu-samples 100`~~ done 2026-10-13 (D-032 addenda; `results/halo_phot/abell2744_all/`:
       one flag, A2744-QSO1, ordinary);
     - ~~`bayes.dat` μ errors~~ — done (D-032 addendum 2026-10-13, `--mu-samples`): model error 0.04–0.07 mag
       median against a 0.32–0.36 mag scatter floor; no new flag. SMACS 6.3 needs a chain (Mahler+22 pins none);
     - ~~MACS0416 halo fluxes~~ — done (D-032 addendum 2026-10-14, `macs0416-canucs`): 5 systems testable, no
       anomaly; the halo-flux line is null in four clusters (MACS1149 needs a Lenstool par); ~~Abell 370~~ done
       2026-10-18 (`results/halo_phot/abell370/RESULT.md`: per-system rms gate opens 27 / 45 systems of
       `abell370-canucs`; 2 images testable, no anomaly, almost no power);
     - observational parity test: `scripts/parity_test.py` (affine fit of image pairs in both handedness families,
       synthetic power per pair; D-032 addenda 2026-10-14): acceptance test failed (87 % vs 95 %; 3 persistent
       flags, all predicted-mirror compact cores; classes depend on the start set). PSF-commuting fit with an
       empirical F200W PSF and a start-stability rule done: worse (76 %, 10 flags), so the PSF is not the cause.
       **Paused** (three acceptance failures). Reopen only with a forward model that is not affine over the fit
       disc (source-plane reconstruction through the published deflection field) or a noise-correlation-calibrated
       Δχ² threshold; no field run before acceptance passes;
   - ~~`arcs` on Abell 2744 and El Gordo~~ — done (D-024 addendum 2026-10-13): null; background 33/35 and 10/12
     aligned, anti counts 0 and 1 against 0.45 and 0.17 from shape noise;
   - ~~`arcs` on MACS0417~~ — done (2026-10-16, `results/arcs/macs0417-relics/RESULT.md`): RELICS v1 fails the
     control (64 % aligned); CANUCS v1 maps 37 / 49 aligned, anti 3 vs 1.52 (P 0.195): the S1063 anti excess does
     not replicate. ~~MACS0647~~ — done (2026-10-16, `results/arcs/macs0647-clash/RESULT.md`): both CLASH Zitrin
     models fail the control (67 %, 50 %), no result. ~~MACS1423~~ — done (2026-10-16,
     `results/arcs/macs1423-canucs/`): CANUCS DR1 maps, screen anti 4 vs 0.95 (P 0.016), but the clump-pair census
     gives P 0.259 (exploratory): not replicated under the clump-aware null. Re-run with the fixed #273 census code
     (2026-10-17, #287): every number unchanged; only 1818 (gate-failed) moves `single` → `undecided`. #279 has its
     own re-run (open). ~~RXJ2129~~ — done (2026-10-16,
     `results/arcs/rxj2129/RESULT.md`): Lenstool and GLAFIC JWST-updated maps both fail the control (8 / 12),
     no result. ~~RXJ2129 on the deeper DJA catalogue~~ done 2026-10-16 (`results/arcs/rxj2129-dja/RESULT.md`: no result, controls fail);
     ~~WHL0137~~ — done (2026-10-17, `results/arcs/whl0137-relics/RESULT.md`): RELICS Lenstool 13 / 18 aligned,
     anti 1 vs 0.15 (P 0.14): does not replicate; ~~RXC J0600~~ — done (2026-10-18,
     `results/arcs/rxc0600-relics/RESULT.md`): RELICS GLAFIC v3 11 / 12 aligned, anti 1 vs 0.20 (P 0.18): does not
     replicate (Lenstool v1 fails the control); ~~SPT0615~~ — done (2026-10-18, `results/arcs/spt0615-relics/RESULT.md`;
     RELICS BPZ photo-z): both RELICS models fail the control (3 / 7, 1 / 4), no result; ~~MACS0308~~ — done
     (2026-10-18, `results/arcs/macs0308-relics/RESULT.md`; RELICS BPZ): GLAFIC v2 17 / 21 aligned, anti 1 vs 0.31
     (P 0.27): does not replicate (Zitrin 15 / 18, same anti); ~~RXC J0018~~ — done (2026-10-20,
     `results/arcs/rxc0018-relics/RESULT.md`; DJA `rj0018p1626`): RELICS Lenstool 3 / 4 aligned, anti 0 vs 0.095: no power
     ("replicates" unreachable, best P 0.091; 21 / 30 strong-shear rows lack a DJA photo-z). ~~Pooled tally~~ done 2026-10-20 (exploratory,
     `results/arcs/pooled-tally/RESULT.md`: six fields 10 anti vs 3.22, P 0.0019; without MACS1423 6 vs 2.27, P 0.029;
     P 0.05 at σ_ε 0.273 / 0.320; S1063 vs a common pooled rate, conditional, P 0.078–0.121, post hoc). ~~Calibrate on
     weak-shear rows~~ done 2026-10-20 (`results/arcs/null-calibration/RESULT.md`: 9 anti vs 11.06 / 9.14, no excess;
     σ̂_ε 0.24 / 0.27, f̂ 0 / 0.06 with wide intervals; frozen label "explains" rests on f̂ only, in substance no power;
     no further calibration from these tables). Next, optional: VENUS (6882) clusters with
     RELICS/CLASH maps (no photo-z yet); Bullet (LAM atlas Lenstool model, JWST 4598) once a photo-z source exists;
   - ~~arc curvature~~ — done (D-024 addenda 5–6): forward-modelled curvature, 5 testable arcs all as predicted (El
     Gordo 2527 explained by its model); critical curves at z_s = 1, 2, 4 still open;
   - shapes for the `arcs` test: S/N ≥ 50 done on the stored tables (`docs/fields/abell370.md`: 3 anti vs 1.33,
     P = 0.15, all intrinsic shapes); F444W vs F200W done (D-024 addendum 4: 81/83 aligned
     stay aligned, F444W anti 3 vs 1.30, same axes, no blend); left: blend rejection in F200W segmentation;
   - exotic screens (D-031, D-034, D-036–D-038, D-040–D-043): SMACS, El Gordo, Sunrise (radial), Abell 2744,
     MACS1149, MACS0717, MACS0416, Abell 370 and Abell S1063 (radial) null. All six HFF clusters are screened. Next:
     - multi-plane checks: `--plane NAME=Z[:SIGMA]` on `validate` / `images` done (D-046 addendum 2026-10-14; MACS0416
       8757 at z 0.268 costs Δχ² ≤ +5.6, no tension). MACS0416 CANUCS off-plane potentials done
       (`offplane_potentials.py`, 2026-10-14: 12 / 222, null; O5 at z 0.114 costs Δχ² +39 at its best σ 150 km/s, unrefitted, so an upper bound on the refitted cost). The same listing is
       done for Abell 2744, El Gordo, Abell 370 and SMACS 0723 (D-046 addenda c–d): null in all five models. Next only if a
       flag comes to depend on an off-plane potential: a local refit (no Lenstool optimiser in the repo yet);
     - orphan pairs (D-048, null): injection-recovery of simulated dark-lens pairs (limit per deg²); CANUCS segmentation
       maps for the same-galaxy test; blending model for pairs < 1″; ~~re-run on Abell 2744~~ — done (D-048 addendum
       2026-10-16, `results/orphan_pairs/abell2744/`): null (P 0.16 / 0.11 under (e) / (f)); efficiencies as in the
       deep fields; no 0.3–3″ positive control; per-pair vetting done (CHANGELOG 2026-10-16,
       `results/orphan_pairs/abell2744/vetting/RESULT.md`): 4 ordinary (1 visual), 7 insufficient (#2, #3, #7, #9, #10,
       #11, #13: E2 not validatable, no public spectra), 0 anomalies; deciding data: spectra of both members (a
       redshift difference would argue against one lensed source; equal redshifts alone would not confirm lensing). Not
       pursued now: a resource decision (owner, 2026-10-16), not an exclusion. Left:
       MACS0417, MACS1423 cluster pointings (published lists now transcribed: Jauzac+19, KALEIDOSCOPE; see the
       image-position item below);
       MACS1149 image list when CANUCS v2 is out;
     - CANUCS models pinned (D-044: `macs0416-canucs`, `abell370-canucs`); next, the 100 MCMC sample maps for μ spreads;
     - `radial`: low-weight veto (relative WHT < 0.5) or aper50 S/N floor; pass `--spike-stars` (Gaia, D-043) on every
       field and re-run earlier fields' radial screens with it (~~the five W1 shear-only fields~~ done 2026-10-17,
       `results/radial_more/RESULT.md`: no flag; RXJ2129/MACS0647 no result; MACS0647 top centre is streak segments;
       ~~W1 streak masks as a veto~~ tried 2026-10-17, not adopted: the masks miss the faint MACS0647 streak, and
       a partial veto thins the null so the centre drops to p 0.03 (`results/radial_streak/RESULT.md`; the twelve-field
       diagnostic changes no flag in the control-passing fields). A streak veto for `radial` needs a faint-streak
       detector checked on that streak, plus its own injections; ~~centre-level axial concentration~~ tried
       2026-10-17, gate K failed (A 0.73 < 0.9, `results/radial_axial/RESULT.md`); next a collinearity rule with a
       new known streak as its gate);
     - ~~re-run MACS0717 radial when a photo-z catalogue exists~~ done 2026-10-17 with ASTRODEEP-FF photo-z (no flag, max 5 lines,
       p 0.69; `results/radial_macs0717/RESULT.md`; photo-z cover ~17 % of rows);
     - `forced_check`: match and search radii scaled to the image-plane rms when it exceeds 1″; classify model copies of
       unpredicted catalogued images automatically (D-041);
     - ~~MACS0417 / MACS1423 published images vs CANUCS v1 maps~~ done 2026-10-21
       (`results/lens_consistency/macs0417_images/RESULT.md`): gate K fails for every model (CANUCS rms 1.12″ / 1.59″, one
       central image each: MACS0417 1b.4, MACS1423 1.5; RELICS v1 9.3″, maps load correctly); post hoc without them
       0.38″ / 0.37″, every other image reproduced, no flag. Next, optional: forced photometry at the CANUCS
       `extra_predicted` positions; re-run when CANUCS v2 (image catalogues) is public;
     - second models as pinned `MapLensModel` entries: Sharon v4cor (MACS1149, done 2026-10-20); ~~MACS0717 Sharon v4cor and Keeton v4~~ done
       2026-10-21 (`results/lens_consistency/macs0717_models/RESULT.md`: Keeton pinned, both gates pass; Sharon G2 fails as
       written, 6.2 off the map, in-footprint use only). ~~the D-041 cross-model readings of the MACS0717 forced rows
       with them~~ done 2026-10-21 (`results/lens_consistency/macs0717_crossmodel/RESULT.md`: gate fails at CATS-fitted z,
       0.55 / 0.59; post hoc on the 10 arcs.txt-z systems 0 to vet; Step 6 stays on scratch readings). Open, low value:
       read each second model's own source z for the 41 z_m_limit systems;
     - ~~MACS1149 `no_reference` rows~~ — done (D-037 addendum 2026-10-14, `m1149_refs.py`): 11 rows, 0 `absent`, none
       unmet in both models;
     - MACS1149: ~~`macs1149-sharon` in `MODELS`~~ done 2026-10-20 (`results/lens_consistency/macs1149_sharon/RESULT.md`: maps gate passes, gold-image gate fails as written on system 9 off the map; use inside the footprint only); ~~D-041 cross-model rule with it on the 16 CATS extra rows~~ done 2026-10-21 (`results/lens_consistency/macs1149_crossmodel/RESULT.md`: gate 38 / 39; one frozen flag, system 16 at its fitted z, D-037's case; tally unchanged); the D-037 rules as code in `forced_check` (open); re-run `images` with CANUCS DR1
       photo-z; ~~fix system 16's z~~ — done (D-037 addendum 2026-10-14: z ≈ 2.5–2.9 from the pair in both models; CANUCS
       5106117 is a possible far image; a spectrum would settle both; ~~published spec-z~~ checked 2026-10-20, none
       (`results/lens_consistency/macs1149_s16_specz/RESULT.md`: Treu+16 and DJA NIRSpec v4.4; 16.2's CANUCS PRISM
       spectrum is ungraded; recheck when DJA grades `1208_5114276`); evaluate the CANUCS
       lens models (see the CANUCS item above);
     - `bayes.dat` position spreads: not needed (D-045: 0.2–0.4″, inside the 1″ radius). UNCOVER v2.0 passes the
       gate (D-041 addendum 2026-10-14). The cross-model rule on every forced row is null (CHANGELOG 2026-10-14);
       sibling reference fluxes for the five both-model `no_reference` extras are done and null
       (`a2744_xmodel_refs.py`, CHANGELOG 2026-10-14). Abell 2744's counter-image line is exhausted with two models;
       revisit only with a third model or deeper data for 33.1 / 62.1.
2. **Cluster fields** (done: #30–#32, D-026). Run the lens-model checks per field:
   - El Gordo: Caminha+2023 multiple images and magnification maps (CDS);
   - Abell 2744: UNCOVER v2.0 maps (validated 2026-10-14, `scripts/a2744_uncover.py`);
   - **Done (D-030):** `validate --model elgordo-caminha23 | abell2744-bergamini23` reproduces Lenstool's χ²
     (82.53 vs 80.22; 146.64 vs 146.60). No image-position anomaly. Next:
     - `images --forced-image`: done for El Gordo (D-034) and Abell 2744 (D-036); 700.1: both models merge a/b at z_spec 2.3151; vetted 2026-10-20 (D-036 addendum c): undetermined, two modelled members nearby (group halo O4 at 4.4″); parked until a local re-fit or spectra;
     - search radii from `bayes.dat` instead of a fixed 1″;
     - ~~Abell 2744 multiplicity residual~~ — done (D-045): model resolution at folds. Next: `bayes.dat` μ errors in
       `fluxratios` (Abell 2744 and El Gordo chains validated; El Gordo's chain uses σ² = a·b, D-045);
     - ~~El Gordo magnification-map check~~ — done (issue #68): `validate` reproduces the z=2 and z=4 maps (median
       ratio 1.00001, parity 100 % at |μ|<10), so El Gordo μ and parities are validated. Its chain is validated too (`forme -10`, D-045).
   - Sunrise: Scofield+2025 image positions done provisionally (no anomaly; deflection unit read as 0.02″ pixels at
     z_ref ≈ 3, D-033 addenda 1–2). RELICS models on the JWST image list done (addendum 3): system 2 fits at
     z ≈ 2.6–3.1 in all four, so its M5-G1 misfit is that model's. ~~Counter-images in the RELICS models~~ — done
     (D-033 addendum 4, `sunrise_counter.py`): 0 flags; the arc's extra images are model-dependent, its counter-image
     undetectable. Reopen with a JWST-era model that fits every knot.
   Field follow-ups (docs/fields/*.md):
   - ~~vet Abell 2744 `5904`, `7987`, `4731`, `264`, `7298`~~ — done (2026-10-14, `a2744_followups.py`): all ordinary
     (two little red dots, a dusty galaxy with clumps, a z_spec-excluded candidate image, a faint galaxy). Optional:
     relabel D-020's "star without host" reason as "point source without host (star or little red dot)";
   - ~~vet El Gordo `2828`, `1601` and Sunrise `1756`, `1907`, `1578`~~ — done (2026-10-14, `field_followups.py`):
     all ordinary (a star's spike, a faint SW-bright galaxy with no `absent` counter-image, two dusty disks, a
     foreground point source with no counter-images). Optional: why DJA v7.5 misses Sunrise 1756;
   - Sunrise `1869` is a two-band detection on a star spike (D-014 spike mask);
   - Sunrise has a rel_weight plateau at 0.45–0.55 (1,909 sources), so the D-011/D-021 0.5 cut splits it on noise;
   - duplicate merged sources (Abell 2744 `3999`/`7179`, 0.25″ apart);
   - extended-veto fallback for sources without DJA r50 (Abell 2744 `7694`);
   - El Gordo module A is a flanking field;
   - reword `acquire._warn_if_reprocessed`.
3. **Two-epoch search** (D-027): null results for SMACS/VENUS and Sunrise o010/o120; Earendel is steady
   (docs/fields/sunrise.md). Third epoch (o052, calibrated with `--controls`): 0 of 57 candidates. Next:
   - `n0022` and `n0150` are detector persistence (D-039). Run `scripts/persistence_check.py` on every future
     single-epoch candidate before vetting. Next: `/vet-candidate n0153` (on sky in all 8 dithers, F150W ×1.9 over
     2.9 yr; host, AGN colours, a pipeline-version check on neighbours);
   - ~~forced photometry on all catalogued sources~~ — done for Sunrise o010/o120 (D-027 addendum 2026-10-15,
     `transient_grid.py`): 0 candidates; −0.75 mag recovered 100 %, −0.44 mag 53 %. Next: the same grid for Sunrise o010/o052 and
     SMACS/VENUS (F150W, F444W), with a spike veto first;
   - ~~MACS0416 grid, validated on Yan et al. 2023~~ — done for PEARLS o211/o212 (D-027 addendum 2026-10-15,
     MACS0416): 0 new; fixes `--recentre-in origin` and mixed-sign hits (SN01 recovered; caustic stars below depth).
     o211/o213 also null (no grid position on the Warhol arc). Arc differencing (`arc_diff.py`, D-027 addendum
     2026-10-15b): o211/o212 recovers 4 of 5 published Ep1→Ep2 arc transients (W1, W2, S1, S2; W3 4.9σ), but thin
     bright arcs still leave a positive ridge. Fixed with pixel kernels (addendum 2026-10-15c, `--kernel-px 4
     --kernel-stamps`): 4/5 Ep1→Ep2 and 2/4 Ep1→Ep3 published transients; unflagged hits = the published ones
     (o212) plus 3 (o213: a rosette wing, the box edge, h108). h108 vetted 2026-10-15: arc-segment matching residual
     (arc coherence in 4 / 5 pairs incl. CANUCS o004; `arc_diff_vet.py`), not a candidate. Arc tiling o211→o212 (28 tiles,
     addendum 2026-10-15f): 0 candidates on the 21 tiles with morph records (persistence, dipoles, core residuals;
     `arc_diff_morph.py`). Morphology veto injection-calibrated (addendum 2026-10-15g): 0.92 of unflagged injected
     hits pass with the strong-lobe rosette rule. Ring and sign_flip verdicts (addendum 2026-10-15h): t−1−4
     46 → 20 passes; all 28 tiles now have morph records. The 9 open passes (6 SW-dominant + t−1−4 h124/h159/h246)
     are vetted: all not candidates (h044 last, by per-dither `_cal` photometry, `cal_dither_phot.py`; CHANGELOG
     2026-10-15). o211→o213 tiles: 28 / 28 done, 27 passes looked at, all ordinary (o211-epoch SW-core class,
     spikes, persistence incl. h002, PSF residuals); h014 closed by per-dither photometry; 0 candidates (CHANGELOG
     2026-10-15). Sunrise o010/o120, o010/o052 and SMACS/VENUS grids re-run under the origin-recentring and
     mixed-sign defaults (D-027 addendum 2026-10-15i): every earlier group but one spike re-found, 5 new rows all
     stars or spikes, SN Atalanta recovered, 0 candidates;
   - ~~Abell 2744 grid~~ — done for DD 2756 o003/o002 (D-027 addendum 2026-10-15d): 22 groups, all ordinary (spikes
     under a 30° roll, persistence, a detector-fixed pixel, edges), 0 candidates. UNCOVER o001 vs 2756 o002 (D-027 addendum
     2026-10-15e): 46 groups, all ordinary (24 spikes, 9 registration dipoles, 3 persistence incl. one in all six bands),
     0 candidates. Next:
     ~~UNCOVER vs ALT F356W single-band~~ — done (D-027 addendum 2026-10-15j): the catalog grid misses PLATE-23a
     (blended into its host), `arc_diff.py` recovers it; 45 groups, 43 ordinary, 2 real changes outside the
     multiple-image region (g2 faded transient, g1 +0.55 mag variable), 0 lensing candidates. Next: an `arc_diff.py`
     tiling of the UNCOVER∩ALT footprint (`--min-bands 1`) to close the host-blend blind spot; vet g1 with later epochs;
     tiling done (CHANGELOG 2026-10-16: 76 survivors, 6 vetted, 0 lensing candidates; h182 closed). Whole-chain
     injections (CHANGELOG 2026-10-16, `arc_diff_inject_chain.py`): efficiency 0.64 at 27.5 mag, 0.19 at 28.0; the
     peer cut drops half the 28-mag stars on smooth light. A local-annulus cut (#265) and that cut plus
     line-coherence and offset-core vetoes (#268) both fail their pre-written sheet rule: not adopted; the tiling
     limit stays ≳ 36 nJy (at 0.64 efficiency on 4 tiles; not a hard completeness limit). PEARLS 1176 o221
     registration (#266, `reg1176/RESULT.md`): a 0.027″ shift against UNCOVER (small affine terms in 3 bands), no
     0.3″ offset; `hit_lightcurve.py --offsets` applies it; h182 stays closed. Registration does not explain the
     on-arc F356W spread. ALT o002 has its own 0.013″ shift (UNCOVER→ALT pairs not re-measured). No further veto
     rules on these four tiles (owner, 2026-10-16); further work goes to new epochs or untested regions;
   - **Sunrise arc `arc_diff`** (D-027 addendum 2026-10-16): o010 → o052 and o010 → o120 (LW only; the SW bands did
     not finish in two attempts, timings not stored) done, 0 candidates on the arc's images (injection recovery
     ≥ 0.83 at 27.5 AB per band; `results/transients/sunrise_arcdiff/RESULT.md`). ~~El Gordo arcs~~ — one-box
     exploratory pilot done (2026-10-16, `results/transients/elgordo_arcdiff/RESULT.md`): PEARLS o241 → VENUS o051,
     F150W + F444W; full-chain brightening efficiency 0.85 / 0.78 / 0.54 at 26.0 / 27.5 / 28.0 AB, fadings not
     calibrated; 4 review entries, explained by one bright source's diffraction spikes (PA_V3 0.84° apart),
     program-to-program scatter or no forced-photometry change; nothing open. Not extended (owner: pilot only);
     **Fading calibration + B1–B3 (2026-10-17, #288, `results/transients/fade_calib/RESULT.md`):** passes on fresh
     held-out seeds after the injection-geometry correction (vanish 26 AB 14 / 14; unchanged 1 / 40 and 4 / 40);
     El Gordo registration fixed (`--reg-tile-arcsec 6`); B1 review 4 → 1 (vetted h024 position), B2 / B3 0 / 0.
     Per-stamp kernel-fit cost: done (D-094, #291; B2 run ×14 screen / ×8–9 injection, same results). Tool choice: keep `arc_diff` (D-095, #293; SFFT B-spline GPU-only, ordinary SFFT somewhat worse on B1/B2). VENUS batch 1 done (#296, `results/transients/batch_venus/RESULT.md`): 20 boxes, 12 calibrated, 5 reviewed, one ordinary transient (PLCK G165 h002, vetted, SN-like; no catalogued host under it, NS_0481 at 3.77″), h003 undetermined. Follow-up #300: h002 identity and type open (TNS not queried, 403; no spectrum); h003 undetermined, no evidence so far beyond matched-control scatter; active follow-up ended (owner, 2026-10-17). Low-cost pending, not blocking: h002 TNS query when reachable (not polled); 8 undeterminable boxes; pairs beyond the cap. VENUS batch 2 done (#304, `results/transients/batch_venus2/RESULT.md`; follow-up #308, `followup/RESULT.md`): 20 boxes, 12 calibrated, 11 of 11 review entries reviewed (0 pending): 3 real changes (ordinary-transient hypotheses, identities open), 2 undetermined, 5 artifacts (incl. MS0451 h011, detector-fixed), 1 no significant change; MACS0417 h002 counter-images testable ≈ 2027-09-14–2027-11-22 (delays re-checked, #309); PLCK G165 2 per-stamp fallback is slow (stamp pixels > max_pix). VENUS batch 3 done (#311, `results/transients/batch_venus3/RESULT.md`): 20 boxes, 10 calibrated, 9 reviewed (0 pending); A370 box 2 h003/h010 rule-level changes unconfirmed, next checks named (matched controls); held calibration-failed boxes now 26. Confirmation round 1 (#323, `results/transients/confirm_r1/DECISION.md`): A370 h003 common flux offset, h010 not supported; MACS0417 h002 = SN Athena (other team; HST 18470 + JWST DD 12781 cover its reappearance) — record our model window vs the outcome when public (HST 18470 visits 1–4 checked: no new point source to 2026-10-07). Same-source feasibility + VENUS batch 4 (#327; D-100, `results/transients/batch_venus4/RESULT.md`): model-free same-source proposal fails its test (none adopted); batch 4 13 of 20 boxes calibrated, 10 reviewed, 8 pending; A370 1 h000 within control scatter, SMACS 1 h005 not decided (control photometry stalled; rerun is the next check); held undeterminable boxes 33. Confirmation round 2 (D-102, `results/transients/confirm_r2/RESULT.md`): SMACS 1 h005 within control scatter; G165 h002 achromaticity test not decided, low power (paused). MACS1423 h000 HST archive (`results/transients/macs1423_h000_hst/RESULT.md`): no change before (−181 d) or after (+180 to +409 d) at the stated limits; no archive step left, identity open. W3 on ZTF DR24 (`results/w3_ztf/RESULT.md`): feasibility trial failed its frozen validation in v1 and bug-fixed v2 (survivors in sample-wide non-detection windows; slow large-ρ edges lost); v3 (counting fixed) still fails; images of 3 of 8 survivors show coverage gaps and quadrant changes, no real drop (`results/w3_ztf/v3/RESULT.md`); v4 on 1,000 new stars with footprint coverage: 0 flags, fails only at ρ 0.1, t_E 100 d (0.405) (`results/w3_ztf/v4/RESULT.md`); v5 exploratory holdout search (486 groups) 0 flags (`results/w3_ztf/v5/RESULT.md`); v6 cross-sky validation (Npix 7134) failed, 10 of 10 cells < 0.5, search group not executed (`results/w3_ztf/v6/RESULT.md`); v7 diagnosis: missing 2022-10 to 2023-10 catalogue epochs (`results/w3_ztf/v7/RESULT.md`); v8: field 672 gap exposures absent from DR24 (HATS = API), image fill fails its control (`results/w3_ztf/v8/RESULT.md`); v9: official psfcat fills the gap for 3 stars after passing the overlap test (`results/w3_ztf/v9/RESULT.md`); v10: psfcat in the chain, same events 711 → 872 of 2,000, 2 of 10 cells ≥ 0.5, validation failed, 539-group search not executed (`results/w3_ztf/v10/RESULT.md`); paused; 527 of 1,128 lost events have < 3 official r-band exposures in the umbra. v11: area chosen by observing conditions (Npix 3839, worst-cell necessary fraction 0.70), new sample passes validation (0.615–0.855), 488 search groups searched, 0 flags (`results/w3_ztf/v11/RESULT.md`); paused. v12: new-field area Npix 2382 fails validation (ρ 0.1 / 0.7 / 30 d at 0.42; sharp-edge test), 2,500 search groups not executed (`results/w3_ztf/v12/RESULT.md`); paused; next step is a method decision for the owner, not more areas. v13: edge-duration rule (from the W3 model) passes validation (cells 0.64–0.87, fade 8.5 %), 2,000 groups searched, 0 vetted survivors, image check finds detector-column artefacts, no confirmed dimming (`results/w3_ztf/v13/RESULT.md`); paused. v14: the 3 newly checked #370 class-4 NDs sit on official DEAD + LOW-RESPONSIVITY mask pixels (656 NDs pending); Npix 8170 passes validation (fade 5.0 % of 1,000); 5,000 groups searched, 1 survivor = near-limit catalogue drop, no confirmed dimming (`results/w3_ztf/v14/RESULT.md`); paused. v15: 3837 searched (4,778 groups, 0 vetted survivors), 183 fails validation (not searched); 11,778 groups searched in total with the v13 rule, no confirmed dimming (`results/w3_ztf/v15/RESULT.md`); paused. v16: 25-partition inventory, frozen 8-area queue, 3 pass, 5,995 groups searched, 2 survivors read as data problems; 17,773 searched in total, no confirmed target event; remaining supply 4 partitions / 17,702 groups / 5.76 GB (`results/w3_ztf/v16/RESULT.md`); paused. v17 survivor evidence check (post hoc): both undetermined; 1624 data problem supported, not proven (official LOW RESPONSIVITY column within 1.70 px of the source in all 9 NDs, ≥ 2.11 px in the detections; column-interpolated ratio ≥ 0.75 in 8 of 9, 6 of 9 exclude a 1-mag drop at 3σ), 1625 consistent with a blend in 4–4.6″ seeing next to a star 3.85 mag brighter, not tested to a decision; no real-dimming evidence (`results/w3_ztf/v17/RESULT.md`); closing either needs a mask-aware PSF fit or difference-image forced photometry (not built). v18 PSF check (official DAOPHOT PSF + photutils, validated on controls first): validation passes in 1 of 9 NDs (1624) and 1 of 3 (1625), both normal (1-mag drop excluded in those exposures; 1625's validated with its own PSF only); the rest unknown; both undetermined; follow-up closed until a seeing-matched independent PSF or difference-image forced photometry exists (`results/w3_ztf/v18/RESULT.md`). v19 audit: v18 validation was mostly a cross-exposure PSF stress test; corrected once (own PSF, phase-matched real controls, real-stamp injections, validated reference): NDs passing the corrected validation: 1624 1 of 9 (206251369, normal), 1625 2 of 3 (normal; 214647008 conditional, dim-level ρ 0.92); v18's 206251129 claim suspended (V1b′ σ defect); both undetermined; closed (`results/w3_ztf/v19/RESULT.md`); archived as "no reliable evidence of a real dimming; some exposures still cannot be judged". v20: frozen queue of 9 from the v16 inventory, 7 started, 2 pass (8126, 2900), 7,724 groups searched, 1 survivor (undetermined; a near-limit catalogue miss is the untested leading explanation); 25,497 searched in total with the v13 rule (17,773 + 7,724, not stored as one number), no confirmed target event; remaining 8124 and 8127 (9,978 searchable, 3.44 GB) (`results/w3_ztf/v20/RESULT.md`); paused. v21: v20 survivor's 4 NDs fail the frozen depth/completeness criteria (limits 0.32–0.60 mag below the g 18.92 baseline, 71–89 % completeness), still undetermined; 8124 and 8127 pass, 9,978 groups searched, 1 survivor undetermined (on an official low-responsivity column; not proven); 35,475 searched in total with the v13 rule (not stored as one number), no confirmed target event; no sampling-qualified v16 partition left unstarted, a new inventory needs the owner (`results/w3_ztf/v21/RESULT.md`); paused. v22: new inventory of 31 partitions (25 qualify, none sharing fields); queue of 8: 1 passes (2344, 1,882 groups searched, 0 survivors), 7 exceed the fade control (10.2–13.8 %; only failing condition in 3); 8124 neighbour also on a low-responsivity column in the shared NDs (8124 still undetermined); 37,357 searched in total with the v13 rule (not stored as one number), no confirmed target event; 17 qualifying partitions (29,698 searchable) remain in this inventory; a method decision on the fade control is the owner's (`results/w3_ztf/v22/RESULT.md`); paused. v23 diagnosis: passing fades rely on unmeasured ("zero-day") edges (129 / 138 in 4991, 76 / 81 in 2344), as do 35–99 recovered W3 per cell; the pre-registered bounds rule removes all fade passes but leaves W3 cells at 43–114 / 200, so its gate failed; no validation or search; v13 stays; next is the owner's method decision on edge information (`results/w3_ztf/v23/RESULT.md`); paused. v24: the undecided edges are mostly without a further covering exposure (64 %) or quality-flagged (27 %); candidate psfcat-fill gaps are 8 % of W3 sides, with an optimistic gain of 0–5 per cell; no download; existing-product restoration closed (`results/w3_ztf/v24/RESULT.md`); paused. ASAS-SN trial: Sky Patrol light-curve service down (3 probes 2026-10-20), 40 stars selected and baselined (63 insufficient W3, 12 insufficient fades), no flux read; re-run the committed plan when the service returns (`results/w3_ztf/asassn/RESULT.md`); paused. No follow-on batch without the owner. Pending, awaiting the owner: the SMACS ERO → VENUS
     whole-box shift (2.11, −1.82) px, not checked;
   - **SMACS 0723 arc `arc_diff`** (D-027 addendum 2026-10-16, system 5 box): ERO → VENUS done, Atalanta recovered,
     no candidate (frozen 1″ rule fails the known case; post hoc 2″ + forced photometry); GO-4043 → DD 12774
     fails its control (Atalanta F115W S ≈ 2.4); h061 undetermined (near-star spike systematics; reopen on a PSF
     model of the star). Next: the rest of the SMACS arcs (more boxes);
     ~~whether AT 2024jfz (1.8″ from image 6.1) has model counter-images~~ done 2026-10-19
     (`results/transients/smacs_at2024jfz/RESULT.md`): yes for z_s ≥ 1; "not testable" in the NIRCam epochs for
     z_s ≤ 2.5 (z 1.71: +31 d, 4.5″), "testable" at z_s 4 only as an unresolved 0.087″ pair in 2024; ~~NIRISS 4368 F150W epoch~~ done 2026-10-19 (`results/transients/smacs_at2024jfz/niriss4368/RESULT.md`: frozen outcome absent at z_s 1.35–1.40; post hoc, single-peaked light curve: absent for z_s 1.35–1.50 and 1.75–2.50, undetermined for 1.55–1.70 next to a static source; the transient is still visible at +116 d); ~~host redshift~~ checked 2026-10-19 (`results/transients/smacs_at2024jfz/host_z/RESULT.md`: frozen outcome untested, no DJA source within 1″, no MUSE z within 3″; post hoc the nearest source is image 6.1 at 1.63″, z_phot 1.60–1.63, inside the undetermined window). ~~PSF-matched difference at the z_s 1.55–1.70 track~~ done 2026-10-19 (`results/transients/smacs_at2024jfz/psfmatch/RESULT.md`: frozen outcome untestable, star core residual; post hoc absent > 0.35″ from the star). Next only a spectrum at the position or a new epoch;
   - **Full-chain known-case check** (2026-10-16, `results/transients/arc_chain/RESULT.md`):
     - MACS0416 Yan+2023: 4 / 4 pass detection, morphology and the arc mask;
     - SMACS Atalanta fails the position stage under both the 1″ rule and the mask;
     - SMACS in-mask full-chain efficiency: 0.50 at 27.5 AB, 0.31 at 28.0 AB.

     Owner chose position as a label (2026-10-16).
   - **Review-list chain** (2026-10-16, `results/transients/review_list/RESULT.md`):
     - known cases 6 / 6 on the list; burden ≤ 4 added entries per box;
     - SMACS whole-box injection efficiency 0.34 at 27.5 AB, below the frozen 0.50 (outside the mask 0.25, lost to
       rosette and dipole vetoes).

     **No new arc box until the off-arc efficiency is understood.** Options for the owner:
     - ~~a no-injection control of the vetoed positions~~ done 2026-10-16 (`results/transients/veto_control/RESULT.md`:
       frozen outcome mixed, 19 / 19 vetoed and 15 / 19 passed off-mask positions carry a pre-existing negative lobe,
       so it does not separate them; next diagnostic: a location-aware lobe measure on re-made injected maps);
     - accept 0.34 and state it per box. This would loosen the frozen R2 threshold, so it needs an owner decision and
       a new plan.

     ~~Open, unvetted by photometry: SMACS h039 and Sunrise o052 h005~~ done 2026-10-19
     (`results/transients/review_list/vet_h039_h005/RESULT.md`: frozen rule inconclusive for both; post hoc within
     matched-control scatter; not candidates).
   - vet `c0049`'s epoch-2 streak, if it recurs elsewhere (a satellite or asteroid trail, or scattered light).

- **Speed** (owner focus, step 2): `find_images` is vectorised. Still open: evaluate published deflection maps
  (UNCOVER v2.0, RELICS, HFF) for fields without a Lenstool `best.par`.

## Then (M1 follow-through, after the M3 items)
4. **Star-stratum sources with NED galaxy matches** (run `20261007T162133Z-2511977c`). Check their classification and photometry.
   - `f090w_1719`: NED `[NDA2023] 00713`, type G, 0.03″.
   - `f200w_2054`: a stellar-locus member 0.23″ from the lensed image NED `[MJR2023] 002.3`; it has spikes, and an
     arc lies next to it.
5. **Vet MIRI's known `G_Lens` object** with `/vet-candidate` and `scripts/vet_evidence.py`.
6. **Spike streaks confirmed by a coincident DJA object** (D-014 "Revisit if"; run `20261007T202419Z-65728e80`): CEERS
   `f200w_2842` (#4) is detected in F200W only, with ellipticity 0.95 and area 376 px. A bright source lies 4.6″
   away, and a DJA object sits 0.04″ off. Test single-band, highly elongated detections near bright sources
   (a spike-geometry check), without dropping lensed arcs.
7. **MIRI confirmation** (D-014 limitation): MIRI passes 103 of 530. Confirm MIRI detections with NIRCam
   counterparts. SMACS and CEERS now both have DJA confirmation.
8. **Control comparability, what remains** (D-022): pipeline-catalog features are jwst 3.0.0 in CEERS and 2.0.1
   in SMACS, until 2736 is reprocessed. CEERS has no stellar locus (too few stars). DJA colours are still to be
   validated against CEERS DR1.0. Then compare the science and control feature distributions.
9. Diffraction-spike mask for multi-band spike detections, only if they appear (D-014 "Revisit if").

## Next (M2)
- Matched-aperture colours for other samples (D-013 is done for SMACS NIRCam):
  - CEERS: done (D-022). Validation against CEERS DR1.0 is still open.
  - MIRI: evaluate fixed-aperture colours.
- DJA photo-z as a ranking feature or stratifier (it is used for vetting since cycle 7).
- Re-fetch program 2736 after MAST reprocesses it with jwst ≥ 3.0 (expected around mid to late October 2026).
  Until then, compare only photometric columns across programs (D-010).
- Image embeddings (Zoobot via `timm`, DINOv2) on cutouts, evaluated against the baseline with
  injection-recovery. Keep the evaluation fields out of any training data.
- Crossmatch and classify: propagate Gaia positions to the JWST epoch, so fast-moving stars and brown
  dwarfs land in the star stratum (D-012). Identify lenses by name or literature, not only by otype.
- Cutouts: north-up panels and an option for shared brightness scaling across panels.
- Second epochs and medium bands: check the galaxy top k against VENUS (6882) and other later programs
  with `scripts/epoch_compare.py` (D-017). Moving sources are stars; variable ones are AGN or transients.
- NIRCam F444W grism spectra of program 4043 (SOURCES) for top-k sources they cover: redshifts and lines. Needs
  a reuse check (grizli or the jwst WFSS pipeline products).

## Owner setup (needs-human)
- Run `claude.exe` interactively in the repo once and accept workspace trust.
- `claude.exe setup-token`, then `gh secret set CLAUDE_CODE_OAUTH_TOKEN`; install the Claude GitHub App;
  test with an `@claude` issue (docs/operations.md §5).
- Cloud routine: live since 2026-10-08 (docs/cloud-routine-prompt.md). Check that the environment's Custom network
  allowlist has every host in docs/operations.md §3. Merge any PR labelled `merge-ready` (cloud runs merge with the
  GitHub MCP tool when they can).
- Zenodo DOI at the first tagged release; optionally a separate Code of Conduct contact.

## Later
- Reuse `classify`'s full-catalog SIMBAD/Gaia matches for the top-k cross-match (re-query only NED).
- M3: published SMACS 0723 lens models (Mahler+2022, RELICS, Caminha+2022); AnomalyMatch for lens finding.
  Evaluate model shear and magnification at candidate positions (open question in
  `docs/candidates/jw02736-o001_t001_nircam_f200w_1096.md`).
- A weekly scheduled link check of the state files (`scripts/check_links.py`).
- Switch Dependabot from pip to uv once a `uv.lock` exists. Keep the version in `pyproject.toml`,
  `__init__.py` and `CITATION.cff` in sync on release.
- NIRISS imaging of program 2736 (F115W, F200W).
