# Changelog

Newest first. Results, failed approaches and the handoff state; not a diary.


## 2026-10-10: W3 MOA-II gb18 re-run under CHAIN_VERSION 2026-10-10.1: 0 survivors, limits tracked
- Hypothesis before running: 0 survivors; gb18-R-9-4-24509 (the D-068 chain survivor, a red giant's slow dimming)
  should now fail inside the chain. A survivor that held up on the contact sheet would have changed the plan.
- `run-field --procs 4` in 1,890 s: all 10 pre-screen chunks reused (83,855 light curves, 1,497 shape passes, 70 off
  shared epochs); 70 fits 400 s, 67 flags; 2,200 injections 1,126 s. **0 survive.** gb18-R-9-4-24509 fails
  `slow_dip_seasons` (ΔBIC 8.6). Contact sheet inspected: box dips, slow dimmings, variables; caustic spikes only
  in the models. CDS XMatch worked from the cloud (Gaia DR3 RP for 39 / 67 flags).
- Limits (`results/w3_moa/limits_gb18.ecsv`, docs/exotic_limits.md): Γ₉₅ ≈ 2.8 × 10⁻⁷–2.3 × 10⁻⁶ per star per year
  (N_s = 1.55 × 10⁷, 78 / 80 subfields); the t_E = 300 d, ρ = 0.1 cell recovered 0 / 200 (no limit there).
  Vetting keeps 83 / 162 flagged W3 injections (51 %); PSPL controls: 0 false W3 calls.
- All eleven D-068 fields are now re-run under the new chain: 0 survivors in each.
- Failed approach (tooling): a `pgrep -f "<pattern>"` wait loop matches its own shell command line and never exits;
  wait on the PID instead.
- **Next:** `combine` (check it handles an infinite-Γ cell), the cloud re-vet of the lenient-reference fields
  (gb7, gb11, gb16, gb19–gb22), then gb13.

## 2026-10-10: W3 MOA-II gb17 re-run under CHAIN_VERSION 2026-10-10.1: 0 survivors, limits tracked
- Hypothesis before running: 0 survivors; gb17-R-6-1-3829 (the D-068 chain survivor, a ~30-d periodic variable)
  should now fail inside the chain. A survivor that held up on the contact sheet would have changed the plan.
- `run-field --procs 4` in 1,797 s: all 13 pre-screen chunks reused (D-078 addendum 3, 0 s instead of ~8 min;
  100,448 light curves, 1,437 shape passes, 85 off shared epochs); 85 fits 479 s (fit chunks still key on the whole
  `Params`), 75 flags; 2,200 injections 1,175 s. **0 survive.** gb17-R-6-1-3829 fails `exotic_chi2_cap`
  (6.67 > 5.46), as the D-068 addendum intended. Contact sheet inspected: box dips, slow dimmings, variables;
  caustic spikes only in the models.
- CDS XMatch worked from the cloud again (Gaia DR3 RP for 53 / 75 flags).
- Limits (`results/w3_moa/limits_gb17.ecsv`, docs/exotic_limits.md): Γ₉₅ ≈ 2.3 × 10⁻⁷–1.9 × 10⁻⁶ per star per year
  (N_s = 1.80 × 10⁷, Nunota et al. 2024, 79 / 80 subfields; the 300-d cells rest on 2 recoveries). Vetting keeps
  86 / 156 flagged W3 injections (55 %); PSPL controls: 0 false W3 calls.
- **Next:** gb18 (the last re-run field), then `combine`; the fit-chunk reuse fix (D-078 audit) would save ~8 min here.

## 2026-10-10: W3 MOA-II gb15 re-run under CHAIN_VERSION 2026-10-10.1: 0 survivors, limits tracked
- Hypothesis before running: 0 survivors as in D-068; a survivor that held up on the contact sheet would have
  changed the plan (vetting record, owner notified).
- `run-field --procs 4` in 2,686 s (pre-screen re-streamed ~810 s, 8 chunks, 56–200 light curves/s as archive
  throughput fell from 78 to 22 MB/s; 83,145 light curves, 4,936 shape passes, 122 off shared epochs; fits 2 chunks,
  118 flags; 2,200 injections 1,182 s). **0 survive.** First failing test: repeated deficit 30, residual deficit 24,
  eclipse dip 18, bracketing 14, fit domain 12, χ² cap 10, other 11. Contact sheet inspected: smooth U/V dimmings,
  eclipse-like boxes, variables; caustic spikes only in the models.
- CDS XMatch worked from the cloud again (Gaia DR3 RP reference for 88 / 118 flags).
- Limits (`results/w3_moa/limits_gb15.ecsv`, docs/exotic_limits.md): Γ₉₅ ≈ 3.2 × 10⁻⁷–3.4 × 10⁻⁶ per star per year
  (N_s = 1.04 × 10⁷, Nunota et al. 2024, only 62 / 80 subfields: the low end of 1.04–1.34 × 10⁷; both t_E = 300 d
  cells rest on 2 recoveries). Vetting keeps 114 / 171 flagged W3 injections (67 %); PSPL controls: 0 false W3 calls.
- The re-streamed gb15 pre-screen rows equal the tracked ones row for row (only header metadata differs), a third
  check of D-078 addendum 3; the tracked chunks are kept unchanged.
- **Next:** gb17, gb18 (pre-screen now reused, D-078 addendum 3), then `combine`.

## 2026-10-10: Acceleration V2 audit, first pass: MOA pre-screen chunks reused across vetting changes; injection-loss map
- Owner brief "Autonomous research acceleration V2": measured audit in `docs/audits/2026-10-10-acceleration-v2.md`
  (what the repository already covers, the gb12 injection-loss funnel, ranked next units).
- **Redundancy removed:** a vetting-only `Params` change invalidated every streamed MOA pre-screen chunk; gb12
  (342 s) and gb7 (487 s) re-streamed rows byte-identical to the tracked ones. Chunk validity now depends only on
  the pre-screen parameters, the tracked-row constants and `PRESCREEN_CODE` (D-078 addendum 3). gb15, gb17, gb18:
  all 31 tracked chunks are reused (0 s instead of ~6–8 min each).
- **Where W3 injections are lost (gb12):** Cut-0 45 % kept, shape cut 23 %, shared-epoch test 75 %, vetting 63 %.
  Faint sources (I_s > 19.5) carry most of the loss (near the noise). Testable design losses: the shared-epoch veto
  (25 % of shape passes) and short / long t_E vetting (`jackknife_nights` at 3 d, `feature_bracketed` at 300 d).
  Sparse sampling does make real peak-dip-peak signals fail as ordinary dips at t_E = 3 d.
- **Next:** gb15, gb17, gb18 with the reuse fix, then `combine`; a shared-epoch masking study on a dev/validation
  injection split with the real flags as the false-positive check.

## 2026-10-10: W3 MOA-II gb12 re-run under CHAIN_VERSION 2026-10-10.1: 0 survivors, limits tracked
- Hypothesis before running: 0 survivors as in D-068; a survivor that held up on the contact sheet would have
  changed the plan (vetting record, owner notified).
- `run-field --procs 4` in 1,528 s (pre-screen re-streamed 342 s at ~215 light curves/s, CPU 99 %; 73,386 light
  curves, 5,532 shape passes, 51 left after the shared-epoch test; 51 fits 203 s, 46 flags; 2,200 injections). **0 survive.**
  First failing test: repeated deficit 17, residual deficit 9, eclipse dip 9, bracketing 4, χ² cap 3, other 4.
  Contact sheet inspected: box dips, quasi-periodic variables, caustic spikes only in the models.
- The batched CDS XMatch (cat2 `vizier:I/355/gaiadr3`, 1″) worked from the cloud this time: Gaia DR3 RP reference
  for 32 / 46 flags, so the source-flux bound used real references, unlike the earlier cloud runs.
- Limits (`results/w3_moa/limits_gb12.ecsv`, docs/exotic_limits.md): Γ₉₅ ≈ 2.3 × 10⁻⁷–1.5 × 10⁻⁶ per star per year
  (N_s = 1.61 × 10⁷, Nunota et al. 2024, 79 / 80 subfields). Vetting keeps 99 / 156 flagged W3 injections (63 %);
  PSPL controls: 0 false W3 calls.
- **Next:** gb15, gb17, gb18, then `combine`; audit which vetting tests remove injected W3 signals
  (37 % loss per field).

## 2026-10-10: W3 MOA-II gb7 re-run under CHAIN_VERSION 2026-10-10.1: 0 survivors, limits tracked
- Hypothesis before running: 0 survivors as in D-068 (gb7-R-8-6-94052 was the known jackknife case); a survivor
  that held up on the contact sheet would have changed the plan (vetting record, CDS-dependent tests where CDS works).
- `run-field --procs 4` in 2,670 s (pre-screen re-streamed 487 s at ~150 light curves/s, CPU 99 %; 75,328 light
  curves, 1,450 shape passes, 44 off shared epochs; 44 fits 289 s, 42 flags; 2,200 injections 1,297 s).
  **0 survive**; gb7-R-8-6-94052 still fails `jackknife_nights`. Contact sheet inspected.
- Limits (`results/w3_moa/limits_gb7.ecsv`, docs/exotic_limits.md): Γ₉₅ ≈ 2.7 × 10⁻⁷–1.4 × 10⁻⁶ per star per
  year (N_s = 1.63 × 10⁷, Nunota et al. 2024, 78 / 80 subfields; the t_E = 300 d, ρ = 0.1 cell rests on 3
  recoveries). Vetting keeps 100 / 169 flagged W3 injections (59 %); pre-screen pass 6–12.5 %. PSPL controls: 0
  false W3 calls. Quiet χ²/dof 95th percentile 6.1.
- CDS XMatch failed again from the cloud ("Too many jobs", recorded in `vetting_gb7.json`): lenient default
  source-flux bound; 0 survivors stands.
- **Next:** gb12, gb15, gb17, gb18, then `combine`.

## 2026-10-10: W3 MOA-II gb11 re-run under CHAIN_VERSION 2026-10-10.1: 0 survivors, limits tracked
- Hypothesis before running: 0 survivors as in D-068; a survivor that held up on the contact sheet would have
  changed the plan (vetting record, CDS-dependent tests where CDS works).
- `run-field --procs 4` in 1,981 s (pre-screen re-streamed 334 s at ~190 light curves/s, CPU 99 %; 64,090 light
  curves, 1,579 shape passes, 44 off shared epochs; 44 fits 238 s, 41 flags; 2,200 injections 973 s). **0 survive**.
  First failing test: eclipse dip 13, repeated deficit 12, residual deficit 6, bracketing 4, robust errors 2,
  χ² cap 1, neighbour 1, smooth dip 1, slow dip seasons 1. Contact sheet inspected (eclipse-like few-day boxes, slow
  dimmings, a deep V-shaped ~15-d dip in gb11-R-2-7-83307 that the box model misfits; caustic spikes only in models).
- Limits (`results/w3_moa/limits_gb11.ecsv`, docs/exotic_limits.md): Γ₉₅ ≈ 4.0 × 10⁻⁷–1.1 × 10⁻⁶ per star per
  year (N_s = 1.09 × 10⁷ from Nunota et al. 2024, 76 / 80 subfields; every cell ≥ 6 recoveries). Vetting keeps
  103 / 163 flagged W3 injections (63 %); pre-screen pass 5.5–11 %. PSPL controls: 0 false W3 calls. Quiet χ²/dof
  95th percentile 6.1.
- CDS XMatch failed again from the cloud ("Too many jobs"; recorded in `vetting_gb11.json`): lenient default
  source-flux bound; 0 survivors stands.
- **Next:** `run-field` for gb7, gb12, gb15, gb17, gb18, then `combine`.

## 2026-10-10: W3 MOA-II gb16 re-run under CHAIN_VERSION 2026-10-10.1: 0 survivors, limits tracked
- Hypothesis before running: 0 survivors as in D-068; a survivor that held up on the contact sheet would have
  changed the plan (vetting record, CDS-dependent tests where CDS works).
- `run-field --procs 4` in 2,148 s (pre-screen re-streamed 381 s at ~170 light curves/s, CPU 99 %; 64,353 light
  curves, 1,427 shape passes, 32 off shared epochs; 32 fits, 31 flags; 2,200 injections 929 s). **0 survive**.
  First failing test: repeated deficit 10, neighbour 6, residual deficit 4, exotic domain 3, eclipse dip 3,
  bracketing 2, feature sampled 1, χ² cap 1, jackknife 1. Contact sheet inspected (slow dimmings, short few-point
  dips, a deep ~20-d dip in gb16-R-9-6-97939 that the box model misfits; caustic spikes only in the models).
- Limits (`results/w3_moa/limits_gb16.ecsv`, docs/exotic_limits.md): Γ₉₅ ≈ 2.9 × 10⁻⁷–2.3 × 10⁻⁶ per star per
  year (N_s = 1.50 × 10⁷ from Nunota et al. 2024, 79 / 80 subfields; both t_E = 300 d cells rest on 2 recoveries).
  Vetting keeps 71 / 122 flagged W3 injections (58 %); pre-screen pass 3–11 %. PSPL controls: 0 false W3 calls.
  Quiet χ²/dof 95th percentile 6.2.
- CDS XMatch failed again from the cloud (non-VOTable reply): lenient default source-flux bound; 0 survivors stands.
  `summary` now keeps this `gaia_rp_xmatch` status in the tracked `vetting_gbN.json` (TASKS follow-up).
- **Next:** `run-field` for gb11, gb7, gb12, gb15, gb17, gb18, then `combine`.

## 2026-10-10: W3 MOA-II gb19 re-run under CHAIN_VERSION 2026-10-10.1: 0 survivors, limits tracked
- Hypothesis before running: 0 survivors as in D-068; a survivor that held up on the contact sheet would have
  changed the plan (vetting record, CDS-dependent tests where CDS works).
- `run-field --procs 4` in 1,596 s (pre-screen re-streamed, rows identical to the tracked tables; 51 fits; 2,200
  injections 686 s). 48 flags, **0 survive**; gb19-R-4-4-31159 now fails `feature_bracketed`. Contact sheet
  inspected (variables, slow dimmings, season steps, eclipse-like boxes; no caustic spikes in the data).
- Limits (`results/w3_moa/limits_gb19.ecsv`, docs/exotic_limits.md): Γ₉₅ ≈ 4.1 × 10⁻⁷–2.9 × 10⁻⁶ per star per
  year (t_E 3–300 d, N_s = 1.21 × 10⁷ from Nunota et al. 2024; the 2.9 × 10⁻⁶ at t_E = 300 d, ρ = 0.01 rests on 2
  recoveries). Vetting keeps 78 / 112 flagged W3 injections (70 %); pre-screen pass 2–8.5 % stays the limiting
  factor. PSPL controls: 0 false W3 calls. Quiet χ²/dof 95th percentile 6.4.
- CDS XMatch failed again from the cloud ("Too many jobs"): lenient default source-flux bound; 0 survivors stands.
- **Next:** `run-field` for gb16, gb11, gb7, gb12, gb15, gb17, gb18, then `combine`.

## 2026-10-10: E1 re-pinned to one GBM snapshot (8cc2…, 4,392 bursts): every E1 conclusion unchanged (D-078 addendum 2)
- Hypothesis before running: only bursts added or revised after the old pins differ, so every E1 channel stays
  null; a family p below 0.01 would have changed the plan.
- Fresh fetch: only `fermigbrst.vot` changed (sha256 8cc2823f…, 424,672 B, 4,392 rows; D-078 lists 4,391 rows for
  this digest). ICECAT-1, GWTC and CHIME Cat 2 match their pins. `data/manifests/e1_gw_events.ecsv` is removed:
  every E1 script (GW ones included) now reads the single `e1_events.ecsv`, so the pins cannot diverge again.
- Like-for-like check (observed counts only): the new snapshot minus bn261008763 and bn261007236 (triggered
  2026-10-07/08, last modified 2026-10-09) reproduces all 75 tracked D-074 observed counts. The other revision
  after 2026-10-09 (bn260930833) changes no cell. No trigger is later than 2026-10-08, so a trigger-date cut is a
  no-op. With the new bursts, 5 GBM–GBM cells gain 1–10 pairs (1d–7d all 17,523 → 17,533).
- Re-runs (same seeds and scramble counts as tracked; 4 cores): `e1_events.py` 456 s, antipodal (`--n 2000`),
  signed lag, CHIME flag, GW directional, GW signed (15–51 s each).

  | Result | old pin | new pin |
  |---|---|---|
  | D-074 five requested channels, family p (jit) | 0.278 | 0.273 |
  | D-074 wide + GW cells, pooled p (jit / jitday) | 0.0058 / 0.0033 | 0.0053 / 0.0029 |
  | D-074 all 75 cells, global p (perm / jit / jitday) | 0.0021 / 0.0036 / 0.0047 | 0.0020 / 0.0043 / 0.0041 |
  | Antipodal global p (addendum 2) | 0.557 | 0.531 |
  | Signed-lag global p (addendum 3) | 0.676 | 0.568 |
  | GW sky-map channels global p (addendum 4) | 0.829 | 0.833 |
  | Signed-lag GW cells global p (addendum 5) | 0.555 | 0.551 |
  | CHIME flag tests (CHIME only) | — | identical |

  The wide + GW driver is still CHIME–CHIME 1 h–1 d wide, which the calibrated rate-modulated null explains
  (D-074 addendum; CHIME-only, not re-run). Changes outside GBM cells are Monte Carlo noise from shifted
  random streams: one RNG per scramble serves all catalogues in order, so two extra GBM bursts reshuffle the
  CHIME and ICECAT draws too (|Δz| ≤ 0.36, |Δp| ≤ 0.01 in counts.ecsv). Injection outputs in limits.ecsv move more
  (40 trials per cell): `ul95_rate_per_anchor` ≤ 24 %, `eff` up to 19 % (CHIME–CHIME 1d–7d wide 0.999 → 0.806),
  `n50_detect` GBM–ICECAT 0s–10s wide 3 → 5.
- Fixed after `/code-review`: the null-ensemble cache name now carries the input digests (an old-snapshot cache
  was silently reused without `--refresh`); `e1_gw_events.ecsv` and its override removed (see above).
- Not done (follow-up): per-catalogue random streams in the nulls (`SeedSequence` keyed by catalogue) would make
  re-pins like-for-like outside the changed catalogue, but they change every tracked number, so they need their
  own PR.
- **Open risk / next (local session):** `main` pins a snapshot that only this ephemeral session held (cloud
  sessions cannot create releases, D-078). HEASARC serves 8cc2… until its next GBM update, so a local session
  should fetch it right away, check the digest, publish it to the derived-data store and wire `e1_events.fetch`
  to it. If the digest has moved on, run `e1_events.py --refresh` there and publish that snapshot instead.

## 2026-10-10: W3 MOA-II gb20 re-run under CHAIN_VERSION 2026-10-10.1: 0 survivors, limits tracked
- `run-field --procs 4` in 2,144 s (pre-screen re-streamed 299 s, rows identical to the tracked tables; 51 fits;
  2,200 injections 1,075 s). 48 flags (as in D-068), **0 survive**; contact sheet inspected (eclipse-like boxes,
  slow dimmings, season steps, variables; no caustic spikes in the data).
- Limits (`results/w3_moa/limits_gb20.ecsv`, docs/exotic_limits.md): Γ₉₅ ≈ 5.0 × 10⁻⁷–6.5 × 10⁻⁶ per star per
  year (t_E 3–300 d, published N_s; the 6.5 × 10⁻⁶ at t_E = 300 d, ρ = 0.1 rests on 1 recovered injection).
  Quiet χ²/dof 95th percentile 6.7 (gb22 5.3, gb21 6.3; D-068 revisits above ~10). Vetting keeps 73 / 124 flagged W3 injections (59 %); the four new tests cost
  8. Pre-screen pass 2.5–11 % stays the limiting factor. PSPL controls: 0 false W3 calls.
- Failed approach (env): CDS XMatch from the cloud failed (gb22/gb21 too), so the source-flux bound used the lenient
  default; 0 survivors stands (docs/exotic_limits.md, note under the gb20 table).
- gb19 was planned in the same run but not started (40-min run budget).
- **Next:** `run-field` for gb19, gb16, gb11, gb7, gb12, gb15, gb17, gb18, then `combine`.

## 2026-10-10: First derived-data release `derived-data-20261010` (GW sky maps); pinned GBM snapshot not recoverable (D-078 addendum)
- Local session on the owner's machine (cloud sessions cannot create releases); details in the D-078 addendum.
  Asset `gw_skymaps_nside32.npz`, 4,337,802 B, sha256 ab0330588b9a380f2f49b0d4dae809ea6cd1ed6ac1d99718e23f66d00c8e1ad5
  (pre-release, not latest), pinned in `data/manifests/derived_data.ecsv`.
- 200 MB rule: the 704 MB of tarballs were streamed into memory and never saved; only reduced products stay in the
  data root (the 4.3 MB npz and four per-tarball caches, 28 MB). The re-reduction reproduces every pin.
- Fresh session (empty data root): 2.7 s wall instead of streaming (row added to the Part 0 speed table below).
- GBM: neither pinned `fermigbrst.vot` is on this machine; nothing downloaded or invented.
- Fixed after `/code-review` (8 findings, 8 addressed): `derived_publish.py` creates the tag at the recorded code commit
  (`--target`; this release's tag sits on c92249c, its code commit is 67fffba); the TASKS re-pin item now names
  the procedure (one snapshot, a copy for the GW manifest, maps first, a like-for-like window).
- **Next (cloud routine):** E1 GBM re-pin (TASKS "Part 0"); TXS SkyLLH benchmark; round 2. Next local session:
  publish the new GBM snapshot and wire `e1_events.fetch` to the durable copy.

## 2026-10-10: W3 MOA-II gb22 and gb21 re-run under CHAIN_VERSION 2026-10-10.1: 0 survivors, first new-chain limits
- `w3_moa.py --field gbN run-field` (~20 min per field on 4 cores). gb22: 30 flags, gb21: 36 flags, **0 survive**;
  contact sheets inspected (dips, season steps, scatter; no caustic spikes in the data). Pre-screen rows identical to
  the tracked tables (only the `Params` header changed).
- Limits (`results/w3_moa/limits_gb22.ecsv`, `limits_gb21.ecsv`; docs/exotic_limits.md "Re-run under
  CHAIN_VERSION 2026-10-10.1"): gb22 Γ₉₅ ≈ 1.4–10 × 10⁻⁶, gb21 5.3 × 10⁻⁷–2.8 × 10⁻⁶ per star per year
  (t_E 3–300 d; gb21 has a published N_s). The efficiency is
  pre-screen-limited (3–12 % pass); vetting keeps 62 % (gb22) / 56 % (gb21) of flagged injections, and the four
  new tests cost ~5 %. PSPL controls: 0 false W3 calls.
- Failed approach (ops): `pkill -f q.sh` in a watcher killed the watcher itself (its own command line matched); kill
  queue processes by PID.
- **Next:** `run-field` for gb20, gb19, gb16, gb11, gb7, gb12, gb15, gb17, gb18 (gb17/gb18 have the most light
  curves, ~35 min each), then `combine`; the pre-screen dominates the efficiency loss, so a better pre-screen is the
  lever for deeper limits.

## 2026-10-10: E-NF1b IceTracks-DR2 short-lag ghost pairs: null; DR2 cannot reach R_g ~ 10⁻³ (D-079)
- Streamed IceTracks-DR2 v3.1 (28 files, ~210 MB, 3 s, sha256 pinned; raw not kept): 1,643,355 events, 4,963.4 d
  good-run union, no duplicates. Known case: TXS 0506+056 2014–15 box 6 on vs 1.11 expected (p = 1.0 × 10⁻³).
- Forecast first (most optimistic: every event astrophysical): reachable R ≥ 0.055 (1 h–1 d), 0.14 (1–7 d),
  0.27 (7–30 d), 0.64 (30–180 d) for every northern cut, so DR2 cannot improve ICECAT-1 beyond 1 h and R_g ~ 10⁻³ is
  out of reach at every lag (TASKS target withdrawn).
- Pre-registered short-lag test (6 cells, 5,000 uptime-aware scrambles; 33 s wall on 4 cores: null 12.4 s,
  injections 3.8 s): **null**, min Bonferroni
  p = 0.48. 95 % CLs limits per northern track: R ≤ 0.002 (log10 E ≥ 4, ≤ 100 s), ≤ 0.005 (≥ 4.5, ≤ 100 s),
  ≤ 0.01–0.02 (100 s–1 h); R_g per astrophysical neutrino is R / f_astro (not estimated).
- Found: the 9 pairs at 0–10 s (log10 E ≥ 4) are all from IC40/IC59, where the northern rate above the cut is
  ×7–18 the IC79/IC86 rate (the energy proxy differs by configuration); the jitter null absorbs it (p = 0.08).
  Post hoc IC86-only: 0 pairs at ≤ 100 s. Any future DR2 energy cut must be per season.
- `/code-review` (10 findings, 9 fixed): good-run union instead of a sum; cache tied to the manifest digest;
  `--mjd-min` also restricts the null's uptime and labels full-release fields; (run, event) group by unique
  index; retry with back-off; injections on the process pool (×4); shared `empirical_p`/`DETECT_P`/floor;
  `wide_pair_counts` reuses `wide_pair_stats`.
- Failed approach: urllib gets 403 from Dataverse (use requests/curl); Dataverse md5 is of the original CSV, so the
  manifest pins the sha256 of the served `.tab`.
- **Next:** NF-H04 needs a quantitative prediction before more compute (TASKS); TXS SkyLLH benchmark stays
  open for other DR2 uses.

## 2026-10-10: Owner brief: Neutrino Frontier round 1 (E-NF1 ghost pairs: null) and Part 0 speed-ups (D-077, D-078)
- **Neutrino Frontier** (new program, `docs/neutrino_frontier/`). N-A (literature-blind) wrote four frameworks;
  N-B (adversarial, checked references) judged three KNOWN-REDUNDANT: stochastic distance (NF-H01), pseudo-Dirac
  partners (NF-H02; N-A's δm² window is already excluded, arXiv:2406.06476), and the ν–γ Shapiro-delay test
  (NF-H03). Only NF-H04, delayed "ghost" copies at unrelated positions, is CONDITIONALLY VIABLE, as phenomenology
  (ladder B). The data audit confirms IceTracks-DR2 (doi:10.7910/DVN/MMIIZA, 1.64 M tracks, CC0) and that every
  audited host is reachable; track catalogs carry no flavor.
- **E-NF1 (NF-H04), ICECAT-1 v4, 340 alerts:** wide pairs in 7 lag bins out to 180 d, counts and signalness
  products (14 cells), 20,000 scrambles. It reproduces D-074 exactly (36 and 138). **Null:** pooled global
  p = 0.058; the largest cell is 1 h–1 d weighted, z = 2.7, Bonferroni p = 0.094. The 36 pairs were inspected:
  they span all 13 years, 6 are same-run and two are alert triplets (ordinary day-scale clustering). 95 % CLs
  limits on R_g per astrophysical alert: ≤ 0.02 (Δt ≤ 10 s), ≤ 0.05 (10 s–1 h), ≤ 0.2 (1 h–1 d), ≤ 0.1
  (1–7 d), ≤ 0.5 (7–30 d); at 30–180 d only R_g = 1 is excluded. Label: statistical null; a parameter limit on an ad-hoc model.
- `/code-review` (10 findings, 9 fixed): manifest written only after a successful release; deterministic pin
  choice; refusal of an underpowered `--n`; CLs limit needs every larger R_g to pass; durable copy accepted by
  content digest; NaN distance maps for unmapped GW events; unused masks removed; ghost Dec and error from one
  event. Not changed: `wide_pair_stats` duplicates `count_channel`'s binning (it adds weights; D-074 counts
  reproduce exactly).
- Failed approaches (E-NF1): classical limits excluded R_g = 0 on a low fluctuation, so CLs replaced them.
  Injections tested against the observation plus ghosts gave limits that were too tight; they now use a null draw
  plus the ghost excess. Without subtracting the random pairs that extra ghost events add, sensitivity at large R_g
  was overstated.
- **Speed-ups (Part 0)**, measured on this 4-core cloud machine:

  | Item | Before | After | Gain |
  |---|---|---|---|
  | S-2 `classify_gw_gw` per scramble (≈ 644 real GW–GW pairs, 30 seeds, identical classes) | 808 ms | 1.79 ms (numba); 28 ms (numpy fallback) | ×452; ×29 |
  | S-2 GW–GW pairs/s | 797 | 360,419 | ×452 |
  | E1 GW `count_all` per scramble (4 channels) | ≈ 0.8 s | 8.1 ms | ≈ ×100 |
  | `GWMaps` build (skip the 109 unmapped events' distance maps) | 150 s (36 maps + 355 uniform) | 27–29 s (282 maps) | ×5 |
  | E1 GW directional, full run (`--n 1000` null + injections + map build, same inputs, 60 cells identical, global p 0.833) | 333 s | 52 s | ×6.4 |
  | S-4 GW sky-map tarballs (704 MB, parallel ranges into memory, per-tarball cache) | > 40 min, unfinished (2026-10-09) | 159 s + 122 s for GWTC-4.1/5.0 (1.7–2.2 MB/s); cached re-run 1 s | — |
  | E-NF1 null (ICECAT, 14 cells; jit pass counted only to 7 d after review) | 1,375 scrambles/s | 1,982 scrambles/s (4 cores) | ×1.4 |
  | S-1 GW sky maps in a fresh session (owner's machine, empty data root, `derived-data-20261010`) | ~5 min to > 40 min of streaming (675 s on the owner's machine) | 2.7 s wall | ≈ ×250 |
- S-1: the store (`derived_store.py`, `derived_publish.py`, durable-first fetch in `e1_gw_skymaps.py`) is built and
  tested. **Release creation from cloud sessions is refused** (HTTP 403, "not permitted for this session type"),
  so the first asset waits for the owner or a local session. Today's re-reduction reproduces the pinned maps
  digest exactly (c5c2e4ff…).
- Found: the live HEASARC GBM table (8cc2…, 4,391 rows) matches neither pinned snapshot (D-074 5b42…, #122
  0cb3…), so both GBM-dependent E1 scripts stop in a fresh session until the snapshot is in the durable store or
  re-pinned. `nf_ghost.py` now fetches only ICECAT-1.
- S-5: the claim rule is now a 10-minute heartbeat and a takeover only after 60 minutes with no commit and no
  heartbeat (docs/cloud-routine-prompt.md, research-cycle skill, docs/operations.md).
- **Next:** publish `derived-data-20261010` (GW maps plus the GBM snapshot) from a session that may create
  releases; E-NF1 on IceTracks-DR2 (stream per-season files, good-run null); TXS 2014–15 SkyLLH benchmark first.

## 2026-10-10: W3 MOA-II: the four D-068 chain gaps closed as one post-hoc change (D-068 addendum, CHAIN_VERSION 2026-10-10.1)
- Hypothesis tested: the chain's survivors gb17-R-6-1-3829 (periodic variable) and gb18-R-9-4-24509 (slow dip with
  season levels) got through four gaps; new tests `exotic_chi2_cap`, `residual_deficit`, `periodic_variable`,
  `slow_dip_seasons` (thresholds in `Params`, ASSUMPTIONs).
- Validation: 6 synthetic W3 events in white noise pass every new test; unit tests for each model.
- Real-data controls (`vet` re-run, local only; tracked summaries wait for `run-field`): gb17 75 flags, gb18 67
  flags, **0 survivors**. gb17-R-6-1-3829 now fails `exotic_chi2_cap` (6.67 > 5.46); gb18-R-9-4-24509 fails
  `slow_dip_seasons` (ΔBIC +8.6, σ 415 d, p 5.2: the vetting note's model, found independently). First failures
  at the new tests: `residual_deficit` 13 + 11, `exotic_chi2_cap` 6 + 3, `slow_dip_seasons` 0 + 2.
- Failed approach: residual tests over every epoch failed a synthetic W3 event fitted at another u0 optimum
  (χ²/dof 1.69 vs 1.06, all on the spike nights): caustic-spike epochs are now left out. `/code-review` fixes:
  bright-source wings not spikes, whole-curve residual scan (window-only scans lost the red-noise normalisation),
  F-test error scale, P ≥ 2.5 d (nightly Nyquist), untestable cases named, no duplicate season refit.
  Not changed: `fit_slow_dip_seasons` duplicates `fit_smooth_dip`'s optimiser loop.
- **Next:** `run-field` for every field under the new chain (injections + limits; ~20 min per field), then `combine`.
  Every 2026-10-09.1 vetting record is stale; no limit may be quoted until then.

## 2026-10-10: D1 kinematic D_s/D_ds on 76 lenses: null, max 2.49σ vs 5.78σ (D-073 addendum 4)
- Hypothesis: one lens's kinematics demand a D_s/D_ds no FLRW model gives its redshifts. Data: the hierArc
  kinematic likelihood pickles in TDCOSMO 2025 (TDCOSMO 8, SLACS KCWI 13, SLACS SDSS 41, SL2S 14), read numpy-only.
- `d1_distance.py kinematic` (~1 min): LOO pull vs the others' offset and intrinsic scatter (τ = 0.15 in ln). Max
  SDSSJ2302−0840 −2.49σ; shuffled-z null median 2.71σ (p = 0.79); upper prior bound ×10 changes nothing.
  Detectable ×0.21 / ×5.6 (median lens; KCWI ×0.30–0.45 / ×1.6–2.9). Plot `results/d1_distance/kinematic_lenses.png` checked: SL2S sits low as a sample
  (δ −0.36 ± 0.15), KCWI vs SDSS differ by 0.2 on the same lenses; sample systematics, no single-lens outlier.
- Failed approach: mean/sd summaries and a grid-edge refusal; the hierArc term has a power-law upper tail
  (error ∝ prediction), so quantiles and a stated grid bound (D_s/D_ds ≤ 40) replace them. `/code-review`: γ_pl grid
  0.1 → 0.025 steps (moved SDSSJ1538 2.56σ → 2.48σ), axes by name, ASSUMPTIONs into `Params.kin_*`, τ-bound flag.
  `_ZGRID` now reaches z = 5 (SL2S sources to z 3.35).
- **Next:** W3 MOA-II chain gaps and re-injection (D-068); D1 is exhausted on public data until per-lens λ_int-free
  kinematics (JWST/KCWI IFU) or new time-delay lenses appear.

## 2026-10-10: D1 handoff: hierArc kinematic pickles are readable with the existing safe unpickler
- RXJ1131's `*_const_processed.pkl` is a `DdtHistKin` dict with numpy-only globals; details and next step in TASKS
  "Now" 0 (D1). Not yet used for any result.

## 2026-10-10: D1 composite model on SDSS1206: null; model choice moves its pull 0.61σ → 0.54σ (D-073 addendum 3)
- Hypothesis: the D-073 null depends on the power-law mass model; a composite (stars + NFW) model could shift one lens
  out of line. Only SDSS1206 has a public composite D_dt chain in TDCOSMO 2025.
- `d1_distance.py tdcosmo --j1206 final_composite` (κ_ext included in the chain): J1206 0.54σ, max |pull| 1.14σ vs
  5.39σ, null. Validation first: TDCOSMO's final power-law chain gives 0.606σ vs this pipeline's 0.608σ.
- Failed approach: the pre-LOS composite pickle cannot be paired per sample (κ_pert 20,000 vs 400,000 rows).
- **Next:** new per-lens D_d via a numpy-only reader for the hierArc `*_const_processed.pkl`; W3 MOA-II chain gaps
  and re-injection (D-068).

## 2026-10-10: W3 MOA-II eleven fields: 582 flags, 0 candidates; limits withdrawn pending the corrected chain (D-068, #95)
- Lands the streamed MOA-II work: gb7, gb11, gb12, gb15–gb22 (712,780 Cut-0 light curves, 582 flags). Nine fields
  have 0 chain survivors; gb17-R-6-1-3829 (periodic variable) and gb18-R-9-4-24509 (red-giant slow dimming) passed
  the chain and are explained only by vetting notes (docs/candidates/). `derived`. No W3 candidate.
- Merged `claude/w3-moa-fields-alt2` (injections vetted against the injected source flux, `CHAIN_VERSION`, refusal
  without MulensModel, truncated-tar and size-pin fixes). Every old-chain `limits_*` and `efficiency_*` table is
  deleted: none may be quoted.
- `/code-review` fixes: the vetting record carries `CHAIN_VERSION` and `limit` refuses a mismatch (the eleven
  tracked records predate the stamp, so every field is re-vetted before a limit); a non-finite night-jackknife refit
  now fails the test (it returned −inf and passed; this could only have kept flags, so 0 survivors stands);
  `fit_step_ramp` without an admissible step time; ragged rows go to the line parser; chunk ranges are checked.
  Not changed (follow-ups): duplicate check over all streamed ids, `prefetch` floor above `--conns`, repeated chunk
  reads, per-field metadata pass, unnamed quiet-carrier thresholds in `is_quiet`.
- **Handoff / next:** add the four chain gaps in D-068 as one post-hoc change, re-inject every field, `combine`;
  gb13 stopped at 25 / 78 fits (run died 2026-10-09 13:02Z); then more Nunota Table 1 fields.

## 2026-10-10: E1 signed-lag GW cells with sky-map classes: null, global p = 0.56 (D-074 addendum 5)
- Hypothesis: GW events lead (or follow) GBM, ICECAT-1 or CHIME events in a given sky-map class (same, wide,
  antipodal) at lags no ordinary path explains: D = N_after − N_before.
- `scripts/e1_gw_signed.py --n 1000` (71 s): 45 cells (3 channels × 5 lags × 3 classes; GW–GW left out because D
  is antisymmetric there). Pooled global p = 0.56, min analytic p × 45 = 0.78. Largest: GW–CHIME 1 h–1 d antipodal,
  D = 10 vs −0.3 ± 4.3 (z = 2.4; it was z = 3.1 with 40 scrambles: the smoke-run null sd was too small).
- Control GW170817 → GRB 170817A gives D = +1 at 0–10 s same; it is left out of the family.
- Sensitivity: 3 one-sided pairs at ≤ 10 s, 3–10 at 10–100 s, 10–30 at 100 s–1 h, mostly > 30 at 1 h–7 d. Per-cell 95 %
  upper limits on extra after-pairs are in `results/e1_events/gw_signed.json`. Cells are tested against the null mean
  (catalogue edges make it non-zero, e.g. GW–CHIME 1–7 d wide −66.5 ± 52.7).
- **Next:** D1 composite-model D_dt chains; E1 IceTracks-DR2 / Swift / Einstein Probe when they can be fetched.

## 2026-10-10: E1 GW channels with sky maps: null, global p = 0.83 (D-074 addendum 4)
- Hypothesis: GW events have same-direction, wide or antipodal partners (GBM, ICECAT-1, CHIME, GW) at lags no
  ordinary path explains. Until now the GW channels were lag-only.
- `scripts/e1_gw_skymaps.py` streamed the four GWTC PE sky-map tarballs (704 MB, about 3.5 min, one process per `--tar` then a merge run; the
  2026-10-09 run got 0.15–1 MB/s and did not finish) into 282 nside-32 maps. The 109 unmapped GWTC-4.1 / 5.0 entries
  have no PE parameters in the GWOSC CSV either.
- `scripts/e1_gw_directional.py --n 1000` (266 s, 4 cores): 60 cells, pooled global p = 0.83, min analytic p × 60 = 1.
  Largest z: GW–ICECAT 1–7 d antipodal (14 vs 8.3 ± 2.9). Control GW170817 × GRB 170817A is `same` at 0–10 s; it
  is left out of the family and reported apart.
- Sensitivity: 3 injected pairs at ≤ 10 s, 3–10 at 10–100 s, 10–30 at 100 s–1 h, mostly > 30 at 1 h–7 d.
- Catalogue snapshot pinned separately (`data/manifests/e1_gw_events.ecsv`): the GBM TAP table grew since the
  D-074 pin, so `e1_events.fetch` refuses the old one. Rerunning older E1 scripts needs `--refresh` or that pin.
- **Next:** signed-lag GW cells with same/wide classes; then D1 composite-model chains.

## 2026-10-09: E1 signed-lag ("which event comes first") channels: null, global p = 0.68 (D-074 addendum 3)
- Hypothesis: one catalogue leads another at lags no ordinary path explains (sign asymmetry D = N_after − N_before).
- `scripts/e1_signed_lag.py`: 45 cells (6 cross channels × 5 lags × class), `jit` null (2,000 scrambles, 52 s).
  Pooled global p = 0.68. The largest |z| is GW170817 → GRB 170817A (D = −1, Skellam p = 0.05), the known ordinary
  ordering.
- Failed approach: a Gaussian tail on one-pair cells (z = −4.4 for GW170817 alone, Bonferroni 6 × 10⁻⁴). Replaced by a
  Skellam tail where the null sd < 1. The empirical p could not be used instead: its floor (1/2001 × 45) makes
  injections undetectable.
- Sensitivity: 3–10 one-sided injected pairs at ≤ 100 s.
- **Next:** E1 GW sky maps (GW cells get directions; redo antipodal and signed lags), then D1 composite-model chains.

## 2026-10-09: E1 antipodal lag channels: null, global p = 0.56 (D-074 addendum 2)
- Hypothesis: a dependent event appears near the antipode of a first event, at 0 s–7 d lags. Ordinary: the D-074
  observer-side common causes.
- `scripts/e1_antipodal.py`: 180° − sep ≤ max(10°, 3σ_comb); 30 cells (6 localized channels × 5 lag bins); `jit` null
  (2,000 scrambles, 19 s on 4 cores). Pooled global p = 0.56. Max z = 1.7 (CHIME–GBM 1 h–1 d, 258 vs 230 ± 16).
- Injection-calibrated: 3 pairs at ≤ 100 s detected at family-wise 3σ in most channels, 3–30 at 100 s–1 h, 10–100 at 1 h–7 d.
- Fixed before the run: a same-catalogue channel must count the injected sample with itself (`pairs_within`
  assumes ta is tb); counting the original against the injected copy undercounted.
- **Next:** E1 GW sky maps (GW–X same/wide/antipodal channels), signed-lag ("which comes first") channels.

## 2026-10-09: E1 CHIME 100 s–1 h wide cell under the calibrated rate-modulated null: z 3.7 → 2.5, family-wise null (D-074 addendum)
- Hypothesis: the 100 s–1 h CHIME–CHIME wide excess (396 vs 330, z = 3.7 `jit`) is lag dependence beyond a smooth
  detection-rate modulation. Ordinary: week-scale rate modulation, daily duty cycle, sub-day outages.
- `scripts/e1_chime_flag.py --cell 100s-1h` (2,000 scrambles, 30 s): 7-day null z = 2.52 (analytic p × 45 = 0.26);
  time of day kept z = 2.42; 3-day z = 0.48 (not primary, D-074 (d)). The 7-day null keeps ≥ 97 % of 300 injected pairs.
- Result: null at the family level; a local 2.5σ residual remains that only hour-scale uptime could decide.
- Failed approach (not run, reasoned): an hour-resolved running mean of the catalogue's own counts would absorb the
  tested pairs, as the in-day null did.
- **Next:** E1 GW sky maps (GW channels get directions), then "which event comes first" / antipodal channels.

## 2026-10-09: Hypothesis round 3: derived small effects in large-N data; no testable survivor; new rounds paused (D-076)
- Brief: a derived small amplitude in a large-N public quantity never analysed for it. **Results:**
  - R3-A (same-exposure pair correlation, about 2.4 × 10⁻⁵): fails as written (axiom A5 gives zero). The Holometer
    disfavours it at about 67 % (reviewer estimate).
  - R3-B (condensate scalar charge in neutron stars): upper half (s ≥ 8 × 10⁻³) excluded at 6.7σ by J1738+0333; self-consistent
    minimum disfavoured at about 2σ; only the no-feedback corner survives. α_NS < 2.2 × 10⁻³ (95 %, reviewer's combination).
  - R3-C (SN Ia rate vs "clock depth"): not excluded but untestable. About 5 × 10³ public SNe Ia give 0.5–2σ, and
    environment systematics are 10–40× the signal (reviewer estimates).
- **Rule:** three rounds (9 worlds) produced survivors only in round 1, and those were null in data. New rounds are
  paused; compute goes to data-limited open searches (D-076).

## 2026-10-09: D1 on TDCOSMO 2025 power-law chains: 8 lenses null with κ_ext; without it RXJ1131 is 5.7σ off (D-073 addendum 2)
- Hypothesis: one lens's D_dt disagrees with the H0 + Ωm of the others (a sightline-specific distance anomaly).
  Ordinary explanations: κ_ext, internal mass sheet / model choice, time-delay systematics. A flag with κ_ext applied
  would have gone to composite models and κ_ext vetting.
- `scripts/d1_distance.py tdcosmo`: power-law D_dt^model chains (6 H0LiCOW lenses + DES0408 + WGD2038), each lens's
  own TDCOSMO κ_ext PDF applied identically, H0 + Ωm LOO; 100 shuffled-z nulls; injections.
- **With κ_ext: max |pull| 1.12σ (RXJ1131), threshold 5.39σ; null.** Reach ×0.41–0.78 / ×1.36–1.85 (WGD2038 ×0.26 /
  ×3.2). LOO H0 72.2–74.4. κ_ext files checked uncorrelated with their chains (|r| < 0.003).
- **Without κ_ext: RXJ1131 −5.67σ, DES0408 +3.9σ**, both removed by the measured κ_ext (positive control: the test
  sees a ~7 % line-of-sight convergence).
- Safety: SDSS1206's pre-LOS file is a pickle; read with a numpy-only unpickler (arbitrary classes refused, tested).
- Data: whole TDCOSMO2025_public checkout (661 MB, cloud) streamed by git, used, deleted; 16 files pinned in
  `data/manifests/d1_distance.ecsv`. Failed approach: `git clone --no-checkout` then `git checkout <sha>` fetches the
  whole tree; check out paths only.
- **Next:** composite-model chains (TDCOSMO SDSS1206 `final_composite_*`, others where released) for a model-choice
  check; per-lens D_d needs a data-only reader for `*_const_processed.pkl`.

## 2026-10-09: Hypothesis round 2: three worlds with derived Lorentz invariance; nothing testable here survives (D-075)
- System A (no literature) was bound by the round-1 failure modes:
  - derive Lorentz invariance;
  - no long-range links;
  - no achromatic dimming;
  - only derived amplitudes.
- **Results:**
  - R2-A (every two records share a future, so horizons re-emit everything): **excluded**. Ergoregion instability;
    GWTC-4.0 null; supporting: GW250114 reflectivity < 0.35 % (preprint).
  - R2-B (gravity as phase inference, so saddle images mirror polarisation): internally inconsistent; B0218+357
    polarisation variations correlate with the wrong sign.
  - R2-C (four axioms; massless lightest neutrino, so Σm_ν = 58.8 meV): not novel (minimal seesaw). 2.2–3.4σ
    tension with DESI DR2 + CMB (reviewer estimates; DESI quotes 3.0σ); tracked externally.
- **Lesson (rule):** deriving Lorentz invariance from a symmetric primitive is easy. The theory's risk then sits in
  one precision-measured sector, where O(1) effects are already excluded. Round 3 should look for derived small
  effects in large-N quantities that nobody has analysed for them.

## 2026-10-09: E1 causal event network: no wide-separation dependence between GRBs, neutrinos, GW events and FRBs (D-074)
- Hypothesis (owner idea 4): events in different directions that depend on each other at lags no ordinary path
  explains. Ordinary explanations: observer-side common causes (uptime, exposure, Sun, follow-up chains, shared
  triggers, duplicates).
- Data (4.7 MB, pinned): Fermi GBM 4,390; ICECAT-1 340; GWTC 391 (no sky positions, so lag only); CHIME/FRB Cat 2
  3,641 sources. Lag × separation pair counts in 75 cells. Nulls keep declination and hour angle; GBM is shifted in
  whole orbits.
- **Five requested channels: family p = 0.28.** All wide/GW cells: pooled family p = 0.0058 (jit) / 0.0033 (jitday).
  The driver, CHIME–CHIME 1 h–1 d wide (z ≈ 3.9; 5.8 per cell in the post-hoc nominal-flag subset), disappears
  under a calibrated rate-modulated null (7-day running mean of daily counts, keeps 87 % of an injected signal):
  z = −0.86. Same-direction CHIME excess: one unflagged six-burst episode (FRB20230825D–I).
- Positive controls: GBM re-triggers recovered. GW170817 × GRB 170817A is the only GBM–GW pair within 10 s, but
  blind lag counting gives only p = 0.058 without GW sky maps.
- Limits (95 %, family-wise 3σ detection, anchors inside the partner's live time) in `results/e1_events/limits.ecsv`.
- **Failed approaches (rules):**
  - Observed and scrambled p computed with different rank formulas lost the trials factor beyond the ensemble.
    Pool the observation with the scrambles.
  - A "keep each event's day" null removes ~85 % of an injected dependent signal, so it is not a vetting test.
  - Resampling times with replacement self-pairs events.
  - The CHIME Cat 2 exposure file (216 MB, downloaded and deleted) has no time axis, so it cannot model uptime.

## 2026-10-09: D1 handoff: TDCOSMO 2025 per-lens files listed
- File names in TDCOSMO2025_public d7f38db recorded in SOURCES.md: D_dt chains for the 6 D-073 lenses plus DES0408 and
  WGD2038; new per-lens D_d only inside hierArc likelihoods. Next: TASKS "Now" 0.
- An equal-weight NE2001/YMW16 mixture cannot flag more than the per-model runs (its tail is their average); a wider
  ISM error only lowers the pulls. Neither can turn the FRB null into a flag, so both wait.

## 2026-10-09: D1 FRBs with YMW16 beside NE2001: still null (D-073 addendum)
- Question: is the D-073 FRB null (and its 3.26σ extreme) an artefact of NE2001? Ordinary explanation for any
  low-side outlier: Galactic-model error at low latitude. A flag under YMW16 alone would have needed vetting.
- `scripts/d1_distance.py frb --ism ymw16`: YMW16 to 30 kpc per RA/Dec, same chain and 5.81σ threshold. Known case
  FRB 20121102A: 287 pc cm⁻³ (published YMW16 value). NE2001 path re-run: identical pulls and flags.
- **Result: no flag under YMW16 or under either model.** Max low pull FRB 20220319D 4.09σ (NE2001: 3.26σ;
  DM_obs 111 below both models, 127 and 211); max high 20190520B 2.37σ; injections 94/94 per side.
  YMW16/NE2001 median 0.84, range 0.36–1.91; > 20 % apart for 54/94 bursts. Outputs `results/d1_distance/*_ymw16*`,
  `frb_ism_compare.ecsv`.
- **Failed route:** `pip install pygedm` fails in the cloud image (its NE2001 extension needs libf2c) and the
  package import fails on current SciPy (`integrate.simps`); built without `ne21c`, `ymw16` imported directly.
- **Next:** TASKS "Now" 0 (per-lens TDCOSMO D_d; a model-marginalised ISM term).

## 2026-10-09: S2 on DES-SN5YR: no flat-kernel term; γ < 0.019 mag per unit T/⟨T⟩ (D-070 addendum)
- Same chain as D-070 (`scripts/s2_flat_kernel.py --sample des`, `results/s2_flat_kernel/fit_des_*.json`) on the
  DES-SN5YR Dovekie Hubble diagram: DES-discovered SNe only, MU/MUERR (BEAMS-renormalised), P(Ia) ≥ 0.5 (ASSUMPTION),
  0.1 < z < 1.3: 1,518 SNe, 1,414 fitted after coverage and 5σ cuts; LS DR9 galaxies (47.2 k) as before.
- **α = 2: γ_F = +0.0009 ± 0.0025 (scramble p = 0.68); trimmed −0.0012 ± 0.0044; colour term +0.0008 ± 0.0018.**
  λ_F = 0.86 (0.81 trimmed); whole-chain recovery 0.57 (limit divided by it). **One-sided 95 %: γ < 0.019 mag per
  unit T/⟨T⟩** (Pantheon+: 0.025). Still above A3's kill threshold 0.005; excludes the upper half of its 0.01–0.04
  fiducial only if the galaxy-to-matter bias is ~1 (ASSUMPTION).
- α = 1: kernel-degenerate again (corr(X_L, X_F) = 0.99), no limit, as in D-070.
- **Positive control:** γ_L = −0.018 ± 0.011 (lensing sign, scramble p = 0.095; λ_L = 0.60). Pantheon+ gave
  −0.017 ± 0.013; naively combined −0.018 ± 0.008 (~2σ), but the samples overlap (DES 3YR spectroscopic SNe are in
  Pantheon+), so this is indicative only. The control is still not a detection with z < 21 LS DR9 galaxies.
- Not done (failed route): the DES Y3 Gold / deeper-galaxy column the D-070 revisit asked for; LS DR9 at z < 21 was
  kept so that the two samples share one calibrated chain.
- **Next:** a deeper galaxy column (DES Y3 Gold or LS DR10 z < 22 with photo-z quality cuts) so the lensing control
  detects; the full DES-SN5YR STAT+SYS covariance (inverse matrices in the repo) instead of diagonal errors.

## 2026-10-09: D1 distance self-consistency: no inconsistent sightline in strong lenses or localized FRBs (D-073)
- Hypothesis (owner idea 1): one sightline whose independent distance measures no ordinary model reconciles.
  Ordinary explanations: mass-sheet degeneracy, line-of-sight convergence, kinematic anisotropy, substructure, wrong
  redshifts; for FRBs, Galactic ISM model error and host DM.
- Lenses: H0LiCOW public posteriors. The H0-free R = D_dt/((1+z_d)D_d) = D_s/D_ds for B1608, RXJ1131, PG1115, J1206
  (prior-predictive ΛCDM / wCDM and leave-one-out), plus leave-one-out D_dt for those and HE0435, WFI2033. **Max pull
  1.29σ** against 5.53σ (18 trials). Shuffled redshifts give a median max pull of 3.8σ (C) / 15σ (D).
- Injection-calibrated reach (|z| crossing the threshold, baseline pull included): R tests ×2–8 up / ×0.2–0.4 down;
  D_dt ×1.3–1.85 / ×0.54–0.78. A ≲ 30 % mismatch on one sightline is not excluded.
- FRBs: 94 localized FRBs (FRBs/FRB repo) against the Macquart predictive distribution (grid convolution, tails to
  ≪ 1e-10). Nothing flagged at the one-sided 5.81σ; the lowest is FRB 20220319D at 3.26σ (DM below its NE2001 DM_ISM,
  a known low-latitude case). Injections flag 94/94 low and 94/94 high.
- **Failed approaches (rules):**
  - A Monte Carlo predictive tail floored at 1/N capped z at 4.26σ below a 5.93σ threshold, so the screen could
    never flag (found by review). Tails must reach past the threshold: test it.
  - Sensitivity factors from a Gaussian extrapolation ignored baseline pulls. Read them from injections.
  - TDCOSMO 2025's SDSS1206 D_d file duplicates H0LiCOW's samples and adds no lens.
- **Next:** see TASKS "Now" 0.

## 2026-10-09: S3 hybrid images: no copies in 6 COSMOGRAIL doubles; sensitivity only r ≈ 0.4–2.4 (D-072)
- Hypothesis (A1 P2): image j carries a faint copy at lag s_i − s_j; the copy into the trailing image can precede the
  leading image. Ordinary mimics: quasar red noise (correlates at all lags), delay errors, microlensing.
- Data: COSMOGRAIL XIX (CDS J/A+A/640/A105), 7 doubles with secure delays. SIS lag window for f ∈ [0, 1]
  (D_ls/D_s = 0.46–0.68). `scripts/s3_hybrid.py`, `results/s3_hybrid/`.
- Known-case gate: the main-term fit recovers the published delay (sign, ≠ 0, within 3σ + 4 d) in 6 of 7. J1620 fails
  with a wrong-sign minimum at τ_BA = −206 d (published +171.5 d) and is not screened.
- **Result:** no copy. Window maxima against same-width off-model windows (|lag| < 10 d excluded from both): lowest
  p = 0.02 (J1226 into-trail, a window on the visible main-image wing) of 12 windows; the windows overlap, so p is
  approximate and the run-wide chance of p ≤ 0.02 is ~0.2. Sensitivity r95 ≈ 0.38 (J1455 into-trail) to 2.4: only
  copies comparable to the main image are excluded, far from B-on-A1's ~10⁻². Not quoted as limits: the
  window-maximum injection efficiency is 0.67–4.7 (> 1 is the copy–main degeneracy; capped at 1 in r95), the
  null windows are not matched in |lag| (windows next to the wing get low p), and for HE0047 the delay gate
  (±14.5 d around 10.4 d) is weak.
- **Failed approaches (rules):** microlensing splines with 120–730-d knots absorb the quasar variability and fit any
  delay (6 of 7 published delays missed); use the delay recovery as the gate for every setting. Quasar variability
  is too slow for copies at lags ≲ 100 d: the copy is collinear with the main image (a synthetic 20-d DRW recovers
  r = 0.1 ± 0.03, so the loss is the data). /code-review: inject through the window maximum, not the fit at the
  injected lag; hold the observed window to the null's wing rule; inject a smoothed copy (the template's noise
  otherwise rides along); fit the main term at the validated delay.
- **Next:** S3 needs a sharp template: SN Refsdal / SN H0pe pre-explosion and post-peak HST/JWST imaging at the
  model lags (B-on-A1: a few × 10⁻²). Low priority; S1 and S2 (DES-SN5YR) first.

## 2026-10-09: S1 burst twins in Fermi GBM: null; < 7.2 × 10⁻³ twin pairs per eligible burst (D-071)
- Hypothesis (D-069 S1, from A2 P1): pairs of GRBs far apart on the sky whose light curves match at any delay.
  Ordinary mimics: single-envelope look-alikes, re-triggers, ordinary lensing (excluded by the position cut).
- Data: HEASARC `fermigbrst` (4,390 bursts) and bcat HDU 2 light curves (4,389; streamed, sha256 recorded).
  `scripts/s1_ingest.py`, `scripts/s1_twins.py`, `results/s1_twins/`.
- 1,414 multi-pulse bursts; 510,294 pairs with inconsistent positions (sep > 3·√(σ₁² + σ₂² + 2·3.7²) deg).
  Pulse-shuffle null ρ* = 0.964 (a pre-screen: the surrogates under-predict the real high-ρ tail ~12×). 6 pairs
  flagged, all bright single-envelope look-alikes failing the χ² twin test (p < 10⁻¹⁰⁰). **0 survivors.**
- Injections with a 15 % per-band mismatch measured from the data: ratio-1 efficiency 29.5 %. **95 %: < 7.2 × 10⁻³
  twin pairs per eligible burst (1.4 % of bursts have a twin) at flux ratio 1, 1.6 × 10⁻² averaged over ratios 1–0.1.**
  Sensitivity: 5.1–9.4 × 10⁻³ for 5–20 % mismatch.
- **Failed approaches (rules):**
  - bcat HDU 1 `PHTCNTS` fill values made 84 % of bursts look multi-pulse.
  - Self-copy injections shared the target's noise (recovery 83 % vs 37 %).
  - Injections whose noise matches the χ² errors exactly pass χ² by construction: inject a measured band mismatch.
  - A ρ threshold chosen so that k = 0 (ρ > 0.90) is post hoc and illustrative, not a limit.
- **Next:** CHIME/FRB Catalog 2 (the proxy resets CANFAR downloads); TTE light curves for short and faint bursts; an
  s ≠ 1 chain with energy-channel mapping; a generative pulse-model null.

## 2026-10-09: S2 flat-kernel SN residuals: no flat-kernel term; γ < 0.025 mag per unit T/⟨T⟩ at α = 2 (D-070)
- Hypothesis (A3 P2b): SN Ia residuals track a flat-kernel foreground column ∫(1+δ)^α dχ, with sign +γ (fainter),
  after the lensing-kernel column is fitted. Ordinary mimics: lensing magnification, grey dust (colour), host-group
  overlap, leverage of a few sightlines, disc coverage, shot noise in the column (regression dilution).
- Data: Pantheon+ (944 unique Hubble-flow SNe, z 0.1–1.3) and LS DR9 photo-z galaxies (countmap selection, dereddened
  z < 21 mag, 2′ discs, 32.3 k galaxies). `scripts/s2_flat_kernel.py`, `results/s2_flat_kernel/`.
- **α = 2 (A3's fiducial), 505 SNe** (shells expecting ≥ 1 galaxy, so z_s ≳ 0.2; coverage and 5σ cuts remove 15):
  γ_F = −0.0003 ± 0.0018 per unit X (scramble p = 0.83); 1 % trimmed +0.0024 ± 0.0043; no colour term
  (−0.0017 ± 0.0012). Dilution λ = 0.92 (0.75 trimmed). Whole-chain injection (γ = 0.01 on full counts, measured on
  half-thinned counts, 20 trials) recovers 1.04 (full) / 0.78 (trimmed) of the input; the limit is divided by 0.78.
  **One-sided 95 %: γ < 0.025 mag per unit galaxy-traced T/⟨T⟩.** A3's fiducial is 0.01–0.04 (B-on-A3 expected a
  ~0.015 reach), and its kill threshold is 0.005. Converting to a matter column (bias b ~ 1–2, ASSUMPTION) weakens the
  limit by up to ~b².
- Lensing column: γ_L = −0.017 ± 0.013 (the lensing sign, p = 0.21). The positive control is not detected, as
  expected at this N (Smith+2014: 1.4σ with 608 SNe).
- **Failed approaches (rules):**
  - The linear (α = 1) flat column is 0.996-correlated with the lensing column in 2′ disc counts. It cannot separate
    the kernels: the chain injection is not recovered, so it sets no limit.
  - N(N−1)/E in shells with E ≪ 1 (z_g ≈ 0.045 in a 2′ disc) gave X ≈ 20–54 from 3–5 galaxies. These sightlines
    drove a spurious "γ < 0.0093" in the first pass. Require E ≥ 1 per shell.
  - Scrambled-column injections recover γ exactly but cannot see regression dilution. Use the Poisson-only λ and the
    thinned-count chain injection. λ on a trimmed sample must trim the simulations too (found by /code-review).
  - Raw `mag_z` with no maskbits/brick_primary/DUP cuts was the first query; use the repository's countmap selection.
- **Next:** DES-SN5YR + DES Y3 Gold (Shah+2024's setup: lensing detected at 6σ, the positive control; deeper
  galaxies give more shells with E ≥ 1); the full Pantheon+ covariance in place of the diagonal errors.

## 2026-10-09: Owner direction: System A/B hypothesis rounds; round 1 (3 worlds, 3 reviews, 3 conditional survivors) (D-069)
- Owner ideas (distance self-consistency, missing light with paired excess, transient connectivity, causal event
  network, Ontological Reset). Distance self-consistency was never tested here; the deficit half of missing light is
  W3/W5.
- System A (no literature): A1 states + transitions; A2 no identity, one quantity "carry"; A3 bridged Markov
  substrate with a fixed pattern of distance disagreement. System B (verified references) failed 14 of 17
  predictions (plus A2's foundations). All three link mechanisms are disordered locality (arXiv:0903.5303) renamed, and achromatic dimming
  is the known cosmic-opacity test.
- **Survivors (conditional):** S1 burst twins (GRB/FRB, wide separation, any delay: untested by 2006.07095 and
  2204.06014); S2 flat-kernel SN residuals (new ~0.015 mag limit reachable); S3 hybrid images (limit only).
- B corrected A1's "blindness theorem": a timing dependence between events ≳ 0.1° apart is provably spacelike at Gpc distances with spectroscopic redshifts (B-est); the
  obstacle is chance coincidences and common causes.
- **Next:** S1 on the Fermi GBM catalogue and CHIME/FRB Catalog 2, with time-scrambled nulls and lensing/duplicates as
  positive controls (docs/hypotheses/round-1/summary.md).
## 2026-10-09: W3 MOA-II six fields null; limits withdrawn after review (D-068, #95)
- Cloud runs streamed gb11, gb16, gb19, gb20, gb21 and gb22 (296,618 Cut-0 light curves; no tar on disk): 234 flags,
  **0 survive** (vetting records for gb20-R-4-0-49379 and gb19-R-4-4-31159; gb11 contact sheet inspected: eclipse
  dips, repeated deficits, unsampled features). New tests: `exotic_in_domain`, `smooth_dip`, `feature_bracketed`,
  `step_ramp`. `derived`.
- **Limits withdrawn:** `/code-review` found that injections were vetted against the lenient 14.2 mag default
  reference flux, while real flags used their own magnitude, so efficiencies were biased high. Injections now use
  I_s; `CHAIN_VERSION` makes `limit`/`combine` refuse old tables. Other fixes: a truncated tar looped forever in
  `read_segment`; a 200 reply skipped the size pin; the Gaia XMatch is retried, and a failure marks flags incomplete.
- **Failed approaches (rules):** a fresh cloud venv without the `mulens` extra skipped the parallax and binary-lens
  tests silently (the stages now refuse to run); `pkill -f` killed the calling shell again (kill by PID).
- **Handoff / next:** a run that died at 03:16 left gb11's injections unfinished; it was taken over at 05:12. Re-run
  `w3_moa.py --field gbN run-field` for the six fields (injections + limit, ~20 min each), then `combine`, in a new
  `[field: W3 MOA-II limits]` PR; then more Nunota Table 1 fields.

## 2026-10-09: W5/W1 Euclid Q1 radial-shear survey: all Deep Fields (60 deg²) null; R calibrated (D-067)
- Hypothesis as D-066 (negative mass shears background galaxies radially). Ordinary radial patterns: PSF-anisotropy
  gradients, blends, tile edges; known clusters must give the opposite (tangential) sign.
- **R = 0.56 ± 0.16** (z_s = 1, scaled by √(χ²/dof); 0.50–0.67 for z_s = 0.8–1.2) from the 4 SZ clusters against
  NFW haloes of their M500: the assumed R = 0.5 is consistent.
- **All 344 Q1 tiles of EDF-F/S/N** fetched by `tileid` (4.78 M rows, 4.60 M of them galaxies; 25 min) and screened at θ_E = 1′, 2′, 4′:
  **0 flags** (field maxima S = 3.3–4.4 against thresholds 4.3–5.1); known clusters 8/9 negative (EDF-S minimum is
  ACT-CL J0405.9-4915). Injection efficiency 0.43–0.58 at 1′ (no limit: below the 0.5 gate), 1.0 at 2′–4′.
  **n₉₅ ≈ 0.049 / 0.051 deg⁻²** at θ_E = 2′ / 4′. With D-063, negative point masses are limited from 2′ to 32′.
- **Failed approaches (rules):** injecting into the observed moments before the PSF deconvolution boosted the
  injected shear ~1.6× over the calibrated R (found by /code-review before merge; the first run's 1′ limit 0.062 deg⁻²
  is withdrawn and the D-066 pilot efficiencies are optimistic): inject where R is calibrated. 16 IRSA threads → 504s,
  no gain; killing a fetch truncated a cache file (writes
  are atomic now); `pgrep -f`/`pkill -f` waiters matched their own shell again (wait on the log file).
- **Next:** per-tile star-ellipticity gradient test (data in `survey.json`); θ_E ≤ 1′ needs fainter or better shapes
  (Euclid DR1); DR1 scale-up when public; θ_E ≈ 1° needs a larger contiguous area than Q1.

## 2026-10-09: J0728+2607 decided with an empirical PSF: no light between the images to F814W ≈ 23 (D-064 addendum)
- The PSF is the median of 6 unsaturated Gaia stars in the same HST F814W cutout. The control lens J2308+3201 is
  detected (S/N 49). J0728's residual is centred on its images, and after the halo correction there is no light
  between them (S/N 1.6). Whole-chain injections are recovered to F814W = 23, ≥ 1.8 mag below any ordinary lens.
  J0130+0725 has 1 usable star, so it has no empirical-PSF result.
- Both D-064 HST pairs therefore have no visible deflector; binary quasars remain the untested ordinary explanation
  (spectra of both images, or time delays, are needed).
- **Failed approach:** the lens-vs-mirror aperture test. The control gives only 3.1σ, because real lens positions
  depart from the SIS flux-ratio rule.
## 2026-10-09: W5/W1 radial-shear screen on Euclid Q1 MER shapes: pilot null, first limits at θ_E = 1–2′ (D-066)
- Hypothesis: a negative-mass lens shears background galaxies radially; ordinary mass gives tangential shear.
- **Convention found wrong in the archive docs:** MER `position_angle` is PA east of north (25 VIS cutouts, median
  |Δ| 0.78°), not THETA_IMAGE "CCW/x". Four SZ clusters in EDF-S then show tangential shear (S = −1.5 to −5.0).
- Pilot: 0.3° discs in EDF-F and EDF-S, θ_E = 30″/1′/2′: **0 flags** (field-max p_random 0.07–0.99); off-grid,
  whole-chain injections (R = 0.5 ASSUMPTION) 0.88–1.0 at 1′, 1.0 at 2′, ≤ 0.18 at 30″. n₉₅ ≈ 8.1 deg⁻² (1′),
  12 deg⁻² (2′).
  The DR10 count floor was 6′.
- **Failed approaches (rules):** unquoted `position_angle` breaks IRSA's ADQL parser (quote it); IRSA returns query
  errors as a VOTable with HTTP 200 (check the body); the EDF-N 0.3° row query hung > 15 min; IBE cutouts come
  gzip-compressed; `pkill -f` on a pattern in your own command kills your shell.
- **Next:** calibrate R on the SZ clusters (shear vs. their M500 NFW prediction); scale to all Q1 (63 deg²; fetch
  in ≤ 0.2° discs in parallel, EDF-N too); star-ellipticity maps for PSF gradients; flag vetting via `/vet-candidate`.

## 2026-10-09: Efficiency rules 8–13 for routine cycles (owner text)
- docs/cloud-routine-prompt.md "EFFICIENCY RULES" gets the owner's rules 8–13: merge main before the final review;
  state-file conflict handling (own CHANGELOG block first, both entries kept); one CI wait per head SHA; verify
  outcomes on GitHub; agents and the owner share the DD-Ching account; stop when the next unit needs a human.
- **Owner action:** paste the updated prompt into the routine (`trig_01PNAmgcfqef8CvhPAY8ggbP`).

## 2026-10-09: D-064 addendum: archival HST of the untestable rejected pairs; J0130+0725 has no lens light to F814W ≈ 23
- Coverage: 2 of the 11 untestable pairs have HST F814W imaging (program 17308), and so do 2 controls.
- Method: a two-PSF fit, halo correction, and the residual flux between the images, with empirical sky-aperture
  errors. Validated on both controls (S/N 42, 16) and on whole-chain injections (recovered to F814W = 23).
- **J0130+0725:** no residual (S/N 0.45). The limit F814W ≈ 23.0 is ≈ 3.5 mag below a typical ordinary lens and
  ≥ 0.7 mag below a 2σ under-luminous one. A binary quasar remains the untested ordinary explanation; it needs
  spectra of both images. Not a candidate.
- **J0728+2607:** inconclusive (residuals at the image cores); it needs an empirical PSF.
- **Failed approaches:**
  - the HAP "combined_skycells" product has a WCS that does not describe its pixels;
  - without the halo correction, PSF mismatch alone gives S/N 18–24;
  - white-noise errors on drizzled pixels overstate S/N by ~3.5;
  - a 60-px cutout cannot hold the pair, because the catalogue position is one image.
- **Next:** spectra or two-epoch flux ratios for J0130+0725 (binary vs lens); an empirical PSF for J0728+2607;
  Euclid DR1 when public for the 9 pairs without HST imaging.
## 2026-10-09: W5 in Euclid Q1: deeper counts cannot open θ_E < 6′; shapes could (D-065)
- Hypothesis: Euclid Q1's deeper counts lower the D-063 floor (blind below θ_E ≈ 6′). Measured (IRSA TAP counts,
  three deep fields): VIS < 24.5 extended galaxies 6.8 × 10⁴ deg⁻², 1.81× DR10, same count slope. With the
  clustering-inflated scatter measured in DR10, the count S/N gains only ×1.0–1.4 (Z scatter at 4′ is 1.6× Poisson
  and clustering does not shrink with depth), on 63 instead of 340 deg². Count screen not built.
- Forecast (model_prediction): a radial-shear test on the same galaxies reaches θ_E ≈ 26″ (1.5–3 θ_E annulus, S/N ≥ 6,
  σ_γ = 0.3, ASSUMPTIONs); the shear sign separates a negative-mass lens from ordinary foreground mass.
- **Failed approach (rule):** IRSA TAP does not index plain RA/Dec ranges (0.25 deg² box > 5 min); use
  `CONTAINS(POINT, CIRCLE)` (40 s for 12 k rows).
- **Next:** Euclid Q1 radial-shear screen: PSF-anisotropy check on stars, synthetic shear injections, trial centres
  on a grid, cross-field null (EDF-N/F/S).

## 2026-10-08: D-056 amendment — the LS pair must be the catalogued pair; one entry per lens (W1/W2)
- `lenscats.pair_match` (moved from D-064's `w12_niq`) now gates the D-056 quasar pair test: a status from an LS image
  pair whose separation differs from the catalogued 2θ_E by > 0.5″ is undecided. Only 115252+004733 changed (LS pair
  4.18″ vs 3.34″ expected; the cutout shows an 18.1 mag lens galaxy with an unrelated faint pair). None of the 15
  pair-decided systems has a catalogued θ_E, so the pipeline checks 0 of 15; a one-off manual check of the 3 also in
  SQLS matches (1.88/2.99/2.03″ vs LS 2.00/3.01/2.01″). For quads (≥ 3 LS images) only a pair wider than 2θ_E + 0.5″
  is a mismatch (a fold/cusp pair is closer than 2θ_E); the check applies only to "deflector"/"none" statuses.
- Same lens listed twice beyond the 3″ merge (MG0414+0534, B2114+022, B2319+052; ~11″ apart): decided systems are
  grouped by designation within 30″ (ASSUMPTION; decimal-degree names give no key) over all covered sensitive
  systems, and a deflector at any copy (decided or not) explains the lens. Each had counted once as
  "deflector" and once as "none".
- New limits (typical, `derived`): quasar 3/15 < 0.52, radio 0/10 < 0.30, all **3/25 < 0.31**; conservative 0/5 < 0.60
  (were 0.48 / 0.23 / 0.27 / 0.50). No new unexplained system; the 3 CHITAH pairs stay open as before.
- **Rules:** a pair test must check that its pair is the catalogued one; name-match decided systems before counting N.

## 2026-10-08: Efficiency rules for routine cycles (owner text)
- docs/cloud-routine-prompt.md gets the owner's "EFFICIENCY RULES" after "MOVE FAST, SAFELY": result first, a review
  stopping rule (fix only result/provenance/reproducibility/guarded-file findings, list the rest, ≤ 3 rounds), verify
  the branch after forked skills (they can leave HEAD detached, which stranded two commits on #98), fail fast on
  data access, calibrate before flagging, a 35-minute time box, one unit per cycle.
- **Owner action:** paste the updated prompt into the routine (`trig_01PNAmgcfqef8CvhPAY8ggbP`).

## 2026-10-08: W5 count deficits in Legacy Surveys DR10: 340.5 deg² null, first W5 limit (D-063, #94)
- Two sessions (a cloud run's 10° pilot, then a worktree worker on the full 20° × 10° regions; see the coordination
  entry below). DR10 Tractor galaxies (r < 23.5) counted per `nest4096` HEALPix pixel on the Data Lab server
  (~170 MB, not a catalogue download); predicted deficit profile from `exotic_sim.count_ratio` with the measured
  number-count slope (`model_prediction`): ~10 % of galaxies inside θ_E missing, ratio 0.02 at θ/θ_E = 0.1.
- Screen: matched filter at θ_E = 2–32′ over desA (RA 20–40°) and desB (RA 50–70°), Dec −30 to −20; cross-region
  null calibrated on vetted peaks (galaxy clustering makes Poisson errors 1.1–5.8× too small). 247,361 peaks →
  **40 flags → 0 survivors** (NGC 1398 sky over-subtraction ×3, bright stars, depth/tile edges; the rest consistent
  with the null). Contact sheet inspected by the worker and the coordinator. `derived`.
- 5,232 injections through screen + vetting: efficiency 0.50–0.73 for θ_E = 8–32′, 0 below 4′. **95 % limit on the
  sky density of W5 lenses n₉₅ ≈ 0.012–0.018 deg⁻² at θ_E = 8–32′** (≈ 3 × 10¹⁰–5 × 10¹¹ M☉ at 1 kpc, 3 × 10¹³–5 ×
  10¹⁴ M☉ at 1 Mpc; geometry ASSUMPTION). Untestable below θ_E ≈ 6′ with DR10 counts.
- **Failed approaches (rules):** Data Lab ADQL rejects sub-selects, CASE, SIGN and GROUP BY on expressions (group by
  `nest4096` only); RA/Dec chunks split edge pixels — sum them, or every chunk border looks like a deficit; an
  unbounded bright-end count slope makes the profile blow up at x → 0 (fix 0.6); a null from the other region's raw
  peaks inherits its artefacts (vet them first); veto radii growing with θ_E removed 90 % of random positions at 32′
  (veto only mimics that can empty ≥ 10 % of the core); HyperLEDA returns sexagesimal unless `_RAJ2000`/`_DEJ2000`
  are requested; a large-galaxy veto is needed (sky over-subtraction around NGC 1398); `pkill -f` killed the calling
  shell; JSON writers must end the file with a newline (pre-commit end-of-file-fixer failed CI twice).
- **Next:** deeper counts (HSC, Euclid) for θ_E < 6′; larger contiguous area for θ_E ≈ 1°; review leftovers: number
  counts divided by the full box area (small bias), per-job rebuild of region state in `inject`.

## 2026-10-08: W1/W2 in rejected lensed-quasar pairs: LS DR10 cannot decide them (D-064)
- Hypothesis: a dark deflector hides among pairs rejected as lenses for lack of a lens galaxy (Lemon 2023 UQP/QSO
  pair; SQLS "no lensing object", "QSO pair", "binary"). 123 rejected and 106 control lenses (≤ 3″) went through the
  D-056 chain.
- **Control efficiency 0 / 5** (1.9–2.6″; lens light blended into the images, fitted as PSFs): an LS "none" at these
  separations carries no information. Rejected: 24 "none". Of these, 10 are colour-mismatched, 2 are catalogued
  binaries and 1 has two redshifts; 11 are untestable. No limit, no candidate.
- **Failed approaches:** several input-format traps (sexagesimal Hennawi coordinates, SQLS two-row pairs, Lemon
  `z2` semantics, time-stamped VizieR headers) and sample-definition bugs were caught in review and fixed; the list is
  in D-064 "Evidence".
- **Next:** HST/Euclid/HSC image models or spectra for the 11 colour-matched pairs; the LaTeX-only NIQ tables.

## 2026-10-08: W3 MOA-II: calibrated baseline test, LF-drawn injections, streaming pipeline; gb22 re-run null (D-068)
- Cloud runs (session that opened #95, taken over at 23:12 UTC after 35 min idle to land it). The variable-baseline
  threshold is now the 95th percentile of the field's quiet χ²/dof (5.31 in gb22; ASSUMPTION), injections are drawn
  from the luminosity function, 200 per cell.
- gb22: same 30 flags as D-062 (event IDs identical; contact sheet inspected in D-062), **0 survivors**.
  **Γ₉₅ ≈ 1.1–3.3 × 10⁻⁶ per star per year for t_E = 10–300 d, 4–7 × 10⁻⁶ at 3 d** (103 / 2,000 recovered,
  0 / 400 PSPL controls; `derived`, N_s ASSUMPTION). Supersedes D-062's numbers.
- Streaming per-field pipeline (`moa_stream`, `w3_moa.py --field`): no tar on disk; gb22 by HTTP in 81 s;
  258 light curves/s on gb21 (test only, no gb21 result yet).
- **Failed approach:** a parent process reading ranges for its workers was OOM-killed on gb21.
- **Handoff / next:** `w3_moa.py --field gbN run-field` for the Nunota et al. 2024 Table 1 fields, smallest first
  (gb21 12.3 GB, gb20 15.2 GB, gb19 16.3 GB …; sizes in docs/exotic_limits.md), then `combine`.

## 2026-10-08: Coordination and dispatch rules for concurrent sessions (owner text)
- Two sessions worked PR #94 (W5) at once: a cloud run started from its 20:39 skeleton, another session pushed the full
  screen at 21:14, and the run's push was rejected (its pilot is parked on `claude/w5-clustered-null`, findings on #94).
  Commit-age heuristics cannot see a session that is coding but has not pushed.
- docs/cloud-routine-prompt.md gets the owner's "COORDINATION AND DISPATCH" section (dispatch first; `claimed` label +
  claim comment with a 10-minute heartbeat, stale after 20; re-check before every push, `-alt` branch on collision;
  "D-069" until merge). It replaces the 15-minute commit-age rule; docs/operations.md and its label table follow.
- **Owner action:** paste the new prompt into the routine (`trig_01PNAmgcfqef8CvhPAY8ggbP`); this file is only the copy.

## 2026-10-08: Governance — owner decisions on scope, parallelism, cloud disk and layout
- Owner direction (2026-10-08, given by the owner in the session that opened this PR) with a one-time
  authorisation to edit CLAUDE.md, including "Owner decisions", in this single PR and to merge it once CI is green
  and `/code-review` findings are fixed: four owner-decision bullets added verbatim (scope: any public dataset; smart
  parallelism and `[field: <unit>]` claim PRs; cloud disk streaming; layout). "Prefer single-thread work" removed
  from Budget. Later edits to "Owner decisions" again need the owner.
- Layout: root CLAUDE.md keeps owner decisions, merge policy, budget and pointers; module rules moved to
  `src/jwst_anomaly/CLAUDE.md` (stage contracts, signature layer), `scripts/CLAUDE.md` (screen / vet / inject /
  limit conventions, parallel topology) and `data/manifests/CLAUDE.md` (data, manifests, cloud disk). New
  directory-scoped skills `scripts/.claude/skills/w3-survey` and `w12-lenscats` carry the recipes and every
  failed-approach rule recorded for those searches; a `w5-counts` skill follows when W5 (#94) lands.
- Inventory before the change: one CLAUDE.md (root, 97 lines) and three root skills (`research-cycle`,
  `reuse-check`, `vet-candidate`), kept at the root; `research-cycle` now points to the new files and the claim rule.

## 2026-10-08: W3 MOA-II gb22 limit corrected: injections re-run after the `s_min` review fix (D-062)
- #92's limit (Γ₉₅ ≈ 0.5–5 × 10⁻⁶, 171 / 600 recovered) came from injections run before the review fix that made
  the pre-screen `s_min` the weaker of the two significances; the merged code already has the fix. Re-run of the whole
  chain on that code: **101 / 600** W3 injections recovered, 0 / 100 PSPL controls; **Γ₉₅ ≈ 1.2–4.4 × 10⁻⁶ per star per
  year for t_E = 10–300 d, 0.6–1.1 × 10⁻⁵ at 3 d** (`derived`; N_s ASSUMPTION unchanged). The real-data result is
  unchanged and reproduced with main's code: 30 flags, 0 survivors. Tracked fit table regenerated.
- Largest efficiency loss: the variable-baseline test (73 of 188 vetted injections; 35 % of quiet carriers have
  χ²/dof > 2 alone). **Rule:** re-run injections after any pre-screen change before quoting a limit.

## 2026-10-08: W3 in MOA-II Cut-0 light curves (gb22): 30 flags, 0 survivors, first W3 rate limit (D-062)
- Cloud runs (two sessions merged on `claude/w3-moa`). Hypothesis: the MOA-II 9-year release publishes every Cut-0
  object (positive *or negative* difference-image detections) before any bump cut, so unlike OGLE (0 / 600) and
  Gaia (0 / 240) its selection should keep W3 dimming events. `gb22.tar` (3.5 GB, smallest field, sha256 pinned).
- Pre-screen on all 18,599 light curves: 1,058 shape passes, **30** after the shared-epoch tests (1,022 removed at
  field-wide shared epochs, 6 at chip-wide ones); 0 / 100 PSPL controls pass. 30 flags (ΔBIC < −10), **0 survive
  vetting** (funnel and the three eclipse-like dips in docs/exotic_limits.md "W3 in MOA-II"). `derived`.
- Injection-recovery through Cut-0 emulation, pre-screen, fit and vetting (600 W3, 100 PSPL): recovered 8–53 % of
  injections per cell. **Γ₉₅ ≈ 0.5–5 × 10⁻⁶ per monitored star per year** for t_E = 3–300 d, strongest at ~30 d
  (N_s = 3.5 × 10⁶ is an ASSUMPTION from Nunota et al.'s N_s-per-object ratio; 1.43× weaker at its low end).
- **Failed approaches:** fitting all 1,058 shape passes (chunk 1/8 flagged 131 / 133: dip-shaped variables at
  shared epochs prefer negative-flux models; superseded by the shared-epoch cut, table kept as
  `fits_gb22_chunk1of8.ecsv.gz`). The first vetting chain left 9 survivors; the contact sheet showed one- or
  two-night drops and flat dips without caustic spikes, so the eclipse model, chip-level shared-epoch test and
  ≥ 3-night rule were added (their cost to W3 is in the injection numbers). Waiting with `pgrep -f`/`pkill -f`
  matches the waiting shell itself; wait on a log line or output file.
- Reproduced in a later cloud run from the tracked fit table (`prescreen`, `merge-chunks --n 1`, `vet`, `sheet`).
- Final `/code-review` fixes: Pool workers get the MOA fit bounds and shared-epoch populations through an
  initializer (spawned workers on Windows inherited neither, so the shared-epoch cut was silently skipped there);
  `fit --chunk` writes its own file instead of replacing the table `vet` reads; `merge-chunks` also checks the
  pre-screen `Params`; the sheet labels a PAR winner's PSPL curve as such.
- **Next:** a field with a published N_s (Nunota et al. 2024 Table 1) to make the limit model-light; the Cut-0
  spurious-detection filters are not emulated, so ε (and the limit) are optimistic by an unknown factor.

## 2026-10-08: W3 in the Gaia DR3 microlensing candidates: the published selection rejects W3 too (D-061)
- Cloud run. Hypothesis: the Gaia DR3 candidates (363; Wyrzykowski et al. 2023) are selected less PSPL-shaped than
  Mróz et al., so they could limit W3. Tested the selection before fitting, as D-057 requires.
- New `gaia_mulens.GaiaDR3Microlensing` adapter (TAP + DataLink epoch photometry, 8 parallel batches of 12 ids,
  ~3 s per source on the server; Table D.1 sample labels from the pinned arXiv source) and `scripts/w3_gaia.py`.
- Emulated Sample A cuts pass 143 / 163 real Sample A events. W3 injections: **2 / 240 selected, 0 / 240 selected
  and flagged**; PSPL controls 17 / 120. W3 dimming fails skewness < 0 and the skew–Abbe cut, as designed for
  brightenings. The fits of all 363 give one flag (4053892503992268288, ΔBIC −40.3). It is an event truncated at the
  window end on a variable baseline (light curve inspected), not a candidate.
- **Failed approaches (rules):** the Extractor cuts (n points, duration > 135 d, max σ > 50) cannot be emulated from
  the paper. The guessed definitions fail 126 / 163 real events, so they are left out. A single-id DataLink request
  answers with bare CSV, not a zip. Truncated chunked replies happen, so retry them.
- **Next:** W3 needs light curves taken before any microlensing selection. Gaia DR3 has epoch photometry only for
  its variable sources (vari_summary; ~11.7 M per the DR3 release, recheck), selected by variability, not shape. Next, check whether W3 survives that
  classifier (`vari_classifier_result`) by injection, and whether a sky-region subset is small enough to screen.
  Other options: KMTNet public seasons; OGLE EWS (owner decision, terms).

## 2026-10-08: W2 deflector test at HST resolution from the Hubble Source Catalog: not decisive (D-060)
- Cloud run. Hypothesis: HST resolution decides the lensed quasars that Legacy Surveys left blended or too close
  (D-056). `scripts/w12_hsc_probe.py`: HSC v3 summary sources within 4″ of each of the 444 galaxy-scale quasar/radio
  systems (MAST catalogs API, 58 s); ≥ 2 point sources (CI < 1.3) are the images, an extended source (CI ≥ 1.5)
  near their centroid and > 0.2″ from an image is the deflector; sources in < 2 HSC images are dropped as likely
  artifacts (MAST's recommendation; #90 review) (ASSUMPTIONs in `Params`).
- 71 of 444 systems have HSC sources in ≥ 2 HSC images (91 with any). Validation on systems with a published lens redshift (a lens galaxy is known
  to exist): 13 deflector, 15 none, 38 undecided → **efficiency 13/28 = 0.46** (`derived`; 15/35 = 0.43 without the
  artifact cut). Misses include quads and doubles (13 of 16 "none" have 2 point images). In four inspected misses (H1413+117, HE1104−1805, SBS0909+532,
  HE2149−2745) the HSC rows within 4″ are only the quasar images: the lens galaxy is absent from the catalogue, not
  mis-typed (likely lost in the quasar PSF; hypothesis, no cutouts inspected). Without a lens redshift: 1 none
  (HS0810+2554), 4 undecided. The 3 open SuGOHI IX pairs (D-056) have no HSC sources.
- **Failed approach (rule):** HST *catalogue* photometry cannot decide a dark deflector in lensed quasars — a
  "none" is more likely a missed lens (0.54) than a dark one. No limit, no candidate; HS0810+2554 is not flagged.
- **Next:** PSF-subtracted HST image modelling (e.g. drizzled frames from MAST, quasar PSF + Sérsic fit) validated
  on the same known-lens set, or HSC PDR3 photometry; until then W2 in wide imaging stays at D-056.

## 2026-10-08: Warp literature check: still nothing testable
- Subagent search (arXiv API 2023–2026, INSPIRE citations of Clough et al. 2024): no imaging or lensing prediction
  for a distant observer; Lentz & Felton 2024 give order-of-magnitude EM fluxes for a bubble 100 lyr away but no
  template that separates it from ordinary transients (found by review); the collapse-burst waveform is not public; an O3 search for superluminal-source GW
  bursts (Kuwahara & Cannon 2023) is already null. D-047 stands; recorded in SOURCES.md and docs/exotic_lensing.md.

## 2026-10-08: W3 OGLE bulge (all 5,790 events): no candidate; the published selection rejects every W3 event (D-057)
- Worktree worker, unbounded π_E (predates D-058; bounding can only add flags). Bulge: 0 fit failures; best ordinary
  PSPL 5,377 / PAR 401 / FSPL 12; ΔBIC(best exotic) 5/50/95 % = −3.6 / 5.7 / 12.0; 127 flags < −10 (`derived`).
- Vetting, cumulative: 127 → 113 refit all ordinary → 80 robust errors → 73 variable baseline → 14 season offsets →
  9 season drifts → 7 binary source/lens, VSX/Gaia, arXiv → 1 `feature_coverage` → 1 jackknife → **0** two unrelated
  events (new `revet` tests: ≥ 3 epochs where the models differ by > 3σ; drop up to 3 influential epochs keeping ≥ 3 in
  the feature; two independent PSPL bumps). BLG519.21.110304's exotic spikes sat on the 2011 event and a 1-day bump in
  2015 (two PSPL bumps better by ΔBIC 24.4). Disk: 6 → 0.
  Contact sheet of the 7 late survivors inspected. The D-059 chunk survivors fail `feature_coverage` here.
- Injections: 600 W3 events (n = 1, ε < 0; t_E 3–300 d, ρ 0.01/0.1) on real bulge cadences + 300 PSPL controls. The
  fitter flags 42–97 %, but **0 / 600 pass the emulated Mróz selection** (controls 15–43 %; cuts failed most: one
  bump, PSPL χ², χ₃₊). The emulation passes 63.9 % of the real selected events (somewhat stricter). No rate limit is
  derivable; ε_W3/ε_PSPL < 0.12–0.33 (95 %). Wall time ≈ 5 CPU h (bulge fit 10.4 ks, vet 6.6 ks, inject 2.0 ks).
- **Failed approaches (rules):** exotic fits started only from the PSPL solution miss the W3 geometry — start with
  caustic spikes on pairs of maxima and on absolute t_E; ΔBIC alone is not a candidate test (season blends and
  feature coverage remove 120 of 127); a jackknife must keep ≥ 3 epochs in the exotic feature or it kills real short
  events; an exotic fit whose spikes sit on two bumps years apart needs the two-unrelated-events test; pin BLAS threads (`OMP_NUM_THREADS=1`) with process pools; a published
  PSPL-selected sample cannot constrain a non-PSPL signal — inject through its selection before fitting it.

## 2026-10-08: W3 OGLE bulge, chunks 2–3/12: 20 flags, no candidate; chunk tables tracked (D-059)
- Cloud run. Chunk 1/12's fit table lived only in the ephemeral session and is lost. `fit --chunk K/N` now also
  writes a deterministic gzipped copy to `results/w3_ogle/` (~0.22 MB per chunk); `merge-chunks --n 12` joins
  chunks 1..12 into the table `vet` reads and marks it the whole sample only when every chunk is present, fitted
  with the current `Params` and holds exactly its own events (none skipped) (D-059). A chunk takes ~9 min on 4 cores (not ~17).
- Chunk 2/12 (483 events, 0 errors): best ordinary PSPL 454 / PAR 28 / FSPL 1; ΔBIC(min exotic) 5/25/50/75/95 % =
  −4.7 / 3.8 / 5.8 / 6.6 / 11.2; 49 below 0, 9 flags below −10 (`derived`). `vet`: 8 fail; **BLG624.18.69573**
  (no EWS name; t_E ≈ 180–240 d, best ordinary PAR) passes every automated test (N1neg ΔBIC −19.0 vs PAR; season
  offsets −15.0, drifts −15.5; binary source / lens −10.1; 0 VSX / Gaia matches). Contact sheet and residuals
  inspected: the N1neg model puts its first caustic spike inside a season gap (no data at t − t0 ∈ [−50, 0] d; the
  models differ by > 5 % over [−89, −7] d); its Δχ² comes from 2 peak points (−9.4) and 16 post-peak points (−6.7).
  An unsampled caustic plus a sparsely sampled peak is not evidence: not a candidate (ASSUMPTION-level judgement).
- Chunk 3/12 (483 events, 0 errors): PSPL 447 / PAR 36 / FSPL 0; ΔBIC(min) 5/50/95 % = −2.6 / 5.7 / 11.1; 39 below
  0, 11 flags; **0 survive** `vet` (season offsets/drifts remove 10, the refit of all ordinary models the 11th; one
  VSX match). Contact sheet inspected; in several flags the exotic and ordinary curves also differ mainly in gaps
  (e.g. BLG597.28.9837 has no peak data).
- **Failed approach:** chaining chunk runs with `while pgrep -f 'chunk 2/12'` — the waiting shell matches its own
  pattern and never starts the next chunk.
- **Next:** a `gap_coverage` vetting test (require data where the exotic and best ordinary models differ, else the
  flag fails); refit chunk 1 and fit chunks 4–12 (two or three per run); `merge-chunks`, `vet`, `sheet`; the
  empirical ΔBIC null and xallarap fit for BLG667.04.62161 (scratch null-simulation design: PSPL best fit plus white
  noise, and plus season-wise circularly shifted residuals); then `inject` / `limit`.

## 2026-10-08: W3 OGLE bulge, chunk 1/12 (483 events): one marginal flag survives automated vetting
- Cloud run. `w3_microlensing.py fit --chunk K/N` fits events K−1, K−1+N, … so sessions fit disjoint, field-balanced
  parts of the bulge sample (measured ~10 s CPU per bulge event, ~4 h for all 5,790 on 4 cores; one chunk of 12 is
  ~17 min wall). Chunk 1/12: 483 events, 0 errors; best ordinary PSPL 449 / PAR 32 / FSPL 2; ΔBIC(min exotic)
  5/25/50/75/95 % = −3.5 / 3.7 / 5.6 / 6.5 / 12.4; 41 below 0, 12 flags below −10 (`derived`).
- `vet`: 11 of 12 fail; **BLG667.04.62161 (OGLE-2015-BLG-1250) passes every automated test**: E2neg ΔBIC −12.1 vs
  PSPL; season offsets −12.7, season drifts −12.4, binary source / binary lens no better than PSPL, 0 outliers,
  baseline χ²/dof 0.87, 0 arXiv records, no VSX / Gaia variable within 1″. Contact sheet inspected: faint source
  (peak flux ~2× baseline), large scatter, a flattened peak and a few low points ~+20…+50 d; the E2neg plateau fits
  those. ASSUMPTION-level judgement: marginal, not a candidate — with 483 trials and 12 flags below −10 from a
  heavy-tailed ΔBIC distribution, one −12 survivor is expected without any exotic lens; untested ordinary
  explanations: xallarap, per-season error underestimation near the peak, blending/difference-imaging systematics
  of a faint source, and a calibrated null (the injection/limit stage).
- **Next:** chunks 2–12 (one per run: `fit --chunk K/12`, then `vet` / `sheet`); for BLG667.04.62161, an
  empirical ΔBIC null from the same chunk's PSPL-simulated light curves (does −12 occur at rate ≥ 1/483?) and an
  xallarap fit before any further attention; then `inject` / `limit`. `limit` refuses vetting of a single chunk
  (#86 Codex), so a `merge-chunks` step (concatenate `fits_bulge2019_chunk*of12.ecsv`, then `vet`) comes first.

## 2026-10-08: W3 OGLE disk re-fit with bounded parallax: still null (D-058)
- Cloud run. ASSUMPTION `Params.pie_max = 5`: PAR fits with |π_E| > 5 are rejected in every fit and vetting refit
  (the unbounded fits reached π_E ~ 30–1,400). Disk (460 events, 376 s): best ordinary PSPL 408 / PAR 52 / FSPL 0;
  36 of 368 PAR fits sit on the bound; ΔBIC(min exotic) 5/50/95 % = −4.0 / 3.6 / 6.9.
- 7 flags (was 6; new GD1217.10.8703 at ΔBIC −10.0, `E2pos`), **0 survive** `vet` (all tests complete). The new flag
  is four post-peak points 0.15–0.35 below baseline that neither model fits; per-season offsets remove it (ΔBIC 0.9).
  Contact sheet inspected.
- Absorption check (simulated: 40 E2pos/N1neg injections, 700 uniform epochs, white noise σ = 0.05–0.3): flags 25
  bounded vs 25 unbounded; ΔBIC shifts ≤ 0.46. On this cadence unbounded parallax does not absorb exotic signals,
  so the bound is a physical prior, not a sensitivity gain. Not tested: seasonal gaps and correlated systematics.
- Environment: install the `mulens` extra (`-e ".[dev,cloud,mulens]"`); without MulensModel `fit` silently skips PAR.
- Follow-up (#84 Codex): checkpoint rows carry a `Params` hash; `fit` refits rows from other Params, so a chunked
  bulge run never mixes bounded and unbounded parallax fits.
- **Next:** bulge `fit` in chunks or locally, then `vet` / `sheet`; `inject` / `limit` with season-drift vetting.

## 2026-10-08: W1/W2 in published lens catalogues: no dark deflector; weak limits (D-056)
- Worktree worker: lenscat (32,838), Euclid Q1 Discovery Engine (2,584) and SuGOHI (3,961) merged into 35,862 systems
  (`src/jwst_anomaly/lenscats.py`, `scripts/w12_lenscats.py`; pinned by sha256). Deflector test in Legacy Surveys
  DR10 Tractor; footprint and depth from the DR10 brick summary, not from detected sources.
- 17,555 galaxy-scale systems in the footprint. Galaxy-finder (17,102) and sub-mm (110) systems cannot show a dark lens
  and give no limit; only lensed-quasar (325) and radio-interferometric (18) systems are tested (pair test: two point
  images, nothing bright enough near the expected deflector). Blended, too-close and faint-galaxy cases are undecided.
- 29 decided (13 with a deflector, 16 without); 13 of the 16 have a literature lens galaxy, and 3 SuGOHI IX CHITAH
  pairs (090434−005328, 091517+040747, 104122−005618) are open only in the typical variant — a lens below the local
  LS depth explains them conservatively. Not candidates. Cutout sheets inspected.
- 95 % limits on the dark-deflector fraction, test completeness assumed (not measured): typical f_dark < 0.48
  (quasar, k = 3, N = 16), < 0.23 (radio, k = 0, N = 13), < 0.27 (all, k = 3, N = 29); conservative < 0.50 (quasar,
  N = 6). Earlier 1.5 × 10⁻⁴ / 8.2 × 10⁻³ / 0.13 are withdrawn (PR #81 reviews). No W1 geometry (no image positions).
- **Failed approaches (rules):** a limit is valid only over systems where the test could have found the signal —
  prove it by injection; a recovery factor from deleting deflectors and re-running the same code is 1 by
  construction; coverage must come from footprint/depth products, never from "a source nearby"; evaluate exclusion
  flags on every system; Data Lab TAP takes no table uploads or q3c (batch box ORs); lenscat types cluster-survey
  entries as "galaxy", has AGEL declination and SPT position errors and rounded positions, and keeps rejected
  candidates; "no lens redshift" ≠ "no lens".

## 2026-10-08: W3 OGLE-IV disk sample: null (D-057)
- Cloud run. `jwst_anomaly.ogle` (Mróz et al. 2019/2020 adapter) and `scripts/w3_microlensing.py`
  (`fit` / `vet` / `sheet` / `inject` / `audit` / `limit` / `manifest` / `summary`). Manifest
  `data/manifests/ogle_mroz.ecsv`. `www.astrouw.edu.pl` reachable from the cloud (2026-10-08).
- **Disk (Mróz 2020, all 460 Table B1 events; 525 s on 4 cores):** best ordinary PSPL 406 / PAR 54 / FSPL 0.
  ΔBIC(min exotic) quantiles 5/50/95 % = −3.9 / 3.6 / 6.9; 6 flags below −10, all `E2pos`. `vet`: 0 survivors;
  every flag loses the exotic preference under per-season baseline offsets and/or drifts (ΔBIC −4.6 … +8.0), and
  BLG568.12.9169 is already −5.1 on refit. arXiv 0 mentions, no VSX / Gaia DR3 variable within 1″. Contact sheet
  inspected: two flags have sparse peak coverage; two have post-peak points below baseline (the Ellis
  demagnification shape) that one season's drift absorbs.
- **Limits:** none yet (injections use bulge light curves); the disk null is a flag count, not a rate limit.
- Review follow-up (#82 Codex): binary-lens α starts were passed in radians to MulensModel (degrees), so only
  0.5–5.8° was searched; fixed and disk re-vetted: still 0 survivors (BL BICs move by ≤ 18). Injection vetting now
  refits PAR with season trends as `vet` does; failed XMatch queries leave a flag unvetted; `limit` refuses a
  zero-event limit unless the bulge vetting is a complete null.
- **Caveats:** unbounded parallax fits reach π_E ~ 30–1,400 (unphysical); they can absorb an exotic signal and cut
  sensitivity. The season-drift test may also absorb real W3 dips: calibrate both with injections.
- Timing: the bulge fit is ~1.1 s/event on 4 cores (~1.8 h for 5,790); `fit` checkpoints to
  `fits_<key>.partial.jsonl`, but cloud disks are ephemeral: fit in chunks per run or locally.
- **Next:** bound π_E (ASSUMPTION, e.g. |π_E| < 5) and re-fit disk; bulge `fit` (local or chunked), `vet`, `sheet`;
  `inject` / `limit` with season-drift vetting inside the injection loop.

## 2026-10-08: W1 shear (aperture-mass) screen: four clusters null; limits 5–8× stronger than radial (D-053)
- Cloud run. `exotic_screens.py shear` builds D-050: PSF-deconvolved catalogue ε, cluster shear removed, spike
  segments vetoed, Schirmer 10″ aperture-mass S/N map, rotation null, B-mode check. `scripts/inject_shear.py` reuses
  D-049's W1 painting.
- **Validation on known signals.** Measured ε along the cluster model's g rises with |g| (MACS0416, Abell 2744), with
  responsivity R = 0.41–0.48, not 1. The screen now removes R g; injected images keep R of their lens shear.
- **Real fields** (Abell 2744, MACS0416, MACS1149, Abell 370; 30.2 arcmin²): **null**. The first run (no spike veto)
  had Abell 370 E 4.51 / B 4.68 and Abell 2744 E 3.92, all p ≤ 0.01 against the rotation null. Diffraction spikes
  were the cause: the veto (77–178 segments per field) brought them to 3.48 / 3.46 and 3.39. B-mode extremes still
  beat the rotation null in three fields, so an E peak must also beat the field's max |S_×| (post hoc,
  conservative). None passes both.
- Injections (200 per field and mass; spike rows never painted): 33 / 197 / 405 of 800 at 2 × 10¹² / 8 × 10¹² / 2 × 10¹³ M☉ (radial 8 / 84 /
  156 of 1,600); 1 at 2 × 10¹¹, 0 at 2 × 10¹⁰. 95 % limits 7.7 × 10³ / 1.3 × 10³ / 6.5 × 10² deg⁻² (radial 6.1 × 10⁴ /
  7.0 × 10³ / 4.0 × 10³). docs/exotic_limits.md.
- **Failed approaches (rules):** subtracting the full model g from catalogue moments (leaves a W1-signed radial
  residual); any radial-alignment statistic without the spike veto; trusting the rotation null without a B-mode
  check; point-mass 1/x² and top-hat filters (lost to Schirmer 10″ in injections).
- Data: the DJA SMACS v7.4 and El Gordo v7.0 photo-z tarballs return 404 (2026-10-08), so these fields are out. MAST
  catalogues and CANUCS DR1 catalogues (~30 MB each) downloaded fine.
- Wall time: 58–168 s per field for 1,000 injections; the base screen takes 6–13 s. Review rounds changed the injection model (spike and
  near-core and unresolved rows not painted, lens change from raw moments); every re-run kept the four real fields
  null. The counts moved from 59/264/459 (first model) to 33/197/405 (final); not painting unresolved rows was the
  largest step (from 50/238/449; docs/exotic_limits.md caveats).
- Final review: painted images now face the spike veto too (409 → 405 at 2 × 10¹³ M☉; limits unchanged at two
  digits). `psf_sigma_px` (now shared with `inject_radial.py`) raises when no S/N > 50 row has a size, where it
  returned NaN. Left as is (maintainability only): the E/B summary is computed in both `cmd_shear` and
  `inject_shear.run_field`; the null recomputes the |e|² noise term per draw.
- **Handoff:** TASKS "JWST focus" 1 W1: lower the B-mode floor (PSF anisotropy, blends), more fields, stacking.
## 2026-10-08: Orphan-pair null (e) fixed, companion-aware null (f): no deep-field excess (D-055)
- Cloud run. Hypothesis: the D-051 flanking-field orphan excess (246 vs 211, P = 0.010) comes from null (e)'s
  colour cell (member i only, non-finite colours in the 0–0.3 bin) or from physical companions the 10–30″
  reference misses.
- `pair_cells` uses both members' colour bins (unordered), a non-finite bin and a ≥ 5-pair cell floor (S/N × size
  fallback); null (f) conditions on 3–6″ pairs. Nine fields re-run (≈ 1 min each, in parallel), orphans unchanged.
  Deep fields: (e) P = 0.13, (f) P = 0.35; clusters P ≥ 0.13. **The excess was a null-model artefact.** Table:
  docs/orphan_pairs.md.
- Injections re-run (6 fields in parallel, ~4 min): no-candidate limits reproduce D-051; background-aware limits
  tighten by ×0.76 (docs/exotic_limits.md). Orphan set unchanged, so D-051's inspected contact sheets stand.
- Not separated: the shares of the cell changes.
- **Handoff:** W2 next is the segmentation-map same_galaxy rule (TASKS "Now" 1); PR #78 (W1 shear) in flight.
## 2026-10-08: Survey-agnostic signature layer; any public dataset in scope (D-054)
- Owner direction: find evidence of traversable wormholes / negative-mass objects or warp-drive spacetimes in any
  public dataset. `jwst_anomaly.signatures` registers W1/W2/W3/W5 (prediction, injection, screens, ordinary
  mimics, limits) and defines the `LightCurveSurvey` / `CatalogueSurvey` adapters.
- Reuse-check: MulensModel (extra `mulens`) for ordinary microlensing fits; OGLE-IV Mróz et al. 2019/2020 samples
  first for W3 (with published efficiencies), then Gaia DR3 `vari_microlensing`, then KMTNet; lens catalogues
  (lenscat, Euclid Q1, SuGOHI) for W1/W2. No published survey light-curve search for negative-mass or Ellis
  events exists, so a W3 limit would be new (owner first). OGLE EWS use waits for the owner (its terms).
- **Handoff:** W3 re-analysis of the Mróz OGLE-IV samples (TASKS "Now" 1).

## 2026-10-08: W3 multi-epoch dimming / inverted-microlensing screen: null, limits (D-052)
- Worktree worker: `scripts/dimming_screen.py` cross-matches per-epoch level-3 catalogues (F200W + F444W), flags
  vanishing, achromatic-dimming and rise-dip-rise sources, vetoes catalogue effects (incl. a `bright_neighbour`
  veto: a source ≥ 100× brighter within 1.5″), and confirms with S3 byte-range forced photometry. Fields:
  NEXUS-Center (8 epochs), MACS0416 (5; F200W baseline 0.35 yr, F444W 3.26 yr), Abell 2744 (3). El Gordo rejected
  (one shared band); JADES deferred.
- 5,177 catalogue flags → 1,441 after catalogue tests → 6 confirmed by forced photometry (MACS0416) → 0 after the
  bright-neighbour veto (all six sit on or beside saturated stars; cutout sheet inspected: the star's wings and
  spikes rotate through the aperture with each epoch's mosaic orientation). No candidate.
- Injection (`exotic_sim.inject_light_curve`, n = 1, ε < 0) with per-copy vetoes: efficiency 0.04–0.24. Headline
  95 % limits (MACS0416 only, the one field with calibrated forced errors, 34 sources): rate < 0.056 per source
  per yr at t_E = 1 yr (umbra fraction τ < 0.18); all fields indicative < 0.015 (τ < 0.049). The per-copy
  single-epoch veto makes full vanishes of two-epoch sources undetectable (main efficiency loss).
- **Failed approaches (rules):** a catalogue non-detection is zero flux only if the source would have been ≥ 10σ;
  vanish only when every testable band vanishes; `is_extended == False` is not a point-source cut (use CI_70_30);
  inject multiplicatively with matched noise; Gaia masks must drop the source's own match before looking for a
  neighbour; forced fluxes need each epoch's pixel solid angle; never credit exposure to sources flagged before
  injection; don't run parallel or multi-target S3 cutout jobs through the proxy (s3fs "bucket does not exist").
- Caveats: forced counts are from the 06:39–07:03 UTC run (later S3 failures), re-calibrated offline; NEXUS and
  Abell 2744 forced errors uncalibrated; PEARLS (MAST lists 2.0.1, headers 3.0.0) + CANUCS processing differ.
- Wall time: screen 4–76 s, injections 20–136 s, forced photometry 3–32 min per field.
- Final /code-review fixes (cloud run): the forced-stage SIMBAD/NED label read non-existent `*_otype` columns
  (always empty; now `best_match_*` of `crossmatch.XMATCH_COLUMNS`); `calibrated` is False when an epoch image was
  unreadable; the dead saturated-star branch of `bright_neighbour` removed (`near_star` covers it). Counts and
  limits unchanged (no survivor reached the cross-match). A failed SIMBAD/NED service is now named in
  `forced.ecsv` meta (`n_<service>` == -1), `inject` also treats `unread_images` as uncalibrated, and
  `check_params` guards the saturated-star / near_star coupling (#77).

## 2026-10-08: Orphan-pair cutout footprints: frame-token polygons, deterministic visit order
- Cloud run. Review follow-ups to #73 that its final squash did not carry: `_in_region` parses `POLYGON ICRS …`
  (a frame token used to match nothing, so every cutout read "outside"); MAST observations are sorted by obs_id
  before "first covering visit" is chosen (archive row order no longer picks the cutout). Counts and limits do
  not depend on either (cutouts only). TASKS: null (e) colour cell asymmetry.
- **Handoff:** TASKS "Now" 1, the D-050 W1 shear screen, is next; it needs the SMACS `_cat.ecsv` and DJA zout.

## 2026-10-08: Dark-deflector (orphan-pair) screen in deep fields, with injection-recovery: null (D-051)
- Worktree worker: `orphan_pairs.py` now reads CANUCS and DJA catalogues through one column layout (D-048 counts
  reproduced bit for bit) and runs on six deep fields: the five CANUCS NIRCam flanking fields and DJA GOODS-N
  (108 arcmin²). `scripts/inject_pairs.py` paints `exotic_sim` point-mass, W2 Ellis and W1 pairs into the real
  catalogues and runs the unchanged screen.
- 355 orphans against 315.4 expected (strictest null, P = 0.015); the excess sits in the flanking fields (246 vs
  211, P = 0.010), GOODS-N matches every null. Read as unmodelled physical companions (lensing would need ~140 dark
  galaxy-mass deflectors per arcmin²); open, not a candidate. 90 top orphans inspected: 50 ordinary (knots,
  satellites, groups, artefacts), 40 faint chance-like pairs. Nothing for `/vet-candidate`.
- Per-deflector efficiency ≤ 0.8 % (point mass, W2) and 0.1–18.6 % (W1, θ_E 0.3–1.5″). Background-aware 95 % limits:
  W1 < 1.3 × 10⁴ deg⁻² at θ_E 1.5″ (4.5 × 10¹¹ M☉); W2 < 2.9 × 10⁵ deg⁻² at 0.7″ (throat ≈ 10 pc). This covers
  θ_E ≤ 1.5″, where `radial` (D-049) is blind. docs/exotic_limits.md, docs/orphan_pairs.md.
- **Failed approaches (rules):** DJA Kron apertures for the same-galaxy rule (2.6× CANUCS; use 3.3 × flux_radius);
  the standalone DJA GOODS-N zout (older catalogue; stream the tarball member); a global null (c) in deep fields
  (under-predicts orphans; use the conditioned null (e)); unretried S3 reads (spurious NoSuchBucket via the proxy).
- Downloads > 200 MB (GOODS-N catalogue 224 MB, photo-z tarball 371 MB streamed) stated in D-051.
- Wall time: 3–20 s per field search; 190–266 s per field of injections.

## 2026-10-08: W1 limits re-run: mass-parametrised, blend-aware, independent nulls (D-049 update)
- #70 merged an intermediate version. This re-run fixes the review findings: one lens mass with θ_E per source
  redshift (photo-z, else z_s = 2), overlapping image pairs painted as one blend, an independent 200-draw null per
  batch of 10 trials, screen grid/null/defaults shared with `cmd_radial` (byte-identical SMACS output), a stated
  S/N ≥ 5 detection floor, measured footprint-border excess (3–16 %).
- Recovered of 1,600 per mass (2 × 10¹⁰ / 2 × 10¹¹ / 2 × 10¹² / 8 × 10¹² / 2 × 10¹³ M☉): 0 / 0 / 8 / 84 / 156.
  Headline 95 % limits (six photo-z fields, 38.0 arcmin²): none below 10¹² M☉; < 6.1 × 10⁴, 7.0 × 10³, 4.0 × 10³
  deg⁻² at 2 × 10¹², 8 × 10¹², 2 × 10¹³ M☉ (all eight fields, optimistic: 3.7 × 10⁴, 3.8 × 10³, 2.1 × 10³).
- **Failed approaches (rules):** one θ_E for every source; two lines for overlapping images; one fixed null for all
  trials; a `nanmin(S/N)` detection floor; recovery at the 3-line peak instead of p_random; catalogue shapes
  without PSF deconvolution; headline limits including fields without photo-z.
- Wall time: 371–993 s per field (1,000 injections, 100 independent nulls; 4 parallel).
- **Handoff:** a W1-specific screen benchmarked with this harness (TASKS).

## 2026-10-08: W1 shear screen reuse-check (D-050)
- `/reuse-check` for the W1-specific screen: build a catalogue aperture-mass map (−M_ap: a negative-mass lens gives
  negative tangential shear) with scipy cKDTree. TreeCorr NG only gives a stacked ⟨N M_ap⟩(R) and has no Windows
  wheels; lenspack works on pixelised maps. Rough S/N 1–2 at θ_E = 1″ in a 5″ aperture, so gains are expected
  mainly at θ_E ≥ 2–3″.
- **Handoff:** implement D-050 in `exotic_screens.py` and benchmark it with `scripts/inject_radial.py` on SMACS
  first (smallest field; needs its MAST `_cat.ecsv` and DJA zout in the cache).

## 2026-10-08: PR #70 merged after review fixes; Abell 2744 radial doc brought to the post-D-034 result
- Cloud run. PR #70 (D-049) was conflicted with main: merged `origin/main`, then `/code-review` on the final diff.
  Main finding: MACS0717 and Abell S1063 have no photo-z, so the injector painted their members and foreground
  galaxies as W1 images (Abell S1063 had the highest efficiency of all fields). The headline limits now use the six
  photo-z fields (38.1 arcmin²): < 7.0 × 10⁴ / 6.8 × 10³ / 3.2 × 10³ deg⁻² at θ_E = 3″ / 6″ / 10″, about 27× weaker
  than Takahashi & Asada; the all-field values stay in `limits.json` as optimistic. `inject_radial` now reads the
  screen defaults from `exotic_screens.radial_defaults()`.
- docs/fields/abell2744.md: radial result with the D-034 spike veto (134 `anti` arcs, max 5 lines, p 0.505; null).
- PR #61 (n0153, `needs-human`) brought up to date with main; it still waits for the owner.
- **Handoff:** the W1-specific screen (TASKS "Now" 1) is next; it fits one full cycle.

## 2026-10-08: W1 injection-recovery through `radial`: the screen is blind to negative-mass lenses (D-049)
- Worktree worker: `scripts/inject_radial.py` paints `exotic_sim` W1 lenses (n = 1, ε < 0) into the real catalogues
  of all eight null `radial` fields and runs the unchanged screen (`exotic_screens.radial_candidates`, a pure
  refactor of `cmd_radial`'s selection). 200 lenses per field and θ_E; 51.2 arcmin² screened.
- Recovered: 0 / 1,600 at θ_E = 0.3″ and 1″; 10 at 3″; 80 at 6″; 181 at 10″. 95 % limits on W1 lens surface
  density (six photo-z fields): none below 3″; < 7.0 × 10⁴ deg⁻² at 3″ (|M| ≈ 1.4–4.3 × 10¹² M☉), < 3.2 × 10³ deg⁻²
  at 10″. MACS0717 and Abell S1063 (no photo-z, members get lensed) are excluded as optimistic. About 27×
  weaker than Takahashi & Asada spread over 0 < z < 1. docs/exotic_limits.md.
- Why: an image reaches e ≥ 0.5 only for β ≲ 2.3 θ_E, and `anti` against the cluster keeps a third, so a lens puts
  1–3 arcs into a screen whose null needs 5–8 lines.
- **Failed approaches (rules):** recovery is the screen's p_random, not the 3-line peak; cache the null draws once
  per field (0.2–1 s per lens instead of 15–150 s); deconvolve the PSF before applying the lens Jacobian.
- Wall time: 322–877 s per field for 1,000 lenses (5 θ_E).
- **Handoff:** a W1-specific screen (collinear radial pairs flanking an empty centre, orientation against the
  candidate centre, local null), benchmarked with this harness.

## 2026-10-08: El Gordo lens model reproduces the published magnification maps (issue #68)
- `validate` now compares Lenstool-par models with published signed μ maps (`mag_map_files`, fetched only by
  `validate`); `map_check` reports `parity_agree` for signed maps. El Gordo (CDS J/A+A/678/A3, z=2 and z=4, |μ|<10,
  37-px sub-grid): median |μ| ratio 1.00001 / 1.00001, p95 relative difference 0.04 % / 0.07 %, parity 100 %;
  χ²pos unchanged (82.53). At full resolution (250,000 points), >20 % mismatches sit only on critical curves or
  within ~1″ of member cores (0.4″ map pixels).
- El Gordo `bayes.dat` `Chi2` explained (D-045 amendment): the sampling run uses `forme -10`, an image-plane χ² with
  σ² = a·b from the image list (19 of 56 images at 1.24″). It reproduces `Chi2` row by row to 0.1 % and the
  ln(Lhood) offset (75.904) exactly. The chain is validated. **Failed approach:** a source-plane χ² with free
  per-family weights (held-out ρ 0.81) fit only partly; the definition came from Lenstool's source.
- **Handoff:** El Gordo μ, parities and the MCMC chain are validated. Next: posterior μ spreads in `fluxratios`.

## 2026-10-08: Exotic-lens predictions first: wormhole / negative-mass searchable, warp not (D-047)
- Owner focus (2026-10-08): search only for signatures of traversable wormholes / negative-mass lenses and warp-drive
  spacetimes, predictions first. Research worker; every citation fetched from arXiv / Crossref.
- Searchable: W1 negative-mass dark lens (radial pair beside an empty centre), W2 Ellis pair without deflector, W3
  inverted microlensing (umbra between caustic spikes), W5 count deficit. Not searchable: the 4 % Ellis gutter, µas
  shifts, and every warp signature (Alcubierre exterior is flat; no published imaging/photometric prediction for a
  distant observer). The warp branch is stopped.
- `jwst_anomaly.exotic_sim`: Kitamura+2013 power-law lens family (either sign of ε), finite-source light curves,
  `inject_images` / `inject_light_curve` (`simulated`) for injection-recovery.
- Review fixes: an exact finite-source integral, checked against inverse ray shooting, puts the negative-mass spike
  peaks at ×7.0 / 2.35 / 1.53 (ρ = 0.01 / 0.1 / 0.3; first version 9.2 / 2.6 / 1.7). Exact demagnification onset
  added (2/(n+1) is KNA13's large-n estimate); KNA13's n = 3 "~10 %" is rounding of their Fig. 2c (13–14 %). NaN
  epochs no longer read as an umbra.
- **Handoff:** injection-recovery for `radial` (W1) and the dark-lens search (W2) to turn nulls into limits; a
  dimming class for the transient screen (W3); counts around `radial` centres (W5).

## 2026-10-08: Multi-plane lens models (D-046)
- `LensModel.split_planes` moves named potentials (e.g. a foreground galaxy fitted as a member) to their own
  redshift; `MultiPlaneLensModel` solves the multi-plane lens equation, and `find_images` / `backtrace_images` /
  `imageplane_residuals` accept it.
- Reproduces the D-042 system-51 result with library code: σ ≤ 70 km/s at z 0.268 leaves exactly 3 images.
- **Handoff:** use `split_planes` in `/vet-candidate` for any extra or missing image near a non-member galaxy.

## 2026-10-08: Orphan image pairs, a blind dark-deflector screen: null (D-048)
- Worktree worker: `scripts/orphan_pairs.py` looks for SED-matched close pairs with no published system and no visible
  galaxy between them in the CANUCS DR1 catalogues of MACS0416, MACS1149 and Abell 370.
- Pair excesses come from same-redshift groups (null (c)); 38 orphans against 34.5 expected. The 35 top orphans are
  knots, group members or chance matches, all at |μ| ≈ 1–2.5. No candidate.
- **Failed approaches:** a midpoint-only lens rule (missed galaxies between members); no per-band S/N cut (25 % false
  SED matches).
- Wall time: 5–10 s per field per search, 30–60 s with cutouts.
- **Handoff:** injection-recovery so the null becomes a limit; segmentation-map same-galaxy test; more clusters.

## 2026-10-08: Lenstool MCMC posteriors; Abell 2744 multiplicity residual explained (D-045)
- `lensmodel.read_lenstool_bayes` / `posterior_par` and `lens_consistency.py posterior` rebuild published models at
  MCMC samples (potfile rescaling, sampled family redshifts). Validated on Abell 2744: best.par is a chain row, and
  two random rows give χ²pos 173.5 / 173.0 against the chain's 178.1 / 174.1.
- Abell 2744 (12 samples plus best.par): 3.2a/b, 34.1a/b and 700.1a/b stay merged in every model. CATS v4.1 splits
  34.1a/b, and 3.2a/b sit on the caustic in both models. That is model resolution at folds; 0 surviving. 4.2c moves
  only 0.2–0.4″ (μ 8.7–10.0), well below its galaxy-scale systematics (D-036).
- **Failed approach:** the El Gordo CDS chain's `Chi2` column (54–77) does not track our χ²pos (93–106), not even in
  rank, so that chain is not validated.
- Wall time: about 70 s per Abell 2744 sample on a 0.25″ grid; 15 min for 13 models.
- **Handoff:** feed the posterior μ spread into `fluxratios` (then recheck SMACS 6.3 once a SMACS chain is pinned) and
  the position spread into `forced_check` search radii.

## 2026-10-08: CANUCS Lenstool models pinned (D-044)
- `macs0416-canucs` reproduces Lenstool's image-plane χ²pos (330.8 against 344.30; rms 0.51″ over 303 images).
  `abell370-canucs` is a source-plane fit (image-plane rms 2.3″).
- They are the independent second model for vetting (the D-042 system-27 and system-51 checks used scratch code).
- `fluxratios` now refuses gated image lists and map models, as `images` does; its `--offset-arcsec` defaults to
  the model's pinned frame offset.
- **Handoff:** a two-plane option in `LensModel`; use `macs0416-canucs` in `images` runs as a direct cross-check.

## 2026-10-08: Abell 370 and Abell S1063 radial screens: null (D-043)
- Worktree worker: CANUCS 1208 for Abell 370 (CANUCS DR1 photo-z; frame offset −0.121″, −0.015″ pinned) and GLIMPSE
  3293 for S1063 (DJA v7.5 photo-z). Only `radial` applies (image lists gated, D-035).
- Abell 370 raised two flags (15 and 8 lines, p ≤ 0.01). Both are diffraction-spike chains from Gaia stars that are
  off the mosaic or saturated and absent from the catalogue. The Gaia-seeded veto is now in the code
  (`radial --spike-stars`, `scripts/gaia_stars.py`); with it, p 0.945. S1063: p 0.435, unchanged by the veto.
- **Failed approach:** spike vetoes seeded from the pipeline catalogue miss saturated and off-mosaic stars, whose
  spikes reach 37″.
- Wall time: under 2 min of pipeline per field, plus cutout vetting.
- **Handoff:** a low-weight veto or an aper50 S/N floor in `radial`; use `--spike-stars` on every field.

## 2026-10-08: MACS0416 counter-images and flux ratios: null; CANUCS model as second model; two-plane check (D-042)
- With the D-040 solver, the MACS0416 image list is open. `images` and `fluxratio` raised 21 flags; 0 survive. The most
  persistent was system 51's fourth image: two independent models predict it, and the 51.1–51.3 photometry says it
  should appear at 9–25σ, but it is not seen. It comes from a z 0.268 foreground galaxy modelled as a member. A
  scratch two-plane model with that galaxy at σ ≤ 60 km/s gives exactly the 3 observed images.
- System 27's two bright `absent` predictions (S/N 247–341) are CATS-only galaxy-scale caustics. The JWST-era CANUCS
  Lenstool model (222 potentials, 111 spectroscopic systems) reproduces system 27 with exactly its 3 images, and the
  cutouts show empty sky there.
- Flux-ratio flags: 45.2 is blended with a bright galaxy 0.5″ away; 38.1 is marginal.
- Wall time: about 3 min of pipeline plus about 5 min of vetting.
- **Handoff:** pin the CANUCS models (MACS0416, MACS1149, Abell 370) as `MODELS` entries for second-model vetting.

## 2026-10-08: Sunrise transient candidates `n0022` and `n0150` are detector persistence (D-039)
- New `scripts/persistence_check.py`: per-exposure photometry on level-2 `_cal` files (S3 byte ranges), plus the
  same detector pixel in earlier exposures on that detector. Validated on a synthetic afterimage and on two real
  sources (`n0153`, Earendel: 0 suspect detections).
- `n0022`: afterimages of a bright galaxy (2.0 % and 0.8 % of its flux) in o010 dithers 3 and 4; the dither geometry
  stacks both on one sky position. `n0150`: in all three epochs an afterimage (0.04–0.07 %) of a saturated star, one
  exposure later; its "motion" is the dither vector. Cutouts inspected (dithers 1–4 at fixed detector pixels).
- `n0153` is on sky in all 8 dithers and F150W is ×1.9 brighter in 2025 at the `_cal` level: still a variable
  candidate.
- **Lesson:** a single-epoch source in a mosaic can be an afterimage that two dithers place on one sky position.
  The archived `_cal` files carry no DQ flag there.
- **Handoff:** `/vet-candidate n0153`; persistence-check new single-epoch candidates before vetting.

## 2026-10-08: MACS0717 screens: null (D-041)
- Worktree worker on VENUS 6882 o029 (10 bands, the only public NIRCam association) with `macs0717-cats` (rms 3.21″).
- 51 flags, 0 surviving. 29 are the model's own copies of catalogued images it does not reproduce (1.6–5.6″ off), 6 are
  CATS-only extra images (Sharon v4cor and Keeton v4 predict none), and 5 have μ more than 2× model-dependent.
  System 65's flux ratio is a 0.6″ catalogue offset. `radial` p ≥ 0.70.
- **Failed approach:** a fixed 1.5″ match radius for a 3.2″-rms model makes most catalogued images "unpredicted",
  and their model copies then flag as absent or confused.
- Wall time: about 7.5 min of pipeline plus about 20 min of vetting.
- **Handoff:** in `forced_check`, rms-scaled match and search radii and automatic copy classification; Sharon v4cor and
  Keeton v4 as pinned `MapLensModel` entries.

## 2026-10-08: `find_images` fold refinement; frame offsets for map models (D-040)
- Cells on a critical curve near the source are subdivided into ≤ 0.02″ sub-cells. MACS0416 system 26 is now solved
  (rms 1.57″ → 0.76″), and its image list is open with offset (0.208, −0.025). No other model changes beyond 0.03″
  rms; Abell 2744 `validate` takes 41 s instead of 29 s.
- **Bug fixed:** `apply_frame_offset` was a no-op for map models (the maps are looked up by sky position), and
  `radial` / `arcs` skipped it for them. No earlier result used a map model with an offset. MACS0416's radial screen
  was re-run in the JWST frame: still null.
- **Handoff:** `images` / `fluxratio` on MACS0416 (CANUCS photo-z); radial on Abell 370 and Abell S1063.

## 2026-10-08: MACS0416: system 26 is a solver-grid artefact; radial null (D-038)
- Worktree worker on CANUCS 1208 (`jw01208-o004_t002`, 8 bands). The CATS system-26 residual (11″) is the 0.25″
  `find_images` grid missing a merging pair near the critical curve; the source lies 0.001–0.005″ from the caustic. On a 0.1″ grid the rms is 0.811″ (1.13×
  quoted). The model stays map-only until the solver refines its grid near high |μ|.
- `radial`: 98 centres against 99.6 random; the 7-line centre (p 0.29) is low-S/N noise segments. With CANUCS DR1
  photo-z (background cut): 77 against 85.0, p 0.225. 0 flags.
- **Failed approach:** jwst 3.0.0 isophotal S/N admits noise segments (66 of 137 anti arcs have aper50 S/N < 3). The
  radial screen needs an aper50 S/N floor.
- Wall time: about 4.5 min (lens and exotic scripts plus vetting stamps; no pipeline `run`).
- **Handoff:** adaptive grid refinement in `find_images`; an aper50 S/N floor in `radial`.

## 2026-10-08: CANUCS DR1 photo-z for MACS1149 (D-037 addendum)
- **Failed approach:** the worker checked only DJA, but docs/landscape.md already listed CANUCS DR1 (PSF-matched
  EAzY photo-z for A370, MACS0416, MACS0417, MACS1149 and MACS1423). Check docs/landscape.md before reporting that a
  field has no photo-z.
- MACS1149 system 16.2: z_phot 2.25 (95 % 0.23–2.33). This disfavours the CATS-fitted 4.419.
- Radial re-run with the background cut: 12 peaks against a null mean of 12.8; null. 6 of 68 matched lensed images
  get a blended z < 0.6.

## 2026-10-08: MACS1149 screens: null (D-037)
- Worktree worker on CANUCS 1208 (8 bands) with `macs1149-cats`. The image-plane rms is 0.673″ (gate passes) and
  the frame offset is under 0.02″.
- 3 flags, 0 surviving: system 16 (fitted z 4.419; at z 2.5 in CATS or 3.0 in Sharon v4cor the image lands on 16.2), system 2 (CATS topology error,
  checked against Sharon v4cor), Refsdal-host knot 1192 (next to the BCG; μ differs >2× between models).
  `fluxratio` 0 flags; `radial` p = 0.91.
- **Failed approach:** forced photometry with an annulus on a BCG core gives negative fluxes; use high-pass.
- No DJA mosaic for MACS1149 (v7).
- Wall time: about 5 min of pipeline plus about 10 min of vetting.
- **Handoff:** `macs1149-sharon` (or the CANUCS models) in `MODELS`; the D-037 rules in `forced_check`; system 16's
  redshift (spectroscopic z, or forced photometry on its far-image track); `images` with CANUCS photo-z.

## 2026-10-08: Abell 2744 screens: null (D-036)
- A worktree worker ran `validate`, `images --forced-image`, `fluxratio` and `radial` on Bergamini+2023b, with
  cutouts of every flag. 15 of 16 flags were ordinary. The survivor, 4.2c, was re-run under the D-034 rules: it
  is `no_reference` (resolved-knot references, BCG halo), not a candidate.
- Radial: 35 peaks against a random mean of 50.4, p ≥ 0.575.
- **Lesson:** a "252σ absent" image can come from a resolved-knot reference. D-034's compact-reference rule now
  catches this.
- Family 4's c images were vetted (4.1c 4–8× underluminous after BCG subtraction; 4.2c undetected). They are
  explained by μ systematics next to member 34423 (3.9–28.7 under ±30 %; CATS 7.3). **Rule:** a μ that moves by
  more than 2× under member perturbation is untestable.
- **Handoff:** HFF field runs (D-035 models), `bayes.dat` position spreads, UNCOVER v2.0 cross-check.

## 2026-10-08: Sunrise third epoch with calibrated significances; two open transient candidates (D-027)
- `transient_forced.py --controls`: noise scale from ordinary sources, applied before thresholding.
- Sunrise o010 against VENUS o052 (2.9 yr; F150W, F444W): scales 1.30/1.18; **0 of 57 catalog candidates pass**;
  Earendel steady (Δm ≤ 0.12 mag, < 1σ).
- Among the 200 controls, `n0022` (gone after 2022-07 in four SW bands), `n0150` (a different position in each
  epoch) and `n0153` (+0.5–0.6 mag in both bands) change. Not vetted; ordinary explanations (supernova, moving
  object, AGN, edge artefact) are untested. docs/fields/sunrise.md has the numbers.
- **Failed approach:** the catalog stage plus the two-band rule misses single-pair, blue transients.
- Follow-up #53: an empty control selection is flagged uncalibrated instead of being skipped silently.
- **Handoff:** `/vet-candidate` for `n0022` and `n0150` (level-2 `_cal` exposures per filter: is the source in
  every dither? Does `n0150` move within one visit?); then a grid of forced photometry (all sources, not only
  catalog candidates) for every epoch pair.

## 2026-10-08: Six HFF clusters as CATS map models (D-035)
- MACS0416, MACS1149, Abell 370, MACS0717, Abell S1063 and Abell 2744 (CATS v4/v4.1) are pinned map models. All six
  reproduce their published z = 2 magnification maps (median 9e-5 to 6e-3).
- Image lists pass the rms gate for MACS1149 and MACS0717. MACS0416 has a `params.txt` but fails the gate (system 26).
  Abell 370, Abell S1063 and Abell 2744 have no `params.txt`, and their placeholder redshifts give 9–12″
  residuals. All four stay map-only.
- **Failed approach:** using CATS `arcs.txt` redshifts as published.
- **Handoff:**
  - per-cluster JWST field runs (configs, catalogs, `radial` on all six, `images` / `fluxratio` on MACS1149 and
    MACS0717);
  - `/vet-candidate` for Abell 2744 4.2c (the field worker's survivor).

## 2026-10-08: El Gordo counter-images and radial screen: no candidate (D-034)
- Forced photometry: 3 `absent` images on the first pass, all ordinary on vetting (cutouts plus numbers):
  - 23: model position error;
  - 6: reference on a galaxy wing; photo-z-consistent counterparts at 2.6–3.1″;
  - 7: HST→JWST frame offset of 0.22″.
- Radial screen: max 4 lines, p = 0.64. A 6-line "centre" was a star's diffraction spikes.
- Four new rules (D-034): frame offset, compact and consistent references, a residual-scaled search radius, and a
  spike-segment veto. SMACS re-run: unchanged, 0 absent.
- Tally: El Gordo screened 17 uncatalogued predicted images and 37 anti arcs; flags 3 + 1; **surviving 0**.
  Wall time: about 50 s for `images`, 21 s for `radial`.
- **Handoff:** Abell 2744 (worker running), then the HFF/RELICS map fields via `MapLensModel` (#47).

## 2026-10-08: Map-based lens models; Sunrise radial screen null (D-033)
- `lensmodel.MapLensModel` evaluates published deflection maps. The first one is `whl0137-relics-lenstool`, which
  reproduces RELICS κ to 3.5e-5 and μ(z = 6.2) to 2.0e-5 (medians).
- Sunrise `exotic_screens radial`: 29 usable anti arcs, 3 centres against a null mean of 2.2, max 3 lines, p 0.885.
  **No candidate.** Wall time 17 s, with cached maps and catalogs.
- Cycle tally:
  - SMACS + Sunrise screened: 60 + 0 images, 34 + 29 arcs;
  - flags: SMACS fluxratio 2 (removed by the two-band compactness rule), radial 0 significant;
  - surviving vetting: 0.
- **Handoff:** more map fields (HFF: Abell 2744, MACS0416, MACS1149, Abell 370; RELICS clusters), El Gordo and
  Abell 2744 `images --forced-image` and screens.

## 2026-10-08: Flux-ratio and colour test of catalogued images: no anomaly in SMACS or El Gordo (D-032)
- New: `lens_consistency.py fluxratios`.
  - It compares each image's DJA `mag_auto` + 2.5 log|μ| and its F150W−F444W colour with the other images of its
    system (leave-one-out).
  - It drops blends, segments larger than 20,000 px and counterparts whose photo-z excludes the system redshift.
- Results (`derived`):
  - SMACS: 8 images in 4 systems flux-tested, rms 0.60 mag.
  - El Gordo (0.3″ after a +0.221″ frame shift): 23 images in 11 systems, rms 0.49 mag.
  - No colour outlier (rms 0.11 and 0.07 mag).
  - The flagged pairs SMACS 6.1/6.3 and El Gordo 18b/18c (and 7b/7c, at 0.5″ without the shift) fail vetting.
    Forced photometry brings SMACS 6 down to 0.61–0.78 mag. The El Gordo pairs are low-S/N or chromatic, so they
    are measurement failures. **No anomaly.**
- **Failed approaches**, each checked on cutouts or SEDs:
  - DJA aperture × `tot_corr` is not a total flux for arcs.
  - SMACS 1.1 looked 1.4 mag too bright, from host-halo light in a 32,864 px segment.
  - El Gordo 9a/9c differed by 1.2 mag in colour, because 9a's DJA counterpart is a z_phot 0.89 interloper.
  - Bare eazy 95 % intervals exclude good images.
  - A median that includes the image itself halves pair differences.
- Limit: DJA misses most arcs inside cluster light (SMACS: 38 of 60 images unmatched within 0.3″).
- **Handoff:**
  - El Gordo `images --forced-image` (TASKS "Now" 2);
  - Abell 2744 `fluxratios` with the 233 MB DJA catalogue (needs a DECISIONS entry);
  - `bayes.dat` μ uncertainties, then recheck SMACS 6;
  - BCG/ICL-subtracted totals, to test the core images.
## 2026-10-08: Exotic-lens screens; SMACS null (D-031)
- New `scripts/exotic_screens.py`:
  - `fluxratio`: two-band forced photometry, luminosity ratio against sibling images, compactness and chromatic
    gates;
  - `radial`: anti-tangential arcs whose axes converge on a dark centre, with a false-alarm rate from randomised
    position angles.
- SMACS:
  - fluxratio: 6 compact images consistent, 0 flags. A first, single-band compactness gate had flagged system 7,
    an ordinary knot-vs-whole-arc mismatch;
  - radial (background sources only): 12 centres against a null mean of 9.0 (p95 15); max 4 lines, p 0.945.
  - **No exotic candidate.**
- **Failed approaches (now rules):**
  - fixed-aperture flux ratios on resolved arcs: surface brightness is conserved, so the ratios scale with 1/|μ|
    (systems 5 and 10);
  - a single-band compactness gate (knots of a clumpy arc pass in F150W);
  - a uniform-angle null for the radial screen (the selected arcs point at the mass centre).
- **Handoff:** fan out per cluster (El Gordo, Abell 2744, Sunrise, then HFF/RELICS). Each runs `validate`,
  `images --forced-image`, `exotic_screens fluxratio` and `radial`, with contact sheets of all flags.

## 2026-10-08: `find_images` 3–5× faster with identical images
- Seeds are pre-filtered with boolean sign tests on the mapped grid corners, and Newton steps run for every seed in
  one `fields_xy` call.
- Benchmark: all catalogued systems at a 0.25″ grid, images identical (max |Δ| 0 arcsec, same counts):
  - SMACS: 23.2 → 6.7 s;
  - El Gordo: 44.3 → 13.0 s;
  - Abell 2744: 118 → 23.5 s (at ±190″).
- **Handoff:** the deflection grid itself (one-time and cached) is now the main cost. Published deflection maps
  (UNCOVER, RELICS, HFF) could replace it for fields without a Lenstool model (TASKS).

## 2026-10-08: Image-plane χ² reproduces Lenstool for SMACS, El Gordo and Abell 2744 (D-030)
- Merged #40 (counter-images, D-029) after its last commit, which GitHub had not attached to the PR, was picked up
  by a follow-up commit.
- `validate` now computes the exact image-plane χ² for every model; El Gordo and Abell 2744 are `MODELS` entries.
  Issue #41's numbers are now reproducible in the repository: χ² 30.87/30.91 (SMACS), 82.53/80.22 (El Gordo),
  146.64/146.60 (Abell 2744). No catalogued image is off by more than 3σ in any field.
- Parser fixes: letter-suffixed image ids, several ids per `z_m_limit`, 6-decimal `_kpc` rounding.
- **Failed approach:** the per-image error column for El Gordo (χ² 52.0); a uniform 0.621″ matches Lenstool.
- CDS reset connections through the proxy in this run; the cache was seeded from a `curl` copy with the same sha256.
- **Handoff:** run `lens_consistency.py images` with forced photometry on El Gordo and Abell 2744 (TASKS "Now" 2);
  then the exotic screens (docs/exotic_lensing.md).

## 2026-10-08: SMACS counter-images: no predicted image is absent (D-029)
- New: `lensmodel.find_images` (an image-plane solver on a cached deflection grid), `lens_consistency.py images`, and
  `--forced-image`, which runs forced photometry on S3 byte-range stamps.
- Results:
  - ICLv2 reproduces all 60 catalogued images within 0.04–0.91″, with 4 demagnified central images.
  - Of the 11 testable uncatalogued images: 3 are recovered (systems 9 and 8, and system 17 marginally: flux ratio
    0.36 on the BCG gradient), 1 is
    confused, 1 is undetectable, 6 have no reference flux, and **0 are absent**.
  - System 8's model z = 11.76 is contradicted by its F090W/F150W detections.
- **Failed approach:** pipeline-catalog flux references near cluster galaxies. They falsely made system 9 `missing`.
- The Mahler κ tarball now comes from raw.githubusercontent.com, because github.com/raw returns 403 in cloud runs.
- **Handoff:** TASKS "Now" 1 is done for counter-images. Next:
  - parser fixes (issue #41), then the same test on El Gordo and Abell 2744;
  - the exotic screens of docs/exotic_lensing.md (demagnified images, radial images around a dark centre).

## 2026-10-08: El Gordo and Abell 2744 lens models, preliminary: no image-position anomaly (issue #41)
- Preliminary validation on published products. The numbers come from scratch code and are not yet reproducible
  in the repository; they need a re-run with #40's solver.
  - `lensmodel.py` reproduces the Caminha+2023 El Gordo magnification maps (median |Δμ|/μ 4e-5).
  - Exact image-plane solves match Lenstool's χ²: El Gordo 82.5 vs 80.22 (rms 0.754″, paper 0.75″); Abell 2744
    (Bergamini+2023b) 146.64 vs 146.60 (rms 0.427″).
- Results:
  - No image-position anomaly in either field.
  - Predicted uncatalogued counter-images: El Gordo 16, inconclusive, because the MUSE Lyα images are too faint in
    continuum. Abell 2744 30, not yet testable, because the pipeline catalogs miss most arcs in the core.
  - docs/fields/elgordo.md, docs/fields/abell2744.md.
- **Failed approach:** source-plane back-trace χ² for image-plane-optimised models (El Gordo 121.6 against 80.22).
- **Bugs found on `main`:** the parser fails on letter-suffixed image IDs and on 6-decimal `_kpc` rounding (#41).
- Also merged #38 after bringing it up to date with `main`.
- **Handoff:** #40 (image solver, local session); then reproduce these numbers in the repository, the parser fixes,
  and the counter-image flux test on DJA photometry (TASKS "Now" 2).

## 2026-10-08: Sunrise two-epoch search, a null result; Earendel steady (D-027 amendment)
- WHL0137 2282 o010 against o120 (same pipeline): 372 catalog candidates → 58 in ≥ 2 bands → 52 after the bright-Gaia
  mask → 0 by recentred forced photometry in all four bands. 2 pass in one or two bands, both ordinary on visual check
  (an epoch-2 streak; galaxy outskirts at the footprint edge). docs/fields/sunrise.md.
- Earendel: |Δm| ≤ 0.18 mag (≤ 1.5σ) in F090W/F115W/F277W/F356W over 164 days.
- Noise calibration on 150 ordinary sources: ERR-based significances are 1.2–1.5× too large.
- New: `scripts/transient_combine.py` (band coincidence and bright-Gaia mask) and
  `transient_forced.py --recentre-arcsec`.
- **Failed approach:** a fixed aperture at a literature position. Earendel seemed to brighten by 0.76 mag at 5.7σ,
  because the aperture sat 0.15″ off the source and the PSF rotates by about 180° between the epochs.
- **Handoff:** lens-model checks per field (#35 merged; El Gordo and Abell 2744 model validation in issue #41).
  Time domain: rescale the thresholds by the control std, and add VENUS 6882 o052 (F150W, F444W) as a third epoch.

## 2026-10-08: Cloud runs use the GitHub MCP tools; conflicting PRs get no CI (D-028)
- The first routine run opened #38 with the session's GitHub MCP tools. The prompt, the research-cycle skill and
  docs/operations.md §3 now prefer them, merge with `mcp__github__merge_pull_request` after checking the policy,
  and fall back to REST plus the `merge-ready` label.
- **Failed approach:** #37, a REST squash-merge allow rule, was closed. The glob `pulls/*/merge` also matches other
  `pulls/...` writes, so the Bash merge API stays denied.
- **Lesson:** a PR that conflicts with `main` gets no `pull_request` CI at all (#38). Merge `origin/main` in first.
- The cloud merge gate checks author, head repository (no forks), branch, labels, every page of files and the CI jobs
  (not the skipped `claude` runs), and pins `expectedHeadSha`. D-028 records the decision.
- **Handoff:** TASKS "Now" 1 (counter-images, then parity and flux ratios; see the lens-model entry below). #38
  needs `origin/main` merged in before its CI and merge.

## 2026-10-08: Lens-model stage validated; SMACS arc orientations agree with ICLv2 (D-024)
- `lensmodel.py` is a Lenstool dPIE port (from PyAutoGalaxy, MIT) that evaluates a published `best.par`. It
  reproduces Mahler+2022 ICLv2:
  - κ map: median |Δκ| 1.6e-5;
  - back-trace χ² of the 60 catalogued images: 31.18, against Lenstool's 30.91.
- `scripts/lens_consistency.py` provides `validate` and `arcs`.
- SMACS arcs: 21 of 25 elongated strong-shear background sources, each tested at its own photo-z range, are
  aligned with the predicted stretch (p = 2.6e-7). The 6 anti candidates are all ordinary:
  - 2 segmentation blends;
  - 1 galaxy at z ≈ 0.77 with an intrinsic shape;
  - 3 noisy low-surface-brightness shapes.
  This is a null result.
- **Lesson:** moment orientations from pipeline segments pick up blends and low-S/N shapes, so filter them before
  calling a source anti-tangential.
- **Handoff:** TASKS "Now" 1, which leaves counter-images (predicted but missing, or observed but not predicted)
  and parity/flux ratios. Then El Gordo and Abell 2744; check their model profiles first.


## 2026-10-08: Hourly cloud routine
- Routine `trig_01PNAmgcfqef8CvhPAY8ggbP` runs `/research-cycle` loops hourly at :07 UTC (Opus 5.5, no connectors)
  with the prompt in docs/cloud-routine-prompt.md.
- First test run: clone, venv, skills and MAST downloads work.
  - GitHub GraphQL is blocked in cloud sessions, so `gh pr` / `gh issue` fail. The REST equivalents and a single
    background CI-wait loop are in docs/operations.md §3, and the research-cycle skill points to them.
  - The merge API is denied by `.claude/settings.json`, so cloud runs stop at merge-ready with the label
    `merge-ready`. The owner or a local session merges.
  - Coordination: runs last about 40 minutes and skip `local-wip` PRs and branches committed to in the last 15
    minutes. An unmerged PR's handoff is read from its branch.
- **Handoff:** unchanged (TASKS "Now (M3)"); the cloud routine continues it.

## 2026-10-08: Two-epoch transient search, a null result; exotic-lens signatures (D-027)
- `scripts/transient_search.py` produces catalog-level candidates; `scripts/transient_forced.py` checks them with
  forced aperture photometry on byte-range cutouts of both epochs.
- SMACS against VENUS: catalog 441 (F444W) and 628 (F150W) → 144 in both bands → 4 by forced photometry (6 with
  WCS-centred apertures) → 0 after a visual check. All are bright stars or a source next to a bright neighbour, where
  fixed-aperture systematics dominate. No credible transient at |Δm| ≥ 0.3 mag, ≥ 5σ.
- **Failed approach:** catalog-only comparison across pipeline versions, which picks up deblending differences.
- docs/exotic_lensing.md records the verified exotic-lens signatures (wormhole demagnification, radial images for
  negative mass), the degeneracies with ordinary lensing, and existing limits. Warp drives have no imaging
  prediction to test.
- **Handoff:** the lens-model stage (TASKS "Now (M3)" 1), then Sunrise's same-pipeline epoch pair.
## 2026-10-08: Three cluster fields; robustness rules from them (D-025, D-026, cycle 17)
- Parallel worktree workers added Abell 2744 (#32), El Gordo (#31) and Sunrise/WHL0137 (#30). PR #29 fixed downloads
  for MAST products whose listed size is stale after the 2026-10-01 reprocessing. Proposals are folded into
  D-026 and SOURCES.
- **D-025** (found on Abell 2744):
  - catalogue "stars" with DJA r50 too large for a point source (bright cluster-galaxy cores in Gaia) go back to the
    galaxy ranking, before the locus calibrates;
  - very red spiky sources are exempt from the D-020 host test.
  - Run `20261007T211229Z-883a5066`: 8 vetoed, the locus recalibrated, the star top 10 is now mostly point sources with spikes, and
    `5904` is kept.
- **Handoff:** TASKS "Now (M3)". The lens-model stage (worker, `claude/lens-model`) is in flight; then per-field
  model checks and the transient forced photometry.


## 2026-10-08: Owner priorities, lensing-violation search first (D-023)
- The owner merged PR #20 (`2804`, inconclusive) and reprioritised: M3 lens-model consistency first; more clusters in
  parallel (Abell 2744, El Gordo, Sunrise; MACS J0416 has no DJA v7 catalog); a two-epoch transient search; no
  pauses; one review per PR.
- **Handoff:** TASKS "Now (M3)" items 1–3. Worker PRs with label `batch-clusters` are in flight.


## 2026-10-08: CEERS control on DJA matched photometry (D-022, cycle 16)
- `ceers_t021` now uses the DJA v7.4 `ceers-full` catalog, a 250.5 MB download with the reason stated in D-022.
  The control field is built the same way as SMACS for colours, confirmation, screening and the top-k metric.
- Run `20261007T202419Z-65728e80`:
  - 3,069 of 7,582 sources match DJA, and 3,156 pass the gate (2,757 before);
  - the top 20 is 85% known objects with 0% cutout-flagged, and 4 sources were screened;
  - the stellar locus does not apply (2 catalogued stars).
- **Found:** `f200w_2842` (#4) is a single-band, ellipticity-0.95 streak, probably a spike from the bright source
  4.6″ away. A coincident DJA object confirmed it (D-014 "Revisit if").
- **Handoff:** TASKS "Spike streaks confirmed by a coincident DJA object", then the star-stratum NED matches.


## 2026-10-08: Low-weight screening of the top-k pool (D-021, cycle 15)
- `stages.cutouts.screen_low_weight` re-applies D-011's 0.5 weight threshold with the cutout's own core weight. It
  runs on every non-star stratum, and `screened.ecsv` now records each source's `reason`.
- Run `20261007T194733Z-ac5f1b59`: the stripe artifact `3034` (#3) is removed: short-wavelength only, not in DJA, core weight 0.468.
  The galaxy top 20 is now 40% known objects, with 1 flagged source (`2915`, kept on purpose) and 4 screened
  (one by D-021, two by D-019, one by D-020).
- **Handoff:** TASKS "Star-stratum sources with NED galaxy matches", then MIRI's `G_Lens` vetting.


## 2026-10-08: Host test for spiky sources (D-020, cycle 14)
- `cutouts.host_ratio` measures annulus light against the peak around the spike peak. D-019 screening now also
  removes spiky sources with no host light (< 0.004), so the two-star segment `940` (0.0017) leaves the top 20.
  `2915` (0.0076, a galaxy nucleus) stays.
- Run `20261007T190843Z-db3c83fe`: the galaxy top 20 is 40% known objects, with 1 spiky source (`2915`, legitimate) and 3 screened.
- **Next contamination:** `3034` (#3) is a diagonal stripe artifact in a low-weight region (cutout flag `low_weight`).
  It passed the D-011 gate.
- **Handoff:** TASKS "Stripe artifact at galaxy rank 3".


## 2026-10-08: Top-k composition metric; `940` re-diagnosed (cycle 13)
- Every stratum now reports its top-k composition: known, lens-related, catalogued star, cutout-flagged,
  spikes and screened. It is stored in `run_record.json` as `samples[].topk` and defined in
  docs/methodology.md. The crossmatch rows now keep `is_lens_related`.
- Baseline, run `20261007T183007Z-4342c5c2`: the galaxy top 20 is 35% known objects, 10% lens-related, 0% stars and 15% flagged.
- **`940` is not saturated.** Its F200W cutout has no no-data pixels. The pipeline merged two stars into one
  5,788 px segment whose centroid lies 0.20″ from the brighter star's peak, so aper50 at the centroid is faint
  (25.7 against isophotal 20.1). It is a single-band detection confirmed by the DJA match (0.18″) to the bright
  star.
- **Failed approaches** (saturation tests):
  - aper50 − isophotal magnitude: cluster galaxies reach 10–21 mag, so `940` (5.6) is not an outlier;
  - a no-data core within 0.3″: `940` has none, and only `1345` does.
- **Handoff:** TASKS "Star-pair segments in the galaxy ranking" (a host-vs-PSF profile test or a centroid–peak
  offset rule) and "Star-stratum sources with NED galaxy matches".


## 2026-10-08: Vetted `f200w_2915`, a catalogued quiescent galaxy with a compact core (cycle 12)
- Run `20261007T174424Z-540be92d` (clean `a27d628`), galaxy rank 21, kept ranked by D-019 as a spiky source with non-stellar colours.
- It is SIMBAD `[VBG2023] SMACS 1060` (Valentino+2023 atlas of colour-selected quiescent galaxies at z > 3). DJA
  z_phot is 2.84 (95% 2.68–3.10), mildly below that selection. Its F090W−F150W break is 2.5 mag. The spikes come from its compact core. D-019 behaved as
  intended.
- **Two epochs:** the position residuals disagree between bands (7 mas F444W, 34 mas F150W), as expected for centroid
  shifts of an extended source, and Δm is F444W only (+0.11; F150W +0.01). So it shows neither motion nor
  variability. `epoch_compare.py` now documents that its significance assumes a point-like target.
- **New resource:** program 4043 has NIRCam F444W grism spectra over part of SMACS 0723 (SOURCES).
- **Handoff:** TASKS "Saturated stars in the galaxy ranking", then the star-stratum NED matches.


## 2026-10-08: Spike screening of the galaxy top k (D-019, cycle 11)
- The galaxy strata get cutouts for 2k. Spike-flagged sources with stellar colours (`1571`, `2242`) are removed
  and backfilled. Spiky sources with non-stellar colours stay ranked and are noted: `940` (a saturated star with
  corrupted colours) and `2915` (a bright nucleus, possibly an AGN).
- **Failed approach:** screening every spiky source removed `2915`. Spikes also mark bright galactic nuclei.
- Run `20261007T170539Z-b369bfa2`: the top 20 now runs to rank 22, with no star+galaxy blends left.
- **Handoff:** TASKS "Vet `2915`" (a bright red nucleus), "Saturated stars in the galaxy ranking" and "Star-stratum
  sources with NED galaxy matches".


## 2026-10-08: Diffraction-spike flag on cutouts (D-018, cycle 10)
- `cutouts.spike_statistic`: hexagonal-harmonic power around the brightest peak near the target. NIRCam cutouts with
  `spike_s6` ≥ 3 get the `spikes` flag, and the report lists them.
- Run `20261007T162133Z-2511977c`: it flags exactly the three PSF-like blends in the galaxy top 20 (`940`, `2242`, `1571`) and all of
  the star-stratum top 10. `940` is really a saturated star: its pipeline aper50 is 5 mag fainter than DJA's
  measurement, and the stellar locus missed it because its colours are corrupted.
- **Failed approach:** centring on the catalog centroid, which misses blends whose star is offset.
- PR #20 (`2804`, needs-human) was revised after review: F410M−F444W is now an observed colour, and the line
  is a hypothesis.
- **Handoff:** TASKS "Spike-flagged sources in the galaxy ranking" (do flagged sources leave the top k?) and
  "Star-stratum sources with NED galaxy matches".

## 2026-10-08: Vetted `f200w_2804`, inconclusive, a little-red-dot-like source (cycle 9b)
- Run `20261007T151316Z-5706ae30` (clean `6fade35`), galaxy rank 6. Uncatalogued, unresolved, F444W 22.9.
- **Second epoch:** no proper motion over 4.0 years (VENUS 6882, 2026-06-05; 1.9 mas = 0.95σ, so
  |μ| ≲ 1.5 mas/yr) and no F444W variability (0.004 mag).
- F444W is 0.725 mag brighter than F410M in the same epoch. That is an observed colour; a line at 4.25–4.98 µm
  and F410M absorption are both hypotheses.
- It meets Kokorev+2024 little-red-dot colour cuts. A brown dwarf is disfavoured on two counts: the 1.5–2.8 µm
  colours are red rather than blue, and there is no motion.
- Verdict: `inconclusive: needs NIRSpec spectroscopy`. Record: `docs/candidates/jw02736-o001_t001_nircam_f200w_2804.md`.
- **Handoff:** the owner decides whether to pursue it (`needs-human`).


## 2026-10-08: Two-epoch vetting tool (D-017, cycle 9a)
- `scripts/epoch_compare.py`: proper motion and variability of one target from two public level-3 catalogs of one
  filter, with a local median frame tie. Program 6882 (VENUS, 2026-06-05) gives SMACS 0723 a 4-year second epoch in
  F150W, F182M, F210M, F300M, F410M and F444W (SOURCES "Second epochs").
- **Handoff:** the vetting record for `f200w_2804` follows in its own PR (inconclusive, `needs-human`).


## 2026-10-08: One-sided stellar-locus size test (D-016, cycle 8)
- Bright calibration stars have inflated r50 (5–8 px under 20.5 mag), so D-015's two-sided band around 2.86 px
  missed the unsaturated stellar sequence at r50 ≈ 2.0 px. The size test is now 0.5–1.2 × r50_psf.
- Run `20261007T144011Z-11342d69`: the locus adds 90 stars, not 16. The star stratum has 109; its top 10 is all PSF-like.
  `f200w_1874` is now a star.
- Point-like sources left in the galaxy top 20 are faint (mag > 24), blends, or `2804`, whose red
  F356W−F444W makes it a brown-dwarf-like candidate to vet. Details are in D-016.
- **Handoff:** TASKS Now 1 (vet `2804`), Now 2 (blended stars).


## 2026-10-08: Vetted four high-ranked SMACS galaxies; DJA photo-z for vetting (cycle 7)
- PR #16 (D-015) review fixes merged. Run `20261007T132538Z-96e911bb` (clean `725eac3`) reproduces the D-015 numbers.
- `scripts/vet_evidence.py --photoz`: the nearest DJA eazy entry, labelled model_prediction (SOURCES "Vetting
  (cycle 7)"). `/vet-candidate` now names the script and the real `candidates` subcommands.
- All four are ordinary; none is an artifact. Records are in `docs/candidates/`.

  | Source | Rank | Verdict |
  |---|---|---|
  | `829` | 3 | red, dusty galaxy at spectroscopic z = 2.74 (NED), a published F150W-dropout "z ≈ 11–20" candidate |
  | `1032` | 8 | lensed arc `[MJR2023] 028.2` (Mahler et al. 2023), a tangential pair with 028.1; blue colours likely from emission lines (hypothesis) |
  | `438` | 11 | uncatalogued; Balmer-break galaxy at z_phot 4.8 (hypothesis, 6-band photo-z) |
  | `1243` | 12 | catalogued red point-like high-z candidate, photometric z ≈ 5.6–5.75 |

- **Lesson:** after D-011/D-014/D-015 the top of the galaxy ranking is dominated by real but rare populations
  (lensed arcs, dusty and high-z galaxies), most already catalogued. The pipeline recovers known lensing
  features. Uncatalogued sources such as `438` are the discovery set.
- **Side finding:** the D-015 stellar locus missed the PSF-like star `f200w_1874` (F200W 21.4 mag, spikes)
  because its r50 (2.10 px) is 27% below r50_psf (TASKS Now 1).
- **Handoff:** TASKS Now 1 (stellar-locus lower bound) and Now 2 (blended stars).


## 2026-10-08: Stellar locus for stars missing from Gaia (D-015, cycle 6)
- `classify.stellar_locus`: a source counts as a star when its DJA detection-image r50 is within 20% of the
  catalogued stars' median and its colours are stellar. The DJA join now carries `r50_pix` and `mag_auto`.
- Run `20261007T130207Z-c8ffc470`: the locus added 16 stars to SMACS, and 3 of the 6 PSF-like sources left the galaxy top 20.
  `829`, `1032` and `438` are ranks 3, 8 and 11. Details are in D-015.
- **Failed approaches:**
  - pipeline CI and `semimajor_sigma` do not separate stars from galaxies, because spikes and
    saturation inflate them;
  - size alone admits compact galaxies, because their colours are too broad.
- **Handoff:** TASKS Now 1 (vet `829`, `1032`, `438`, `1243`). Three blended stars remain in the galaxy top
  20 (TASKS Now 2).


## 2026-10-07: Rank only confirmed detections (D-014, cycle 5)
- Quality gate: best-band S/N floor (`low_snr`) and single-band confirmation (`single_band`, confirmed by a
  DJA match). Too few survivors now fall back to the D-011 tests instead of ranking ungated.
- S/N is now inverted exactly from the pipeline's `abmag_err` everywhere: the features S/N ≥ 3 gate,
  the quality gate and the scripts. The linear approximation overstated S/N at low S/N.
- Run `20261007T124944Z-bbad4ab3`: spike and stripe detections in the SMACS galaxy top 20 fell from 5–6 to 1, and near-noise
  sources from about 6 to 0. `[YML2023] F150DB-C-4` and `[MJR2023] 028.2` are ranks 4 and 9. Evidence and
  counts are in D-014.
- **Failed approach:** an S/N floor in the reference band only removed red dropouts, including `829`
  (D-014).
- **Handoff:** TASKS Now 1 (stellar locus: about 5–6 faint PSF-like stars remain), Now 2 (vet `829`,
  `1032`, `438`, `1243`).

## 2026-10-07: Size-robust colours from DJA matched apertures (D-013, cycle 4)
- New `photometry.py`. The DJA v7.4 catalog is verified against its sha256 before it is cached, and is
  joined one-to-one within 0.2″. Colours now come from 0.5″ matched apertures.
- `features`: `ref_mag` and the morphology features keep the pipeline aperture's S/N gate, so sources
  without a DJA match keep their morphology. With matched photometry, a band counts as detected when DJA
  measured it at S/N ≥ 3. DAOFind sharpness/roundness apply only when 0 < CI_50_30 ≤ 1.8.
- **The size bias is gone, measured on the same 1,423 rows** (run `20261007T122054Z-c52935ec`, F200W−F277W, F200W aper50
  S/N > 10). Spearman between log area and colour:
  - aper50: +0.181 (p = 6.7e-12);
  - isophotal: −0.312 (p = 1.5e-33);
  - DJA 0.5″: −0.012 (p = 0.66).

  Re-measure with `python scripts/feature_size_bias.py --run-dir <outputs>/runs/20261007T122054Z-c52935ec --sample
  smacs0723_nircam --compare dja05`.
- Join coverage: 2,729 of 5,254 sources matched one-to-one. Nearest-neighbour matching gave 2,920, of
  which 191 were fragments sharing one DJA object.
- **Ranking:** the sources vetted in #13 dropped (`2925` to rank 1297, `2559` to 482, `1096` to 296). Two
  catalogued lens-related objects are in the galaxy top 20 (observed cross-matches):
  - rank 7 `f277w_829`, SIMBAD `[YML2023] F150DB-C-4` (also NED G_Lens);
  - rank 13 `f200w_1032`, NED `SMACS J0723-73:[MJR2023] 028.2` (G_Lens).
- **Top-20 purity is still mixed** (visual, unvetted):
  - 4 faint PSF-like stars missing from Gaia;
  - 5–6 detections on diffraction spikes or parallel stripes (a new artifact class the D-011 gate misses);
  - about 6 faint, compact or near-noise sources, one of them hot-pixel-like;
  - 3 arc-like sources;
  - 1 catalogued elongated galaxy.
- Failed approaches, both fixed after `/code-review` (15 findings, all addressed):
  - nearest-neighbour joining (fragments shared photometry);
  - gating morphology on DJA S/N, which erased the morphology of 44% of sources.
  The first run, `20261007T120528Z-378af5a1`, used both.
- **Handoff:** TASKS Now 1 (mask detections on diffraction spikes and stripes), Now 2 (S/N floor for
  ranking), Now 3 (stellar locus), Now 4 (vet the catalogued lensed images and arcs).

## 2026-10-07: First vetting: three arc-like galaxy candidates (cycle 3)
- New `scripts/vet_evidence.py` gathers vetting evidence:
  - per-band catalog rows;
  - radius and tangential alignment from both the catalog orientation and image moments (they agree
    to within 1–7°);
  - nearest multiple image in Mahler+2022 `arcs.dat`;
  - nearest star;
  - SIMBAD/NED/Gaia cross-match at 1″ and 3″;
  - six-band S3 cutouts.
- Vetted SMACS NIRCam `f200w_2925`, `f200w_2559` and `f200w_1096` (records in `docs/candidates/`; store status
  `explained`):
  - none is among the 62 Mahler+2022 constraint images (nearest 24–45″); this does not show they are
    singly imaged, which needs the model's critical curves (M3);
  - 2925 and 2559 are not tangentially aligned (45–60°); 1096 is at 19°, consistent with weak shear
    (`hypothesis`);
  - NED knows 2925 (z 1.98) and 1096 (z 1.36, with catalogued clumps) as background galaxies.
- **Systematic finding (derived):** the baseline's colour and shape features are size-biased.
  - aper50 F200W−F277W reddens with isophotal area (median +0.77 → +0.99, Spearman 0.21,
    p ≈ 1e-16), while isophotal colour does the opposite (+0.71 → +0.02). EE apertures differ in angular
    size between SW and LW.
  - Reproduce with `python scripts/feature_size_bias.py --run-dir <outputs>/runs/20261007T040938Z-01527ace
    --sample smacs0723_nircam`: sources detected in F200W and F277W with F200W aper50 S/N > 10, F200W
    isophotal-area bins [0, 50, 200, 1000, ∞) px.
  - DAOFind sharpness reaches +145σ on clumpy extended galaxies.

  This explains these three candidates (verdict `catalog effect`) and sets the next priority (TASKS Now 1).
- **Handoff:** TASKS Now 1, size-robust features (DJA matched-aperture colours; sharpness only for compact
  sources).

## 2026-10-07: Stars ranked as their own stratum (D-012)
- New `classify.classify_sources`. One bulk Gaia DR3 + SIMBAD cross-match over the whole catalog
  (about 8 s per sample) labels stars. The runner ranks them as `<sample>-stars` with their own top k;
  galaxies are ranked without them.
- Failed approach: the pipeline's `is_extended` flag labels every bright, saturated star in the top 20
  as extended.
- Real data, run `20261007T040938Z-01527ace`:
  - Stars found: SMACS 52, MIRI 12, CEERS 7. Of these, 40, 8 and 7 passed the gate and were ranked as
    `-stars`.
  - Galaxy-stratum top 20s: 0 cross-matched stars (SMACS had 7), 1/60 image-quality flags.
- Visual check (unvetted): the SMACS galaxy top 20 is now mostly interacting, clumpy or elongated
  galaxies, including three arc-like sources (`f200w_2925`, `f200w_2559`, `f200w_1096`). About 4
  faint PSF-like sources without Gaia counterparts remain.
- `/code-review` found 15 issues, all fixed. The main ones:
  - a classification leaving fewer than 2 galaxies could abort the run;
  - stars were dropped entirely when there were fewer than `min_stars`;
  - the offline suite made live CDS calls;
  - the rule used any match instead of the nearest one.
- **Handoff:** TASKS Now 1 (vet the arc-like sources against published lens models), then Now 2
  (stellar locus for faint stars).

## 2026-10-07: Quality gate before ranking (D-011)
- New `quality.assess_sources` with `cutouts.sample_weight_map`. A coarse WHT map (1″ cells,
  concurrent S3 byte ranges, about 10 s per 1.8 GB mosaic) gives each source's relative weight and
  edge distance. A CI_50_30 test
  catches sources sharper than the PSF. Ranking now sees only sources that pass.
- Real data, run `20261007T033832Z-cfcf6032`: top-20 sources with a cutout image-quality flag fell from
  20/60 to 2/60.
  - SMACS NIRCam: 6→1. MIRI: 10→0. CEERS: 4→1.
  - Gated out: 28% of SMACS NIRCam, 55% of MIRI, 28% of CEERS.
- `/code-review` found 15 issues in the first version, all fixed:
  - edge distance was snapped to cells, so `edge` never fired for LW or MIRI;
  - a gate passing too few sources could abort the run;
  - the gate ran even when unconfigured;
  - noise was labelled as artifacts.
- Visual check (unvetted): the SMACS NIRCam top 20 is now about 12 bright stars plus 8 galaxies, which
  include interacting pairs and the two arc-like sources (`f200w_2925`, `f200w_2559`). Low-weight
  noise and streaks are gone. Stars dominate, so star/galaxy separation is next (TASKS Now 1).
- The agent self-merge policy was adopted (owner decision; CLAUDE.md "Merge policy and version
  control"), and PRs #1–#10 were merged.
- **Handoff:** TASKS Now 1 (star/galaxy separation).

## 2026-10-07: M0 bootstrap batch landed; first real-data run (M1 slice)
- Nine parallel units became PRs #1–#9: archive, catalog, rank, cutouts, crossmatch, runner/CLI, OSS,
  agent harness, landscape. The integration branch merges all of them and adds a contact sheet to the
  runner. 297 offline tests pass.
- First end-to-end run on real public data: run `20261007T020124Z-4bfabaaa`, config `reference_sample_v1`,
  about 3 minutes on the owner's laptop.
  - Merged sources: SMACS 0723 NIRCam (6 bands) 5,254; MIRI (4 bands) 530; CEERS t021 (7 bands) 7,582.
  - The top 20 per sample get cutouts through S3 byte-range reads (no full i2d download) and a
    SIMBAD/NED/Gaia cross-check.
- Results. These come from the integrating agent's visual classification of the contact sheets
  plus the cross-match; nothing is vetted.
  - SMACS NIRCam top 20: mostly bright stars with diffraction spikes (7 confirmed by Gaia/SIMBAD,
    about 11 PSF-like), 2 image artifacts (a streak and low-weight noise), and galaxies or blends.
    Ranks 6 and 19 are elongated and arc-like; they are hypotheses to vet against published lens models.
  - MIRI top 20: 10/20 carry image-quality flags (edge or low weight). One known lens-related object
    (NED `G_Lens`) was recovered; it is a known object, not a discovery.
  - CEERS control: with F200W only, about 16/20 were artifacts (edges, hot-pixel-like compact sources,
    linear streaks). With all 7 bands, about 5/20 are artifacts and the rest are structured galaxies, so
    cross-band features matter.
  - Interpretation: the baseline mostly surfaces stars and instrument or processing artifacts, as
    expected. Quality gating and star/galaxy separation come next (TASKS Now 1–2).
- Injection-recovery (unit 3, NIRCam): only LOF clearly beats random on every injection type. The
  combined score reaches precision@50 of 0.06–0.20. Nothing is tuned toward one injection type until
  vetted labels exist (D-004).
- Data caveats: program 2736 catalogs were made with jwst 2.0.1 and CEERS with jwst 3.0.0, and deblending
  changed in between (D-010). The CEERS `s_ra/s_dec` is the MIRI target; use `s_region` instead.
- Process lessons:
  - `setup-uv` has no floating major tag, so pin the full version.
  - Workers collided in the shared scratchpad; each now gets its own subdirectory (CLAUDE.md).
  - Unit 5's `gh pr create` was blocked, so the coordinator opened #9 with an inline body;
    `--body-file` is denied by project settings.
- Guardrails:
  - `main` is protected: a PR and 4 CI checks are required, admins included.
  - A ruleset blocks force-push on every branch.
  - Private vulnerability reporting and Dependabot security updates are on.
  - Direct commits to `main` happened only in Step 0, before protection existed.
- **Handoff:** the owner merges #1–#9 in any order, then this integration PR (updated by merging `main`
  into it). The next cycle starts at TASKS "Now" item 1.

## 2026-10-07: M0 bootstrap started
- Created repository skeleton: interfaces (`schema.py`, stage stubs with fixed signatures), state files,
  CI, BSD-3-Clause license, verbatim agent charter (`docs/agent-charter.md`) imported by `CLAUDE.md`.
- Verified on MAST: program 2736 level-3 NIRCam/MIRI imaging is public; catalogs ~3 MB vs NIRCam i2d
  ~1.8 GB, hence the catalog-first slice (D-001). The S3 mirror key pattern is confirmed.
- Owner decisions: public repo; every change via PR and the owner merges; chat in Chinese.
