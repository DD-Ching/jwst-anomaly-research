# Plan: does the RELICS GLAFIC v3 model treat RXC J0600 581 as a cluster member? (frozen 2026-11-08)

Context: D-024 addendum 18 found, post hoc, a galaxy-scale κ peak (5.24) 0.2″ from 581 in the GLAFIC v3 best fit.
No GLAFIC member list is public (RELICS `models/glafic/v3/` holds maps only; headers carry no input list, checked
2026-11-08). This plan tests the peak as a membership proxy on galaxies with spectroscopic redshifts.

- **Hypothesis H0 (ordinary):** GLAFIC puts galaxy-scale haloes on the cluster members it selected; a peak at 581 then
  means the model treats 581 as a member, so its GLAFIC axis prediction at 581 is not a test of a background arc.
- **Alternative:** peaks also sit on background galaxies or appear at random, so a peak says nothing about membership.
- **What changes the plan:** gate passes and 581 has a peak → 581's GLAFIC row is a model-input member, the open
  question for 581 is membership alone (spectroscopy), and Lenstool is the only model treating it as background.
  Gate fails → the κ route is closed for membership (`undetermined`).

## Inputs (read in the scratchpad, deleted after use; nothing stored)
- GLAFIC v3 best-fit κ (`hlsp_relics_model_model_rxc0600m20_glafic_v3_kappa.fits`, 2400², 0.1″); secondary:
  Lenstool v1 κ (1000², 0.2″). RELICS HLSP, DOI 10.17909/T9SP45. sha256 recorded.
- NED `objdir` (TAP) within 2′ of (90.0340, −20.1358) with z not null. Spectroscopic = `zflag` in {SLS, S1L, SST,
  IST} (no `?`); only `prefphytype` G or G_Lens. Deviation (labelled): this NED table was fetched once before this
  plan to count spectroscopic rows (78 G / G_Lens SLS + S1L); no κ value had been read.
- 581 position: stored GLAFIC run (`../rxc0600-relics/rxc0600-relics-glafic/arcs.ecsv`).

## Frozen rules (ASSUMPTIONs)
- Peak statistic at a position: in a 2.0″ × 2.0″ box centred on it, the κ maximum lies ≤ 0.4″ from the position and
  max / box median ≥ 1.5 → `peak`; else `no_peak`. Positions within 1.0″ of the map edge are skipped.
- Member: |z − 0.43| ≤ 0.03 (1 + 0.43). Background: z > 0.43 + 0.03 (1 + 0.43). Foreground (z lower) is not used.
- Gate G (per map): n_members ≥ 8 and n_background ≥ 8; member peak recall ≥ 0.70; background peak rate ≤ 0.10.
- Label for 581 (GLAFIC, the deciding map): gate passes and `peak` → `GLAFIC member (model input)`; gate passes and
  `no_peak` → `GLAFIC background (model input)`; gate fails → `undetermined`. Lenstool is reported with the same rules,
  secondary. Labels are about model inputs, never about 581's true redshift.
- Known before this plan (addendum 18, post hoc): 581's GLAFIC peak numbers. The new information is the validation.
