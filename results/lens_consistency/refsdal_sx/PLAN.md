# SN Refsdal (MACS J1149.5+2223): measured SX delay and magnification ratio vs the blind predictions: plan

Frozen 2026-10-26 05:16Z, before the model table (Treu et al. 2016, ApJ 817, 60, arXiv:1510.05750, the blind
SX predictions compiled before the reappearance) and the measurement (Kelly et al. 2023, Science 380, eabh1322,
arXiv:2305.06367) are read in this session.

## Question
Do the measured SX − S1 time delay and SX / S1 magnification ratio of SN Refsdal agree with each blind model in
Treu et al. 2016? Unlike SN Requiem, the delay is measured to a few per cent, so this is the one cluster-scale
case where a blind delay prediction can fail at > 3σ on its own.

## Ordinary explanations of a disagreement
Model systematics (the published spread between models; H0 and cosmology fixed at each model's fiducial value),
microlensing of S1 (flux ratio, not delay), substructure near SX or S1, measurement error.

## Frozen rule (ASSUMPTION; same structure as D-034 addenda bn–bp)
- Delay: z_m = (Δt_obs − Δt_m) / sqrt(σ_obs² + σ_m²), linear, σ as half the 16–84 (or ±1σ) width of each
  published interval; the models' fiducial H0 is used as published (no rescaling; a rescaling would be post hoc).
- Magnification ratio: z_m in ln space, ln(r_obs / r_m) / sqrt(σ_ln,obs² + σ_ln,m²), σ_ln = half the ln width.
- A model is flagged when |z| > 3 for either quantity. Outcome F (vet first) only if every model with a published
  prediction is flagged; N if at least one model is within |z| ≤ 3 on both; U if no model gives both quantities.
- Report per-model z, the number of models with |z| > 3 (as model ranking, not an anomaly), and the reach: the
  smallest delay offset (days) and ratio factor exp(3 σ_ln,obs) a flag needs.
- Only models whose numbers were published before the reappearance (blind) enter the outcome; later or revised
  models are listed separately, post hoc.

## Inputs
arXiv e-prints only (sha256 recorded); script `scripts/refsdal_sx_check.py` with rule functions and tests written
before it is run on the transcribed tables.
