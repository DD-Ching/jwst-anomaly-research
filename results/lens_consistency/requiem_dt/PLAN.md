# SN Requiem (MACS J0138.0−2155): measured time delays vs the seven blind lens models: plan

Frozen 2026-10-26 03:26Z, before Rodney et al. 2021 is read in this session. The model table (Suyu et al.,
arXiv:2509.12319, Table "tab:snrequiem") has not been read yet either; only the SN Encore table was (D-034 addendum bo).

## Question
Do the measured delays of SN Requiem images 2b and 2c relative to 2a (Rodney et al. 2021, Nature Astronomy,
arXiv:2106.08935) agree with the seven blind models? The primary statistic is the H0-free ratio
R = Δt_2c,2a / Δt_2b,2a. Individual delays depend on H0 and are reported at the models' fiducial H0 70 only.

## Ordinary explanations of a disagreement
Measurement error (sparse 2016 HST light curves), microlensing (chromatic, affecting colour-curve delays), model
systematics (seven-model spread), substructure near one image.

## Frozen rule (ASSUMPTION; same structure as addenda bn / bo)
- z_m = (ln|R_obs| − ln|R_m|) / sqrt(σ_obs² + σ_m²), ln space, σ as half the 16–84 width in ln. σ_obs from the two
  measured delays in quadrature (no published covariance); σ_m likewise from the model's two delay intervals.
  If a delay is consistent with zero (interval crossing 0), the ratio is not defined in ln: outcome U for the ratio,
  and the per-delay test below decides.
- Per-delay test (H0 70): z on Δt itself, linear, symmetrised half-width errors in quadrature.
- Flag: |z| > 3 for every model (ratio or a delay). Outcome N, F (vet first) or U.
- Power: report the smallest ratio factor exp(3 σ_obs) a flag needs.

## Inputs
arXiv e-prints only, sha256 recorded; script `scripts/requiem_dt_check.py` with rule functions and tests written
before it is run.
