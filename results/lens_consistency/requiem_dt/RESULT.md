# SN Requiem (MACS J0138.0−2155): measured time delays vs the seven blind lens models: result

Plan `PLAN.md` (frozen before either table was read). Published tables only, no images (no contact sheet).
Numbers: `summary.json`, written by `scripts/requiem_dt_check.py` (seeded draws). Labels: observed (Rodney et al.
2021 Table 1 Δt_obs: single-epoch HST colours plus a lens-model magnification prior from their own model, not one
of the seven), model_prediction (Suyu et al. Table "tab:snrequiem", the seven blind models of addendum bo, H0 70),
derived (z, ratio), assumption (|z| > 3).

## Inputs as used
- Image identity by position: Rodney SN1 / SN2 / SN3 = Suyu 2b / 2a / 2c (Table 1 coordinates against Suyu et al.'s
  offsets from the BCG: within 0.15″; positions not stored).
- Deviation (input form): Rodney et al. publish SN2 and SN3 relative to SN1, not to 2a. The two published delays are
  drawn independently (split normal, no covariance published) and combined: Δt_2b,2a −114.1 +31.0 / −27.8 d,
  Δt_2c,2a −129.3 +35.0 / −33.4 d. Their shared SN2 term makes the ratio tighter than either delay:
  R = Δt_2c,2a / Δt_2b,2a = 1.148, σ_ln 0.146.

## Frozen outcome: **N** (no flag)
- Ratio (H0-free): models 1.12 (Lenstool I) to 2.37 (MrMARTIAN); z +0.02 to −1.97. Every model except Lenstool I
  predicts a larger ratio than measured (z ≤ −0.54); a flag needs a factor ≥ 1.55.
- Delays at H0 70: Δt_2b,2a z −0.65 (Lenstool I) to −2.86 (MrMARTIAN): the measured delay is longer in every
  model, as Suyu et al. note (model delays about half the measured one). Δt_2c,2a z −0.61 to −1.57.

## Reading (no exotic interpretation)
No model is 3σ from the measured delays or their ratio. The consistent sign (measured 2b–2a delay longer, ratio
smaller, all seven models) is a 2σ-level pattern (largest |z| 2.86), not a flag. The delays rest on one epoch per
image and on the SN Ia assumption, and on Rodney et al.'s own magnification prior. Not a candidate.

## Deviations
- The input form above (delays relative to SN1, combined by draws).
- Rodney et al.'s lens-model-independent colour-curve ages exist only as a figure (Extended Data Fig. 5); figure
  readings are not used (D-083 T5).

## What would decide it
Image 2d reappearance (predicted 3,711–4,516 d after 2a in the seven models; Suyu et al. table, not stored) or a
re-analysis of the 2016 light curves with published delay covariance. Rerun `scripts/requiem_dt_check.py`.

## Cost
Two arXiv e-prints (≈ 25 MB, not stored; scratch, deleted after use), < 2 s compute.
