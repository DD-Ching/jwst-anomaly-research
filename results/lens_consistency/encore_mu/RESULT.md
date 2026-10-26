# SN Encore (MACS J0138.0−2155): standard-candle magnifications vs the seven blind lens models: result

Plan `PLAN.md` (frozen before any table was read). Published tables only, no images (so there is no contact sheet;
the tables were read in the e-prints). Numbers: `summary.json` (frozen rules) and `posthoc.json`, both written by
`scripts/encore_mu_check.py`. Labels: observed (Pierel et al., arXiv:2509.12301: SN Ia standard-candle μ of images
1a and 1b, BayeSN, host E(B−V) < 0.1, and Δt_1b,1a), simulated (their micro- and millilensing terms),
model_prediction (Suyu et al., arXiv:2509.12319, Table "tab:snencore", the seven blind models), derived (z, offsets,
λ, H0_MST), assumption (thresholds in the plan).

## Inputs as used
- Observed μ: the final table (1a 21.8 +9.3 / −8.6, 1b 32.4 +11.1 / −11.0). Its intervals already contain the
  published micro- and millilensing terms, so the plan's σ_ml is not added a second time (reading of the plan's
  "published term" branch, stated here). σ_ln 0.428 (1a), 0.355 (1b).
- Seven blind models: glafic, GLEE, Lenstool I, Lenstool II, MrMARTIAN, WSLAP+, Zitrin-analytic (the set Pierel et
  al. Table "tab:h0" combine). The GLEE-baseline and the two post-blind columns are not used. |μ| is compared (1a is
  a saddle point, negative in every model).
- No term shared by both images is published separately, so the ratio test adds every error in quadrature (the
  plan's conservative branch).

## Frozen outcome: **N** (no image flag, no ratio flag)
- Image z: −1.75 (WSLAP+, 1a) to +1.36 (Zitrin-analytic, 1b); the rule needs |z| > 3 in all seven.
- Ratio |μ_1b/μ_1a|: observed 1.486; models 0.60 (WSLAP+) to 1.76 (glafic); z −0.29 to +1.33.
- **Power (derived):** with a zero model error a flag needs |μ_obs/μ_model| beyond 3.62 (1a) or 2.90 (1b). The
  test can only catch gross disagreements; that is the micro/millilensing budget (millilensing ≈ 29 % per image at
  f_sub 5 %, simulated; Pierel et al. text, not stored).

## Offset (the H0pe question)
- Per-model offset 2.5 log10(|μ_model| / μ_obs), mean over 1a, 1b: −0.405 (Zitrin-analytic) to +0.472 mag (WSLAP+);
  pooled +0.086 mag; 5 of 7 positive. Frozen reading: **same sense as H0pe** (pooled > 0 and ≥ 5 of 7). The pooled
  value is 0.086 against H0pe's 0.633 mag (bn), and no per-model offset exceeds about one σ_obs (0.355–0.428 in
  ln μ, i.e. 0.39–0.46 mag; conversion not stored), so the sign agreement is not a measured shared
  over-prediction (not tested further).

## Mass sheet
- H0_MST (one λ for μ and Δt_1b,1a, fiducial H0 70): glafic 66.8, GLEE 73.2, Lenstool II 65.2, MrMARTIAN 75.4 (in
  64–80); Zitrin-analytic 58.7, Lenstool I 116.6, WSLAP+ 244.8 (outside). Frozen: **4 of 7** in range.
- Post hoc 1σ intervals (`posthoc.json`; observed μ errors fully correlated across images): 6 of 7 reach 64–80; WSLAP+
  (167.8–357.2) does not; its Δt_1b,1a −112 ± 32 d is 2.2σ from the observed −39.8 d (not stored as a z;
  the H0 paper gives it a weight of 10⁻⁴⁶³, Pierel et al. Table "tab:h0", not stored).

## Post hoc checks (`posthoc.json`)
- Statistical μ errors with the two simulated terms added in quadrature instead of the final intervals: image z
  −1.77 to +1.36, ratio z −0.29 to +1.34; outcome unchanged.
- The directly fitted ratio 1.5 ± 0.5 (SNTD, statistical only): z −0.42 to +1.75; no flag.

## Reading (no exotic interpretation)
SN Encore's two image magnifications and their ratio agree with all seven blind models within errors dominated by
simulated micro- and millilensing. The H0pe shared over-prediction (bn) is neither confirmed nor excluded here:
the errors are 0.355–0.428 in ln μ per image against an H0pe offset of 0.633 mag (0.58 in ln, not stored). The
docs/exotic_lensing.md "demagnified images" signature is not tested beyond a factor of about 3. Not a candidate.

## Deviations
- The rule tests (`tests/test_encore_mu_check.py`) were written after the first run (scripts/CLAUDE.md asks for
  them before); the rule functions were in the script at the first run and are unchanged by the tests.
- The plan's frozen-time stamp was corrected from 03:25Z to 03:10Z (typo; the plan was committed at 03:10Z,
  before the e-prints were downloaded).

## What would decide it
Template images of the host (approved HST/JWST follow-up, per Pierel et al.) shrink the photometric term; the
reappearance of image 1d (predicted Δt 2,994–3,544 d after 1a in the seven blind models; Suyu et al. table, not stored) gives a second delay. Rerun
`scripts/encore_mu_check.py` with the new tables.

## Cost
Two arXiv e-prints (≈ 40 MB, not stored; scratch, deleted after use), < 1 s compute.
