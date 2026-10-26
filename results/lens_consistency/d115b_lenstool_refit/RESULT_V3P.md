# D115b stage V3′: small-step gate V3 (result, 2026-10-26)

Plan `PLAN_V3P.md` (committed 16:13Z, fd23864e); tool `scripts/d115b_lenstool_refit.py --stage v3p`, `choose_start`
and `gate_v3p` tested in `tests/test_d115b_lenstool_refit.py` (17badd05), both before the run. Numbers:
`result_v3p.json` (gate), `posthoc_v3p.json` (post hoc, `--stage v3p-posthoc`). All χ² are `model_prediction` of
variants of the pinned `macs0416-canucs` model (D-044) with 8515 at its `best.par` σ; no D115b number follows.

## Result
- **Gate V3′: `uninformative`.** Every ladder start (0.001, 0.002, 0.003, 0.005, 0.0075, 0.01, 0.015, 0.02 of each
  range, alternating sign over the 28 O1–O6 parameters) is infeasible, so no optimiser run took place.
- Per the frozen plan, the hard-constraint local route (Nelder–Mead on Lenstool `inverse 0` χ² with lost or merged
  images rejected) is closed for this model.

## Post hoc (after the gate; `posthoc_v3p.json`)
- Tool check: `best.par` through the tool gives 344.42, feasible.
- **The χ² is discontinuous at `best.par` itself.** Single + offsets shrunk from 10⁻³ to 3 × 10⁻⁵ of the range change
  family D52b's χ² by a fixed +12.1 (D52c +14.2, D34 −10.0 in some), not smoothly:
  - O2 angle_pos: +12.12 already at 3 × 10⁻⁵; O3 v_disp: +12.13 at 3 × 10⁻⁵;
  - O1 ellipticity: 0 at 10⁻⁴, +12.11 at 3 × 10⁻⁴; O5 cut_radius and O6 v_disp: 0 up to 3 × 10⁻⁴, +12.1 at 10⁻³.
  - Reading (not tested): a D52b image prediction switches branch, so `best.par` sits on a branch edge. A 1 %
    tolerance (3.44) is smaller than the jump, and a simplex step of any size crosses it.
- Single ±10⁻³ offsets: 33 / 56 feasible (3 % offsets: 7 / 56); only 12 / 56 stay within 1 of 344.42.
- The all-parameter start at 10⁻⁴ (below the frozen ladder) is feasible, χ² 357.66 (the D52b branch).

## What this changes
- Local re-fits from `best.par` with a hard image-loss rule and a χ² tolerance below ≈ 12 cannot be validated on this
  model: the published optimum is on a discontinuity. The next stage needs a different handling (proposed, not
  frozen): a continuous objective (e.g. Lenstool's source-plane χ², or image-plane χ² with D52b's image matching
  fixed), validated first. D115b stays `inconclusive`.

## Deviations
- The micromamba binary came from micro.mamba.pm piped into tar (worked here; refused in an earlier run, D-109 addendum 2).
- archive.stsci.edu now serves a chain ending at "emSign Root CA - G1", absent from this container's system store; the
  first two launches failed on TLS before any evaluation, and the run used `SSL_CERT_FILE` = the proxy bundle +
  certifi's roots (scratchpad file).
- The post hoc writer (`posthoc_v3p`) was added after the gate was read.
