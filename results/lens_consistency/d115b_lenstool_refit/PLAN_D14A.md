# D115b Lenstool: D14a's family-row χ² term (plan, frozen 2026-10-26 12:40Z, before any run of this round)

Continues `RESULT_VPRIME_F.md` (post hoc: at f ≤ 0.90 D14a's family row is ≈ 3748 above its image rows; cause not
diagnosed). Tool `scripts/d115b_lenstool_refit.py --stage d14a`, rule `rule_p` tested in
`tests/test_d115b_lenstool_refit.py`, both committed before the run. Lenstool 8.15.5 (conda-forge, D-109).

## What is already known before this plan (stated so the test is not mistaken for a blind one)
- **Code reading** (Lenstool 8.15.5 source, `src/o_chires.c` `chi2_img`, tarball sha256 = the conda-forge recipe's
  `b6abd0c8…6901`): in the image-plane χ² Lenstool adds a "chi2 proximity bug fix" penalty (30″)² / σ² to the family
  χ² (`chipi`), not to any image row, whenever two predicted images of one family coincide within 10⁻³″ in both x
  and y; it prints `chipenalty <total> <rest>` and a stderr line `WARNING: chi2 penalty for image <id>`.
- With σ = `sigposArcsec` 0.49, (30 / 0.49)² = 3748.44. The stored `chires_f0.90.dat.gz` already has the line
  `chipenalty 3748.44 355.70` (seen while writing this plan).

**Hypothesis H_P:** the D14a excess at every f ≤ 0.90 is this penalty: two of D14a's predicted images (expected: the
close pair D14a.1 / D14a.2, about 1″ apart) are back-traced onto the same image. Ordinary alternative: a genuine
image-position term, or a second effect only at some f.

## Run
Each f ∈ {1.00, 0.90, 0.85, 0.80, 0.75, 0.70}: the stage F run file (unchanged), stderr kept. Per f
(`model_prediction`): χ²_pos, D14a family row and image-row sum, `chipenalty` total P(f) (0 if the line is absent),
image ids in the penalty warnings, Δχ²_np(f) = [χ²_pos(f) − P(f)] − [χ²_pos(1.00) − P(1.00)].

**Rule P (frozen, `rule_p`):** H_P `confirmed` if, at every f, P(f) = k · 3748.44 for an integer k ≥ 0 (|·| ≤ 0.02),
and the D14a family row − its image-row sum equals P_D14a(f) within 0.05 (rounding); the D14a share counts the
penalty warnings that name a D14a image, times 3748.44. Otherwise `rejected`, and the residual is reported.

## Reported, not gated
- Rule F's arithmetic on Δχ²_np(f_L*) with the stored f_L* = 0.80. It still scores lost images (D115a.1 / D115a.2)
  at χ² 0, so it is **not** a cost of removing the D115b image; it only says what remains once the penalty is out.

## What changes the plan
- `confirmed` → the F2 plan must handle the penalty (it acts as a 3748 cliff in any optimiser run that merges
  D14a.1 / D14a.2) and must count lost images; the next step is that frozen F2 plan.
- `rejected` → D14a's term is diagnosed further before any F2 work.
