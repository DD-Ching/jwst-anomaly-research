# D115b Lenstool: D14a's family-row χ² term is Lenstool's image-proximity penalty (result, 2026-10-26)

Plan `PLAN_D14A.md` (frozen 12:21Z before any run of this round; it states what was already known from reading
Lenstool's source and the stored `chires_f0.90.dat.gz`, so rule P is a check, not a blind test). Tool
`scripts/d115b_lenstool_refit.py --stage d14a`, rule `rule_p` tested in `tests/test_d115b_lenstool_refit.py`, both
committed before the run. Numbers: `result_d14a.json` (rule P), `posthoc_d14a.json` (post hoc). All χ² and positions
are `model_prediction` of the pinned `macs0416-canucs` model (D-044) with member 8515's σ scaled by f; the image list
is `observed`.

## Result
- **Rule P: `confirmed`** (`rule_P_failures` empty). At every f ≤ 0.90 Lenstool printed `chipenalty 3748.44` (one
  unit, (30″ / 0.49″)² = 3748.4382, `penalty_unit`), and its penalty warning names a D14a image (D14a.1 at f 0.90 and
  0.75, D14a.2 at 0.85, 0.80, 0.70). D14a's family row minus its image rows equals one unit within 0.05 at every f.
  At f 1.00 there is no penalty.
- **The penalty is a code rule, not a position residual.** In `src/o_chires.c` (`chi2_img`, Lenstool 8.15.5 source,
  sha256 = the conda-forge recipe's) Lenstool adds (30″)² / σ² to a family when two of its predicted images coincide
  within 10⁻³″ in x and y. It is added to the family row only, so `chimul` jumps by 3748.44 at once.
- **Profile without the penalty** (`dchi2_no_penalty`): +11.28 (f 0.90), +7.28 (0.85), +10.33 (0.80), +15.04 (0.75),
  +21.15 (0.70). It still scores the lost D115a images at χ² 0, so it is **not** the cost of removing the D115b image
  (plan: reported, not gated).

## Post hoc (after rule P; `posthoc_d14a.json`)
- **D14a.1 is traced onto D14a.2's image.** From f 0.90 down, the predicted positions of D14a.1 and D14a.2 agree to
  0.006–0.008″ (limited by chires' 2-decimal offsets); at f 1.00 they are 1.259″ apart. So in Lenstool's
  fixed-parameter χ² the model no longer reproduces D14a.1 at its own position once 8515's σ drops by 10 %. Its own
  row prices that only by the distance to D14a.2's image (χ² 8.88 at f 0.90, 0.31 at f 1.00), with no warning.
- This is a second lost-image effect besides D115a's, and it starts at f 0.90, above Lenstool's f_L* 0.80 for the
  D115b flag image. Not tested here: whether a re-fit with the cluster halos free keeps D14a.1 at lower σ(8515).

## Reading
D14a's ≈ 3748 is Lenstool's image-proximity penalty, triggered because D14a.1 and D14a.2 back-trace to one image.
Removing it leaves Δχ² 7.28–21.15 over the grid, but that still does not count lost images, so the stage F reading
(not a usable cost; D115b stays `inconclusive`) is unchanged. An F2 optimiser run will meet the penalty as a 3748
step wherever D14a.1 merges into D14a.2, so its plan must state how lost and merged images are counted.

## Deviations
- The plan's first time stamp ("12:40Z") was a typo, corrected to 12:21Z in the tool commit, before the run.
- The first three invocations stopped before any Lenstool run: archive.stsci.edu needs the proxy bundle and
  certifi concatenated in `SSL_CERT_FILE` (as in earlier cloud runs). Only the fourth invocation ran.
- The post hoc writer `posthoc_d14a` and its test were added after rule P was read.
