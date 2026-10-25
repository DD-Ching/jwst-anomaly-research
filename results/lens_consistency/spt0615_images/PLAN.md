# SPT-CL J0615−5746 published images vs the RELICS Lenstool v1 and GLAFIC v2 maps — frozen plan

Frozen and pushed before any SPT0615 map is read in this run (2026-10-25). Rules are those of
`../slice_relics_macs0358/PLAN.md` (D-034 addenda be, bh); the differences are listed. Script
`scripts/spt0615_images.py`; rule and image-list tests in `tests/test_spt0615_images.py`.

## Question and hypotheses
- **Q1 (amplitude; resolves addendum bk's "intermediate"):** addendum bk found the SPT0615 Lenstool v1 map **smaller**
  than GLAFIC v2 over |μ_ref(z 2)| ≥ 10 (R 0.772). Which map, if either, reproduces the published images at their
  spectroscopic redshifts? Rule 7 class per usable model; Q1 answer as in addendum bh (`H_sys`: Lenstool v1
  `lenstool_like` and GLAFIC `near_1`; `H_none`: both `near_1`; else `undecided`).
  - Ordinary expectation: both `near_1` (the HST-era maps were fitted to these images).
- **Q2 (positions):** does any published image miss every usable model by > 3″ (reading B: spectroscopic families at
  z_spec; reading A: the three less-secure system 3 images)?
  - Ordinary expectation: no miss in B; misses in A are possible (the paper's model 1 excludes them).
- **What would change the plan:** a usable model whose spectroscopic families need a factor ≠ 1 → that map is not used
  without the addendum-bd check (practical rule, not an anomaly). A joint miss in B → `redshift_residual_input`,
  checked first against the paper's own rms (Table 2) and the map / model-version match, before any interpretation.

## Inputs (observed)
Paterno-Mahler et al. 2018 (ApJ 863, 154; arXiv:1805.09834, e-print sha256 83c2bc6c…6cdb26, `spt0615_arxiv.tex`
sha256 25b66b8d…c807), Table 2 ("Properties of Lensed Galaxies"), transcribed by hand into `spt_images.ecsv`
(sexagesimal → degrees in the script, `TABLE2`). z_lens 0.972.
- **Family 1** (1.1, 1.2, 1.3; z_spec 1.358): clump 1 only; clumps 10–12 of the same source are not used (as SLICE
  clump 1 in addendum bh).
- **Family 2** (2.1, 2.2): no z_spec; table z 2.43 = the paper's model 1 value (model prediction; used only by rule 7).
- **Family 3** (3.1, 3.2, 3.5; z_spec 4.013): the secure images ("used in models 1 and 3").
- **Family 3b** (3.1–3.6, six images; reading A; `held_out`): the paper's model 1 omits 3.3, 3.4, 3.6 (less secure;
  3.6 photo-z 0.47 [0.15, 4.12]). Reading A is **exploratory, not held-out** (trial rule T5): whether the RELICS v1/v2
  maps used them is unknown (no constraint list in the HLSP directories). Read only for 3.3, 3.4, 3.6.
- **Known weaknesses, stated in advance:**
  - Families 1–3 are the models' own constraints, so K0 and reading B test reproduction, not prediction.
  - Which of the paper's four models the HLSP Lenstool v1 map is, is not stated in the HLSP; ASSUMPTION model 1.
  - A foreground structure at z ≈ 0.4 is not in single-plane maps (ASSUMPTION single plane at z_lens 0.972).

## Rules (ASSUMPTIONs; per model unless stated)
1–6 as addendum bh (fitted z on the 60-step grid; **gate K0** on families 1–3 at fitted z: rms ≤ 2.0″ and matched
fraction ≤ 3″ ≥ 0.8; **frame F** ≤ 0.5″ per axis; **gate P** κ and μ(z 2) ratios in [0.95, 1.05]; usable = P and K0
and F; **reading B** families 1 and 3 at z_spec, **reading A** 3.3, 3.4, 3.6 in family 3b at its fitted z; joint
flag = missed by every usable model; extra / fewer predicted images listed).
7. **Amplitude class** as addendum bh, with one difference (deviation from bh, fixed here): the median is over
   families 1–3 only (3b repeats family 3).

## Not done here
No imaging or photometry. A post hoc listing (labelled) only after the frozen run.
