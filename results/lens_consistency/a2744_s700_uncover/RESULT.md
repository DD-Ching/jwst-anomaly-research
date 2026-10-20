# Abell 2744 system 700.1 in the UNCOVER v2.0 maps (result)

Plan: `PLAN.md` (frozen and committed before the run). Script: `scripts/a2744_s700.py --model uncover2`. Numbers:
`summary.json` (`known_cases` per radius with `failed_detail` / `n_failed`, `outcome`, `uncover_counterparts`,
`known_excluded`, `target_systems_fit`, `target_pairs`). Computed positions are `model_prediction`; catalogued
positions and redshifts are `observed` (catalogue values); pair geometry is `derived`.

**Frozen outcome: "tool not validated"; no z scan was run.** On UNCOVER's own list (60 systems with z > 0 and ≥ 2
images, each at its own z, without the target's source 78.1–78.3) the one-to-one rule reproduces 31 / 60 (0.517) at
1.0″ and 44 / 60 (0.733) at 1.5″, both below the 0.8 gate (ASSUMPTION). At 1.5″ the 16 failures are 0 missing, 6
shared (every catalogued image has a predicted image within 1.5″, but two share it) and 10 offset (some catalogued
image's nearest predicted image is 1.77–2.31″ away). So the rule is too strict for a map model with a 0.84″
own-list rms (D-041), and a shared predicted image also occurs in known UNCOVER systems.

**Post hoc findings (not part of the frozen rule).**
- *Identity.* 700.1a/b/c are UNCOVER's 78.21 / 78.22 / 78.23 (0.02″ each; `uncover_counterparts`). The UNCOVER table
  lists system 78 (knots 78.1–78.3) at z 2.3151, `Redshift_type` spec, reference "Price+24", JWST/NIRSpec. The paper
  and spectrum were not checked here.
- *UNCOVER's own fit of the target* (`target_systems_fit`; system 78 is one of its constraints, so this is a fit,
  not a prediction): at z 2.3151 every knot gets 3 predicted images, but the a and b images of each knot share one
  (`nearest_pred` [0, 0, 2] for 78.1 and 78.2; [0, 0] for 78.3). 78.2: a 0.56″, b 1.34″ from the shared image, c
  0.99″ from its own.
- *Bergamini+23 at the spectroscopic z*: at z 2.30–2.325 it fails as at every z (rms 10.86–10.90″;
  `../a2744_s700_zscan/summary.json`).
- *Pair geometry* (`target_pairs`): the a/b pairs of the three knots are nested, with separations 0.50″ (78.3),
  1.21″ (78.1) and 1.74″ (78.2), and midpoints agreeing within 0.05″ (x 0.54–0.59″, y −0.16 to −0.21″ from 78.21).

**Reading.** The redshift is not free: 2.3151 per the UNCOVER table. At that z, neither model gives a and b
separate images, although both use the system as a constraint (Bergamini+23 with its z free, best fit 1.217;
UNCOVER at 2.3151). The nested, mirror-symmetric knot pairs are what a fold arc with a critical curve at the common midpoint
gives. They are also what one edge-on galaxy with symmetric knots gives. So 700.1 is undetermined: either a
critical curve that both models misplace, which an unmodelled local perturber could cause (ordinary, untested),
or a/b are one image. It is not a candidate. Under #412's frozen plan, B there plus UNCOVER also failing makes it
worth `/vet-candidate`. The UNCOVER failure here is post hoc and the rule failed its gate, so this only motivates
that vetting step; it does not satisfy the plan's condition.

**Deviations (labelled).** Post hoc, after the first run: the per-system `failed_detail` / `n_failed` fields, the
`uncover_counterparts` look-up, the exclusion of the target's source (78.1–78.3) from the known cases (63 → 60
systems; with it the fractions were 0.492 / 0.698, so the outcome does not change), `target_systems_fit` and
`target_pairs`. The plan's timestamp line was corrected to the commit time (21:12Z) before results were written.

**Limitations.** One rule (one-to-one, fixed radii). UNCOVER's fit to 78 is not independent. The redshift comes
from the UNCOVER table, not from the paper or the spectrum. No cluster-member or perturber check was made near the
arc.
