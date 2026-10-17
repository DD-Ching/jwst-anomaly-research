# W1 calibrated streak veto (plan, frozen before the screen and injection runs)

Cloud routine, 2026-10-17 (claim PR #307). TASKS "JWST focus" 1 (W1) and the handoff of
`results/w1_shear_rxj_m0647/RESULT.md`: a streak veto run through injections, then MACS1423 back into the limit.

**Hypothesis (H).** A negative-mass (W1) lens gives a radial aperture-mass peak (D-053 sign convention).
**Ordinary explanation tested here.** Linear image features (detector stripes, scattered-light or satellite trails,
long spike streaks) put spurious or contaminated shapes in a line through an aperture and mimic the W1 sign. Three
fields show this by eye or by a post hoc test: Abell S1063 (top peaks on parallel streaks), MACS1423 (the D-053 hit is
carried by 15 sources on a vertical stripe, `results/w1_shear_more/vet_macs1423.json`), MACS0647 (top centres on two
crossing streaks). The post hoc V1 test only looked at image columns near one peak and was not run through
injections, so it cannot be used for a limit (scripts/CLAUDE.md).

**Veto (frozen; every number an ASSUMPTION).** `scripts/w1_streak_mask.py`, one mask per field, from pixels only:
1. Stream the SCI extension of the field's F200W `_i2d.fits` (the mosaic of the screen's catalogue; S3, not stored),
   block-median 4 × 4 (zero fill ignored).
2. Normalise: (image − median) / (1.4826 × MAD); off-mosaic pixels → 0; pixels > 5σ and a 3-pixel dilation → 0
   (source mask); mean-rebin 2 × 2 (8 × 8 native pixels in all).
3. STScI `acstools.findsat_mrt.TrailFinder` (median Radon transform; Stark et al. 2022) with its defaults except
   `threshold=5`; keep trails with status 2 (passed the persistence check). Mask = union of their trail masks.
4. A shape row is vetoed when its position lies within 1″ of a masked pixel (the mask dilated by 1″).

The veto removes rows from the shape table before the screen (`inject_shear.py --streak-masks`), so the real
screen and every injection see the same reduced catalogue; painted lensed images that land in the mask are vetoed
too. Injections therefore carry the veto's cost (scripts/CLAUDE.md: a vetting test that removes signal is
calibrated).

**Known-case checks (computable; `w1_streak_mask.known_case`).**
- K1: the MACS1423 mask covers the 2026-10-17 hit position (`screen.s_max_radec` of
  `results/w1_shear_more/inject_summary_macs1423.json`), with the 1″ dilation.
- K2: the MACS1423 screen with the veto does not pass the D-053 rule.

**Adoption rule (computable; `w1_streak_mask.adopt`).** The veto is adopted for the W1 limit only if K1 and K2
both hold. If adopted, the headline limit is the `--combine-only` combination of all twelve fields (the eleven of
`results/w1_shear_rxj_m0647/` plus MACS1423) run with the veto, leaving out any field that passes the D-053 rule
under the veto (such a field goes to `scripts/w1_peak_vet.py` / `/vet-candidate`, as before). If not adopted, the
eleven-field limits stay and the masks are reported as a failed approach.

**Reported regardless.** Per field: the number of status-2 trails, the share of used shape rows vetoed, S_max,
p_rot, b_max with and without the veto, and the injection efficiencies at the five D-053 masses (seed 50, 200 per
mass). The contact sheet of each mask over its downsampled mosaic is looked at before any conclusion.

**What changes the plan.** A field with more than 20 % of used rows vetoed (ASSUMPTION) is still combined but
named in RESULT. A field whose mosaic cannot be streamed is left out of the veto run and named.
